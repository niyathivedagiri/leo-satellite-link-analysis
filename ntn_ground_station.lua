--[[ =====================================================================
   5G NTN LEO SATELLITE TRACKING TERMINAL  -  CoppeliaSim simulation script
   =====================================================================
   A ground-station dish follows a LEO satellite across the sky, keeps a
   5G NTN radio link to it, and uses the link quality to correct its own
   pointing.

   HOW TO USE
   1. Build the scene by hand (see the guide), or skip that: if the script
      cannot find an object called "PanJoint", it builds the whole scene itself.
   2. Add -> Script -> Simulation script -> Non-threaded -> Lua.
      Paste this whole file into it. Press Play.

   SECTIONS
     0. Settings              5. State machine
     1. Small helpers         6. Scene (find or build)
     2. Satellite orbit       7. Visuals
     3. Radio link            8. Telemetry panel
     4. Antenna control       9. Main loop
======================================================================== ]]

do local ok, m = pcall(require, 'sim'); if ok then sim = m end end
local okUI, simUI = pcall(require, 'simUI')

-------------------------------------------------------------------------
-- 0. SETTINGS
-------------------------------------------------------------------------
local SET = {
    satName = 'NTN-LEO-01', altitudeKm = 600, inclinationDeg = 25, raanDeg = 300, startDeg = 60,
    stationLat = 12.97, stationLon = 77.59,                  -- Bengaluru

    freqGHz = 20, bandwidthMHz = 400,                        -- Ka-band 5G NTN downlink
    eirpDbm = 60, dishGainDbi = 39.7, noiseDbm = -88.4,      -- noise + interference over 400 MHz
    beamwidthDeg = 1.75,                                     -- 0.6 m dish at 20 GHz

    -- Phase 1 (validation): no hidden errors and no step-track, so pointing error = tracking only.
    -- Phase 4 demo: set biasElDeg = 1.0 to see a 1-degree pointing error (beam misses, LINK DEGRADED).
    -- Old demo values were biasAzDeg = 0.6, biasElDeg = -0.4, linkAssist = true.
    biasAzDeg = 0, biasElDeg = 0,                            -- dish installation error (truth only)
    linkAssist = false,                                      -- step-track correction at start
    windDeg = 1.5,                                           -- push from the "Wind gust" button
    slewDegPerSec = 6,                                       -- motor speed limit
    probeDeg = 0.2,                                          -- step-track probe size

    minElevation = 10,                                       -- connect only above this
    detectDb = -7, connectDb = -2, lostDb = -5,              -- SINR thresholds
    passSpeed = 5,                                           -- real seconds per sim second in a pass
    skyRadius = 6.5,                                         -- size of the sky dome in the scene
    beamZoom = 5,                                            -- beam drawn 5x wider so errors are visible
    panSign = 1, tiltSign = 1,                               -- set to -1 if your hand-built joint turns the wrong way
}

-------------------------------------------------------------------------
-- 1. SMALL HELPERS (vectors, angles, positions)
-------------------------------------------------------------------------
local rad, deg = math.rad, math.deg
local function add(a, b) return {a[1] + b[1], a[2] + b[2], a[3] + b[3]} end
local function sub(a, b) return {a[1] - b[1], a[2] - b[2], a[3] - b[3]} end
local function mul(a, s) return {a[1] * s, a[2] * s, a[3] * s} end
local function dot(a, b) return a[1] * b[1] + a[2] * b[2] + a[3] * b[3] end
local function cross(a, b) return {a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1]} end
local function unit(a) local l = math.sqrt(dot(a, a)); return l > 1e-9 and mul(a, 1 / l) or {0, 0, 1} end
local function clamp(x, lo, hi) return math.max(lo, math.min(hi, x)) end
local function wrap180(a) return (a + 180) % 360 - 180 end

-- direction in the scene for an azimuth / elevation (X = East, Y = North, Z = Up)
local function dirOf(az, el)
    return {math.sin(rad(az)) * math.cos(rad(el)), math.cos(rad(az)) * math.cos(rad(el)), math.sin(rad(el))}
end

-- orientation (quaternion) that makes an object's Z axis point along 'dir'
local function lookAlong(dir, hint)
    local z = unit(dir)
    local x = hint and sub(hint, mul(z, dot(hint, z))) or cross({0, 0, 1}, z)   -- sideways axis
    if dot(x, x) < 1e-6 then x = cross({0, 1, 0}, z) end
    x = unit(x)
    local y = cross(z, x)
    local t = x[1] + y[2] + z[3]
    if t > 0 then local s = math.sqrt(t + 1) * 2
        return {(y[3] - z[2]) / s, (z[1] - x[3]) / s, (x[2] - y[1]) / s, s / 4}
    elseif x[1] > y[2] and x[1] > z[3] then local s = math.sqrt(1 + x[1] - y[2] - z[3]) * 2
        return {s / 4, (y[1] + x[2]) / s, (z[1] + x[3]) / s, (y[3] - z[2]) / s}
    elseif y[2] > z[3] then local s = math.sqrt(1 + y[2] - x[1] - z[3]) * 2
        return {(y[1] + x[2]) / s, s / 4, (z[2] + y[3]) / s, (z[1] - x[3]) / s}
    else local s = math.sqrt(1 + z[3] - x[1] - y[2]) * 2
        return {(z[1] + x[3]) / s, (z[2] + y[3]) / s, s / 4, (x[2] - y[1]) / s} end
end

local function place(h, pos, quat)
    local pose = {pos[1], pos[2], pos[3], quat[1], quat[2], quat[3], quat[4]}
    if not pcall(sim.setObjectPose, h, pose, sim.handle_world) then sim.setObjectPose(h, sim.handle_world, pose) end
end

local function paint(h, rgb, glow)
    sim.setShapeColor(h, nil, sim.colorcomponent_ambient_diffuse, rgb)
    sim.setShapeColor(h, nil, sim.colorcomponent_emission, glow or {0, 0, 0})
end

-------------------------------------------------------------------------
-- 2. SATELLITE ORBIT: where is the satellite, seen from the station?
-------------------------------------------------------------------------
local EARTH_R, MU, EARTH_SPIN, LIGHT = 6371, 398600.4, 7.2921e-5, 299792.458
local orbitR = EARTH_R + SET.altitudeKm
local orbitRate = math.sqrt(MU / orbitR ^ 3)            -- radians per second (Kepler's 3rd law)
local stn = {EARTH_R * math.cos(rad(SET.stationLat)) * math.cos(rad(SET.stationLon)),
             EARTH_R * math.cos(rad(SET.stationLat)) * math.sin(rad(SET.stationLon)),
             EARTH_R * math.sin(rad(SET.stationLat))}

-- returns azimuth, elevation (degrees) and range (km) at time t (seconds)
local function lookAngles(t)
    local u = rad(SET.startDeg) + orbitRate * t                     -- position along the orbit
    local px, py = orbitR * math.cos(u), orbitR * math.sin(u)
    local i, O = rad(SET.inclinationDeg), rad(SET.raanDeg)
    local x = px * math.cos(O) - py * math.cos(i) * math.sin(O)       -- tilt and turn the orbit
    local y = px * math.sin(O) + py * math.cos(i) * math.cos(O)
    local z = py * math.sin(i)
    local spin = EARTH_SPIN * t                                       -- the Earth turns underneath
    local X, Y = math.cos(spin) * x + math.sin(spin) * y, -math.sin(spin) * x + math.cos(spin) * y
    local d = {X - stn[1], Y - stn[2], z - stn[3]}                    -- station -> satellite
    local la, lo = rad(SET.stationLat), rad(SET.stationLon)
    local east = -math.sin(lo) * d[1] + math.cos(lo) * d[2]
    local north = -math.sin(la) * math.cos(lo) * d[1] - math.sin(la) * math.sin(lo) * d[2] + math.cos(la) * d[3]
    local up = math.cos(la) * math.cos(lo) * d[1] + math.cos(la) * math.sin(lo) * d[2] + math.sin(la) * d[3]
    local range = math.sqrt(dot(d, d))
    return deg(math.atan(east, north)) % 360, deg(math.asin(up / range)), range
end

-- next pass: step forward until the satellite rises, then follow it until it sets
local function predictPass(t)
    while select(2, lookAngles(t)) > 0 do t = t + 10 end            -- skip a pass in progress
    while select(2, lookAngles(t)) <= 0 do t = t + 10 end
    local pass = {rise = t, path = {}, maxEl = 0}
    pass.riseAz = lookAngles(t)
    while true do
        local az, el = lookAngles(t)
        if el <= 0 then break end
        pass.path[#pass.path + 1] = {t = t, az = az, el = el}
        if el > pass.maxEl then pass.maxEl, pass.topAz = el, az end
        t = t + 5
    end
    return pass
end

-------------------------------------------------------------------------
-- 3. RADIO LINK: how good is the 5G connection right now?
-------------------------------------------------------------------------
local function radioLink(range, el, pointingErr)
    local L = {}
    L.pathLoss = 92.45 + 20 * math.log(range, 10) + 20 * math.log(SET.freqGHz, 10)   -- free-space loss
    L.airLoss = 0.5 / math.sin(rad(math.max(el, 5)))                                 -- more air near the horizon
    L.pointLoss = math.min(12.04 * (pointingErr / SET.beamwidthDeg) ^ 2, 30)         -- Gaussian beam: 3 dB at half beamwidth
    L.power = SET.eirpDbm + SET.dishGainDbi - L.pathLoss - L.airLoss - L.pointLoss - 0.5
    if el <= 0 then L.power = -200 end                                               -- Earth is in the way
    L.rsrp = L.power - 35.0                                                          -- per resource element
    L.sinr = L.power - SET.noiseDbm
    local eff = L.sinr < -10 and 0 or math.min(0.6 * math.log(1 + 10 ^ (L.sinr / 10), 2), 4.4)
    L.throughput = eff * SET.bandwidthMHz                                            -- Mbps
    L.delay = range / LIGHT * 1000                                                   -- ms, one way
    return L
end

-------------------------------------------------------------------------
-- 4. ANTENNA CONTROL: predicted direction + learned correction (step-track)
-------------------------------------------------------------------------
local ant = {az = 0, el = 90, corrX = 0, corrY = 0, probe = 0, samples = {}, assist = SET.linkAssist, wind = 0}
local PROBES = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}}         -- right, left, up, down

-- step-track: try 4 small offsets, remember the signal at each, then move toward the strongest
local function stepTrack(lastRsrp)
    if ant.probe > 0 then ant.samples[ant.probe] = lastRsrp end
    local s = ant.samples
    if s[1] and s[2] and s[3] and s[4] then
        local k = 0.4 * SET.beamwidthDeg ^ 2 / 24 / (2 * SET.probeDeg)
        ant.corrX = ant.corrX + clamp(k * (s[1] - s[2]), -0.08, 0.08)
        ant.corrY = ant.corrY + clamp(k * (s[3] - s[4]), -0.08, 0.08)
        ant.samples = {}
    end
    ant.probe = ant.probe % 4 + 1
    return PROBES[ant.probe][1] * SET.probeDeg, PROBES[ant.probe][2] * SET.probeDeg
end

-- move the dish toward a target, never faster than the motors allow
local function moveDish(targetAz, targetEl, dt)
    local step = SET.slewDegPerSec * dt
    ant.az = (ant.az + clamp(wrap180(targetAz - ant.az), -step, step)) % 360
    ant.el = ant.el + clamp(targetEl - ant.el, -step, step)
end

-------------------------------------------------------------------------
-- 5. STATE MACHINE
-------------------------------------------------------------------------
local LABEL = {WAIT = 'WAIT', PREPOINT = 'PRE-POINT', ACQUIRE = 'ACQUIRE', TRACK = 'TRACK',
               CONNECTED = 'CONNECTED', DEGRADED = 'LINK DEGRADED', LOST = 'LINK LOST', PARK = 'PARK'}
local COLOR = {WAIT = {0.45, 0.48, 0.55}, PREPOINT = {0.55, 0.65, 0.95}, ACQUIRE = {0.95, 0.85, 0.35},
               TRACK = {1.0, 0.70, 0.15}, CONNECTED = {0.15, 0.95, 0.55}, DEGRADED = {1.0, 0.50, 0.10},
               LOST = {0.95, 0.20, 0.20}, PARK = {0.35, 0.55, 1.0}}
local onSky = {ACQUIRE = true, TRACK = true, CONNECTED = true, DEGRADED = true, LOST = true}

local function nextState(s, x)   -- x = {el, rising, sinr, shortfall, toRise, parked, age}
    local setting = x.el <= 0 or (not x.rising and x.el < SET.minElevation)
    local settled = x.age >= 2                       -- link states last at least 2 s (no flicker)
    if s == 'WAIT' and x.toRise <= 60 then return 'PREPOINT' end
    if s == 'PREPOINT' and x.el > 0 then return 'ACQUIRE' end
    if s == 'PARK' and x.parked then return 'WAIT' end
    if onSky[s] and setting then return 'PARK' end
    if not settled then return s end
    if s == 'ACQUIRE' and x.sinr >= SET.detectDb then return 'TRACK' end
    if s == 'TRACK' and x.sinr >= SET.connectDb and x.el >= SET.minElevation then return 'CONNECTED' end
    if (s == 'CONNECTED' or s == 'DEGRADED') and x.sinr < SET.lostDb then return 'LOST' end
    if s == 'TRACK' and x.sinr < SET.detectDb and x.el >= SET.minElevation then return 'LOST' end
    if s == 'CONNECTED' and x.shortfall > 3 then return 'DEGRADED' end
    if s == 'DEGRADED' and x.shortfall < 1.5 then return 'CONNECTED' end
    if s == 'LOST' and x.sinr >= SET.detectDb then return 'TRACK' end
    return s
end

-------------------------------------------------------------------------
-- 6. SCENE: use the hand-built objects, or build them if they are missing
-------------------------------------------------------------------------
local H = {}            -- handles of scene objects
local made = {}         -- objects this script created (removed again at Stop)
local PIVOT, FEED = {0, 0, 2.30}, 0.74

local function find(name)
    local ok, h = pcall(sim.getObject, '/' .. name, {noError = true})
    return (ok and h and h >= 0) and h or nil
end

local function shape(kind, size, pos, rgb, parent, axis, glow)
    local t = ({box = sim.primitiveshape_cuboid, cyl = sim.primitiveshape_cylinder, ball = sim.primitiveshape_spheroid})[kind]
    local h = sim.createPrimitiveShape(t, size, 0)
    made[#made + 1] = h
    local q = {0, 0, 0, 1}
    if axis == 'x' then q = {0, 0.7071068, 0, 0.7071068} elseif axis == 'y' then q = {-0.7071068, 0, 0, 0.7071068} end
    place(h, pos, q)
    paint(h, rgb, glow)
    if parent then sim.setObjectParent(h, parent, true) end
    return h
end

local function rod(a, b, r, rgb, parent)        -- cylinder from point a to point b
    local h = sim.createPrimitiveShape(sim.primitiveshape_cylinder, {2 * r, 2 * r, math.sqrt(dot(sub(b, a), sub(b, a)))}, 0)
    made[#made + 1] = h
    place(h, mul(add(a, b), 0.5), lookAlong(sub(b, a), {1, 0, 0}))
    paint(h, rgb)
    if parent then sim.setObjectParent(h, parent, true) end
end

-- a round surface made by spinning a profile {{radius, distance}, ...} around the dish axis (+Y)
local function spun(profile, rgb, parent)
    local v, idx, n = {}, {}, 48
    for _, p in ipairs(profile) do
        for k = 0, n - 1 do
            local a = 2 * math.pi * k / n
            v[#v + 1] = PIVOT[1] + p[1] * math.cos(a); v[#v + 1] = PIVOT[2] + p[2]; v[#v + 1] = PIVOT[3] + p[1] * math.sin(a)
        end
    end
    for i = 0, #profile - 2 do
        for k = 0, n - 1 do
            local a, b = i * n + k, i * n + (k + 1) % n
            idx[#idx + 1] = a; idx[#idx + 1] = a + n; idx[#idx + 1] = b
            idx[#idx + 1] = b; idx[#idx + 1] = a + n; idx[#idx + 1] = b + n
        end
    end
    local ok, h = pcall(sim.createShape, 0, 0.5, v, idx)
    if not ok then h = sim.createMeshShape(0, 0.5, v, idx) end
    made[#made + 1] = h
    paint(h, rgb)
    if parent then sim.setObjectParent(h, parent, true) end
    return h
end

local function buildScene()
    local grey, steel, dark, white = {0.62, 0.62, 0.60}, {0.78, 0.80, 0.83}, {0.28, 0.30, 0.33}, {0.95, 0.96, 0.97}
    -- site
    shape('cyl', {28, 28, 0.02}, {0, 0, 0}, {0.20, 0.22, 0.23})                         -- gravel apron
    shape('box', {4.2, 4.2, 0.16}, {0, 0, 0.08}, grey)                                  -- concrete pad
    shape('box', {1.5, 1.0, 1.15}, {2.85, -2.1, 0.74}, {0.88, 0.89, 0.86})              -- equipment shelter
    shape('box', {1.62, 1.12, 0.06}, {2.85, -2.1, 1.345}, dark)                         -- shelter roof
    shape('box', {0.55, 0.02, 0.85}, {2.75, -2.61, 0.62}, dark)                         -- shelter door
    -- pedestal + PAN joint (azimuth)
    shape('cyl', {0.9, 0.9, 0.08}, {0, 0, 0.20}, dark)
    shape('cyl', {0.40, 0.40, 1.16}, {0, 0, 0.82}, steel)
    H.pan = sim.createJoint(sim.joint_revolute or sim.joint_revolute_subtype, sim.jointmode_kinematic or sim.jointmode_passive, 0, {0.2, 0.04})
    made[#made + 1] = H.pan
    place(H.pan, {0, 0, 1.40}, {0, 0, 0, 1})
    -- azimuth structure (turns with the pan joint)
    shape('cyl', {0.66, 0.66, 0.16}, {0, 0, 1.48}, dark, H.pan)                         -- azimuth bearing
    shape('box', {1.62, 0.44, 0.08}, {0, 0, 1.60}, steel, H.pan)                        -- turret
    shape('box', {0.22, 0.20, 0.22}, {0, -0.34, 1.50}, dark, H.pan)                     -- azimuth motor
    shape('box', {0.10, 0.30, 0.72}, {0.76, 0, 1.99}, steel, H.pan)                     -- yoke arm east
    shape('box', {0.10, 0.30, 0.72}, {-0.76, 0, 1.99}, steel, H.pan)                    -- yoke arm west
    shape('cyl', {0.22, 0.22, 0.16}, {0.82, 0, 2.30}, dark, H.pan, 'x')                 -- elevation bearings
    shape('cyl', {0.22, 0.22, 0.16}, {-0.82, 0, 2.30}, dark, H.pan, 'x')
    shape('box', {0.16, 0.22, 0.24}, {0.90, 0, 2.00}, dark, H.pan)                      -- elevation motor
    H.led = shape('ball', {0.07, 0.07, 0.07}, {0.76, -0.17, 2.25}, COLOR.WAIT, H.pan)
    -- TILT joint (elevation), axis East-West
    H.tilt = sim.createJoint(sim.joint_revolute or sim.joint_revolute_subtype, sim.jointmode_kinematic or sim.jointmode_passive, 0, {0.2, 0.04})
    made[#made + 1] = H.tilt
    place(H.tilt, PIVOT, {0, 0.7071068, 0, 0.7071068})
    sim.setObjectParent(H.tilt, H.pan, true)
    -- dish (built facing North), feed horn at the focus, struts, counterweight
    rod({-0.74, 0, 2.30}, {0.74, 0, 2.30}, 0.05, dark, H.tilt)                          -- elevation axle
    local dish = {}
    for j = 0, 12 do local r = 0.62 * j / 12; dish[#dish + 1] = {r, 0.14 + r * r / 2.4} end
    spun(dish, white, H.tilt)                                                           -- parabolic reflector
    spun({{0.17, 0.04}, {0.57, 0.26}}, {0.86, 0.87, 0.88}, H.tilt)                     -- back structure
    spun({{0.065, 0.72}, {0.03, 0.85}}, dark, H.tilt)                                   -- feed horn
    shape('box', {0.09, 0.14, 0.09}, {0, 0.92, 2.30}, {0.86, 0.87, 0.88}, H.tilt)       -- LNB (receiver)
    shape('box', {0.40, 0.14, 0.26}, {0, -0.26, 2.25}, dark, H.tilt)                    -- counterweight
    for _, a in ipairs({45, 135, 225, 315}) do
        local c, s = math.cos(rad(a)), math.sin(rad(a))
        rod({0.6 * c, 0.29, 2.30 + 0.6 * s}, {0.05 * c, 0.80, 2.30 + 0.05 * s}, 0.011, {0.55, 0.57, 0.60}, H.tilt)
    end
    -- satellite (built at the origin, moved every step): body, solar wings, ground antenna
    local o = {0, 0, -100}
    H.sat = shape('box', {0.34, 0.34, 0.46}, o, {0.83, 0.64, 0.20}, nil, nil, {0.06, 0.04, 0})
    for _, sg in ipairs({1, -1}) do
        rod(add(o, {sg * 0.17, 0, 0.05}), add(o, {sg * 0.32, 0, 0.05}), 0.018, {0.55, 0.57, 0.60}, H.sat)
        for k = 0, 2 do shape('box', {0.29, 0.44, 0.015}, add(o, {sg * (0.48 + 0.31 * k), 0, 0.05}), {0.07, 0.12, 0.30}, H.sat, nil, {0.01, 0.03, 0.10}) end
    end
    shape('box', {0.26, 0.26, 0.02}, add(o, {0, 0, -0.24}), dark, H.sat)               -- antenna panel (faces Earth)
    H.beacon = shape('ball', {0.06, 0.06, 0.06}, add(o, {0, 0, 0.27}), {0.85, 0.12, 0.10}, H.sat)
end

local function setupScene()
    H.pan, H.tilt, H.sat = find('PanJoint'), find('TiltJoint'), find('Satellite')
    if H.pan and H.tilt and H.sat then                                    -- hand-built scene
        H.led, H.beacon = find('StatusLED'), find('Beacon')
        local ok, p = pcall(sim.getObjectPosition, H.tilt, sim.handle_world)
        if ok then PIVOT = p end
    else
        buildScene()
    end
    -- Phase 1 check: a marker that rides on the dish. With both joints at 0 the dish faces
    -- North (+Y), so the marker's own Y axis IS the dish boresight. Every step we read where
    -- the scene really points it and compare with where the maths says it points.
    sim.setJointPosition(H.pan, 0)
    sim.setJointPosition(H.tilt, 0)
    H.bore = sim.createDummy(0.03)
    made[#made + 1] = H.bore
    place(H.bore, PIVOT, {0, 0, 0, 1})
    sim.setObjectParent(H.bore, H.tilt, true)
    -- two small camera views: through the dish, and a close-up of the antenna
    H.views = {}
    local function camera(eye, look, fov, parent, x, y, title)
        local cam = sim.createVisionSensor(2 + 4, {320, 240, 0, 0}, {0.9, 80, rad(fov), 0.05, 0, 0, 0, 0, 0, 0, 0})
        made[#made + 1] = cam
        place(cam, eye, lookAlong(sub(look, eye)))
        if parent then sim.setObjectParent(cam, parent, true) end
        local ok, v = pcall(sim.floatingViewAdd, x, y, 0.24, 0.27, 0)
        if ok and v >= 0 then H.views[#H.views + 1] = v; pcall(sim.adjustView, v, cam, 0, title) end
    end
    camera(add(PIVOT, {0, 0.16, 0}), add(PIVOT, {0, 5, 0}), 14, H.tilt, 0.86, 0.46, 'Through the dish')
    camera({3.6, -3.4, 2.9}, add(PIVOT, {0, 0, -0.25}), 50, nil, 0.86, 0.17, 'Antenna close-up')
    local ok, cam = pcall(sim.getObject, '/DefaultCamera')
    if ok and cam >= 0 then H.cam, H.camPose = cam, sim.getObjectPose(cam, sim.handle_world) end
end

-------------------------------------------------------------------------
-- 7. VISUALS: sky, path, trail, beam, data packets, lights, camera
-------------------------------------------------------------------------
local D = {}
local function draw(kind, size, rgb, max)
    local h = (sim.addDrawingObject or sim.createDrawingObject)(kind + (sim.drawing_emissioncolor or 0), size, 0, -1, max or 1000, rgb)
    D[#D + 1] = h
    return h
end
local function line(h, a, b) sim.addDrawingObjectItem(h, {a[1], a[2], a[3], b[1], b[2], b[3]}) end
local function dot3(h, p) sim.addDrawingObjectItem(h, p) end
local function clear(h) sim.addDrawingObjectItem(h, nil) end
local function sky(az, el, r) return add(PIVOT, mul(dirOf(az, el), r or SET.skyRadius)) end

local BEAM = {ACQUIRE = {{0.55, 0.50, 0.22}, 1}, TRACK = {{0.70, 0.50, 0.12}, 1},
              CONNECTED = {{0.15, 0.95, 0.55}, 2}, DEGRADED = {{1.0, 0.50, 0.10}, 1}}

local function setupVisuals()
    pcall(function()                                                  -- dusk sky background
        H.bg = {sim.getArrayParam(sim.arrayparam_background_color1), sim.getArrayParam(sim.arrayparam_background_color2)}
        sim.setArrayParam(sim.arrayparam_background_color1, {0.10, 0.14, 0.22})
        sim.setArrayParam(sim.arrayparam_background_color2, {0.01, 0.02, 0.06})
    end)
    local stars = draw(sim.drawing_points, 2, {0.85, 0.88, 1.0})
    for k = 1, 300 do dot3(stars, sky((k * 137.5) % 360, deg(math.asin(0.08 + 0.92 * ((k * 0.618) % 1))), 40)) end
    local grid = draw(sim.drawing_lines, 1, {0.13, 0.20, 0.32}, 2000)  -- faint sky grid
    for az = 0, 330, 30 do for el = 0, 75, 5 do line(grid, sky(az, el), sky(az, el + 5)) end end
    for _, el in ipairs({30, 60}) do for az = 0, 355, 5 do line(grid, sky(az, el), sky(az + 5, el)) end end
    local horizon = draw(sim.drawing_lines, 2, {0.55, 0.62, 0.75})
    local range = draw(sim.drawing_lines, 2, {0.95, 0.60, 0.15})      -- 10 deg: communication range
    for az = 0, 355, 5 do
        line(horizon, sky(az, 0), sky(az + 5, 0))
        if az % 10 == 0 then line(range, sky(az, SET.minElevation), sky(az + 5, SET.minElevation)) end
    end
    D.future = draw(sim.drawing_spherepts, 0.045, {0.55, 0.70, 0.95})                    -- predicted path
    D.trail = draw(sim.drawing_lines + sim.drawing_cyclic, 3, {0.30, 0.85, 1.0}, 600)    -- where it has been
    D.sight = draw(sim.drawing_lines, 1, {0.90, 0.92, 0.95})                             -- line of sight
    D.packets = draw(sim.drawing_spherepts, 0.07, {0.6, 1.0, 0.75})
    D.beam = {}
    for s, b in pairs(BEAM) do D.beam[s] = draw(sim.drawing_lines, b[2], b[1], 300) end
end

local function frameCamera(pass)                 -- stand behind the station, facing where the pass peaks
    if not H.cam then return end
    local f = dirOf(pass.topAz or 0, 0)
    local eye = add(PIVOT, {-f[1] * 10.5, -f[2] * 10.5, 3.4})
    place(H.cam, eye, lookAlong(sub(add(PIVOT, {f[1] * 2, f[2] * 2, 1.9}), eye)))
end

local lastSatPos, clock = nil, 0
local function updateVisuals(s, dt)
    clock = clock + dt
    if H.led then paint(H.led, COLOR[s.state], COLOR[s.state]) end
    local satPos = sky(s.az, s.el)

    -- satellite: on the sky dome when above the horizon, antenna facing the station
    if s.el > 0 then
        local aheadAz, aheadEl = lookAngles(s.t + 5)
        local ahead = sub(sky(aheadAz, aheadEl), satPos)
        local up = dirOf(s.az, s.el)
        place(H.sat, satPos, lookAlong(up, cross(up, ahead)))
        if lastSatPos then line(D.trail, lastSatPos, satPos) end
        lastSatPos = satPos
        if H.beacon then paint(H.beacon, {0.85, 0.12, 0.10}, clock % 1 < 0.2 and {0.85, 0.12, 0.10} or {0, 0, 0}) end
    else
        place(H.sat, {0, 0, -100}, {0, 0, 0, 1})
        lastSatPos = nil
    end

    -- predicted path: only the part still ahead of the satellite
    clear(D.future)
    if s.pass and s.state ~= 'PARK' then
        for _, p in ipairs(s.pass.path) do if p.t > s.t then dot3(D.future, sky(p.az, p.el)) end end
    end

    -- line of sight, beam and data packets
    clear(D.sight); clear(D.packets)
    for _, h in pairs(D.beam) do clear(h) end
    local feed = add(PIVOT, mul(s.dishDir, FEED))
    if s.el > 0 and onSky[s.state] then
        local seg = sub(satPos, feed)
        for k = 0, 23, 2 do line(D.sight, add(feed, mul(seg, k / 24)), add(feed, mul(seg, (k + 1) / 24))) end
    end
    local h = D.beam[s.state]
    if s.el > 0 and h then
        -- beam points where the dish REALLY points; the error is drawn beamZoom times bigger
        local dir = unit(add(s.satDir, mul(sub(s.dishDir, s.satDir), SET.beamZoom)))
        local len = SET.skyRadius - FEED
        local r = len * math.tan(rad(SET.beamZoom * SET.beamwidthDeg / 2))
        local ex = unit(cross({0, 0, 1}, dir)); local ey = cross(dir, ex)
        local function edge(f, a) return add(add(feed, mul(dir, len * f)), add(mul(ex, r * f * math.cos(a)), mul(ey, r * f * math.sin(a)))) end
        for k = 0, 15 do line(h, feed, edge(1, k * math.pi / 8)) end
        for _, f in ipairs({0.35, 0.7, 1}) do for k = 0, 31 do line(h, edge(f, k * math.pi / 16), edge(f, (k + 1) * math.pi / 16)) end end
        if s.state == 'CONNECTED' or s.state == 'DEGRADED' then
            local n = clamp(math.floor(s.link.throughput / 120), 2, 10)
            for i = 0, n - 1 do dot3(D.packets, add(satPos, mul(sub(feed, satPos), (clock * (0.25 + s.link.throughput / 2000) + i / n) % 1))) end
        end
    end
    if s.state == 'PREPOINT' and s.changed then frameCamera(s.pass) end
end

-------------------------------------------------------------------------
-- 8. TELEMETRY PANEL (small window, top-left)
-------------------------------------------------------------------------
local ui
local ROWS = {'SATELLITE', 'AZ / EL', 'RANGE', 'DELAY', 'DOPPLER', 'POINTING ERROR', 'RSRP', 'SINR', 'THROUGHPUT', 'CLOCK'}

local function setupPanel()
    if not okUI then return end
    local rows = ''
    for i, name in ipairs(ROWS) do
        rows = rows .. string.format('<label text="%s" style="color: rgb(148,163,184); font-size: 11px;"/><label id="%d" text="-" style="font-family: Menlo, Consolas, monospace;"/>', name, 10 + i)
    end
    ui = simUI.create([[
<ui title="5G NTN TERMINAL - TELEMETRY" closeable="false" resizable="false" placement="relative" position="12,12" layout="vbox" activate="false"
    style="QWidget { background-color: rgb(13,19,31); color: rgb(226,232,240); font-size: 12px; } QPushButton { background-color: rgb(30,41,59); border: 1px solid rgb(71,85,105); padding: 3px 7px; }">
  <label id="1" text="WAIT"/>
  <group layout="form" flat="true">]] .. rows .. [[</group>
  <group layout="hbox" flat="true">
    <button id="30" text="Link assist: ]] .. (ant.assist and 'ON' or 'OFF') .. [[" on-click="onButton"/>
    <button id="31" text="Wind gust" on-click="onButton"/>
    <button id="32" text="Speed: x5" on-click="onButton"/>
  </group>
</ui>]])
end

local function updatePanel(s)
    if not ui then return end
    local c = COLOR[s.state]
    local hex = string.format('#%02x%02x%02x', math.floor(c[1] * 255), math.floor(c[2] * 255), math.floor(c[3] * 255))
    local L, up = s.link, s.el > 0
    local values = {
        string.format('%s  (%d km)', SET.satName, SET.altitudeKm),
        up and string.format('%6.1f° / %5.1f°', s.az, s.el) or 'below horizon',
        string.format('%6.0f km', s.range),
        string.format('%5.2f ms', L.delay),
        string.format('%+7.1f kHz', s.doppler),
        onSky[s.state] and string.format('%5.2f°', s.err) or '-',
        up and string.format('%7.1f dBm', L.rsrp) or '-',
        up and string.format('%5.1f dB', L.sinr) or '-',
        string.format('%6.0f Mbps', s.thr),
        string.format('T+%02d:%02d:%02d  x%d', s.t // 3600, s.t // 60 % 60, math.floor(s.t % 60), math.floor(s.speed)),
    }
    simUI.setLabelText(ui, 1, string.format('<span style="font-size: 17px; font-weight: 700; color: %s;">&#9679; %s</span>', hex, LABEL[s.state]))
    for i, v in ipairs(values) do simUI.setLabelText(ui, 10 + i, v) end
end

local speeds, speedIx, passSpeed = {2, 5, 10, 25}, 2, SET.passSpeed
function onButton(_, id)
    if id == 30 then
        ant.assist = not ant.assist
        if not ant.assist then ant.corrX, ant.corrY, ant.samples, ant.probe = 0, 0, {}, 0 end
        simUI.setButtonText(ui, 30, 'Link assist: ' .. (ant.assist and 'ON' or 'OFF'))
    elseif id == 31 then
        ant.wind = ant.wind > 0 and 0 or SET.windDeg
        simUI.setButtonText(ui, 31, ant.wind > 0 and 'Wind gust: ON' or 'Wind gust')
    elseif id == 32 then
        speedIx = speedIx % #speeds + 1
        passSpeed = speeds[speedIx]
        simUI.setButtonText(ui, 32, 'Speed: x' .. passSpeed)
    end
end

-------------------------------------------------------------------------
-- 9. MAIN LOOP
-------------------------------------------------------------------------
local t, state, pass, lastEl, lastRsrp, log, age = 0, 'WAIT', nil, -90, -150, nil, 0

function sysCall_init()
    setupScene()
    setupVisuals()
    setupPanel()
    pass = predictPass(0)
    frameCamera(pass)
    local home = os.getenv('HOME') or os.getenv('USERPROFILE') or '.'
    log = io.open(home .. '/ntn_log.csv', 'w')
    if log then
        log:write('time_s,state,sat_az_deg,sat_el_deg,range_km,range_rate_kms,dish_az_deg,dish_el_deg,'
               .. 'pointing_error_deg,scene_check_deg,doppler_khz,delay_ms,prx_dbm,rsrp_dbm,sinr_db,throughput_mbps\n')
        print('Logging to ' .. home .. '/ntn_log.csv')
    end
    -- the settings Python needs to re-check the run independently
    local f = io.open(home .. '/ntn_settings.json', 'w')
    if f then
        f:write(string.format('{"altitudeKm": %g, "inclinationDeg": %g, "raanDeg": %g, "startDeg": %g, '
            .. '"stationLat": %g, "stationLon": %g, "freqGHz": %g, "bandwidthMHz": %g, "eirpDbm": %g, '
            .. '"dishGainDbi": %g, "beamwidthDeg": %g, "slewDegPerSec": %g, "biasAzDeg": %g, "biasElDeg": %g, '
            .. '"linkAssist": %s, "minElevation": %g, "earthRadiusKm": %.10g, "mu": %.10g, "earthSpin": %.10g}\n',
            SET.altitudeKm, SET.inclinationDeg, SET.raanDeg, SET.startDeg, SET.stationLat, SET.stationLon,
            SET.freqGHz, SET.bandwidthMHz, SET.eirpDbm, SET.dishGainDbi, SET.beamwidthDeg, SET.slewDegPerSec,
            SET.biasAzDeg, SET.biasElDeg, tostring(SET.linkAssist), SET.minElevation, EARTH_R, MU, EARTH_SPIN))
        f:close()
    end
    print(string.format('NTN terminal ready. Next pass in %.1f min, max elevation %.0f deg', pass.rise / 60, pass.maxEl))
end

function sysCall_actuation()
    local dt = sim.getSimulationTimeStep()
    -- time: fast-forward while waiting, normal speed during a pass
    local speed = passSpeed
    if state == 'WAIT' then speed = clamp((pass.rise - 60 - t) / (20 * dt), passSpeed, 600) end
    local h = dt * speed
    t = t + h

    -- 1. where is the satellite?
    local az, el, range = lookAngles(t)
    local _, _, rangeBefore = lookAngles(t - 0.5)
    local _, _, rangeAfter = lookAngles(t + 0.5)
    local rangeRate = rangeAfter - rangeBefore                                  -- km/s (central difference)
    local doppler = -rangeRate / LIGHT * SET.freqGHz * 1e6                      -- kHz

    -- 2. where should the dish point?
    local tAz, tEl, px, py = ant.az, 90, 0, 0                                     -- default: park straight up
    if state == 'PREPOINT' then tAz, tEl = pass.riseAz, 1
    elseif onSky[state] then
        tAz, tEl = az, math.max(el, 0)                                            -- prediction
        if ant.assist and state ~= 'ACQUIRE' then px, py = stepTrack(lastRsrp) end
    end
    if px == 0 and py == 0 then ant.probe, ant.samples = 0, {} end
    if tEl < 90 then
        local c = math.max(math.cos(rad(tEl)), 0.17)
        tAz, tEl = tAz + (ant.corrX + px) / c, tEl + ant.corrY + py               -- add learned correction + probe
    end
    moveDish(tAz, tEl, h)

    -- 3. where does the dish REALLY point? (installation error + wind, unknown to the controller)
    local realAz = ant.az + (SET.biasAzDeg + ant.wind / math.max(math.cos(rad(ant.el)), 0.17))
    local realEl = ant.el + SET.biasElDeg
    local dishDir, satDir = dirOf(realAz, realEl), dirOf(az, el)
    local err = deg(math.acos(clamp(dot(dishDir, satDir), -1, 1)))
    sim.setJointPosition(H.pan, -SET.panSign * rad(realAz))
    sim.setJointPosition(H.tilt, SET.tiltSign * rad(realEl))
    -- Phase 1 check: where does the 3D dish really point? (column 2 of the marker's matrix = its Y axis)
    local m = sim.getObjectMatrix(H.bore, sim.handle_world)
    local sceneDir = unit({m[2], m[6], m[10]})
    local sceneErr = deg(math.acos(clamp(dot(sceneDir, dishDir), -1, 1)))     -- should be ~0

    -- 4. radio link, and how far it is below what a perfectly pointed dish would get
    local L = radioLink(range, el, err)
    local shortfall = radioLink(range, el, 0).sinr - L.sinr
    lastRsrp = L.rsrp

    -- 5. state machine
    local old = state
    state = nextState(state, {el = el, rising = el > lastEl, sinr = L.sinr, shortfall = shortfall,
                              toRise = pass.rise - t, parked = ant.el > 89.5, age = age})
    lastEl = el
    age = (state == old) and age + h or 0
    if state ~= old then
        print(string.format('[T+%.0f s] %s -> %s', t, LABEL[old], LABEL[state]))
        if state == 'WAIT' then pass = predictPass(t); lastSatPos = nil; clear(D.trail) end
    end

    -- 6. show it
    local thr = (state == 'CONNECTED' or state == 'DEGRADED') and L.throughput or 0
    local s = {t = t, state = state, changed = state ~= old, az = az, el = el, range = range, satDir = satDir,
               dishDir = dishDir, err = err, doppler = doppler, link = L, thr = thr, pass = pass, speed = speed}
    updateVisuals(s, dt)
    updatePanel(s)
    if log and state ~= 'WAIT' then
        log:write(string.format('%.4f,%s,%.5f,%.5f,%.4f,%.6f,%.5f,%.5f,%.6f,%.6f,%.4f,%.5f,%.4f,%.4f,%.4f,%.3f\n',
            t, state, az, el, range, rangeRate, realAz, realEl, err, sceneErr, doppler, L.delay,
            L.power, L.rsrp, L.sinr, thr))
    end
end

function sysCall_cleanup()
    if log then log:close() end
    if ui then simUI.destroy(ui) end
    for _, v in ipairs(H.views or {}) do pcall(sim.floatingViewRemove, v) end
    for _, d in ipairs(D) do pcall(sim.removeDrawingObject, d) end
    if H.cam and H.camPose then pcall(sim.setObjectPose, H.cam, H.camPose, sim.handle_world) end
    if H.bg then pcall(function()
        sim.setArrayParam(sim.arrayparam_background_color1, H.bg[1])
        sim.setArrayParam(sim.arrayparam_background_color2, H.bg[2]) end) end
    for i = #made, 1, -1 do pcall(sim.removeObjects, {made[i]}) end
end
