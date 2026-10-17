-- Fake CoppeliaSim 'sim' + 'simUI' modules that CHECK every call (types, sizes, NaN, unit quaternions).
local nextH, calls = 1000, {}
local function H() nextH = nextH + 1; return nextH end
local function n(name) calls[name] = (calls[name] or 0) + 1 end
local function finite(x, what) assert(type(x) == 'number' and x == x and x ~= math.huge and x ~= -math.huge, 'bad number in ' .. what) end
MOCK_CALLS, MOCK_LABELS = calls, {}
JOINTPOS, JORDER = {}, {}

local sim = setmetatable({
    handle_world = -1, primitiveshape_cuboid = 3, primitiveshape_cylinder = 2, primitiveshape_spheroid = 1,
    colorcomponent_ambient_diffuse = 0, colorcomponent_specular = 2, colorcomponent_emission = 3,
    colorcomponent_transparency = 4, joint_revolute = 10, jointmode_kinematic = 1,
    drawing_points = 0, drawing_lines = 1, drawing_spherepts = 7, drawing_cyclic = 0x2000, drawing_emissioncolor = 0x80,
    arrayparam_background_color1 = 1, arrayparam_background_color2 = 2,
    getSimulationTimeStep = function() return 0.05 end,
    createPrimitiveShape = function(t, s, o) assert(type(t) == 'number' and #s == 3); for _, x in ipairs(s) do finite(x, 'size'); assert(x > 0, 'size<=0') end; n('prim'); return H() end,
    createShape = function(o, a, v, i)
        assert(#v % 3 == 0 and #i % 3 == 0 and #i > 0)
        for _, x in ipairs(v) do finite(x, 'vertex') end
        for _, k in ipairs(i) do assert(math.type(k) == 'integer' and k >= 0 and k < #v / 3, 'bad index') end
        n('mesh'); return H() end,
    createJoint = function(t, m, o, s) assert(t == 10 and m == 1 and o == 0); n('joint'); return H() end,
    createDummy = function() return H() end,
    createVisionSensor = function(o, ip, fp) assert(#ip == 4 and #fp == 11); n('cam'); return H() end,
    setObjectPose = function(h, pose, rel)
        assert(type(pose) == 'table' and #pose == 7, 'pose must have 7 numbers')
        for k = 1, 7 do finite(pose[k], 'pose') end
        local q = math.sqrt(pose[4] ^ 2 + pose[5] ^ 2 + pose[6] ^ 2 + pose[7] ^ 2)
        assert(math.abs(q - 1) < 1e-6, 'quaternion not unit: ' .. q); n('pose') end,
    getObjectPose = function() return {0, 0, 0, 0, 0, 0, 1} end,
    getObject = function(p) assert(type(p) == 'string'); return H() end,
    setObjectParent = function(h, p, keep) assert(type(h) == 'number' and type(p) == 'number' and keep == true) end,
    setObjectAlias = function(h, s) assert(type(s) == 'string') end,
    setShapeColor = function(h, name, comp, rgb)
        assert(name == nil); assert(#rgb == (comp == 4 and 1 or 3), 'color size')
        for _, x in ipairs(rgb) do finite(x, 'color'); assert(x >= 0 and x <= 1.0001, 'color out of 0..1') end end,
    setBoolProperty = function(h, k, v) assert(type(k) == 'string' and type(v) == 'boolean') end,
    setJointPosition = function(h, a)
        finite(a, 'joint'); n('joint_set')
        if JOINTPOS[h] == nil then JORDER[#JORDER + 1] = h end
        JOINTPOS[h] = a end,
    -- the marker on the dish: world matrix = Rz(pan) * Rx(tilt), the joint chain the script assumes
    -- (pan about world Z, tilt about East). Plumbing test only; the real check runs in CoppeliaSim.
    getObjectMatrix = function(h, rel)
        n('matrix')
        local p, t = JOINTPOS[JORDER[1]] or 0, JOINTPOS[JORDER[2]] or 0
        local cp, sp, c, s = math.cos(p), math.sin(p), math.cos(t), math.sin(t)
        return {cp, -sp * c, sp * s, 0,  sp, cp * c, -cp * s, 0,  0, s, c, 0} end,
    addDrawingObject = function(t, s, d, p, m, c) assert(#c == 3); n('draw'); return H() end,
    addDrawingObjectItem = function(h, d) if d then for _, x in ipairs(d) do finite(x, 'drawing') end end end,
    removeDrawingObject = function() end,
    removeObjects = function(t) assert(type(t) == 'table') end,
    floatingViewAdd = function(x, y, w, hh, o) n('view'); return H() end,
    adjustView = function(v, c, o, label) assert(type(label) == 'string') end,
    floatingViewRemove = function() end,
    getArrayParam = function() return {0.5, 0.5, 0.5} end,
    setArrayParam = function(p, v) assert(#v == 3) end,
    getStringParam = function() return '.' end,
    stringparam_scene_path = 1,
}, {__index = function(_, k) error('sim.' .. k .. ' is not mocked') end})

package.preload['sim'] = function() return sim end
package.preload['simUI'] = function()
    return {create = function(xml) MOCK_UI_XML = xml; return 'ui1' end,
            setLabelText = function(u, id, t) assert(type(t) == 'string'); MOCK_LABELS[id] = t end,
            setButtonText = function(u, id, t) MOCK_LABELS['b' .. id] = t end,
            destroy = function() end}
end
return sim
