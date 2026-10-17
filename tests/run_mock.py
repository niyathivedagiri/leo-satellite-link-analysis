"""Run ntn_ground_station.lua outside CoppeliaSim against a fake 'sim' API.

    python tests/run_mock.py                 # smoke test: both scene modes, buttons pressed
    python tests/run_mock.py --log OUTDIR    # clean validation run (no buttons), writes
                                             # OUTDIR/ntn_log.csv + OUTDIR/ntn_settings.json
Needs: pip install lupa
"""
import os
import sys

import lupa.lua54 as lupa

BUTTONS = ((1500, 31), (1800, 31), (2200, 30), (2600, 31), (3000, 32))


def run(mode, steps=12000, buttons=BUTTONS, home='/tmp'):
    L = lupa.LuaRuntime(unpack_returned_tuples=True)
    L.execute("package.path='./?.lua;'..package.path")
    btn = ' '.join(f"if i == {i} then onButton(nil, {b}) end" for i, b in buttons)
    L.execute(f"""
    dofile('tests/mock_sim.lua')
    local s = require('sim')
    local MODE = '{mode}'
    rawset(s, 'getObject', function(p, opt)
        if p == '/DefaultCamera' then return 999 end
        if MODE == 'handbuilt' and (p == '/PanJoint' or p == '/TiltJoint' or p == '/Satellite' or p == '/StatusLED') then return 500 + #p end
        if opt and opt.noError then return -1 end
        error('object not found ' .. p)
    end)
    rawset(s, 'getObjectPosition', function(h, rel) return {{0, 0, 2.3}} end)
    os.getenv = function(k) if k == 'HOME' then return '{home}' end return nil end
    OUT = {{}}
    local rp = print
    print = function(x) OUT[#OUT+1] = x end
    sim = nil
    dofile('ntn_ground_station.lua')
    sysCall_init()
    for i = 1, {steps} do
        sysCall_actuation()
        {btn}
    end
    sysCall_cleanup()
    print = rp
    """)
    g = L.globals()
    return list(g.OUT.values()), dict(g.MOCK_CALLS), dict(g.MOCK_LABELS)


if __name__ == '__main__':
    if len(sys.argv) > 2 and sys.argv[1] == '--log':
        out_dir = os.path.abspath(sys.argv[2])
        os.makedirs(out_dir, exist_ok=True)
        out, _, _ = run('auto', steps=12000, buttons=(), home=out_dir)
        for o in out:
            print('  ', o)
        print('log written to', out_dir)
    else:
        for mode in ('auto', 'handbuilt'):
            out, calls, labels = run(mode)
            print('==', mode, calls)
            for o in out[:30]:
                print('  ', o)
