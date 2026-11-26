#!/usr/bin/env python3
# fix_rib_layers.py: rib zones of the inter-layer gap foam must not contain more foam than CAD MEC-301/304.
# MEC-301 in the rib zones (module windows) has only the end blocks (+ small posts at module junctions); the 302/303
# outer strips are full height. Old v2d slabs GapLowRib/GapUpRib = octagon minus 36 windows -> continuous 12 mm inner ribs
# that do not exist in CAD. New: octagon minus ONE square |x|,|y| <= 603.25+114.8 = 718.05 mm (the outer edge of the
# window grid, unchanged) -> only the outer band remains. Small posts are dropped (slightly less foam, accepted by user).
# Same solid names, heights and placements; only the two tessellated solids and their vertex defines are replaced.
import sys, re, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import make_foam_gdml_v2 as m
src, dst = sys.argv[1], sys.argv[2]
t = open(src).read()
outer = np.array(m.octagon(), float)
if m.signed_area(outer) < 0: outer = outer[::-1]
H = m.E_IN          # 723.9 mm = inner edge of the MEC-301 end blocks and of the 302/303 outer strips (CAD); a first version used 718.05
hole = np.array([(-H, -H), (-H, H), (H, H), (H, -H)], float)
if m.signed_area(hole) > 0: hole = hole[::-1]
for nm, h in (("GapLowRib", m.GAP_LO[1] - m.GAP_LO[0]), ("GapUpRib", m.GAP_UP[1] - m.GAP_UP[0])):
    V, F, area = m.flat_prism([outer, hole], h)
    d, s = m.tess_xml("PassiveTrackerFoam%sSolid" % nm, V, F, "Foam%s_v" % nm)
    pre = "Foam%s_v" % nm
    lines = t.split("\n")
    idx = [i for i, l in enumerate(lines) if ('<position name="%s' % pre) in l]
    assert idx and idx[-1] - idx[0] + 1 == len(idx), "defines not contiguous"
    lines[idx[0]:idx[-1] + 1] = d
    t = "\n".join(lines)
    t, n = re.subn(r'    <tessellated aunit="deg" lunit="mm" name="PassiveTrackerFoam%sSolid">.*?</tessellated>' % nm, "\n".join(s), t, flags=re.S)
    assert n == 1
    print(f"{nm}: h={h:.2f} mm, area {area/100:.1f} cm2, volume per gap {area*h/1000:.1f} cm3, {len(V)} vertices (old defines {len(idx)})")
open(dst, "w").write(t)
