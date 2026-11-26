#!/usr/bin/env python3
# fix_strips.py: rotate the 8 Si(Li) strip shapes of a v25.09-type gdml to the reco/hotfix/v300 layout.
#  1) BoxForIntersection: swap x <-> y   (narrow side along local x, long side 89.6 along local y)
#  2) strip intersection translation (0, y) -> (-y, 0)   => shape H (strip 0) at local x=-36.58 ... A at +36.58
#  3) every strip physvol (vid 2xxxxxxxx, k = vid%100 in 0..7) references "HGFEDCBA"[k]  (odd layers were reversed)
import re, sys
src, dst = sys.argv[1], sys.argv[2]
t = open(src).read()
n = collections = None
# 1) boxes
def box(m):
    return f'{m.group(1)}x="{m.group(3)}" y="{m.group(2)}"{m.group(4)}'
t, nb = re.subn(r'(<box lunit="mm" name="BoxForIntersection[^"]*" )x="([^"]+)" y="([^"]+)"( z="[^"]+"/>)', box, t)
# 2) intersection translations
def ipos(m):
    y = float(m.group(3)); assert float(m.group(2)) == 0.0
    return f'{m.group(1)}x="{-y!r}" y="0"{m.group(4)}'
t, ni = re.subn(r'(<position name="ActiveSiLiDetectorCylStripSecSolid[^"]*_pos" unit="mm" )x="([^"]+)" y="([^"]+)"( z="[^"]+"/>)', ipos, t)
# 3) volumerefs of strip physvols
lv = dict(re.findall(r'<volume name="(LVOLForIntersec([A-H])0x[0-9a-f]+)">', t)[i][::-1] for i in range(8))
order = "HGFEDCBA"
changed = 0; total = 0
def pv(m):
    global changed, total
    vid = int(m.group(2)); k = vid % 100; assert 0 <= k < 8, vid
    total += 1
    new = lv[order[k]]
    if new != m.group(4): changed += 1
    return f'{m.group(1)}{new}{m.group(5)}'
t = re.sub(r'(<physvol copynumber="(\d+)" name="(ActiveSiLiDetCylStripVolId[^"]*)">\s*<volumeref ref=")(LVOLForIntersec[A-H]0x[0-9a-f]+)("/>)', pv, t)
open(dst, "w").write(t)
print(f"boxes swapped {nb}, intersection translations {ni}, strip physvols {total}, volumeref changed {changed}")
