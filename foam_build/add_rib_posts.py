#!/usr/bin/env python3
# add_rib_posts.py: put the MEC-301 rib-zone foam (posts at module junctions + end-block fillets) back.
# For each gap L = 0..8 and each rib layer (Low / Up):
#   cut = slice of the 6 v1 MEC-301 rows of THAT gap at the rib layer's mid height
#   rib = outer band (octagon minus |x|,|y|<=723.9)  UNION  (cut  MINUS  real obstacles at that height + 0.3 mm)
#   obstacles = sections of the GDML PassiveTracker{Top,Bot}Window and PassiveTrackerModuleFrame meshes at the same z
# All even gaps are identical to L0 and all odd gaps to L1 (checked: symmetric difference 0 mm2), so only 4 solids are
# built (Even/Odd x Low/Up from gaps L0/L1) and the 18 rib physvols reference them; physvol positions unchanged.
import sys, re, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_foam_gdml_v2 as m
from gdml_mesh import G, slice_tri
from shapely.geometry import Polygon, LineString, Point, box
from shapely.geometry.polygon import orient
from shapely.ops import polygonize, unary_union
src, dst = sys.argv[1], sys.argv[2]
v1 = G(os.path.join(HERE, "input", "gaps-assembly_foam_v1_mec301.gdml"))
cur = open(src).read()
octa = Polygon(m.octagon())
H = m.E_IN
band = octa.difference(box(-H, -H, H, H))
grid = octa.difference(unary_union([box(xc - m.WIN_HW, yc - m.WIN_HW, xc + m.WIN_HW, yc + m.WIN_HW) for xc in m.CELLS for yc in m.CELLS]))

def section_fast(T, z):
    d = T[:, :, 2] - z
    sgn = d > 0
    keep = ~(sgn.all(1) | (~sgn).all(1))
    T, d = T[keep], d[keep]
    segs = []
    for i, j in ((0, 1), (1, 2), (2, 0)):
        c = (d[:, i] > 0) != (d[:, j] > 0)
        f = np.where(c, d[:, i] / np.where(c, d[:, i] - d[:, j], 1), 0)
        segs.append((T[:, i] + f[:, None] * (T[:, j] - T[:, i]))[:, :2]); segs.append(c)
    P = [segs[0], segs[2], segs[4]]; C = np.stack([segs[1], segs[3], segs[5]], 1)
    out = []
    for k in range(len(T)):
        q = [P[e][k] for e in range(3) if C[k, e]]
        if len(q) == 2: out.append([q[0], q[1]])
    return np.array(out)

OBS = None
def obstacles_h(zc, h, margin=0.3):
    # union over 7 heights spanning the layer (windows are wider near their base)
    return unary_union([obstacles(z, margin) for z in zc + 0.5 * h * np.linspace(-0.998, 0.998, 7)])

def obstacles(z, margin=0.3):
    global OBS
    if OBS is None:
        g = G(src); OBS = []
        for cn, vol, pz in re.findall(r'<physvol copynumber="(\d+)" name="[^"]*">\s*<volumeref ref="([^"]+)"/>\s*<position[^>]*z="([^"]+)"', g.t):
            if ("Window" in vol or "Frame" in vol) and vol.startswith("Passive") and "Foam" not in vol:
                OBS.append((float(pz), g, cn))
    polys = []
    for pz, g, cn in OBS:
        if abs(pz - z) > 60: continue
        T = g.placed(cn)[1]
        if T[:, :, 2].min() > z or T[:, :, 2].max() < z: continue
        sec = section(T, z, fast=True)
        if not sec.is_empty: polys.append(sec)
    if not polys: return Polygon()
    u = unary_union(polys)
    # windows are thin hollow shells in the GDML: fill each outline (no foam inside a window cavity)
    filled = unary_union([Polygon(g.exterior) for g in (u.geoms if hasattr(u, "geoms") else [u])])
    return filled.buffer(margin)

def section(T, z, fast=False):
    s = np.round(section_fast(T, z) if fast else slice_tri(T, 2, z), 4)
    s = [x for x in s if not np.allclose(x[0], x[1])]
    A = np.array([x[0] for x in s]); B = np.array([x[1] for x in s])
    def inside(p):     # even-odd ray cast (+x direction) against all section segments
        y0, y1 = A[:, 1], B[:, 1]
        c = (y0 > p.y) != (y1 > p.y)
        xi = A[c, 0] + (p.y - y0[c]) * (B[c, 0] - A[c, 0]) / (y1[c] - y0[c])
        return np.count_nonzero(xi > p.x) % 2 == 1
    faces = [f for f in polygonize([LineString(x) for x in s]) if f.area > 0.5]
    return unary_union([f for f in faces if inside(f.representative_point())])

defines, solids, volumes = [], [], []
log = []
Z_LO0, Z_MID0, Z_MID1, Z_UP1 = 0.863, 11.7, 83.6, 96.013   # relative to the MEC-301 bottom; mid = CAD 301/304 solid middle 11.7-83.6 (thermal_step/foam_cad_summary.txt)
def set_z(txt, cn, z):
    q = re.search(r'(<physvol copynumber="%d" name="[^"]*">\s*<volumeref ref="[^"]+"/>\s*<position[^>]*z=")([^"]+)(")' % cn, txt)
    return txt[:q.start(2)] + "%.6f" % z + txt[q.end(2):]
ZB = {}
for LL in range(9):
    ZB[LL] = min(float(z) for z in re.findall(r'<physvol copynumber="50%d10[1-6]000" name="[^"]*">\s*<volumeref ref="[^"]+"/>\s*<position[^>]*z="([^"]+)"' % LL, v1.t))
for L in range(2):
    par = ("Even", "Odd")[L]
    rows = np.concatenate([v1.placed(500000000 + L * 1000000 + 101000 + 1000 * (n - 1))[1] for n in range(1, 7)])
    for tag, cn_off, za, zb_ in (("LowRib", 101000, Z_LO0, Z_MID0), ("UpRib", 103000, Z_MID1, Z_UP1)):
        h = zb_ - za
        cns = [500000000 + LL * 1000000 + cn_off for LL in range(L, 9, 2)]
        cn = cns[0]
        pm = re.search(r'<physvol copynumber="%d" name="[^"]*">\s*<volumeref ref="([^"]+)"/>\s*<position[^>]*z="([^"]+)"' % cn, cur)
        zc = ZB[L] + 0.5 * (za + zb_)
        for c2 in cns:
            cur = set_z(cur, c2, ZB[(c2 - 500000000) // 1000000] + 0.5 * (za + zb_))
        cut = section(rows, zc, fast=True)
        obs = obstacles_h(zc, h)
        posts = cut.intersection(box(-H, -H, H, H)).difference(obs).intersection(octa)
        posts = posts.buffer(-0.05, join_style=2).buffer(0.05, join_style=2).simplify(0.02)   # drop slivers < 0.1 mm
        rib = unary_union([band, posts]).buffer(0)
        polys = [rib] if rib.geom_type == "Polygon" else list(rib.geoms)
        polys = [orient(p, 1.0) for p in polys if p.area > 5.0]
        V, F, off = [], [], 0
        for p in polys:
            loops = [np.array(p.exterior.coords[:-1])] + [np.array(r.coords[:-1]) for r in p.interiors]
            try:
                v, f, _ = m.flat_prism(loops, h)
            except AssertionError:
                print("  dropped non-closable piece area %.1f mm2 at" % p.area, p.representative_point()); continue
            V.append(v); F.append(f + off); off += len(v)
        V = np.concatenate(V); F = np.concatenate(F)
        nm = "PassiveTrackerFoamGap%s%s" % (tag, par)
        d, s = m.tess_xml(nm + "Solid", V, F, "FoamGap%s%s_v" % (tag, par))
        defines += d; solids += s
        volumes += ['    <volume name="%s">' % nm, '      <materialref ref="%s"/>' % m.MAT, '      <solidref ref="%sSolid"/>' % nm, "    </volume>"]
        for c2 in cns:
            q = re.search(r'<physvol copynumber="%d" name="[^"]*">\s*<volumeref ref="([^"]+)"/>' % c2, cur)
            cur = cur[:q.start(1)] + nm + cur[q.end(1):]
        log.append("gaps %s %-6s z=%.2f h=%.2f: band %.0f cm2 + posts %.0f cm2 (%d pieces), cut before clip %.0f cm2, %d facets"
                   % (",".join("L%d" % ((c2 - 500000000) // 1000000) for c2 in cns), tag, zc, h, band.area / 100, posts.area / 100, len(polys) - 1, cut.intersection(box(-H, -H, H, H)).area / 100, len(F)))
# ---- middle slab = CAD 301/304 solid middle (zb+11.7 .. zb+83.6), Ø101.6 holes as before; shared by the 9 gaps
hmid = Z_MID1 - Z_MID0
circ = []
nseg = 72; rc = m.HOLE_R / np.cos(np.pi / nseg)
for xc in m.CELLS:
    for yc in m.CELLS:
        for dx in (-m.DET_OFF, m.DET_OFF):
            for dy in (-m.DET_OFF, m.DET_OFF):
                a = -np.linspace(0, 2 * np.pi, nseg, endpoint=False) + np.pi / nseg
                circ.append(np.c_[xc + dx + rc * np.cos(a), yc + dy + rc * np.sin(a)])
outer = np.array(m.octagon(), float)
if m.signed_area(outer) < 0: outer = outer[::-1]
V, F, _ = m.flat_prism([outer] + circ, hmid)
d, s2 = m.tess_xml("PassiveTrackerFoamGapMid301Solid", V, F, "FoamGapMid301_v")
defines += d; solids += s2
volumes += ['    <volume name="PassiveTrackerFoamGapMid301">', '      <materialref ref="%s"/>' % m.MAT, '      <solidref ref="PassiveTrackerFoamGapMid301Solid"/>', "    </volume>"]
for LL in range(9):
    c2 = 500000000 + LL * 1000000 + 102000
    q = re.search(r'<physvol copynumber="%d" name="[^"]*">\s*<volumeref ref="([^"]+)"/>' % c2, cur)
    cur = cur[:q.start(1)] + "PassiveTrackerFoamGapMid301" + cur[q.end(1):]
    cur = set_z(cur, c2, ZB[LL] + 0.5 * (Z_MID0 + Z_MID1))
log.append("middle slab: h=%.2f mm (was %.2f), z = zb+%.2f..zb+%.2f (L0: %.2f..%.2f)" % (hmid, 100.0 - m.GAP_LO[1] - m.GAP_UP[1], Z_MID0, Z_MID1, ZB[0] + Z_MID0, ZB[0] + Z_MID1))
# ---- 401 rib layer above L0 (cn 500001000): CAD MEC-401 = solid with Ø126.5 pockets 0-11.7 mm from the bottom
#      (thermal_step/foam_cad_summary.txt) -> octagon minus 144 circles Ø126.5 at the detector centres, minus real
#      windows/frames at that height (+0.3 mm). Own solid; it used the shared GapLowRib solid before.
cn = 500001000
pm = re.search(r'<physvol copynumber="%d" name="[^"]*">\s*<volumeref ref="([^"]+)"/>\s*<position[^>]*z="([^"]+)"' % cn, cur)
zc = float(pm.group(2)); h = m.GAP_LO[1] - m.GAP_LO[0]
circ = unary_union([Point(xc + dx, yc + dy).buffer(126.5 / 2, 128) for xc in m.CELLS for yc in m.CELLS for dx in (-m.DET_OFF, m.DET_OFF) for dy in (-m.DET_OFF, m.DET_OFF)])
cad401 = octa.difference(circ)
obs = obstacles_h(zc, h)
rib = cad401.difference(obs).buffer(-0.05, join_style=2).buffer(0.05, join_style=2).simplify(0.02).buffer(0)
polys = [orient(p, 1.0) for p in ([rib] if rib.geom_type == "Polygon" else list(rib.geoms)) if p.area > 5.0]
V, F, off = [], [], 0
for p in polys:
    loops = [np.array(p.exterior.coords[:-1])] + [np.array(r.coords[:-1]) for r in p.interiors]
    try:
        v, f, _ = m.flat_prism(loops, h)
    except AssertionError:
        print("  dropped non-closable 401 piece area %.1f mm2" % p.area); continue
    V.append(v); F.append(f + off); off += len(v)
V = np.concatenate(V); F = np.concatenate(F)
nm = "PassiveTrackerFoamL0AboveRib401"
d, s2 = m.tess_xml(nm + "Solid", V, F, "Foam401Rib_v")
defines += d; solids += s2
volumes += ['    <volume name="%s">' % nm, '      <materialref ref="%s"/>' % m.MAT, '      <solidref ref="%sSolid"/>' % nm, "    </volume>"]
cur = cur[:pm.start(1)] + nm + cur[pm.end(1):]
log.append("401 rib above L0 z=%.2f h=%.2f: CAD (octagon - 144 x Ø126.5) %.0f cm2, after removing windows/frames %.0f cm2 (%d pieces), %d facets"
           % (zc, h, cad401.area / 100, rib.area / 100, len(polys), len(F)))
cur = cur.replace("</define>", "\n".join(defines) + "\n  </define>", 1)
cur = cur.replace("</solids>", "\n".join(solids) + "\n  </solids>", 1)
k = cur.index("</volume>", cur.index('<volume name="PassiveTrackerFoamGapUpRib">')) + len("</volume>")
cur = cur[:k] + "\n" + "\n".join(volumes) + cur[k:]
open(dst, "w").write(cur)
print("\n".join(log))
