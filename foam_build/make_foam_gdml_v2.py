"""v2: full simplified CAD stack (thermal-sim STEP), material NGX 250 XPS.
Changes w.r.t. v1 (below):
  * material FoamNGX250XPS = polystyrene C8H8, 0.0248 g/cm3 (Owens Corning FOAMULAR NGX 250, EPD min 1.55 lb/ft3)
  * MEC-401 rows rotated 90 deg about z (CAD: 401 rows perpendicular to the 305/306 rows below)
  * NEW 402/403 outer strips beside 401: |x| 723.9..800.1, z 986.78..1037.58     (cn 500007000 / 500008000)
  * NEW top box = MEC-404 + MEC-501/502 merged (user decision): z 1037.58..1082.68, solid (cn 500009000)
  * NEW 201 window layer z 70.45..82.80: foam ribs around 36 square detector windows 229.6 mm (GDML MEC-105 wall max 114.46)
  * NEW bottom box = MEC-201/202/203 (mech. layer 1, under L9): z 0.55..70.45, solid (cn 509109000);
    CAD top 86.78, but 70.45..86.78 holds the L9 bottom windows (z>=72.1) + frames in GDML -> stopped at the
    CAD Ø126.5 pocket floor (70.45).
  all new pieces: 1600.2 x 1600.2 square, corners cut at |x|+|y| <= 1561.4 mm (corner TOF paddles at >= 1565)

v1: Write gaps-assembly_foam.gdml = gaps-assembly.gdml + tracker Ethafoam (78 pieces).

Pieces (decided with the user):
  * 9 inter-layer gaps (under L0 .. under L8), each: 6 x MEC-301 rows (mesh from
    SimpleDet trk-styrofoam-sec1-lay2.gdml, unchanged) + 2 outer strips (MEC-302/303,
    76.2 mm wide, full 100 mm high, corners cut at |u|+|w| <= 1561.4 mm, 1 mm from
    the cube-corner TOF paddles 160-163).
  * 6 x MEC-401 cover rows above L0 (5 STEP slabs, 'Simplified version for sim').
  * nothing below L9.
Material Ethafoam (C2H4, 0.035 g/cm3, as GAPSMaterials.cc). All volumes passive (no SensDet).
Copy number (9 digits)  5 0 L h 0 n 0 0 0
  L = detector layer (L0 = top), h = 0 foam above layer L (only the 401 cover, L = 0),
  h = 1 foam below layer L (gap L / L+1),  n = 1..6 rows (negative -> positive along the
  stacking axis), 7 = outer strip on the negative side, 8 = on the positive side.
Placement: every row uses the world transform of the module at the row end (the CAD
'module #6'), with the CAD chain  MEC-301 local -> STEP mockup -> module local -> world.

Usage: python3 -I make_foam_gdml.py in.gdml trk-styrofoam-sec1-lay2.gdml slabs_mockup.npz out.gdml pieces.txt
"""
import collections
import re
import sys

import numpy as np
from scipy.spatial import Delaunay
from shapely.geometry import Point, Polygon
from shapely.prepared import prep

MODULE_ORIGIN = np.array([-196.596, -52.387, 0.0])     # module #6 origin in the STEP mockup
M301 = np.array([[0, 0, -1], [0, 1, 0], [1, 0, 0]], float)   # MEC-301 local -> mockup
T301 = np.array([-800.1, -150.8, 0.0])
E_IN, E_OUT, CUT = 723.9, 800.1, 1561.4                  # outer strip |w| range, corner cut
CELL = np.array([-603.25, -361.95, -120.65, 120.65, 361.95, 603.25])
RZ90 = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], float)
KEY401 = "GAPS-TKR-MEC-401 Center Section, Layer 11 Bot, Tracker_Simplied version for sim"
C_ID, H_ID = "C0x55bc03689590", "H0x55bc036898b0"


def R(ax, a):
    c, s = np.cos(np.radians(a)), np.sin(np.radians(a))
    return {"x": np.array([[1, 0, 0], [0, c, -s], [0, s, c]]), "y": np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]]),
            "z": np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])}[ax]


def gdml_matrix(x, y, z):
    """world = M @ local for a GDML physvol rotation (G4GDML uses (Rz Ry Rx)^-1)."""
    return R("x", -x) @ R("y", -y) @ R("z", -z)


def gdml_angles(M):
    B = M.T                                           # B = Rz(z) Ry(y) Rx(x)
    y = -np.degrees(np.arcsin(np.clip(B[2, 0], -1, 1)))
    if abs(B[2, 0]) < 1 - 1e-9:
        x = np.degrees(np.arctan2(B[2, 1], B[2, 2]))
        z = np.degrees(np.arctan2(B[1, 0], B[0, 0]))
    else:                                             # gimbal lock (y = -+90): take x = 0
        x = 0.0
        z = np.degrees(np.arctan2(-B[0, 1], B[1, 1]))
    assert np.allclose(gdml_matrix(x, y, z), M, atol=1e-7), "angle decomposition failed"
    return x, y, z


# ------------------------------------------------------------ inputs ----
def module_frames(txt):
    out = []
    pat = re.compile(r'<physvol copynumber="\d+" name="PassiveModuleFrameL(\d)R(\d)M(\d)[^"]*">\s*'
                     r'<volumeref ref="[^"]+"/>\s*<position[^>]*x="([^"]+)" y="([^"]+)" z="([^"]+)"/>\s*'
                     r'(?:<rotation[^>]*x="([^"]+)" y="([^"]+)" z="([^"]+)"/>)?')
    for m in pat.finditer(txt):
        a = [float(v) if v else 0.0 for v in m.group(7, 8, 9)]
        out.append(dict(L=int(m.group(1)), pos=np.array([float(v) for v in m.group(4, 5, 6)]),
                        W=gdml_matrix(*a)))
    assert len(out) == 360, len(out)
    return out


def load_mec301(path):
    t = open(path).read()
    pos = {m.group(1): (float(m.group(2)), float(m.group(3)), float(m.group(4))) for m in
           re.finditer(r'<position name="([^"]+)" unit="mm" x="([^"]+)" y="([^"]+)" z="([^"]+)"', t)}
    names = list(pos)
    idx = {n: i for i, n in enumerate(names)}
    tris = np.array([[idx[m.group(1)], idx[m.group(2)], idx[m.group(3)]] for m in
                     re.finditer(r'<triangular vertex1="([^"]+)" vertex2="([^"]+)" vertex3="([^"]+)"', t)])
    return np.array([pos[n] for n in names]), tris


def check_mesh(V, T, label):
    e = collections.Counter()
    for a, b, c in T:
        for u, v in ((a, b), (b, c), (c, a)):
            e[(u, v)] += 1
    bad = sum(1 for (u, v), n in e.items() if n != 1 or e.get((v, u), 0) != 1)
    vol = np.einsum("ij,ij->i", V[T[:, 0]], np.cross(V[T[:, 1]], V[T[:, 2]])).sum() / 6
    print("  %-22s %6d vertices %6d facets  open/non-manifold edges %d  volume %.1f cm3"
          % (label, len(V), len(T), bad, vol / 1000))
    assert bad == 0, label + " mesh is not closed"
    return vol


# ------------------------------------------------- MEC-401 from slabs ----
def loops_from_chords(ch):
    p = np.round(ch, 2)
    nxt = {}
    for s in p:
        nxt.setdefault((s[0], s[1]), []).append((s[2], s[3]))
    loops, used = [], set()
    for start in list(nxt):
        for first in nxt[start]:
            if (start, first) in used:
                continue
            loop, a, b = [start], start, first
            while True:
                used.add((a, b))
                if b == start:
                    break
                loop.append(b)
                cand = [c for c in nxt[b] if (b, c) not in used]
                a, b = b, cand[0]
            loops.append(np.array(loop))
    return loops


def signed_area(L):
    x, z = L[:, 0], L[:, 1]
    return 0.5 * np.sum(x * np.roll(z, -1) - np.roll(x, -1) * z)


def triangulate(loops, region):
    """conforming Delaunay of all loop vertices; boundary edges missing -> split, repeat."""
    loops = [l.copy() for l in loops]
    for it in range(30):
        P = np.vstack(loops)
        tri = Delaunay(P)
        E = set()
        for s in tri.simplices:
            for u, v in ((s[0], s[1]), (s[1], s[2]), (s[2], s[0])):
                E.add((min(u, v), max(u, v)))
        off, missing, new = 0, 0, []
        for l in loops:
            n = len(l)
            out = []
            for i in range(n):
                u, v = off + i, off + (i + 1) % n
                out.append(l[i])
                if (min(u, v), max(u, v)) not in E:
                    out.append(0.5 * (l[i] + l[(i + 1) % n]))
                    missing += 1
            new.append(np.array(out))
            off += n
        if missing == 0:
            break
        loops = new
    else:
        raise RuntimeError("triangulation did not conform")
    cen = P[tri.simplices].mean(1)
    keep = np.array([region.contains(Point(c)) for c in cen])
    return loops, P, tri.simplices[keep]


def slab_prism(loops, y0, y1):
    outer = [l for l in loops if signed_area(l) > 0]
    holes = [l for l in loops if signed_area(l) < 0]
    polys = []
    for o in outer:
        po = Polygon(o)
        polys.append(Polygon(o, [h for h in holes if po.contains(Point(h[0]))]))
    region = prep(polys[0] if len(polys) == 1 else __import__("shapely.ops").ops.unary_union(polys))
    area = sum(p.area for p in polys)
    loops, P, tris = triangulate(loops, region)
    n = len(P)
    V = np.vstack([np.c_[P[:, 0], np.full(n, y0), P[:, 1]], np.c_[P[:, 0], np.full(n, y1), P[:, 1]]])
    F = []
    # in (x, z) a positively oriented triangle has normal -y (x cross z = -y)
    for t in tris:
        if signed_area(P[t]) < 0:
            t = t[::-1]
        F.append([t[0], t[1], t[2]])                 # bottom cap, normal -y
        F.append([t[0] + n, t[2] + n, t[1] + n])     # top cap, normal +y
    off = 0
    for l in loops:                                  # outer loops CCW, holes CW in (x, z)
        m = len(l)
        for i in range(m):
            a, b = off + i, off + (i + 1) % m
            F.append([b, a, a + n])
            F.append([b, a + n, b + n])
        off += m
    F = np.array(F)
    vol = np.einsum("ij,ij->i", V[F[:, 0]], np.cross(V[F[:, 1]], V[F[:, 2]])).sum() / 6
    if vol < 0:                                      # fix global orientation if needed
        F = F[:, ::-1]
    return V, F, area


# ------------------------------------------------------------ output ----
def fmt(v):
    return "%.6f" % v


def tess_xml(name, V, F, vprefix):
    d = ['    <position name="%s%d" unit="mm" x="%s" y="%s" z="%s"/>' % (vprefix, i, fmt(p[0]), fmt(p[1]), fmt(p[2]))
         for i, p in enumerate(V)]
    s = ['    <tessellated aunit="deg" lunit="mm" name="%s">' % name]
    s += ['      <triangular vertex1="%s%d" vertex2="%s%d" vertex3="%s%d" type="ABSOLUTE"/>'
          % (vprefix, a, vprefix, b, vprefix, c) for a, b, c in F]
    s.append("    </tessellated>")
    return d, s


EXTRA_DEFINES = []


def xtru_xml(name, poly, zh):
    """v2f: written as a closed tessellated prism (TGeoXtru triggered loading a second ROOT libGeom in GAPSSim)"""
    P = np.array(poly, float)
    if signed_area(P) < 0:
        P = P[::-1]
    V, F, area = flat_prism([P], 2 * zh)
    d, x = tess_xml(name, V, F, name + "_v")
    EXTRA_DEFINES.extend(d)
    return x


def xtru_xml_old(name, poly, zh):
    a = signed_area(np.array(poly))
    if a > 0:
        poly = poly[::-1]                            # G4ExtrudedSolid: clockwise vertices
    s = ['    <xtru lunit="mm" name="%s">' % name]
    s += ['      <twoDimVertex x="%s" y="%s"/>' % (fmt(x), fmt(y)) for x, y in poly]
    s += ['      <section zOrder="0" zPosition="%s" xOffset="0" yOffset="0" scalingFactor="1"/>' % fmt(-zh),
          '      <section zOrder="1" zPosition="%s" xOffset="0" yOffset="0" scalingFactor="1"/>' % fmt(zh),
          "    </xtru>"]
    return s


def strip_poly(axis, s):
    """outer strip on side s (+-1) of stacking axis 'x' or 'y', world (x, y), corners cut."""
    wi, wo = s * E_IN, s * E_OUT
    c = CUT - E_OUT                                  # |u| where the cut meets the outer edge
    uw = [(-E_OUT, wi), (E_OUT, wi), (E_OUT, s * (CUT - E_OUT)), (c, wo), (-c, wo), (-E_OUT, s * (CUT - E_OUT))]
    # if the cut starts inside the strip the inner corner stays: |u|=E_OUT at |w|=CUT-E_OUT >= E_IN
    assert CUT - E_OUT >= E_IN
    return [(w, u) if axis == "x" else (u, w) for u, w in uw]


def physvol(name, vol, cn, pos, M):
    x, y, z = gdml_angles(M)
    return ['      <physvol copynumber="%d" name="%s">' % (cn, name),
            '        <volumeref ref="%s"/>' % vol,
            '        <position name="%s_pos" unit="mm" x="%s" y="%s" z="%s"/>' % (name, fmt(pos[0]), fmt(pos[1]), fmt(pos[2])),
            '        <rotation name="%s_rot" unit="deg" x="%s" y="%s" z="%s"/>' % (name, fmt(x), fmt(y), fmt(z)),
            "      </physvol>"]


MAT = "FoamNGX250XPS"


def octagon():
    h, c = E_OUT, CUT - E_OUT
    return [(-c, -h), (c, -h), (h, -c), (h, c), (c, h), (-c, h), (-h, c), (-h, -c)]


def add_vol(name, solid_xml, volumes, solids):
    solids += solid_xml
    volumes += ['    <volume name="%s">' % name, '      <materialref ref="%s"/>' % MAT,
                '      <solidref ref="%sSolid"/>' % name, "    </volume>"]


# ---- v2b simplified gap (301/304 + 302/303), relative to the frame centres (mm), from GDML window meshes:
#  top window MEC-103: zc+1.59..zc+11.11 (full-width plate to ~zc+2.3, wall half-width <= 114.45 above)
#  bottom window MEC-105: zc'-13.14..zc'-1.59 (plate from zc'-2.3, wall <= 114.46 below)
GAP_LO = (2.45, 11.40)      # lower rib zone above the lower frame centre
GAP_UP = (2.40, 13.40)      # upper rib zone below the upper frame centre
WIN_HW = 114.8              # square window half width in the rib zones
HOLE_R = 50.8               # Ø101.6 through holes (one per Si(Li) detector)
CELLS = [-603.25, -361.95, -120.65, 120.65, 361.95, 603.25]
DET_OFF = 57.65             # detector centres = module centre +- 57.65 (GDML Si(Li) placements)


def flat_prism(loops, h):
    """polygon loops (outer CCW + holes CW) in (x, y) -> closed mesh extruded along z, centred (-h/2..h/2)"""
    loops2 = [np.array(l, float) for l in loops]
    V, F, area = slab_prism(loops2, -h / 2, h / 2)    # slab_prism extrudes along its 2nd axis: (P0, y, P1)
    V = V[:, [0, 2, 1]]                              # -> (x, y, z); axis swap flips handedness
    F = F[:, ::-1]
    vol = np.einsum("ij,ij->i", V[F[:, 0]], np.cross(V[F[:, 1]], V[F[:, 2]])).sum() / 6
    if vol < 0:
        F = F[:, ::-1]
    check_mesh(V, F, "gap prism h=%.2f" % h)
    return V, F, area


def gap_solids(defines, solids, volumes):
    oct_ = [np.array(p) for p in octagon()]
    outer = np.array(oct_)
    if signed_area(outer) < 0:
        outer = outer[::-1]
    sq = []
    for xc in CELLS:
        for yc in CELLS:
            q = np.array([(xc - WIN_HW, yc - WIN_HW), (xc - WIN_HW, yc + WIN_HW), (xc + WIN_HW, yc + WIN_HW), (xc + WIN_HW, yc - WIN_HW)])
            sq.append(q if signed_area(q) < 0 else q[::-1])
    nseg = 72
    rc = HOLE_R / np.cos(np.pi / nseg)               # circumscribed polygon: hole never smaller than Ø101.6
    circ = []
    for xc in CELLS:
        for yc in CELLS:
            for dx in (-DET_OFF, DET_OFF):
                for dy in (-DET_OFF, DET_OFF):
                    a = -np.linspace(0, 2 * np.pi, nseg, endpoint=False) + np.pi / nseg
                    circ.append(np.c_[xc + dx + rc * np.cos(a), yc + dy + rc * np.sin(a)])
    hlo = GAP_LO[1] - GAP_LO[0]
    hup = GAP_UP[1] - GAP_UP[0]
    hmid = 100.0 - GAP_LO[1] - GAP_UP[1]
    for nm, loops, h in (("GapLowRib", [outer] + sq, hlo), ("GapUpRib", [outer] + sq, hup), ("GapMid", [outer] + circ, hmid)):
        V, F, area = flat_prism(loops, h)
        print("      %-10s h=%.2f mm  area %.1f cm2  (%d vertices)" % (nm, h, area / 100, len(V)))
        d, x = tess_xml("PassiveTrackerFoam%sSolid" % nm, V, F, "Foam%s_v" % nm)
        defines += d; solids += x
        volumes += ['    <volume name="PassiveTrackerFoam%s">' % nm, '      <materialref ref="%s"/>' % MAT,
                    '      <solidref ref="PassiveTrackerFoam%sSolid"/>' % nm, "    </volume>"]


def new_pieces(L, stack, solids, volumes, pvs, report):
    if L == 0:
        z0, z1 = 1037.58, 1082.68
        nm = "PassiveTrackerFoamTopBox"
        add_vol(nm, xtru_xml(nm + "Solid", octagon(), (z1 - z0) / 2), volumes, solids)
        pvs += physvol("PassiveTrackerFoamL0AboveTopBox", nm, 500009000, np.array([0, 0, (z0 + z1) / 2]), np.eye(3))
        report.append("500009000  L0 above  MEC-404+501/502 merged box, z %.2f..%.2f" % (z0, z1))
    if L == 9:
        z0, z1 = 0.55, 70.45
        nm = "PassiveTrackerFoamBottomBox"
        add_vol(nm, xtru_xml(nm + "Solid", octagon(), (z1 - z0) / 2), volumes, solids)
        pvs += physvol("PassiveTrackerFoamL9BelowBox", nm, 509109000, np.array([0, 0, (z0 + z1) / 2]), np.eye(3))
        report.append("509109000  L9 below  MEC-201/202/203 box, z %.2f..%.2f" % (z0, z1))
        # window layer of 201: foam ribs around 36 square windows (half 114.0 mm; the L9 MEC-105 bottom
        # window mesh spans |x|,|y| <= 113.5 for z 72.05..82.9, its 0.7 mm plate to 120.5 from z 82.91)
        z0, z1, hw = 70.45, 82.80, 114.8   # window wall widens to 114.46 at z 82.9 (mesh slice)
        zc, zh = (z0 + z1) / 2, (z1 - z0) / 2
        cen = [-603.25, -361.95, -120.65, 120.65, 361.95, 603.25]
        rail = 120.65 - hw                                  # half width of inner ribs
        edge = 603.25 + hw                                  # outer edge of the window grid
        def box(nm, dx, dy):
            add_vol(nm, ['    <box lunit="mm" name="%sSolid" x="%s" y="%s" z="%s"/>' % (nm, fmt(dx), fmt(dy), fmt(2 * zh))], volumes, solids)
        box("PassiveTrackerFoamBotWinRail", 2 * E_OUT, 2 * rail)
        box("PassiveTrackerFoamBotWinSegIn", 2 * rail, 2 * hw)
        box("PassiveTrackerFoamBotWinSegOut", E_OUT - edge, 2 * hw)
        k = 0
        for yc in (-482.6, -241.3, 0.0, 241.3, 482.6):      # 5 inner rails along x
            k += 1
            pvs += physvol("PassiveTrackerFoamL9BelowWinRail%d" % k, "PassiveTrackerFoamBotWinRail", 509110000 + k,
                           np.array([0, yc, zc]), np.eye(3))
        for s_, tag in ((-1, "Neg"), (1, "Pos")):           # 2 outer rails |y| edge..800.1, corners cut
            nm = "PassiveTrackerFoamBotWinEdgeY%s" % tag
            pts = [(u, s_ * edge) if abs(w - s_ * E_IN) < 1e-6 else (u, w) for u, w in strip_poly("y", s_)]
            add_vol(nm, xtru_xml(nm + "Solid", pts, zh), volumes, solids)
            pvs += physvol("PassiveTrackerFoamL9BelowWinEdge%s" % tag, nm, 509120000 + (s_ > 0), np.array([0, 0, zc]), np.eye(3))
        k = 0
        for yc in cen:                                      # segments between windows, inside each window row
            for xc in (-482.6, -241.3, 0.0, 241.3, 482.6):
                k += 1
                pvs += physvol("PassiveTrackerFoamL9BelowWinSeg%d" % k, "PassiveTrackerFoamBotWinSegIn", 509130000 + k,
                               np.array([xc, yc, zc]), np.eye(3))
            for xs in (-1, 1):
                k += 1
                pvs += physvol("PassiveTrackerFoamL9BelowWinSeg%d" % k, "PassiveTrackerFoamBotWinSegOut", 509130000 + k,
                               np.array([xs * (edge + E_OUT) / 2, yc, zc]), np.eye(3))
        report.append("5091xxxxx  L9 below  201 window layer z %.2f..%.2f: 36 square windows %.1f mm, ribs %.1f mm"
                      " (5 rails + 2 edge strips + %d segments)" % (z0, z1, 2 * hw, 2 * rail, k))


def main():
    gin, sty, npz, gout = sys.argv[1:5]
    txt = open(gin).read()
    frames = module_frames(txt)
    print("module frames:", len(frames))

    print("meshes:")
    V301, T301m = load_mec301(sty)
    v301 = check_mesh(V301, T301m, "MEC-301")
    if v301 < 0:
        T301m = T301m[:, ::-1]
    slabs = np.load(npz)
    import json
    summ = json.load(open(npz.replace(".npz", "_summary.json")))[KEY401]["slabs"]
    parts401 = []
    for i, s in enumerate(summ):
        loops = loops_from_chords(slabs[KEY401 + "||%d" % i])
        V, F, area = slab_prism(loops, s["y_lo"], s["y_hi"])
        vol = check_mesh(V, F, "MEC-401 slab %d" % i)
        print("      slab %d y %.2f..%.2f  area %.1f cm2 (STEP %.1f)  V/(A h) = %.5f"
              % (i, s["y_lo"], s["y_hi"], area / 100, s["area_mm2"] / 100, vol / (area * (s["y_hi"] - s["y_lo"]))))
        parts401.append((V, F))

    defines, solids, volumes, pvs = [], [], [], []
    gap_solids(defines, solids, volumes)              # v2b: simplified 301/304 gaps (replaces MEC-301 mesh + strips)
    parts401 = []                                     # v2e: 401 CAD mesh no longer written
    for i, (V, F) in enumerate(parts401):
        d, s = tess_xml("PassiveTrackerFoamMEC401Slab%dSolid" % i, V, F, "FoamMEC401s%d_v" % i)
        defines += d; solids += s
    prev = "PassiveTrackerFoamMEC401Slab0Solid"
    for i in range(1, len(parts401)):
        name = "PassiveTrackerFoamMEC401Solid" if i == len(parts401) - 1 else "PassiveTrackerFoamMEC401U%d" % i
        solids += ['    <union name="%s">' % name, '      <first ref="%s"/>' % prev,
                   '      <second ref="PassiveTrackerFoamMEC401Slab%dSolid"/>' % i, "    </union>"]
        prev = name
    vols = {}
    zframe = {L: float(np.mean([f["pos"][2] for f in frames if f["L"] == L])) for L in range(10)}
    for k, sol in vols.items():
        volumes += ['    <volume name="PassiveTrackerFoam%s">' % k, '      <materialref ref="%s"/>' % MAT + '',
                    '      <solidref ref="%s"/>' % sol, "    </volume>"]

    report = []
    for L in range(10):
        mods = [f for f in frames if f["L"] == L]
        W = mods[0]["W"]
        assert all(np.allclose(f["W"], W) for f in mods), "modules of a layer not parallel"
        ux = W[:, 0]                                  # module local x in world = along the row
        stack = "x" if abs(ux[1]) > 0.99 else "y"     # rows run along u; they stack along the other axis
        groups = collections.defaultdict(list)
        for f in mods:
            groups[round(float(f["pos"]["xy".index(stack)]), 0)].append(f)
        assert len(groups) == 6 and all(len(g) == 6 for g in groups.values())
        rows = sorted(groups.items())
        layer_pieces = []
        for n, (wc, g) in enumerate(rows, 1):
            m6 = max(g, key=lambda f: float(f["pos"] @ ux))      # row end module (CAD module #6)
            Rw, pw = m6["W"], m6["pos"]
            # check: the 6 module centres land on the MEC-301 cell centres
            loc = np.array([M301.T @ (Rw.T @ (f["pos"] - pw) + MODULE_ORIGIN - T301) for f in g])
            err = max(np.abs(loc[:, 0]).max(), np.abs(np.sort(loc[:, 2]) - CELL).max())
            A301, t301 = Rw @ M301, Rw @ (T301 - MODULE_ORIGIN) + pw
            A401, t401 = Rw, Rw @ (-MODULE_ORIGIN) + pw
            A401, t401 = RZ90 @ A401, RZ90 @ t401          # v2: 401 rows perpendicular to L0
            layer_pieces.append((n, A301, t301, A401, t401, err, wc))
        z301 = [(A @ np.array([0, y, 0]) + t)[2] for y in (0.0, 100.0) for (_, A, t, *_r) in layer_pieces[:1]]
        if L == 9:
            new_pieces(L, stack, solids, volumes, pvs, report)
        if L <= 8:                                    # foam below layer L = gap between frames L+1 (lower) and L (upper)
            zlo, zup = zframe[L + 1], zframe[L]
            for n, nm, z0, z1 in ((1, "GapLowRib", zlo + GAP_LO[0], zlo + GAP_LO[1]),
                                  (2, "GapMid", zlo + GAP_LO[1], zup - GAP_UP[1]),
                                  (3, "GapUpRib", zup - GAP_UP[1], zup - GAP_UP[0])):
                cn = 500000000 + L * 1000000 + 100000 + n * 1000
                pvs += physvol("PassiveTrackerFoamL%dBelow%s" % (L, nm), "PassiveTrackerFoam" + nm, cn,
                               np.array([0, 0, 0.5 * (z0 + z1)]), np.eye(3))
                report.append("%d  L%d below  %-10s z %.2f..%.2f" % (cn, L, nm, z0, z1))
        if L == 0:                                    # 401 cover above L0
            # v2e: simplified 401 + 402/403 (CAD mesh rotated 90 deg overlapped the L0 top windows by 0.03 mm)
            zc0 = zframe[0]
            z0, z1 = zc0 + GAP_LO[0], zc0 + GAP_LO[1]
            pvs += physvol("PassiveTrackerFoamL0AboveRib", "PassiveTrackerFoamGapLowRib", 500001000, np.array([0, 0, 0.5 * (z0 + z1)]), np.eye(3))
            report.append("500001000  L0 above  401 rib layer (36 square windows), z %.2f..%.2f" % (z0, z1))
            z0, z1 = zc0 + GAP_LO[1], 1037.58
            nm = "PassiveTrackerFoamCover401"
            add_vol(nm, xtru_xml(nm + "Solid", octagon(), (z1 - z0) / 2), volumes, solids)
            pvs += physvol("PassiveTrackerFoamL0AboveCover", nm, 500002000, np.array([0, 0, 0.5 * (z0 + z1)]), np.eye(3))
            report.append("500002000  L0 above  401+402/403 solid cover, z %.2f..%.2f" % (z0, z1))
            new_pieces(L, stack, solids, volumes, pvs, report)

    mat = ['    <material name="%s" state="solid">' % MAT, '      <D unit="g/cm3" value="0.0248"/>',
           '      <fraction n="%.6f" ref="%s"/>' % (12.0107 / (12.0107 + 1.00794), C_ID),
           '      <fraction n="%.6f" ref="%s"/>' % (1.00794 / (12.0107 + 1.00794), H_ID), "    </material>"]
    out = txt
    defines += EXTRA_DEFINES
    for tag, block in (("  </define>", defines), ("  </materials>", mat), ("  </solids>", solids)):
        assert out.count(tag) == 1, tag
        out = out.replace(tag, "\n".join(block) + "\n" + tag)
    wv = '    <volume name="WorldVolId00x55bc036a2390">'
    assert out.count(wv) == 1
    out = out.replace(wv, "\n".join(volumes) + "\n" + wv)
    end = "    </volume>\n  </structure>"
    assert out.count(end) == 1
    out = out.replace(end, "\n".join(pvs) + "\n" + end)
    open(gout, "w").write(out)
    print("\n%d foam physvols written to %s" % (len(pvs) // 5, gout))
    open(sys.argv[5], "w").write("\n".join(report) + "\n")
    print("\n".join(report))


if __name__ == "__main__":
    main()
