"""gdml_mesh.py: minimal GDML reader for tessellated solids (used by add_rib_posts.py).
  G(file).placed(copynumber) -> (volume name, world triangles (N,3,3) in mm)
  slice_tri(T, axis, c, keep) -> line segments of the triangles cut by the plane coordinate[axis] = c
Rotation convention: world = Rx(-x) Ry(-y) Rz(-z) of the GDML rotation angles."""
import re
import numpy as np


def R(ax, a):
    c, s = np.cos(np.radians(a)), np.sin(np.radians(a))
    return {"x": np.array([[1,0,0],[0,c,-s],[0,s,c]]), "y": np.array([[c,0,s],[0,1,0],[-s,0,c]]), "z": np.array([[c,-s,0],[s,c,0],[0,0,1]])}[ax]

class G:
    def __init__(s, fn):
        s.t = open(fn).read()
        s.pos = {n: (float(x), float(y), float(z)) for n, x, y, z in re.findall(r'<position name="([^"]+)" unit="mm" x="([^"]+)" y="([^"]+)" z="([^"]+)"', s.t)}
    def tess(s, solid):
        m = re.search(r'<tessellated [^>]*name="%s">(.*?)</tessellated>' % re.escape(solid), s.t, re.S)
        tri = re.findall(r'<triangular vertex1="([^"]+)" vertex2="([^"]+)" vertex3="([^"]+)"', m.group(1))
        return np.array([[s.pos[a], s.pos[b], s.pos[c]] for a, b, c in tri])          # (N,3,3)
    def solid_of(s, vol):
        return re.search(r'<volume name="%s">\s*<materialref[^>]*/>\s*<solidref ref="([^"]+)"' % re.escape(vol), s.t).group(1)
    def placed(s, cn):
        m = re.search(r'<physvol copynumber="%s" name="[^"]*">\s*<volumeref ref="([^"]+)"/>\s*<position[^>]*x="([^"]+)" y="([^"]+)" z="([^"]+)"/>\s*(?:<rotation[^>]*x="([^"]+)" y="([^"]+)" z="([^"]+)"/>)?' % cn, s.t)
        vol = m.group(1); p = np.array(list(map(float, m.group(2, 3, 4))))
        rx, ry, rz = (float(v) if v else 0.0 for v in m.group(5, 6, 7))
        M = R("x", -rx) @ R("y", -ry) @ R("z", -rz)
        T = s.tess(s.solid_of(vol))
        return vol, T @ M.T + p

def slice_tri(T, axis, c, keep=(0, 1)):
    d = T[:, :, axis] - c
    segs = []
    for tri, dd in zip(T, d):
        pts = []
        for i, j in ((0, 1), (1, 2), (2, 0)):
            if (dd[i] > 0) != (dd[j] > 0):
                f = dd[i] / (dd[i] - dd[j]); pts.append(tri[i] + f * (tri[j] - tri[i]))
        if len(pts) == 2: segs.append([pts[0][list(keep)], pts[1][list(keep)]])
    return np.array(segs)
