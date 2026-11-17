"""Production Juris headstock furniture, shared by all three models (C-75, C-75A, C-75T).

Parts are stored in the headstock-local frame of the shared head_c75 outline:
  x  across the head, s along the face from the nut line (same s as head_c75.HEAD_OUTLINE / TUNER_POS),
  n  normal to the face, 0 = top of the neck wood (the face veneer sits on top of it; 0.6 mm on the C-75 / C-75T,
     1.6 mm on the C-75A, so face parts are lifted by the difference).
Parts:
  badge, badge_fill  the INSIGNIA: a 27 mm octagon of 1.2 mm (.048") 304 stainless, laser cut, gutters fibre-laser
                     marked black, screwed on proud of the face (like a Gretsch headstock plate) with two black
                     #0 pan-head screws in the horizontal gutters on the centreline (BADGE_SCREWS). No pocket.
                     badge_fill carries the black marking and the two screw heads.
  truss_cover, truss_screws, truss_engrave   rosewood cover, JURIS engraved and filled
  asset_tag          riveted plate on the back of the headstock
Source: the published C-75 model (insignia_c75.npz, extracted from it in its own frame)."""
import os, numpy as np, manifold3d as m3
from scipy.spatial import ConvexHull
from shapely.geometry import Polygon
from shapely.geometry.polygon import orient
from shapely.ops import unary_union

_D = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'insignia_c75.npz'))
FACE_PARTS = ('badge', 'badge_fill', 'truss_cover', 'truss_screws', 'truss_engrave')
BACK_PARTS = ('asset_tag',)
TRUSS_PARTS = ('truss_cover', 'truss_screws', 'truss_engrave')
C75_FACE_T = 0.6                                    # face veneer the parts were drawn on

_bv = _D['badge_v']
_h = ConvexHull(_bv[:, :2])
RECESS_POLY = _bv[_h.vertices, :2]                  # (x, s), CCW
RECESS_FLOOR = float(_bv[:, 2].min())               # n of the recess floor in the C-75 frame
RECESS_DEPTH = C75_FACE_T - RECESS_FLOOR            # below the finished face
ASSET_TAG_POLY = np.array([(-12, 95), (12, 95), (14, 97), (14, 105), (12, 107), (-12, 107), (-14, 105), (-14, 97)], float)   # 28 x 12 plate, 2 mm corner chamfers
ASSET_TAG_RIVETS = [(-11.0, 101.0), (11.0, 101.0)]        # rivet centres (heads 2.6 mm)
BADGE_POLY = RECESS_POLY                            # outline of the badge (x, s)
BADGE_S = 131.0                                     # badge centre along the face
BADGE_T = 1.2                                       # .048 in stainless
BADGE_SCREWS = [(0.0, BADGE_S + 5.07), (0.0, BADGE_S - 5.07)]   # centres of the two horizontal gutters, on the centreline
BADGE_HOLE, PILOT_HOLE, SCREW_HEAD = 1.8, 1.1, 2.8   # clearance hole in the badge, pilot in the face, head diameter
TRUSS_COVER_S0 = float(_D['truss_cover_v'][:, 1].min())   # 4.5: the cover's nut end


def affine_from(f):
    """3x4 matrix of an affine map f(x, s, n) -> body xyz"""
    o = np.asarray(f(0.0, 0.0, 0.0), float)
    ex, ey, ez = (np.asarray(f(*e), float) - o for e in ((1, 0, 0), (0, 1, 0), (0, 0, 1)))
    return np.c_[ex, ey, ez, o]


def _man(v, f):
    return m3.Manifold(m3.Mesh(vert_properties=np.ascontiguousarray(v, np.float32), tri_verts=np.ascontiguousarray(f, np.uint32)))


def gutters():
    """the insignia's gutters (x, s), as a shapely geometry"""
    import trimesh
    tm = trimesh.Trimesh(_D['badge_fill_v'], _D['badge_fill_f'])
    return unary_union([Polygon(tm.vertices[f][:, :2]) for f, n in zip(tm.faces, tm.face_normals) if n[2] > 0.9]).buffer(0)

def _prism_local(poly, n0, n1, M, holes=()):
    import build as B
    poly = orient(poly, 1.0)
    m = B.prism(np.array(poly.exterior.coords)[:-1], n0, n1)
    for ring in poly.interiors:
        r = np.array(ring.coords)[:-1]
        if not Polygon(r).exterior.is_ccw: r = r[::-1]
        m = m - B.prism(r, n0 - 1, n1 + 1)
    for x, s, d in holes:
        m = m - m3.Manifold.cylinder(n1 - n0 + 2, d / 2, d / 2, 24).translate([x, s, n0 - 1])
    me = m.to_mesh(); v = np.asarray(me.vert_properties)[:, :3] @ M[:, :3].T + M[:, 3]
    tri = np.asarray(me.tri_verts)
    if np.linalg.det(M[:, :3]) < 0: tri = tri[:, ::-1]
    return _man(v, tri)

def badge_parts(M, face_t=C75_FACE_T):
    """the screwed-on badge: plate on the face, black marking, two black pan-head screws"""
    n0 = face_t; n1 = face_t + BADGE_T
    plate = _prism_local(Polygon(BADGE_POLY), n0, n1, M, [(x, s, BADGE_HOLE) for x, s in BADGE_SCREWS])
    g = gutters()
    mark = [_prism_local(p, n1 - 0.01, n1 + 0.02, M) for p in getattr(g, 'geoms', [g]) if p.area > 0.01]
    heads = []
    for x, s in BADGE_SCREWS:
        prof = [(0, 0), (SCREW_HEAD / 2, 0), (SCREW_HEAD / 2, 0.35), (SCREW_HEAD / 2 * 0.8, 0.7), (0, 0.85)]
        h = m3.Manifold.revolve(m3.CrossSection([prof]), 24).translate([x, s, n1])
        me = h.to_mesh(); v = np.asarray(me.vert_properties)[:, :3] @ M[:, :3].T + M[:, 3]
        tri = np.asarray(me.tri_verts)
        if np.linalg.det(M[:, :3]) < 0: tri = tri[:, ::-1]
        heads.append(_man(v, tri))
    fill = mark[0]
    for m in mark[1:] + heads: fill = fill + m
    return plate, fill


def parts(M, face_t=C75_FACE_T, truss_s0=None):
    """manifolds of the headstock parts, placed by the 3x4 map M. face_t: this model's face veneer;
    truss_s0: move the truss cover so its nut end sits at this s (None = as on the C-75)"""
    out = {}
    for k in ('truss_cover', 'truss_screws', 'truss_engrave') + BACK_PARTS:
        v = _D[k + '_v'].copy()
        if k in FACE_PARTS: v[:, 2] += face_t - C75_FACE_T
        if k in TRUSS_PARTS and truss_s0 is not None: v[:, 1] += truss_s0 - TRUSS_COVER_S0
        out[k] = _man(v @ M[:, :3].T + M[:, 3], _D[k + '_f'])
    out['badge'], out['badge_fill'] = badge_parts(M, face_t)
    return out


def recess(M, face_t=C75_FACE_T):
    """no pocket any more: the badge screws on proud of the face (kept so callers need no change)"""
    return m3.Manifold()


if __name__ == '__main__':
    print('recess', RECESS_POLY.round(2).tolist(), 'depth', round(RECESS_DEPTH, 2), 'floor', round(RECESS_FLOOR, 2))
