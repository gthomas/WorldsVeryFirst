"""Production Juris headstock furniture, shared by all three models (C-75, C-75A, C-75T).

Parts are stored in the headstock-local frame of the shared head_c75 outline:
  x  across the head, s along the face from the nut line (same s as head_c75.HEAD_OUTLINE / TUNER_POS),
  n  normal to the face, 0 = top of the neck wood (the face veneer sits on top of it; 0.6 mm on the C-75 / C-75T,
     1.6 mm on the C-75A, so face parts are lifted by the difference).
Parts:
  badge, badge_fill  the INSIGNIA: one etched nickel-silver plate, 27 mm octagon, gutters filled black,
                     set flush in an octagonal recess (RECESS_POLY, RECESS_DEPTH below the face)
  truss_cover, truss_screws, truss_engrave   rosewood cover, JURIS engraved and filled
  asset_tag          riveted plate on the back of the headstock
Source: the published C-75 model (insignia_c75.npz, extracted from it in its own frame)."""
import os, numpy as np, manifold3d as m3
from scipy.spatial import ConvexHull

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
TRUSS_COVER_S0 = float(_D['truss_cover_v'][:, 1].min())   # 4.5: the cover's nut end


def affine_from(f):
    """3x4 matrix of an affine map f(x, s, n) -> body xyz"""
    o = np.asarray(f(0.0, 0.0, 0.0), float)
    ex, ey, ez = (np.asarray(f(*e), float) - o for e in ((1, 0, 0), (0, 1, 0), (0, 0, 1)))
    return np.c_[ex, ey, ez, o]


def _man(v, f):
    return m3.Manifold(m3.Mesh(vert_properties=np.ascontiguousarray(v, np.float32), tri_verts=np.ascontiguousarray(f, np.uint32)))


def parts(M, face_t=C75_FACE_T, truss_s0=None):
    """manifolds of the headstock parts, placed by the 3x4 map M. face_t: this model's face veneer;
    truss_s0: move the truss cover so its nut end sits at this s (None = as on the C-75)"""
    out = {}
    for k in FACE_PARTS + BACK_PARTS:
        v = _D[k + '_v'].copy()
        if k in FACE_PARTS: v[:, 2] += face_t - C75_FACE_T
        if k in TRUSS_PARTS and truss_s0 is not None: v[:, 1] += truss_s0 - TRUSS_COVER_S0
        out[k] = _man(v @ M[:, :3].T + M[:, 3], _D[k + '_f'])
    return out


def recess(M, face_t=C75_FACE_T):
    """cutter for the badge pocket: from the recess floor up through the face"""
    n0 = face_t - RECESS_DEPTH; n1 = face_t + 5.0
    pts = [(x, s, n) for x, s in RECESS_POLY for n in (n0, n1)]
    P = np.array(pts) @ M[:, :3].T + M[:, 3]
    return m3.Manifold.hull_points(P.tolist())


if __name__ == '__main__':
    print('recess', RECESS_POLY.round(2).tolist(), 'depth', round(RECESS_DEPTH, 2), 'floor', round(RECESS_FLOOR, 2))
