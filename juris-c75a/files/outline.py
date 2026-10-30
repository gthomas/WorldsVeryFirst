"""Juris C-75A body outline + facet layout. Units: mm. x+ = treble side, y+ = toward neck, y=0 at tail."""
import numpy as np

IN = 25.4

# Edge lines traced from Juris_C_75_Body_2D.pdf (inches, calibrated on the 18.000 / 13.500 dims)
_E = {
    'b0': ((-3.258, 18.011), (-1.420, 17.101)),   # bass horn, inner (neck side)
    'b1': ((-3.475, 18.011), (-5.227, 15.224)),   # bass horn, outer
    'b2': ((-5.227, 15.224), (-3.949, 11.792)),   # bass upper bout -> waist
    'b3': ((-3.916, 11.430), (-6.670, 5.487)),    # waist -> lower bout
    'b4': ((-6.762, 4.848), (-4.674, 1.199)),     # lower bout -> tail corner
    'b5': ((-4.193, 0.817), (-0.636, -0.079)),    # tail V
    't5': ((1.077, -0.013), (4.173, 0.784)),
    't4': ((4.496, 1.034), (6.762, 4.868)),
    't3': ((6.683, 5.435), (4.094, 11.034)),
    't2': ((4.121, 11.561), (5.122, 14.084)),
    't1': ((5.109, 14.282), (4.648, 15.204)),     # cutaway horn outer
    't0': ((4.404, 15.389), (2.447, 14.736)),     # cutaway floor
}

# lower body is meant to be symmetric; average the traced halves
def _sym(bk, tk):
    (a, b), (c, d) = _E[bk], _E[tk]
    # mirror treble onto bass and average (treble points are listed in opposite direction)
    pa = np.mean([a, (-d[0], d[1])], axis=0); pb = np.mean([b, (-c[0], c[1])], axis=0)
    _E[bk] = (tuple(pa), tuple(pb)); _E[tk] = ((-pb[0], pb[1]), (-pa[0], pa[1]))
for k in ('3', '4', '5'):
    _sym('b' + k, 't' + k)


def _line(k):
    p, q = np.array(_E[k][0]) * IN, np.array(_E[k][1]) * IN
    return p, q - p

def _isect(l1, l2):
    (p, d), (q, e) = l1, l2
    A = np.array([d, -e]).T
    t = np.linalg.solve(A, q - p)
    return p + t[0] * d

def _at_x(l, x):
    p, d = l; t = (x - p[0]) / d[0]; return p + t * d


NECK_HALF = 28.0          # half-width of neck/heel tongue at body (2.2")
TAIL_HALF = 16.0          # small flat at tail for tailpiece hinge / endpin


def sharp_outline(neck_end_y_in=17.6):
    L = {k: _line(k) for k in _E}
    xN = NECK_HALF
    t_tail = _at_x(L['t5'], TAIL_HALF); b_tail = _at_x(L['b5'], -TAIL_HALF)
    ytail = (t_tail[1] + b_tail[1]) / 2
    pts = [
        ('tailR', (TAIL_HALF, ytail)),
        ('tTail', _isect(L['t5'], L['t4'])),
        ('tLower', _isect(L['t4'], L['t3'])),
        ('tWaist', _isect(L['t3'], L['t2'])),
        ('tUpper', _isect(L['t2'], L['t1'])),
        ('tHorn', _isect(L['t1'], L['t0'])),
        ('tCut', _at_x(L['t0'], xN)),
        ('neckT', (xN, neck_end_y_in * IN)),
        ('neckB', (-xN, neck_end_y_in * IN)),
        ('bCut', _at_x(L['b0'], -xN)),
        ('bHorn', _isect(L['b0'], L['b1'])),
        ('bUpper', _isect(L['b1'], L['b2'])),
        ('bWaist', _isect(L['b2'], L['b3'])),
        ('bLower', _isect(L['b3'], L['b4'])),
        ('bTail', _isect(L['b4'], L['b5'])),
        ('tailL', (-TAIL_HALF, ytail)),
    ]
    names = [n for n, _ in pts]
    P = np.array([p for _, p in pts], float)
    P[:, 1] -= ytail          # tail flat at y = 0
    return names, P


# corner radii (mm). convex corners must stay >= chamfer+offset used in build
RADII = dict(tailR=40, tailL=40, tTail=14, bTail=14, tLower=14, bLower=14, tWaist=16, bWaist=16,
             tUpper=14, bUpper=14, tHorn=12, bHorn=12, tCut=19, bCut=8, neckT=0, neckB=0)   # square corners where the heel lands


def fillet(names, P, radii=RADII, seg_deg=6.0):
    out = []
    n = len(P)
    for i in range(n):
        V, A, B = P[i], P[i - 1], P[(i + 1) % n]
        r = radii.get(names[i], 0)
        u1 = (A - V) / np.linalg.norm(A - V); u2 = (B - V) / np.linalg.norm(B - V)
        cosang = np.clip(np.dot(u1, u2), -1, 1); ang = np.arccos(cosang)
        if r <= 0 or ang > np.pi - 1e-3:
            out.append(V); continue
        d = r / np.tan(ang / 2)
        d = min(d, 0.45 * np.linalg.norm(A - V), 0.45 * np.linalg.norm(B - V))
        r_eff = d * np.tan(ang / 2)
        T1, T2 = V + u1 * d, V + u2 * d
        bis = (u1 + u2) / np.linalg.norm(u1 + u2)
        C = V + bis * (r_eff / np.sin(ang / 2))
        a1 = np.arctan2(*(T1 - C)[::-1]); a2 = np.arctan2(*(T2 - C)[::-1])
        da = (a2 - a1 + np.pi) % (2 * np.pi) - np.pi
        k = max(2, int(abs(np.degrees(da)) / seg_deg))
        for t in np.linspace(0, 1, k + 1):
            a = a1 + da * t
            out.append(C + r_eff * np.array([np.cos(a), np.sin(a)]))
    return np.array(out)


def signed_area(P):
    x, y = P[:, 0], P[:, 1]
    return 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)


def offset_polygon(P, d):
    """Per-vertex miter offset of a dense smooth CCW polygon; d>0 = inward. Keeps vertex count (for lofting)."""
    n = len(P); out = np.zeros_like(P)
    for i in range(n):
        a, v, b = P[i - 1], P[i], P[(i + 1) % n]
        e1 = v - a; e2 = b - v
        n1 = np.array([-e1[1], e1[0]]) / np.linalg.norm(e1)   # left normal = inward for CCW
        n2 = np.array([-e2[1], e2[0]]) / np.linalg.norm(e2)
        m = n1 + n2; m /= np.linalg.norm(m)
        out[i] = v + m * d / max(np.dot(m, n1), 0.3)
    return out


def hexagon(cx, cy, half_w, tip, side, tip_neck=None):
    """Bridge plateau: pointed at tail/neck, long sides parallel to strings. CCW from tail tip."""
    tn = tip if tip_neck is None else tip_neck
    return np.array([(cx, cy - tip), (cx + half_w, cy - side), (cx + half_w, cy + side),
                     (cx, cy + tn), (cx - half_w, cy + side), (cx - half_w, cy - side)], float)


def ring_triangulation(outer, inner, start_outer=0, start_inner=0):
    """Loft-style triangulation of the ring between two CCW polygons. Returns list of
    triangles as (('o',i),('o',j),('i',k)) style index tuples."""
    no, ni = len(outer), len(inner)
    i, j = 0, 0
    tris = []
    O = lambda k: outer[(start_outer + k) % no]
    I = lambda k: inner[(start_inner + k) % ni]
    while i < no or j < ni:
        adv_o = i < no and (j >= ni or np.linalg.norm(O(i + 1) - I(j)) <= np.linalg.norm(I(j + 1) - O(i)))
        if adv_o:
            tris.append(((('o', (start_outer + i) % no)), ('o', (start_outer + i + 1) % no), ('i', (start_inner + j) % ni)))
            i += 1
        else:
            tris.append((('o', (start_outer + i) % no), ('i', (start_inner + j + 1) % ni), ('i', (start_inner + j) % ni)))
            j += 1
    return tris
