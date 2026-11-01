"""Body geometry from Chuck's model (WyreCustomGuitars, Juris.dxf), shared by the C-75 and the C-75T.
Data comes from chuck/extract.py -> chuck_body.npz, already in the generators' body frame (tail at y = 0,
z = 0 at the back, x+ treble, mm).

The bevels are Chuck's: flat facets that fall 5/8" (15.9 mm) from the flat of the top (and of the back) to the
edge, leaving a 3/4" square band round the side. They run all round except across the middle of the tail, which
is square; toward the tail corners the flat comes out to meet the edge, so the bevel narrows to nothing there.
`bevel()` sweeps that surface as a band: at each station on the outline the cutter runs from the nearest point
of the flat (on the surface) out past the edge at the drop Chuck's model has there."""
import math, os
import numpy as np
import manifold3d as m3
from shapely.geometry import Polygon, LinearRing, Point
from shapely.ops import nearest_points
import build as B

_D = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'chuck_body.npz'))
T = float(_D['T']); DY = float(_D['DY'])
_O = np.asarray(_D['outl'], float)
if np.allclose(_O[0], _O[-1]): _O = _O[:-1]
OUTL = _O                                          # CCW, same start as the profile samples
RING = LinearRing(OUTL); L = RING.length
TOP_FLAT = Polygon(_D['top']).buffer(0); BACK_FLAT = Polygon(_D['back']).buffer(0)
# Chuck's flats are notched where his mortise is; fill that (the mortise or pocket is cut there anyway) so the bevel
# beside the neck end runs to the real edge rather than to the notch
def _fill_notch(F):
    from shapely.geometry import box
    return F.union(box(-29.0, 355.0 , 29.0, 460.0).intersection(Polygon(_O))).buffer(0.3, join_style=1).buffer(-0.3, join_style=1)
TOP_FLAT = _fill_notch(TOP_FLAT); BACK_FLAT = _fill_notch(BACK_FLAT)

def _fill(a):
    """circular linear interpolation over the gaps in a sampled profile"""
    a = np.array(a, float); n = len(a); idx = np.arange(n); ok = ~np.isnan(a)
    return np.interp(idx, np.r_[idx[ok] - n, idx[ok], idx[ok] + n], np.r_[a[ok], a[ok], a[ok]])
_SS = np.asarray(_D['ss'], float)
DROP_TOP = np.clip(_fill(-np.asarray(_D['zt'])), 0, 15.875)
DROP_BACK = np.clip(_fill(np.asarray(_D['zb']) + T), 0, 15.875)
def drop(s, side): return float(np.interp(s % L, _SS, DROP_TOP if side > 0 else DROP_BACK))

CHAMBER = np.asarray(_D['chamber'], float)
ROUTE_NECK = np.asarray(_D['route_neck'], float); ROUTE_BRIDGE = np.asarray(_D['route_bridge'], float)
COVER = np.asarray(_D['cover'], float)
POTS = np.asarray(_D['pots'], float); TOGGLE = np.asarray(_D['toggle'], float)[0]
JACK = np.asarray(_D['jack'], float)[0]; JACK_N = np.asarray(_D['jack_n'], float); JACK_Z = float(_D['jack_z'][0])
PU_NECK = np.asarray(_D['pu_neck_c'], float)[0]; PU_BRIDGE = np.asarray(_D['pu_bridge_c'], float)[0]
BRIDGE_POSTS = np.asarray(_D['bridge_posts'], float)
MORTISE = np.asarray(_D['mortise'], float); MORTISE_Z = float(_D['mortise_z'][0])
Y_END = float(max(p[1] for p in OUTL if abs(p[0]) < 30))          # body end at the neck
Y_SADDLE = float(BRIDGE_POSTS[0][1])                               # thumbwheel post line
CHAMBER_FLOOR = T - 42.86; CAP = 5.08; ROUTE_D = 17.78

class Body:
    """A bevelled body on `outl` with the flats `top`/`back` (shapely polygons). The drop at each station comes
    from the nearest point of Chuck's edge, scaled by how this body's bevel width compares with Chuck's there, so
    the slope carries over to an outline that differs (the C-75T's wider neck tongue)."""
    def __init__(self, outl, top, back, T_body, max_slope=None, rul_top=None, rul_back=None):
        self.OUTL = np.asarray(outl, float); self.RING = LinearRing(self.OUTL); self.L = self.RING.length
        self.TOP, self.BACK, self.T, self.max_slope = top, back, T_body, max_slope
        self.RUL = {+1: self._prep(rul_top, top), -1: self._prep(rul_back, back)}
    def _prep(self, rul, flat):
        """rulings sorted round the outline as (s_out, s_in) pairs, plus zero-width ones where the flat meets the edge
        (the square tail), so stations can be lofted between them the way Chuck's surfaces are"""
        if rul is None or len(rul) == 0: return None
        F = LinearRing(np.array(flat.exterior.coords)[:-1]); F = F if F.is_ccw else LinearRing(np.array(F.coords)[::-1])
        R = self.RING; pairs = []
        for xi, yi, xo, yo in rul:
            pairs.append((R.project(Point(xo, yo)), F.project(Point(xi, yi))))
        # zero-width rulings at the ends of the square runs
        ss = np.arange(0.0, self.L, 1.0)
        w = np.array([flat.exterior.distance(R.interpolate(x)) for x in ss]); sq = w < 0.3
        for i in range(len(ss)):
            if sq[i] != sq[i - 1]:
                x = ss[i] if sq[i] else ss[i - 1]; p = R.interpolate(x); pairs.append((x, F.project(p)))
        pairs.sort()
        return F, np.array(pairs)
    def _inner(self, s, side, flat, p):
        r = self.RUL[side]
        if r is None: return np.array(nearest_points(flat.exterior, Point(p))[0].coords[0])
        F, P = r; Lo, Lf = self.L, F.length; n = len(P)
        k = int(np.searchsorted(P[:, 0], s, side='right')) - 1
        a, b = P[k % n], P[(k + 1) % n]
        so = (b[0] - a[0]) % Lo; si = (b[1] - a[1]) % Lf
        if si > Lf / 2: si -= Lf
        f = ((s - a[0]) % Lo) / so if so > 1e-9 else 0.0
        return np.array(F.interpolate((a[1] + f * si) % Lf).coords[0])
    def _stations(self, side, step):
        flat = self.TOP if side > 0 else self.BACK; R, Lb = self.RING, self.L; out = []
        for s in np.arange(0.0, Lb, step):
            p = np.array(R.interpolate(s).coords[0]); a = np.array(R.interpolate((s - 0.8) % Lb).coords[0]); b = np.array(R.interpolate((s + 0.8) % Lb).coords[0])
            t = (b - a) / np.linalg.norm(b - a); n_in = np.array([-t[1], t[0]])
            q = self._inner(s, side, flat, p); w = float(np.linalg.norm(q - p))
            sc = RING.project(Point(p)); wc = float(np.interp(sc, _SS, W_TOP if side > 0 else W_BACK)); dc = drop(sc, side)
            d = 0.0 if wc < 0.3 else dc * min(1.0, w / wc)
            if self.max_slope: d = min(d, w * self.max_slope)
            out.append((s, p, n_in, q, w, d))
        return out
    def bevel(self, z_top, side, step=1.5):
        """cutter for the top (side = +1, z_top = T) or back (side = -1, z_top = 0) bevel"""
        pts = []; jit = np.random.default_rng(7 if side > 0 else 11)
        for s, p, n_in, q, w, d in self._stations(side, step):
            if w < 1.0 or d < 1.0: pts.append(None); continue          # stop just short of zero width: no knife-edge contacts
            r = (p - q) / w                                              # ruling direction, outward
            A = q + jit.normal(0, 0.004, 2) - r * 0.05; Bo = p + r * 3.0; zl = z_top - side * d * (w + 3.0) / w
            pts.append([(A[0], A[1], z_top + side * 0.001), (Bo[0], Bo[1], zl), (A[0], A[1], z_top + side * 60), (Bo[0], Bo[1], z_top + side * 60)])
        n = len(pts)
        return B.union([B.hull(pts[i] + pts[(i + 1) % n]) for i in range(n) if pts[i] is not None and pts[(i + 1) % n] is not None])
    def roundover(self, z_top, side, r=1.6, K=8, step=1.5):
        """display-only sanded edge where the bevel (or the flat, across the square tail) meets the side"""
        prof = []
        for s, p, n_in, q, w, d in self._stations(side, step):
            if w < 0.3 or d < 0.3: d = 0.0; u1 = np.array([n_in[0], n_in[1], 0.0])
            else:
                rr = (q - p) / w; u1 = np.array([rr[0] * w, rr[1] * w, side * d]); u1 /= np.linalg.norm(u1)
            E = np.array([p[0], p[1], z_top - side * d]); u2 = np.array([0, 0, -side * 1.0])
            al = math.acos(float(np.clip(u1 @ u2, -1, 1))); tl = r / math.tan(al / 2)
            m = (u1 + u2) / np.linalg.norm(u1 + u2); C = E + m * r / math.sin(al / 2)
            v1 = E + u1 * tl - C; v2 = E + u2 * tl - C
            arc = [C + (v1 * (1 - k / K) + v2 * (k / K)) / np.linalg.norm(v1 * (1 - k / K) + v2 * (k / K)) * r for k in range(K + 1)]
            arc[0] = arc[0] - (C - arc[0]) * 0.004; arc[-1] = arc[-1] - (C - arc[-1]) * 0.004
            prof.append((E - m * 0.4, arc))
        pieces = []; n = len(prof)
        for i in range(n):
            (e0, a0), (e1, a1) = prof[i], prof[(i + 1) % n]
            for k in range(K):
                pieces.append(B.hull([tuple(e0), tuple(a0[k]), tuple(a0[k + 1]), tuple(e1), tuple(a1[k]), tuple(a1[k + 1])]))
        return B.union(pieces)
    def blank(self, round_r=0.0):
        """the bevelled body blank (no routes); round_r > 0 adds the display roundover"""
        body = B.prism(self.OUTL, 0.0, self.T) - self.bevel(self.T, +1) - self.bevel(0.0, -1)
        if round_r > 0: body = body - self.roundover(self.T, +1, round_r) - self.roundover(0.0, -1, round_r)
        return body

W_TOP = np.array([TOP_FLAT.exterior.distance(RING.interpolate(s)) for s in _SS])
W_BACK = np.array([BACK_FLAT.exterior.distance(RING.interpolate(s)) for s in _SS])
RUL_TOP = np.asarray(_D['rul_top'], float) if 'rul_top' in _D.files else None
RUL_BACK = np.asarray(_D['rul_back'], float) if 'rul_back' in _D.files else None
CHUCK = Body(OUTL, TOP_FLAT, BACK_FLAT, T, rul_top=RUL_TOP, rul_back=RUL_BACK)
def blank(round_r=0.0): return CHUCK.blank(round_r)

def widened(dx=7.5, T_body=44.45):
    """Chuck's body with the treble side of the neck tongue and the cutaway eased outward by `dx` (blended to
    nothing toward the horn and down the body), so a standard 2.2" bolt-on pocket keeps a wall on the treble side;
    the bass side of the pocket already falls inside the horn. Outline and flats move together, so every bevel
    keeps Chuck's width and slope."""
    def sm(t): t = np.clip(t, 0, 1); return t * t * (3 - 2 * t)
    def f(P):
        P = np.array(P, float); x, y = P[:, 0], P[:, 1]
        m = (x > 0) * sm((x - 0.0) / 6.0) * sm((95.0 - x) / 50.0) * sm((y - 330.0) / 40.0)
        P[:, 0] = x + dx * m; return P
    o = Polygon(f(OUTL)).buffer(0)
    top = Polygon(f(np.array(TOP_FLAT.exterior.coords))).buffer(0); back = Polygon(f(np.array(BACK_FLAT.exterior.coords))).buffer(0)
    P = np.array(o.exterior.coords)[:-1]
    if not LinearRing(P).is_ccw: P = P[::-1]
    sr = lambda Rr: np.c_[f(Rr[:, :2]), f(Rr[:, 2:])]
    b = Body(P, top, back, T_body, rul_top=sr(RUL_TOP), rul_back=sr(RUL_BACK)); b.shift = f
    return b
