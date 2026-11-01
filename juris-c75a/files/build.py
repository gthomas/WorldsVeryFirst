"""Juris C-75A — faceted acoustic archtop on the C-75 outline.
Parametric generator. Units: mm. Body frame: x+ treble, y+ toward headstock (y=0 tail), z+ up (z=0 back edge).
Run:  python build.py   -> out/ (GLB, STL per part, DXF templates, report.json)
"""
import json, math, os
import numpy as np
import manifold3d as m3
import mapbox_earcut as earcut
from outline import (sharp_outline, fillet, offset_polygon, hexagon, signed_area, IN)

# ----------------------------------------------------------------------------- parameters
P = dict(
    scale=628.65,            # 24.75"
    bridge_y=190.5,          # saddle line, from tail (7.5")
    neck_angle=5.0,
    head_angle=13.0,
    # machining assumptions: 3-axis router, 1/4" (6.35) end mill smallest for pockets/profiles
    tool_r=3.5,              # minimum inside-corner radius on any pocket or cut-out
    # rim: one ring from a standard 1-3/4" mahogany body blank, all cut from the top face in one setup
    # (pocket inside the wall down to the back's glue ledge, cut through inside the ledge, profile)
    rim_single=True,
    plate_edge=3.2,          # top & back thickness at the edge
    rim_h=44.45,             # 1-3/4" blank -> 50.85 (2.0") edge depth, the same as the C-75
    rim_wall=8.0,            # also the top's glue land
    lining_w=12.0, lining_h=8.0,   # integral glue ledges, routed in
    # top (Port Orford cedar): faceted arch
    crease_r=25.0,           # radius blended into the inside (re-entrant side) of every top crease
    top_arch=11.0, plat_half_w=70.5,   # plateau fits the 135 mm Gretsch flat-bottom base (feet land on the ring beam)
    plat_tip=85.0, plat_tip_neck=70.0, plat_side=45.0,
    top_plateau_t=3.5, top_fall_t_at_plateau=2.8,   # lightened (study/run.py v4)
    beam_w=5.0, beam_h=6.0,
    arch_inset=22.0, arch_neck_y=340.0,   # faceted arch zone, flat skirt outside it (sits on the linings)
    slot_area=545.0, collar_h=12.0, collar_w=4.0, xbar_h=5.0,   # deeper collars (was 7.5) keep the air resonance ~127 Hz in the 1-3/4" box
    bar_x=46.0, bar_w=6.0, bar_h=5.0, bar_taper=35.0, bar_end_inset=14.0,   # two parallel tone bars
    # back: shallower faceted arch, uniform thickness
    back_arch=22.0,          # ridge depth below the back edge (raised from 12 to win back air with the 1-3/4" rim)
    arch_upper='follow', arch_upper_pts=((102.0, 359.0), (4.0, 330.0), (-52.0, 386.0), (-90.0, 380.0)),   # upper-bout arch corners: treble shoulder, cutaway, bass by the neck block, bass shoulder
    ridge_y=(90.0, 280.0),   # centre ridge; the back arch uses the same eight-sided footprint as the top's arch
    chamfer=3.5,
    # neck
    fb_above_edge=15.0,      # fingerboard top (centre) above top edge at the joint
    nut_w=42.86, fb_end_w=56.4, fb_radius=304.8, fb_t=7.3, fret_h=1.4, frets=22,
    neck_d1=21.8, neck_d12=24.4,  # '59 roundback total depth incl. fingerboard
    action12=2.0,
    pickup=dict(string_gap=4.0),
    tenon_w=38.0, tenon_len=42.0, tenon_depth=26.0, tenon_fit=0.1,
    truss_w=6.35, truss_d=9.5, truss_len=381.0,   # 15" dual-action rod, 1/4" channel
    fret_slot_w=0.6, fret_slot_d=3.0,
    tuner_hole=8.73,     # 11/32" for Grover 133 Vintage Deluxe press-in bushings (8.8 mm)
    dowel_d=6.0, dowel_depth=10.0,
)

P.update(json.loads(os.environ.get('JURIS_P', '{}')))   # parameter overrides for studies
OUT = os.path.join(os.path.dirname(__file__), 'out'); os.makedirs(OUT, exist_ok=True)
DEG = math.pi / 180

# ----------------------------------------------------------------------------- helpers
def mesh_manifold(V, F):
    return m3.Manifold(m3.Mesh(np.asarray(V, np.float32), np.asarray(F, np.uint32)))

def cap_tris(ring2d, offset, flip):
    idx = earcut.triangulate_float64(np.asarray(ring2d, np.float64), np.array([len(ring2d)], np.uint32)).reshape(-1, 3)
    if flip: idx = idx[:, ::-1]
    return idx + offset

def loft(rings, caps2d):
    """rings: list of (N,3) loops with matching vertex order; caps2d: 2D outlines of first/last ring."""
    n = len(rings[0]); k = len(rings)
    V = np.vstack(rings); S = []
    for s in range(k - 1):
        a, b = s * n, (s + 1) * n
        for i in range(n):
            j = (i + 1) % n
            S += [(a + i, a + j, b + j), (a + i, b + j, b + i)]
    S = np.array(S)
    def oriented_cap(c2d, off, away):
        idx = earcut.triangulate_float64(np.asarray(c2d, np.float64), np.array([len(c2d)], np.uint32)).reshape(-1, 3) + off
        nsum = np.zeros(3)
        for t in idx:
            nsum += np.cross(V[t[1]] - V[t[0]], V[t[2]] - V[t[0]])
        return idx if np.dot(nsum, away) > 0 else idx[:, ::-1]
    c0, c1 = rings[0].mean(0), rings[1].mean(0); cl, cp = rings[-1].mean(0), rings[-2].mean(0)
    C0 = oriented_cap(caps2d[0], 0, c0 - c1); C1 = oriented_cap(caps2d[1], (k - 1) * n, cl - cp)
    for side in (S, S[:, ::-1]):
        man = mesh_manifold(V, np.vstack([side, C0, C1]))
        if man.status() == m3.Error.NoError and not man.is_empty() and man.volume() > 0:
            return man
    raise RuntimeError('loft failed')

def prism(poly2d, z0, z1):
    cs = m3.CrossSection([np.asarray(poly2d, float)])
    return m3.Manifold.extrude(cs, z1 - z0).translate([0, 0, z0])

def hull(pts):
    return m3.Manifold.hull_points(np.asarray(pts, float))

def union(ms):
    return m3.Manifold.batch_boolean(ms, m3.OpType.Add)

def tri_prism(tri3, zfar):
    """convex prism between a 3D triangle and its projection at z=zfar"""
    pts = [p for p in tri3] + [(p[0], p[1], zfar) for p in tri3]
    return hull(pts)

def open_cs(cs, r=None):
    """round every convex corner of a pocket/hole region to the tool radius"""
    r = P['tool_r'] if r is None else r
    return cs.offset(-r, m3.JoinType.Round, circular_segments=96).offset(r, m3.JoinType.Round, circular_segments=96)

def open_poly(poly, r=None):
    cs = open_cs(m3.CrossSection([ccw(np.asarray(poly, float))]), r)
    return ccw(np.array(max(cs.to_polygons(), key=lambda q: abs(signed_area(np.array(q))))))

def ccw(Pts):
    return Pts if signed_area(Pts) > 0 else Pts[::-1]

# ----------------------------------------------------------------------------- outline & facets
names, SHARP = sharp_outline()
OUTL = fillet(names, SHARP)
OUTL = ccw(OUTL)
OUTL_CS = m3.CrossSection([OUTL])

# plateau assignment of each sharp outline vertex (indices into hexagon, CCW from tail tip)
ASSIGN = dict(tailR=[0], tTail=[0, 1], tLower=[1, 2], tWaist=[2], tUpper=[2, 3], tHorn=[3], tCut=[3],
              neckT=[3], neckB=[3], bCut=[3], bHorn=[3], bUpper=[3, 4], bWaist=[4], bLower=[4, 5],
              bTail=[5, 0], tailL=[0])

def facet_triangles(hexa, onames=None, assign=None):
    onames = onames or names; assign = assign or ASSIGN
    tris = []   # each: list of 3 (kind, idx)
    n = len(onames)
    for i, nm in enumerate(onames):
        A = assign[nm]
        for a, b in zip(A[:-1], A[1:]):
            tris.append([('o', i), ('p', a), ('p', b)])
        nxt = (i + 1) % n
        shared = A[-1]; assert shared == assign[onames[nxt]][0], (nm, onames[nxt])
        tris.append([('o', i), ('o', nxt), ('p', shared)])
    return tris

def _arch_polygon():
    """Simplified, symmetric arch zone: an octagon inset from the outline. Outside it the top is a flat
    skirt at edge level that lands on the rim and its glue ledges."""
    from shapely.geometry import Polygon, Point
    ins = Polygon(SHARP).buffer(-P['arch_inset'], join_style=2)
    C = np.array(ins.exterior.coords)[:-1]
    near = lambda nm: C[np.argmin(np.linalg.norm(C - SHARP[names.index(nm)], axis=1))]
    pts = {}
    for k in ('Tail', 'Lower', 'Waist', 'Upper'):
        t, b = near('t' + k), near('b' + k)
        x, y = (abs(t[0]) + abs(b[0])) / 2, (t[1] + b[1]) / 2
        pts['t' + k] = (x, y); pts['b' + k] = (-x, y)
    ys = [p[1] for p in C if abs(p[0]) < 1.0] or [min(C[:, 1])]
    y_tail = min(y for y in np.linspace(0, 60, 601) if ins.contains(Point(0, y + 0.5)))
    if P.get('arch_upper') == 'follow':      # carry the arch into the upper bout, following the outline (asymmetric there)
        (tu, tc, bn, bu) = [tuple(map(float, q)) for q in P['arch_upper_pts']]
        an = ['aTail', 'tTail', 'tLower', 'tWaist', 'tUpper', 'tCut', 'bNeck', 'bUpper', 'bWaist', 'bLower', 'bTail']
        A = np.array([(0.0, y_tail), pts['tTail'], pts['tLower'], pts['tWaist'], tu, tc, bn, bu,
                      pts['bWaist'], pts['bLower'], pts['bTail']], float)
        return an, A
    if P.get('arch_upper'):      # carry the arch on into the upper bout: ten sides
        an = ['aTail', 'tTail', 'tLower', 'tWaist', 'tUpper', 'aNeck', 'bUpper', 'bWaist', 'bLower', 'bTail']
        A = np.array([(0.0, y_tail)] + [pts['tTail'], pts['tLower'], pts['tWaist'], pts['tUpper']] + [(0.0, P['arch_neck_y'])]
                     + [pts['bUpper'], pts['bWaist'], pts['bLower'], pts['bTail']], float)
        return an, A
    an = ['aTail', 'tTail', 'tLower', 'tWaist', 'aNeck', 'bWaist', 'bLower', 'bTail']
    A = np.array([(0.0, y_tail)] + [pts['tTail'], pts['tLower'], pts['tWaist']] + [(0.0, P['arch_neck_y'])]
                 + [pts['bWaist'], pts['bLower'], pts['bTail']], float)
    return an, A

ARCH_NAMES, ARCH = _arch_polygon()
ARCH_ASSIGN = dict(aTail=[0], tTail=[0, 1], tLower=[1, 2], tWaist=[2, 3], aNeck=[3], bWaist=[3, 4],
                   bLower=[4, 5], bTail=[5, 0])
if P.get('arch_upper') == 'follow':
    ARCH_ASSIGN = dict(aTail=[0], tTail=[0, 1], tLower=[1, 2], tWaist=[2, 3], tUpper=[3], tCut=[3], bNeck=[3], bUpper=[3],
                       bWaist=[3, 4], bLower=[4, 5], bTail=[5, 0])
elif P.get('arch_upper'):
    ARCH_ASSIGN = dict(aTail=[0], tTail=[0, 1], tLower=[1, 2], tWaist=[2, 3], tUpper=[3], aNeck=[3], bUpper=[3],
                       bWaist=[3, 4], bLower=[4, 5], bTail=[5, 0])     # the waist is concave, so the upper bout fans from the plateau's neck tip

TOP_HEX = hexagon(0, P['bridge_y'], P['plat_half_w'], P['plat_tip'], P['plat_side'], P['plat_tip_neck'])
# back: single centre ridge. Footprint hexagon (flat skirt outside it), ridge line along the centreline.
_y1, _y2 = P['ridge_y']
BACK_FOOT = ARCH.copy()          # same footprint as the top's arch zone: broad, and it reads as one shape
BACK_HEX = np.array([(0, _y1), (0, _y2)], float)          # the ridge (degenerate plateau)
TOP_TRIS = facet_triangles(TOP_HEX, ARCH_NAMES, ARCH_ASSIGN)
def _back_tris(n):
    """ten- or eight-sided footprint, ridge p0 (tail end) to p1 (neck end). Tail and lower-bout corners
    fan to p0, the rest to p1; one long facet each side joins the lower bout to the far ridge end."""
    h = n // 2                       # index of the neck point
    T = [[('o', 0), ('o', 1), ('p', 0)], [('o', 1), ('o', 2), ('p', 0)], [('o', 2), ('p', 1), ('p', 0)]]
    T += [[('o', i), ('o', i + 1), ('p', 1)] for i in range(2, h)]
    T += [[('o', n - 1), ('o', 0), ('p', 0)], [('o', n - 2), ('o', n - 1), ('p', 0)], [('o', n - 2), ('p', 0), ('p', 1)]]
    T += [[('o', i - 1), ('o', i), ('p', 1)] for i in range(n - 2, h, -1)]
    return T
BACK_TRIS = _back_tris(len(ARCH))

def tri_pts(tri, hexa, z_outline, z_plat, outer=None):
    outer = SHARP if outer is None else outer
    return [(*outer[i], z_outline) if k == 'o' else (*hexa[i], z_plat) for k, i in tri]

# ----------------------------------------------------------------------------- Z levels
Z_BACK_IN = P['plate_edge']                       # top of back at edge = bottom of rim
Z_RIM_TOP = Z_BACK_IN + P['rim_h']                # 47.6
Z_EDGE = Z_RIM_TOP + P['plate_edge']              # 50.8
Z_PLAT = Z_EDGE + P['top_arch']                   # plateau top surface
Z_PLAT_IN = Z_PLAT - P['top_plateau_t']
FALL_IN_AT_PLAT = Z_PLAT - P['top_fall_t_at_plateau']

def top_z_outer(a): return Z_EDGE + a
def top_z_inner(a):
    t = P['plate_edge'] + (P['top_fall_t_at_plateau'] - P['plate_edge']) * a / P['top_arch']
    return Z_EDGE + a - t

BIG = 400.0
BINDING_W = 6.35     # undyed band around the top edge (faux binding), measured in plan from the outline
NECK_HALF_X = 28.0
HEEL_W = 2 * NECK_HALF_X - 3.0     # heel sits 1.5 mm inside the body's sides
HEEL_HOLD = 9.0
OUT_GROW = m3.CrossSection([OUTL]).offset(3.0, m3.JoinType.Round)

from outline import offset_polygon as _offp
SHARP_EXT = _offp(SHARP, -6.0)       # sharp outline pushed out 6 mm; flat apron covers concave fillets

def _terrain(tris, hexa, z_edge, z_plat_edge, z_plateau, z_far, outer=None):
    """Single closed mesh: faceted surface (+ flat apron out to SHARP_EXT) extruded to z_far."""
    import trimesh
    outer = SHARP if outer is None else outer
    nS, nH = len(outer), len(hexa); nE = len(SHARP_EXT)
    V = [(*p, z_edge) for p in outer]                 # 0..nS-1
    V += [(*p, z_plat_edge) for p in hexa]            # nS..nS+nH-1  (fall side of plateau edge)
    V += [(*p, z_plateau) for p in hexa]              # plateau level
    V += [(*p, z_edge) for p in SHARP_EXT]            # apron
    V += [(*p, z_far) for p in SHARP_EXT]             # far cap
    iS = lambda i: i; iHf = lambda i: nS + i; iHp = lambda i: nS + nH + i
    iE = lambda i: nS + 2 * nH + i; iF = lambda i: nS + 2 * nH + nE + i
    F = []
    for t in tris:
        F.append([iS(i) if k == 'o' else iHf(i) for k, i in t])
    for i in range(1, nH - 1): F.append([iHp(0), iHp(i), iHp(i + 1)])
    if abs(z_plateau - z_plat_edge) > 1e-9:
        for i in range(nH):
            j = (i + 1) % nH; F += [[iHf(i), iHf(j), iHp(j)], [iHf(i), iHp(j), iHp(i)]]
    # flat apron: SHARP_EXT with the facet zone as a hole
    ring_pts = np.vstack([SHARP_EXT, outer])
    idx = earcut.triangulate_float64(ring_pts.astype(np.float64), np.array([nE, nE + nS], np.uint32)).reshape(-1, 3)
    for t in idx:
        F.append([iE(k) if k < nE else iS(k - nE) for k in t])
    for i in range(nE):
        j = (i + 1) % nE
        F += [[iE(i), iF(i), iF(j)], [iE(i), iF(j), iE(j)]]
    F += cap_tris(SHARP_EXT, nS + 2 * nH + nE, False).tolist()
    tm = trimesh.Trimesh(np.array(V), np.array(F), process=False)
    tm.merge_vertices()
    trimesh.repair.fix_normals(tm)
    if tm.volume < 0: tm.invert()
    man = mesh_manifold(tm.vertices, tm.faces)
    assert man.status() == m3.Error.NoError and man.volume() > 0, man.status()
    return man

def terrain_below(tris, hexa, z_edge, z_plat_edge, z_plateau, slab_z=None):
    outer = ARCH if tris is TOP_TRIS else (BACK_FOOT if tris is BACK_TRIS else SHARP)
    return _terrain(tris, hexa, z_edge, z_plat_edge, z_plateau, -BIG, outer)

def terrain_above(tris, hexa, z_edge, z_plat_edge, z_plateau, slab_z=None):
    outer = ARCH if tris is TOP_TRIS else (BACK_FOOT if tris is BACK_TRIS else SHARP)
    return _terrain(tris, hexa, z_edge, z_plat_edge, z_plateau, BIG, outer)

# chamfer keep-solids (45 deg) for the top and bottom outer edges
def resample(Pts, step):
    Q = np.vstack([Pts, Pts[:1]]); seg = np.linalg.norm(np.diff(Q, axis=0), axis=1)
    s = np.r_[0, np.cumsum(seg)]; t = np.arange(0.3, s[-1], step)
    return np.c_[np.interp(t, s, Q[:, 0]), np.interp(t, s, Q[:, 1])]

def chamfer_keep(z_edge, up=True):
    """45-degree chamfer keep-solid. The cone starts 1 mm outside the outline (same 45-degree plane)
    so none of its faces coincide with the part walls."""
    c = P["chamfer"]; H = 2.5; g = 1.0
    RS = resample(OUTL, 1.37)            # vertex spacing deliberately unlike the part walls
    # keep the exact square corners at the neck end so the mitred offset stays clean
    Q = np.vstack([OUTL, OUTL[:1]]); seg = np.linalg.norm(np.diff(Q, axis=0), axis=1); sQ = np.r_[0, np.cumsum(seg)]
    seg2 = np.linalg.norm(np.diff(np.vstack([RS, RS[:1]]), axis=0), axis=1)
    for nm in ('neckT', 'neckB'):
        c_ = SHARP[names.index(nm)]
        k = int(np.argmin(np.linalg.norm(OUTL - c_, axis=1)))
        if np.linalg.norm(OUTL[k] - c_) > 1e-6: continue
        dists = np.linalg.norm(RS - c_, axis=1)
        keep = dists > c + H + 0.5        # points this close to a square corner would fold over when offset
        # insert the corner in perimeter order
        sR = np.r_[0, np.cumsum(seg2)][:-1]
        idx = np.searchsorted(sR, sQ[k])
        RS = np.insert(RS, idx, c_, axis=0); keep = np.insert(keep, idx, True)
        RS = RS[keep]
        seg2 = np.linalg.norm(np.diff(np.vstack([RS, RS[:1]]), axis=0), axis=1)
    # chamfer size per vertex: full size everywhere, fading out over FADE mm toward the neck end face
    a_, b_ = SHARP[names.index('neckB')], SHARP[names.index('neckT')]
    ab = b_ - a_
    def dist_end(p):
        t = np.clip(np.dot(p - a_, ab) / np.dot(ab, ab), 0, 1); return np.linalg.norm(p - (a_ + t * ab))
    FADE = 28.0
    f = np.array([np.clip((dist_end(p) - 1.0) / FADE, 0, 1) for p in RS]); f = f * f * (3 - 2 * f)
    ci = -0.3 + (c + 0.3) * f                      # -0.3 = plane passes just outside the corner: no cut
    outer = offset_polygon(RS, -g)
    inner = offset_polygon_var(RS, ci + H)
    zc = lambda P_, h: np.c_[P_, h if np.ndim(h) else np.full(len(P_), h)]
    if up:
        return loft([zc(outer, -BIG), zc(outer, z_edge - ci - g), zc(inner, z_edge + H), zc(inner, BIG)], [outer, inner])
    else:
        return loft([zc(inner, -BIG), zc(inner, z_edge - H), zc(outer, z_edge + ci + g), zc(outer, BIG)], [inner, outer])

def offset_polygon_var(Pts, d):
    n = len(Pts); out = np.zeros_like(Pts)
    for i in range(n):
        a, v, b = Pts[i - 1], Pts[i], Pts[(i + 1) % n]
        e1 = v - a; e2 = b - v
        n1 = np.array([-e1[1], e1[0]]) / np.linalg.norm(e1); n2 = np.array([-e2[1], e2[0]]) / np.linalg.norm(e2)
        m = n1 + n2; m /= np.linalg.norm(m)
        out[i] = v + m * d[i] / max(np.dot(m, n1), 0.3)
    return out

# the edge chamfer fades to nothing approaching the neck end face, so the heel lands on a square face
KEEP_TOP = chamfer_keep(Z_EDGE, True)
KEEP_BOT = chamfer_keep(0.0, False)
BODY_PRISM = m3.Manifold.extrude(OUTL_CS, 2 * BIG).translate([0, 0, -BIG])

# ----------------------------------------------------------------------------- soundholes
def line_through(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float); return a, b - a

def inset_line(a, b, d, toward):
    """line a->b shifted by d toward point 'toward'"""
    a, b = np.asarray(a, float), np.asarray(b, float); e = b - a
    n = np.array([-e[1], e[0]]) / np.linalg.norm(e)
    if np.dot(np.asarray(toward) - a, n) < 0: n = -n
    return a + n * d, e

def y_at_x(line, x):
    p, d = line; t = (x - p[0]) / d[0]; return p[1] + t * d[1]

def soundhole(side):
    """Mirror-image slots in the side facets beside the bridge plateau (plateau long side to the
    lower-bout corner of the arch), long sides parallel to the plateau edge, like f-holes flanking
    the bridge. Kept 9 mm inside every fold line."""
    from shapely.geometry import Polygon
    sg = 1 if side == 'treble' else -1
    LT, UT = TOP_HEX[1], TOP_HEX[2]
    tL = np.array(ARCH[ARCH_NAMES.index('tLower')])
    facet = Polygon([LT, UT, tL])
    inset = facet.buffer(-9.0, join_style=2)
    e = (UT - LT) / np.linalg.norm(UT - LT); n = np.array([e[1], -e[0]])
    if np.dot(tL - LT, n) < 0: n = -n
    band = Polygon([LT - e * 1e3 + n * 9, UT + e * 1e3 + n * 9, UT + e * 1e3 + n * 27, LT - e * 1e3 + n * 27])
    reg = inset & band
    if reg.area > P['slot_area']:            # trim from the lower end to the target area
        lo, hi = 0.0, 300.0
        for _ in range(40):
            mid = (lo + hi) / 2
            cut = Polygon([UT - e * mid + n * 1e3, UT - e * mid - n * 1e3, UT + e * 1e3 - n * 1e3, UT + e * 1e3 + n * 1e3])
            if (reg & cut).area > P['slot_area']: hi = mid
            else: lo = mid
        reg = reg & cut
    pts = np.array(reg.exterior.coords)[:-1]
    pts[:, 0] *= sg
    return ccw(pts)

SLOTS = [open_poly(soundhole('bass')), open_poly(soundhole('treble'))]

# ----------------------------------------------------------------------------- TOP
def bar_segments():
    """two parallel tone bars under the long sides of the plateau, running the length of the arch zone"""
    from shapely.geometry import Polygon, LineString
    arch = Polygon(ARCH).buffer(-P['bar_end_inset'])
    segs = []
    for sx in (-1, 1):
        x = sx * P['bar_x']
        g = LineString([(x, -100), (x, 600)]).intersection(arch)
        segs.append((np.array(g.coords[0]), np.array(g.coords[-1])))
    return segs

def top_inner_at(x, y):
    a = top_height_at(x, y) - Z_EDGE
    return top_z_inner(a)

def _scalloped_bar(a, b):
    L = np.linalg.norm(b - a); d = (b - a) / L; n = np.array([-d[1], d[0]])
    h, w, T = P['bar_h'], P['bar_w'], P['bar_taper']
    fp = np.array([a + n * w / 2, b + n * w / 2, b - n * w / 2, a - n * w / 2])
    up = terrain_below(TOP_TRIS, TOP_HEX, Z_RIM_TOP + 0.5, FALL_IN_AT_PLAT + 0.5, Z_PLAT_IN + 0.5)
    dn = terrain_below(TOP_TRIS, TOP_HEX, Z_RIM_TOP - h, FALL_IN_AT_PLAT - h, Z_PLAT_IN - h)
    bar = (prism(ccw(fp), -BIG, BIG) ^ up) - dn
    zi = lambda s_: top_inner_at(*(a + d * s_))
    prof = np.array([(0, zi(0) + 0.5), (T, zi(T) - h - 1), (T, -BIG), (L - T, -BIG), (L - T, zi(L - T) - h - 1),
                     (L, zi(L) + 0.5), (L + 1, BIG), (-1, BIG)])
    to3 = lambda t: np.array([(*(a + d * s_ + n * t), z) for s_, z in prof])
    return bar ^ loft([to3(-w), to3(w)], [prof, prof])

def cross_bar_segment():
    """upper-bout cross bar: butts the neck block's lower face and runs out to the bass upper bout,
    stopping 2 mm short of the glue ledge"""
    from shapely.geometry import Polygon, LineString
    S = np.array([-34.0, 389.0])
    d = np.array([-84.0, -46.0]); d /= np.linalg.norm(d)
    ledge = Polygon(OUTL).buffer(-(P['rim_wall'] + P['lining_w'] + 2.0))
    g = LineString([S, S + d * 300]).intersection(ledge)
    E = np.array(g.coords[-1]) if g.geom_type == 'LineString' else np.array(list(g.geoms)[0].coords[-1])
    return S, E

def build_cross_bar():
    a, b = cross_bar_segment()
    L = np.linalg.norm(b - a); d = (b - a) / L; n = np.array([-d[1], d[0]])
    h, w, T = P['xbar_h'], P['bar_w'], P['bar_taper']
    zt = Z_PLAT + 5.0
    prof = np.array([(0, zt), (0, Z_RIM_TOP - h), (L - T, Z_RIM_TOP - h), (L, Z_RIM_TOP - 0.01), (L, zt)])
    to3 = lambda t: np.array([(*(a + d * s_ + n * t), z) for s_, z in prof])
    env = loft([to3(-w / 2), to3(w / 2)], [prof, prof])
    # follow the underside of the top: h deep below it (where the arch reaches the upper bout it rises with it)
    up = 0.5
    raised = terrain_below(TOP_TRIS, TOP_HEX, Z_RIM_TOP + up, FALL_IN_AT_PLAT + up, Z_PLAT_IN + up, Z_RIM_TOP + up)
    deeper = terrain_below(TOP_TRIS, TOP_HEX, Z_RIM_TOP - h, FALL_IN_AT_PLAT - h, Z_PLAT_IN - h, Z_RIM_TOP - h)
    return env ^ (raised - deeper)

def build_bars():
    return union([_scalloped_bar(a, b) for a, b in bar_segments()] + [build_cross_bar()])

def collar_cs():
    """plan shape of the soundhole collars: a band round each slot, with any gap to the ring beam
    narrower than the 1/4 in cutter filled in (closing by the tool radius + a little)"""
    from shapely.geometry import Polygon
    beam_in = Polygon(offset_polygon(ccw(TOP_HEX), P['beam_w'])); beam = Polygon(TOP_HEX) - beam_in
    out = []
    for sl in SLOTS:
        s = Polygon(sl); band = s.buffer(P['collar_w'], join_style=1)
        r = 3.6
        u = (band | beam.buffer(-0.5, join_style=2)).buffer(r, join_style=1).buffer(-r, join_style=1)   # edges kept off the beam faces
        # run 1 mm into the beam rather than stopping on its face, so no faces coincide
        reg = (u - beam_in.buffer(1.0, join_style=2)).buffer(0) & s.buffer(P['collar_w'] + 12.0, join_style=1)
        reg = reg - s
        out.append(reg)
    return out

def build_collars(below_inner):
    """a lip hanging under each soundhole, collar_h below the inside surface. It lengthens the port,
    which lowers the air resonance without a deeper rim or smaller holes. Cut with the inside (setup B)."""
    h = P['collar_h']
    deeper = terrain_below(TOP_TRIS, TOP_HEX, Z_RIM_TOP - h, FALL_IN_AT_PLAT - h, Z_PLAT_IN - h, Z_RIM_TOP - h)
    up = 1.0     # overlap into the plate so the union has no coincident faces
    raised = terrain_below(TOP_TRIS, TOP_HEX, Z_RIM_TOP + up, FALL_IN_AT_PLAT + up, Z_PLAT_IN + up, Z_RIM_TOP + up)
    sheet = raised - deeper
    sols = []
    from shapely.geometry import Polygon, MultiPolygon
    for reg in collar_cs():
        polys = list(reg.geoms) if isinstance(reg, MultiPolygon) else [reg]
        for pg in polys:
            ext = ccw(np.array(pg.exterior.coords)[:-1])
            sol = prism(ext, -BIG, BIG)
            for hole in pg.interiors:
                sol = sol - prism(ccw(np.array(hole.coords)[:-1]), -BIG, BIG)
            sols.append(sol)
    return union(sols) ^ sheet

# ---- crease radii: roll a ball of radius crease_r over each faceted surface and fill only the re-entrant side
# of every fold (closing on the outer surface's valleys, opening on the inner surface's ridges). Convex sides stay
# crisp: the exact faceted solid is kept wherever the ball doesn't add material.
CREASE_H = 1.0
def _hf_sample(man, X, Y):
    """heights of the upward faces of a terrain solid on a grid"""
    mh = man.to_mesh(); V = np.array(mh.vert_properties)[:, :3]; Fc = np.array(mh.tri_verts)
    Z = np.full(X.shape, np.nan)
    for f in Fc:
        a, b, c = V[f]
        n = np.cross(b - a, c - a)
        if n[2] <= 1e-9: continue
        x0, x1 = min(a[0], b[0], c[0]) - 1e-6, max(a[0], b[0], c[0]) + 1e-6
        y0, y1 = min(a[1], b[1], c[1]) - 1e-6, max(a[1], b[1], c[1]) + 1e-6
        m_ = (X >= x0) & (X <= x1) & (Y >= y0) & (Y <= y1)
        if not m_.any(): continue
        px, py = X[m_], Y[m_]
        d = (b[1] - c[1]) * (a[0] - c[0]) + (c[0] - b[0]) * (a[1] - c[1])
        l1 = ((b[1] - c[1]) * (px - c[0]) + (c[0] - b[0]) * (py - c[1])) / d
        l2 = ((c[1] - a[1]) * (px - c[0]) + (a[0] - c[0]) * (py - c[1])) / d
        l3 = 1 - l1 - l2
        ins = (l1 >= -1e-9) & (l2 >= -1e-9) & (l3 >= -1e-9)
        z = l1 * a[2] + l2 * b[2] + l3 * c[2]
        cur = Z[m_]
        cur[ins] = np.where(np.isnan(cur[ins]), z[ins], np.maximum(cur[ins], z[ins]))
        Z[m_] = cur
    return Z

def _hf_solid(X, Y, Z, z_far):
    """closed solid between a height grid and a flat cap at z_far"""
    ny, nx = Z.shape; N = nx * ny
    V = np.vstack([np.c_[X.ravel(), Y.ravel(), Z.ravel()], np.c_[X.ravel(), Y.ravel(), np.full(N, z_far)]])
    idx = np.arange(N).reshape(ny, nx)
    a, b, c, d = idx[:-1, :-1].ravel(), idx[:-1, 1:].ravel(), idx[1:, 1:].ravel(), idx[1:, :-1].ravel()
    top = np.r_[np.c_[a, b, c], np.c_[a, c, d]]
    bot = np.r_[np.c_[a, c, b], np.c_[a, d, c]] + N
    ring = np.r_[idx[0, :], idx[1:, -1], idx[-1, -2::-1], idx[-2:0:-1, 0]]
    r2 = np.roll(ring, -1)
    side = np.r_[np.c_[ring, ring + N, r2], np.c_[r2, ring + N, r2 + N]]
    Fc = np.r_[top, bot, side]
    if z_far > Z.max(): Fc = Fc[:, ::-1]
    man = mesh_manifold(V, Fc)
    if man.volume() < 0 or man.is_empty():
        man = mesh_manifold(V, Fc[:, ::-1])
    return man

def _crease_grid():
    from shapely.geometry import Polygon as _P
    x0, y0, x1, y1 = _P(ARCH).buffer(14.0).bounds
    xs = np.arange(x0, x1 + CREASE_H, CREASE_H); ys = np.arange(y0, y1 + CREASE_H, CREASE_H)
    return np.meshgrid(xs, ys)

def _ball(R, slope_deg=28.0):
    r = int(math.ceil((R * math.sin(math.radians(slope_deg)) + 1.0) / CREASE_H))
    d = np.arange(-r, r + 1) * CREASE_H
    D = np.hypot(*np.meshgrid(d, d))
    fp = D <= r * CREASE_H
    return fp, np.where(fp, np.sqrt(np.maximum(R * R - D * D, 0.0)) - R, 0.0)

def round_creases(outer, inner, R, eps=0.02):
    """outer: solid below the outer surface; inner: void below the inner surface. Returns both with every
    re-entrant crease filled to radius R (adds wood only)."""
    from scipy import ndimage as nd
    X, Y = _crease_grid(); fp, st = _ball(R)
    Zo = _hf_sample(outer, X, Y); Zi = _hf_sample(inner, X, Y)
    Zo = np.where(np.isnan(Zo), np.nanmin(Zo), Zo); Zi = np.where(np.isnan(Zi), np.nanmin(Zi), Zi)   # outside the body: flat skirt level
    Zc = nd.grey_erosion(nd.grey_dilation(Zo, footprint=fp, structure=st, mode='nearest'), footprint=fp, structure=st, mode='nearest')
    Zp = nd.grey_dilation(nd.grey_erosion(Zi, footprint=fp, structure=st, mode='nearest'), footprint=fp, structure=st, mode='nearest')
    if P.get('crease_pure'):        # grid surfaces throughout the blend zone (clean input for the FEA mesher)
        rect = prism(np.array([(X.min(), Y.min()), (X.max(), Y.min()), (X.max(), Y.max()), (X.min(), Y.max())]), -BIG, BIG)
        outer2 = (outer - rect) + _hf_solid(X, Y, Zc, -BIG).simplify(0.005)
        inner2 = (inner - rect) + _hf_solid(X, Y, Zp, -BIG).simplify(0.005)
    else:
        outer2 = outer + _hf_solid(X, Y, Zc - eps, -BIG).simplify(0.005)
        inner2 = inner - _hf_solid(X, Y, Zp + eps, BIG).simplify(0.005)
    global CREASE_INFO
    CREASE_INFO = dict(radius_mm=R, outer_added_max_mm=float(np.nanmax(Zc - Zo)), inner_added_max_mm=float(np.nanmax(Zi - Zp)))
    return outer2, inner2
CREASE_INFO = {}

def build_top(chamfer=True):
    below_outer = terrain_below(TOP_TRIS, TOP_HEX, Z_EDGE, Z_PLAT, Z_PLAT, Z_EDGE)
    below_inner = terrain_below(TOP_TRIS, TOP_HEX, Z_RIM_TOP, FALL_IN_AT_PLAT, Z_PLAT_IN, Z_RIM_TOP)
    if P.get('crease_r'): below_outer, below_inner = round_creases(below_outer, below_inner, P['crease_r'])
    shell = (BODY_PRISM ^ below_outer) - below_inner
    shell = shell ^ prism(OUTL, Z_RIM_TOP - 0.001, BIG)
    beam_in = open_poly(offset_polygon(ccw(TOP_HEX), P['beam_w']))
    beam = prism(TOP_HEX, Z_PLAT_IN - P['beam_h'] + 0.0, Z_PLAT_IN + 0.5) - prism(beam_in, -BIG, BIG)
    top = shell + beam + build_bars() + build_collars(below_inner)
    for sl in SLOTS:
        top = top - prism(sl, -BIG, BIG)
    top = top - prism(tenon_poly(P['tenon_fit']), -BIG, BIG)      # notch the tenon drops through
    wx, wy, wd = WIRE_HOLE
    top = top - m3.Manifold.cylinder(2 * BIG, wd / 2, wd / 2, 32).translate([wx, wy, -BIG])
    return top ^ KEEP_TOP if chamfer else top

# ----------------------------------------------------------------------------- BACK
def build_back(chamfer=True):
    a = P['back_arch']
    above_outer = terrain_above(BACK_TRIS, BACK_HEX, 0.0, -a, -a, 0.0)
    above_inner = terrain_above(BACK_TRIS, BACK_HEX, Z_BACK_IN, Z_BACK_IN - a, Z_BACK_IN - a, Z_BACK_IN)
    back = (BODY_PRISM ^ above_outer) - above_inner
    back = back ^ prism(OUTL, -BIG, Z_BACK_IN + 0.001)
    return back ^ KEEP_BOT if chamfer else back

# ----------------------------------------------------------------------------- RIM
Y_NECK_END = SHARP[names.index('neckT')][1]

NECK_BLOCK = np.array([(-44, 392), (44, 392), (44, 470), (-44, 470)], float)
TAIL_BLOCK = np.array([(-34, -5), (34, -5), (34, 24), (-34, 24)], float)
Z_SPLIT = Z_BACK_IN + P['rim_h'] / 2
DOWELS = [(0.0, 399.8), (-34.0, 405.0), (-20.0, 10.0), (20.0, 10.0)]

def rim_regions():
    blocks = (m3.CrossSection([NECK_BLOCK]) + m3.CrossSection([TAIL_BLOCK])) ^ OUTL_CS
    through = open_cs(OUTL_CS.offset(-(P['rim_wall'] + P['lining_w']), m3.JoinType.Round, circular_segments=96) - blocks)
    pad = m3.CrossSection([ctrl_box(CTRL_PAD_L, -30.0, 1.0)]) + m3.CrossSection([jack_box(JACK_PAD_L, -30.0, 1.0)])
    # wall left full ledge-thick behind the control plate, the battery box and the jack
    recess = open_cs(OUTL_CS.offset(-P['rim_wall'], m3.JoinType.Round, circular_segments=96) - blocks - pad)
    return through, recess

def build_rim_halves(machining=True):
    """lower and upper rim halves. Each is a plain 2.5D part machined from its split face:
    pocket the wall recess down to the glue ledge, then cut through inside the ledge, then profile."""
    through, recess = rim_regions()
    hh = P['rim_h'] / 2; lh = P['lining_h']; zs = Z_SPLIT
    ring = m3.Manifold.extrude(OUTL_CS - through, hh)
    pocket = m3.Manifold.extrude(recess, hh - lh)
    dow = lambda z0: union([m3.Manifold.cylinder(P['dowel_depth'], P['dowel_d'] / 2 + 0.05, P['dowel_d'] / 2 + 0.05, 32)
                            .translate([x, y, z0]) for x, y in DOWELS])
    lower = ring.translate([0, 0, Z_BACK_IN]) - pocket.translate([0, 0, Z_BACK_IN + lh]) - dow(zs - P['dowel_depth'])
    upper = ring.translate([0, 0, zs]) - pocket.translate([0, 0, zs]) - dow(zs) - MORTISE
    # control window: a notch through the wall in each half, cut from the split face (3-axis)
    holes = union([m3.Manifold.cylinder(14, 1.0, 1.0, 16).rotate([0, 90, 0]).rotate([0, 0, math.degrees(math.atan2(CTRL_N[1], CTRL_N[0]))])
                   .translate([*(CTRL_C + CTRL_E * se + CTRL_N * -12.0), zs + sz]) for se, sz in CTRL_SCREWS])
    lower = lower - ctrl_cutters() - jack_cutter(); upper = upper - ctrl_cutters() - jack_cutter()          # plate screw pilots are drilled by hand into the rim face
    return lower, upper

def build_rim_single():
    """one ring from a standard 1-3/4" blank, machined entirely from the top face: pocket the wall recess down to
    the back's glue ledge (open at the top, so the top glues to the full wall), cut through inside the ledge, cut
    the open-top neck mortise, profile. The control-plate window is routed through the wall from outside with a
    template after machining."""
    through, recess = rim_regions()
    lh = P['lining_h']; z0 = Z_BACK_IN; h = P['rim_h']
    ring = m3.Manifold.extrude(OUTL_CS - through, h).translate([0, 0, z0])
    pocket = m3.Manifold.extrude(recess, h - lh + 1.0).translate([0, 0, z0 + lh])
    return ring - pocket - MORTISE - ctrl_cutters() - jack_cutter()

def build_rim():
    if P.get('rim_single'):
        rim = build_rim_single()
        return rim ^ KEEP_TOP ^ KEEP_BOT
    lower, upper = build_rim_halves()
    rim = lower + upper
    # endpin bore: drilled by hand after glue-up (the jack lives on the control plate)
    rim = rim - m3.Manifold.cylinder(40, 6.0, 6.0, 48).rotate([90, 0, 0]).translate([0, 20, Z_BACK_IN + P['rim_h'] / 2 - 4])
    return rim ^ KEEP_TOP ^ KEEP_BOT

# ----------------------------------------------------------------------------- NECK FRAME
TH = P['neck_angle'] * DEG
FJ = Z_EDGE + P['fb_above_edge']           # fb top centreline at the joint (u = 0)

def n2b(u, v, w):
    """neck frame (u along neck from joint toward head, v=x, w up from fb plane) -> body"""
    y = Y_NECK_END + u * math.cos(TH) + w * math.sin(TH)
    z = FJ - u * math.sin(TH) + w * math.cos(TH)
    return np.array([v, y, z])

def b2n(y, z):
    dy, dz = y - Y_NECK_END, z - FJ
    u = dy * math.cos(TH) - dz * math.sin(TH)
    w = dy * math.sin(TH) + dz * math.cos(TH)
    return u, w

def NB(pts):
    return np.array([n2b(*p) for p in pts])

W_NUT_STR = 1.9                  # string height over fb plane at nut (fret + 0.5)
W_SAD = 2 * (P['fret_h'] + P['action12']) - W_NUT_STR       # from 12th-fret action
U_SAD = ((P['bridge_y'] - Y_NECK_END) - W_SAD * math.sin(TH)) / math.cos(TH)
U_NUT = U_SAD + math.sqrt(P['scale'] ** 2 - (W_SAD - W_NUT_STR) ** 2)
SADDLE = n2b(U_SAD, 0, W_SAD)
NUT_PT = n2b(U_NUT, 0, W_NUT_STR)

def fret_s(n): return P['scale'] * (1 - 2 ** (-n / 12))
FRET_U = [U_NUT - fret_s(n) for n in range(1, P['frets'] + 1)]
U_FB_END = FRET_U[-1] - 6.0
def fb_width(u): return P['nut_w'] + (U_NUT - u) * (P['fb_end_w'] - P['nut_w']) / (U_NUT - U_FB_END)
def crown(v): R = P['fb_radius']; return -(R - math.sqrt(R * R - v * v))

def section_loop(u, W, d, n_exp, w_top, npts=40):
    """'59 roundback style section below w_top: superellipse half + flat top. CCW in (v,w) seen from +u."""
    pts = []
    for t in np.linspace(0, math.pi, npts):
        c, s = math.cos(t), math.sin(t)
        v = (W / 2) * math.copysign(abs(c) ** (2 / n_exp), c)
        w = w_top - d * abs(s) ** (2 / n_exp)
        pts.append((u, v, w))
    return pts       # from +v, under, to -v (top edge closes it)

def build_neck():
    fbb = -P['fb_t']
    def depth(u):
        s = U_NUT - u
        d = P['neck_d1'] + (P['neck_d12'] - P['neck_d1']) * (s - fret_s(1)) / (fret_s(12) - fret_s(1))
        d = d - P['fb_t']
        heel_start = 70.0
        u = max(u - HEEL_HOLD, 0.0)          # hold full heel depth until the heel bottom has cleared the body face
        if u < heel_start:
            target = FJ / math.cos(TH) - P['fb_t'] + 3.0   # a little proud: trimmed flush with the back plane
            k = (1 - u / heel_start) ** 2.2
            d = d + (target - d) * k
        return d
    stations = list(np.linspace(-14.0, 0.0, 5))[:-1] + list(np.linspace(0.0, 70.0, 36)) + list(np.linspace(73, U_NUT - 25, 70))
    rings = []
    for u in stations:
        uc = max(u, 0.0)
        W = fb_width(uc)
        ub = max(uc - HEEL_HOLD, 0.0)
        if ub < 70: W = W + max(0.0, HEEL_W - W) * (1 - ub / 70) ** 2
        nexp = 2.25 + 0.9 * max(0.0, 1 - ub / 70) ** 2      # keep the '59 round carve all the way into the heel
        d = depth(uc) + (6.0 if u < 0 else 0.0) * min(1.0, -u / 3.0)
        rings.append(NB(section_loop(u, W, d, nexp, fbb, npts=96)))
    neck = loft(rings, [np.array([(p[1], p[2]) for p in [b2nv(q) for q in rings[0]]]),
                        np.array([(p[1], p[2]) for p in [b2nv(q) for q in rings[-1]]])])
    neck = (neck ^ heel_keep()) + build_tenon() + build_tongue()
    neck = neck + build_headstock_wood(holes=False)
    neck = neck - truss_channel()
    return neck

TONGUE_Y0 = 392.0       # tongue runs from the heel face back to the end of the neck block
def build_tongue():
    """Neck tongue under the fingerboard extension: fills the space between the top and the board
    over the neck block, glued to both. The fingerboard overhangs past it toward the pickup."""
    pts = []
    for y in (TONGUE_Y0, Y_NECK_END + 4.0):
        u, _ = b2n(y, Z_EDGE + 8.0)
        hw = fb_width(u) / 2 - 0.5
        for x in (-hw, hw):
            pts += [(x, y, Z_EDGE), (x, y, fb_bottom_z(y))]
    return hull(pts)

def heel_keep():
    """Below the rim top the heel face is the body's vertical end face (y = Y_NECK_END) and the heel
    bottom is flush with the back; above the rim top the face is square to the fingerboard plane, so
    nothing on the neck overhangs the tenon (3-axis friendly in the neck's own frame)."""
    Yn = Y_NECK_END
    lower = m3.Manifold.cube([400, 900, Z_RIM_TOP], False).translate([-200, Yn, 0]) + \
            m3.Manifold.cube([400, 900, Z_RIM_TOP + 300], False).translate([-200, Yn + 150, -300])
    uj, _ = b2n(Yn, Z_RIM_TOP)
    upper = hull([n2b(u, v, w) for u in (uj, uj + 1200) for v in (-200, 200) for w in (-600, 600)])
    upper = upper ^ m3.Manifold.cube([400, 2000, 400], False).translate([-200, -500, Z_RIM_TOP])
    return lower + upper

def nhull(pts):
    """hull of points given in the neck frame (u, v, w)"""
    return hull([n2b(*p) for p in pts])

def truss_channel():
    """1/4\" channel, 9.5 mm below the fingerboard plane; runs out through the headstock face under the
    truss rod cover so the adjuster is reachable."""
    w, d, L = P['truss_w'], P['truss_d'], P['truss_len']
    u1 = U_H0 + 34.0 * math.cos(HA); u0 = u1 - L
    r = w / 2; wb = -P['fb_t'] - d
    pts = []
    for uc in (u0 + r, u1 - r):
        for t in np.linspace(0, 2 * math.pi, 33)[:-1]:
            for ww in (wb, wb + 70):
                pts.append((uc + r * math.cos(t), r * math.sin(t), ww))
    return nhull(pts)

def b2nv(p):
    u, w = b2n(p[1], p[2]); return (u, p[0], w)

def tenon_poly(grow=0.0):
    """tenon plan (body x,y): runs from the heel face into the neck block, far corners at tool radius"""
    w = P['tenon_w'] / 2 + grow; y1 = Y_NECK_END + 2.0; y0 = Y_NECK_END - P['tenon_len'] - grow
    r = P['tool_r']
    pts = [(w, y1), (-w, y1)]
    for cx, a0 in [(-w + r, 180), (w - r, 270)]:
        for t in np.linspace(a0, a0 + 90, 10):
            pts.append((cx + r * math.cos(math.radians(t)), y0 + r + r * math.sin(math.radians(t))))
    return ccw(np.array(pts))

def build_tenon():
    return prism(tenon_poly(), Z_RIM_TOP - P['tenon_depth'], Z_EDGE + 0.5)   # passes up through a notch in the top into the tongue

def mortise():
    """open-top pocket in the upper rim half's neck block; floor in the upper half, walls vertical"""
    return prism(tenon_poly(P['tenon_fit']), Z_RIM_TOP - P['tenon_depth'], Z_RIM_TOP + 5) + \
           prism(np.array([(-P['tenon_w'] / 2 - P['tenon_fit'], Y_NECK_END - 1), (P['tenon_w'] / 2 + P['tenon_fit'], Y_NECK_END - 1),
                           (P['tenon_w'] / 2 + P['tenon_fit'], Y_NECK_END + 5), (-P['tenon_w'] / 2 - P['tenon_fit'], Y_NECK_END + 5)]),
                 Z_RIM_TOP - P['tenon_depth'], Z_RIM_TOP + 5)

# headstock frame: origin at nut back face on the shaft top, face tilted back head_angle
HA = P['head_angle'] * DEG
U_H0 = U_NUT + 5.0
def h2b(v, s, n):
    u = U_H0 + s * math.cos(HA) + n * math.sin(HA)
    w = -P['fb_t'] + 0.0 - s * math.sin(HA) + n * math.cos(HA)
    return n2b(u, v, w)

# Juris headstock: straight sides flaring from the nut, gabled top, bevelled face edge (reference render)
import head_c75 as _H   # the built C-75's headstock: waisted sides, gabled top (shared by all three models)
HEAD_POLY = _H.HEAD_OUTLINE.copy()
HEAD_T = 14.0
HEAD_PLATE_S0 = 55.0     # plate solid starts here; the volute loft forms the headstock below this
HEAD_BEVEL = 0.0   # flat face, square edges: takes standard tuners
_HR = np.array(_H._RIGHT[:-1])
def _head_edge_v(s): return float(np.interp(s, _HR[:, 1], _HR[:, 0]))
TUNER_POS = list(_H.TUNER_POS)
TUNER_LIST = [(sv * v, s, sv) for v, s in TUNER_POS for sv in (-1, 1)]

def head_slab(poly, n0, n1):
    pts = []
    for v, s in poly:
        pts.append(h2b(v, s, n0)); pts.append(h2b(v, s, n1))
    return pts

def build_headstock_wood(holes=True):
    # plate as a loft of the outline at two normal offsets (non-convex safe)
    from shapely.geometry import Polygon, box
    poly = ccw(np.array(Polygon(HEAD_POLY).intersection(box(-200, HEAD_PLATE_S0, 200, 400)).exterior.coords)[:-1])
    ring = lambda pp, n: np.array([h2b(v, s, n) for v, s in pp])
    plate = loft([ring(poly, -HEAD_T), ring(poly, 0.0)], [poly, poly])   # base of the headstock comes from the volute
    tr = build_volute()
    if not holes:
        return plate + tr
    return (plate + tr) - tuner_holes()

def tuner_holes():
    r = P['tuner_hole'] / 2
    return union([m3.Manifold.cylinder(60, r, r, 32).translate([0, 0, -30]).transform(head_mat(v, s))
                  for v, s, _ in TUNER_LIST])

def build_volute():
    """Smooth lofted transition from the neck shaft into the angled headstock, with a volute bump on the
    back just behind the nut. Every section stays inside the neck/headstock outline (no side wings)."""
    fbt = P['fb_t']; sH, cH = math.sin(HA), math.cos(HA)
    def w_neck_back(u):
        s_ = U_NUT - u
        d = P['neck_d1'] + (P['neck_d12'] - P['neck_d1']) * (s_ - fret_s(1)) / (fret_s(12) - fret_s(1))
        return -fbt - (d - fbt)
    def w_head_back(u):
        s_ = (u - U_H0 + HEAD_T * sH) / cH
        return -fbt - s_ * sH - HEAD_T * cH
    def smin(a, b, k):
        return -k * math.log(math.exp(-a / k) + math.exp(-b / k))
    def w_face(u):
        return -fbt if u <= U_H0 else -fbt - (u - U_H0) / cH * sH
    rings, caps = [], []
    for u in np.linspace(U_NUT - 35, U_H0 + 62, 110):
        top = w_face(u) - 0.02
        k = 0.05 + 5.0 * math.exp(-((u - U_H0) / 22.0) ** 2)          # fillet only around the break
        tail = float(np.clip((u - (U_H0 + 18)) / 20.0, 0, 1))          # finish 0.4 mm inside the headstock back
        bot = (smin(w_neck_back(u), w_head_back(u), k) - 4.5 * math.exp(-((u - (U_H0 + 8)) / 16.0) ** 2)
               + 0.02 * tail * tail * (3 - 2 * tail))
        if u <= U_H0:
            half = fb_width(min(u, U_NUT)) / 2
        else:
            half = _head_edge_v((u - U_H0) / cH)
        W = 2 * half - 0.04                                  # matches the shaft, then the headstock sides
        t = float(np.clip((u - (U_H0 + 4)) / 26.0, 0, 1)); t = t * t * (3 - 2 * t)
        loop = section_loop(u, W, top - bot, 2.25 + 45.0 * t ** 3, top, npts=96)   # round neck -> square headstock section
        rings.append(NB(loop))
    caps = [np.array([(q[1], q[2]) for q in map(b2nv, r)]) for r in (rings[0], rings[-1])]
    return loft(rings, caps)

def head_mat(v, s):
    o = h2b(v, s, 0); ex = h2b(v + 1, s, 0) - o; ey = h2b(v, s + 1, 0) - o; ez = h2b(v, s, 1) - o
    return np.array([ex, ey, ez, o]).T.astype(float)  # 3x4

def build_faceplate():
    poly = ccw(HEAD_POLY)
    ring = lambda pp, n: np.array([h2b(v, s, n) for v, s in pp])
    fp = loft([ring(poly, 0.0), ring(poly, 1.6)], [poly, poly])
    return fp - tuner_holes() - _face_part(TRUSS_ACCESS, -1.0, 3.0)

TRUSS_ACCESS = np.array([(-5.0, -1.0), (5.0, -1.0), (5.0, 33.0), (-5.0, 33.0)])

def _face_part(poly2d, n0, n1):
    poly2d = ccw(np.asarray(poly2d, float))
    ring = lambda n: np.array([h2b(v, s, n) for v, s in poly2d])
    return loft([ring(n0), ring(n1)], [poly2d, poly2d])

def build_truss_cover():
    return _face_part(TRUSS_COVER, 1.6, 2.7)

TRUSS_COVER = np.array([(-8.0, 1.0), (8.0, 1.0), (8.0, 31.0), (0.0, 47.0), (-8.0, 31.0)])

def build_badge():
    """riveted asset tag: brushed plate with a QR-style field (pattern is decorative)"""
    plate = _face_part([(-13.5, 140), (13.5, 140), (13.5, 178.5), (-13.5, 178.5)], 1.5, 2.3)
    rivets = union([m3.Manifold.cylinder(0.6, 1.1, 1.1, 12).transform(head_mat(v, s)).translate([0, 0, 0])
                    .translate((h2b(v, s, 2.3) - h2b(v, s, 0)).tolist()) for v in (-11.5, 11.5) for s in (142, 176.5)])
    rng = np.random.default_rng(75); N = 21; m = 1.05; v0, s0 = -N * m / 2, 143.0
    grid = rng.random((N, N)) < 0.5
    for (r0, c0) in [(0, 0), (0, N - 7), (N - 7, 0)]:
        grid[r0:r0 + 7, c0:c0 + 7] = False
        for i in range(7):
            for j in range(7):
                if i in (0, 6) or j in (0, 6) or (2 <= i <= 4 and 2 <= j <= 4): grid[r0 + i, c0 + j] = True
    mods = []
    for i in range(N):
        j = 0
        while j < N:
            if grid[i, j]:
                k = j
                while k + 1 < N and grid[i, k + 1]: k += 1
                v_a, v_b = v0 + j * m, v0 + (k + 1) * m; s_a, s_b = s0 + i * m, s0 + (i + 1) * m
                mods.append(_face_part([(v_a, s_a), (v_b, s_a), (v_b, s_b), (v_a, s_b)], 2.25, 2.45))
                j = k + 1
            else:
                j += 1
    return plate + rivets, union(mods)

def build_fingerboard():
    rings = []
    for u in np.linspace(U_FB_END, U_NUT, 8):
        W = fb_width(u); loop = []
        for v in np.linspace(W / 2, -W / 2, 48): loop.append((u, v, crown(v)))
        loop += [(u, -W / 2, -P['fb_t']), (u, W / 2, -P['fb_t'])]
        rings.append(NB(loop))
    caps = [np.array([(q[1], q[2]) for q in map(b2nv, r)]) for r in (rings[0], rings[-1])]
    return loft(rings, caps) - fret_slots()

def fret_slots():
    """0.6 mm slots, 3 mm deep at the board edge (deeper at the centre), flat-bottomed, full width"""
    out = []
    hw = P['fret_slot_w'] / 2
    for u in FRET_U:
        W = fb_width(u)
        wb = crown(W / 2) - P['fret_slot_d']
        out.append(nhull([(u + du, v, ww) for du in (-hw, hw) for v in (-W / 2 - 2, W / 2 + 2) for ww in (wb, 5.0)]))
    return union(out)

def build_frets():
    out = []
    for u in FRET_U:
        W = fb_width(u) - 0.6; loop_hi, loop_lo = [], []
        for v in np.linspace(W / 2, -W / 2, 20):
            loop_hi.append((v, crown(v) + P['fret_h'])); loop_lo.append((v, crown(v) - 0.3))
        loop = loop_hi + loop_lo[::-1]
        r0 = NB([(u - 1.2, v, w) for v, w in loop]); r1 = NB([(u + 1.2, v, w) for v, w in loop])
        out.append(loft([r0, r1], [np.array(loop), np.array(loop)]))
    return union(out)

def build_nut():
    W = P['nut_w']; pts = []
    for u in (U_NUT, U_NUT + 5.0):
        for v in (-W / 2, W / 2):
            pts += [(u, v, -P['fb_t']), (u, v, W_NUT_STR - 0.3)]
    return hull(NB(pts))

def build_side_dots():
    out = []
    for n in (3, 5, 7, 9, 12, 15, 17, 19, 21):
        u = (U_NUT - fret_s(n - 1) + U_NUT - fret_s(n)) / 2 if n > 1 else 0
        W = fb_width(u)
        offs = [-3.5, 3.5] if n == 12 else [0]
        for o in offs:
            c = m3.Manifold.cylinder(2.0, 1.1, 1.1, 16).rotate([0, 90, 0])   # axis along x
            p = n2b(u + o, -W / 2 - 0.6, -2.4)
            out.append(c.translate(p))
    return union(out)

# ----------------------------------------------------------------------------- HARDWARE
def build_tuners():
    parts = []
    for v_, s, sv in TUNER_LIST:
        if True:
            M = head_mat(v_, s)
            post = m3.Manifold.cylinder(17, 3.1, 3.1, 24).translate([0, 0, 1.6])
            bush = m3.Manifold.cylinder(3.5, 5.0, 5.0, 24).translate([0, 0, 1.6])
            house = m3.Manifold.cube([13, 20, 12], True).translate([0, 0, -HEAD_T - 6])
            shaft = m3.Manifold.cylinder(14, 1.6, 1.6, 12).rotate([0, 90 * sv, 0]).translate([0, 0, -HEAD_T - 6])
            button = hull([(sv * (14 + a), 9 * math.cos(t), 9 * math.sin(t) - HEAD_T - 6)
                           for a in (0, 5) for t in np.linspace(0, 2 * math.pi, 9)[:-1] + math.pi / 8])
            parts.append((post + bush + house + shaft + button).transform(M))
    return union(parts)

def top_height_at(x, y):
    """outer top surface z at (x,y) via ray cast on the facet terrain (analytic)"""
    pt = np.array([x, y])
    if point_in_poly(pt, TOP_HEX): return Z_PLAT
    for t in TOP_TRIS:
        tri = np.array(tri_pts(t, TOP_HEX, Z_EDGE, Z_PLAT, ARCH))
        l = bary(pt, tri[:, :2])
        if l is not None: return float(l @ tri[:, 2])
    return Z_EDGE

def bary(p, T):
    a, b, c = T; v0, v1, v2 = b - a, c - a, p - a
    d = v0[0] * v1[1] - v1[0] * v0[1]
    if abs(d) < 1e-9: return None
    l1 = (v2[0] * v1[1] - v1[0] * v2[1]) / d; l2 = (v0[0] * v2[1] - v2[0] * v0[1]) / d
    l0 = 1 - l1 - l2
    if min(l0, l1, l2) < -1e-9: return None
    return np.array([l0, l1, l2])

def point_in_poly(p, poly):
    x, y = p; inside = False; n = len(poly)
    for i in range(n):
        x1, y1 = poly[i]; x2, y2 = poly[(i + 1) % n]
        if (y1 > y) != (y2 > y) and x < x1 + (y - y1) * (x2 - x1) / (y2 - y1): inside = not inside
    return inside

def string_z_at(y):
    t = (y - SADDLE[1]) / (NUT_PT[1] - SADDLE[1]); return SADDLE[2] + t * (NUT_PT[2] - SADDLE[2])

BRIDGE_FIT = 0.0
BASE = dict(L=135.0, W=14.95, H=11.75, post=37.0)   # Gretsch 008-0636-200 flat-bottom rosewood base, posts at 74 mm
def build_bridge():
    """Stand-in for the purchased bridge: Gretsch flat-bottom base on the plateau, thumbwheels on
    74 mm posts, and a compensated wooden archtop top with a radiused saddle crown."""
    h = SADDLE[2] - Z_PLAT
    y0 = P['bridge_y']; L, W, H = BASE['L'], BASE['W'], BASE['H']
    # base: full width at the centre, tapering to narrow wing tips
    base = hull([(sx * L / 2, y0 + sy * 3.0, Z_PLAT + dz) for sx in (-1, 1) for sy in (-1, 1) for dz in (0.0, 5.0)] +
                [(sx * 30.0, y0 + sy * W / 2, Z_PLAT + dz) for sx in (-1, 1) for sy in (-1, 1) for dz in (0.0, H - 1.5)] +
                [(sx * 12.0, y0 + sy * (W / 2 - 2), Z_PLAT + H) for sx in (-1, 1) for sy in (-1, 1)])
    feet_relief = hull([(sx * 22.0, y0 + sy * 10, Z_PLAT - 1) for sx in (-1, 1) for sy in (-1, 1)] +
                       [(sx * 18.0, y0 + sy * 10, Z_PLAT + 2.5) for sx in (-1, 1) for sy in (-1, 1)])
    base = base - feet_relief
    zt = Z_PLAT + H + 3.0                       # thumbwheels set the gap
    wheels = union([m3.Manifold.cylinder(3.0, 10.0, 10.0, 32).translate([sx * BASE['post'], y0, Z_PLAT + H]) for sx in (-1, 1)])
    posts = union([m3.Manifold.cylinder(h - H, 2.0, 2.0, 16).translate([sx * BASE['post'], y0, Z_PLAT + H]) for sx in (-1, 1)])
    comp = lambda x: -(2.5 - x / 35.0 * 1.0)   # treble ~1.5 mm, bass ~3.5 mm behind the scale line
    top = hull([(sx * 58.0, y0 + sy * 6.0 + comp(sx * 35) * 0.4, zt) for sx in (-1, 1) for sy in (-1, 1)] +
               [(sx * 52.0, y0 + sy * 4.5 + comp(sx * 35) * 0.7, zt + (SADDLE[2] - zt) * 0.6) for sx in (-1, 1) for sy in (-1, 1)] +
               [(sx * 38.0, y0 + sy * 1.4 + comp(sx * 35), SADDLE[2] + 3) for sx in (-1, 1) for sy in (-1, 1)])
    R_ = P['fb_radius']
    top = top ^ m3.Manifold.cylinder(200, R_, R_, 720).rotate([90, 0, 0]).translate([0, y0 + 100, SADDLE[2] - R_])
    return union([base, wheels, posts, top]), h

TAIL_ANCHOR_Y = 50.0
def build_tailpiece():
    # hinge on the tail face, strap over the edge, trapeze plate to the string anchor bar
    zt = top_height_at(0, TAIL_ANCHOR_Y)
    anchor = np.array([0, TAIL_ANCHOR_Y, zt + 4.5])
    hinge = hull([(sx * 25.5, y, z) for sx in (-1, 1) for y in (-3.0, -1.0) for z in (Z_BACK_IN + 10, Z_BACK_IN + 46)])   # standard 51 mm plate
    strap = hull([(sx * 11, y, z) for sx in (-1, 1) for y in (-7.0, -1.5) for z in (Z_BACK_IN + 30, Z_EDGE + 7.5)])
    strap = strap + hull([(sx * 11, y, z) for sx in (-1, 1) for y in (-7.0, 7.0) for z in (Z_EDGE + 4.8, Z_EDGE + 7.5)])
    p0 = np.array([0, 6.0, Z_EDGE + 6.0])
    arms = []
    for sx in (-1, 1):
        a = [(sx * 9, p0[1], p0[2]), (sx * 14, p0[1], p0[2]), (sx * 30, anchor[1], anchor[2]), (sx * 36, anchor[1], anchor[2])]
        arms.append(hull([(x, y, z + dz) for x, y, z in a for dz in (-1.5, 1.5)]))
    bar = hull([(sx * 37, anchor[1] + dy, anchor[2] + dz) for sx in (-1, 1) for dy in (-2, 10) for dz in (-2.5, 2.5)])
    # diamond badge in the middle of the trapeze
    mid = (p0 + anchor) / 2
    dia = hull([(0, mid[1] - 18, mid[2]), (0, mid[1] + 18, mid[2] + 3), (9, mid[1], mid[2] + 1.5), (-9, mid[1], mid[2] + 1.5)]
               + [(0, mid[1] - 18, mid[2] + 2.4), (0, mid[1] + 18, mid[2] + 5.4)])
    return union([hinge, strap] + arms + [bar, dia]), anchor

def rrect(w, l, r, n=6):
    pts = []
    for cx, cy, a0 in [(w / 2 - r, l / 2 - r, 0), (-w / 2 + r, l / 2 - r, 90), (-w / 2 + r, -l / 2 + r, 180), (w / 2 - r, -l / 2 + r, 270)]:
        for t in np.linspace(a0, a0 + 90, n):
            pts.append((cx + r * math.cos(math.radians(t)), cy + r * math.sin(math.radians(t))))
    return np.array(pts)

def string_slope():
    return (NUT_PT[2] - SADDLE[2]) / (NUT_PT[1] - SADDLE[1])

EMG91 = dict(L=86.3, W=38.1, H=10.9, R=5.0)     # EMG 91 floating humbucker (EMG install sheet; 0.43 in max height)
def build_pickup():
    """EMG 91 active floating humbucker, hung from the fingerboard extension on a bracket.
    Black housing, no exposed poles. Sits parallel to the strings, 4 mm under them."""
    pk = P['pickup']
    L, W, H, R = EMG91['L'], EMG91['W'], EMG91['H'], EMG91['R']
    fb_end_y = n2b(U_FB_END, 0, 0)[1]
    yc = fb_end_y - W / 2 - 2.0
    ztop = string_z_at(yc) - pk['string_gap']
    zbot = ztop - H
    inner = rrect(L, W, R)
    lo = [(x, y, 0.0) for x, y in inner]; mid = [(x, y, H - 1.2) for x, y in inner]
    hi = [(x * (L - 2.4) / L, y * (W - 2.4) / W, H) for x, y in inner]
    body = hull(lo + mid + hi)
    ang = math.degrees(math.atan(string_slope()))
    place = lambda man: man.rotate([ang, 0, 0]).translate([0, yc, zbot])
    xs = np.linspace(-L / 2, L / 2, 13); ys = np.linspace(yc - W / 2, yc + W / 2, 9)
    sl = string_slope()
    clr = min(zbot + sl * (y - yc) - top_height_at(x, y) for x in xs for y in ys)
    topmax = max(top_height_at(x, y) for x in xs for y in ys)
    info = dict(model='EMG 91, fingerboard-bracket mount', y=float(yc), z_top=float(ztop), z_bottom=float(zbot), height=H,
                max_height_for_3mm_clearance=float(ztop - 3.0 - topmax), clearance_to_top=float(clr),
                string_to_top=float(string_z_at(yc) - topmax), length=L, width=W)
    return place(body), None, None, info

def STR_SPREAD_AT(y):
    t = (y - SADDLE[1]) / (NUT_PT[1] - SADDLE[1]); return 52.4 + t * (35.0 - 52.4)

# ----------------------------------------------------------------------------- CONTROLS (rim access plate)
def ctrl_frame():
    a, b = SHARP[names.index('tTail')], SHARP[names.index('tLower')]
    e = (b - a) / np.linalg.norm(b - a); n = np.array([e[1], -e[0]])     # n points out of the body (CCW outline)
    return (a + b) / 2, n, e
CTRL_C, CTRL_N, CTRL_E = ctrl_frame()
CTRL_PLATE = (36.0, 36.0, 2.0)           # along the rim, tall, thick (black fibre or aluminium): stacked pot only
CTRL_PLATE_S = 36.0                      # plate centre along the face from CTRL_C (toward the waist)
CTRL_WIN_D = 26.0                        # round window through the wall behind it (the stacked pot body is 24 mm)
CTRL_PAD_L = 116.0
CTRL_SCREWS = [(CTRL_PLATE_S + se, sz) for se in (-13.5, 13.5) for sz in (-13.5, 13.5)]
CTRL_HOLES = [('concentric', CTRL_PLATE_S, 9.5)]   # stacked volume/tone pot. Position along the face, hole dia
# Gotoh BB-04 single 9 V battery box (StewMac #3577 is the same format), toward the tail on the same face
BB = dict(s=-21.0, flange=(70.0, 29.0, 1.6), screws=(61.5, 20.7), body=(56.0, 24.3, 33.0), cut=(57.0, 25.0))
def ctrl_box(L, n0, n1):
    """plan rectangle centred on the control face: +-L/2 along the rim, n0..n1 along the outward normal"""
    c, n, e = CTRL_C, CTRL_N, CTRL_E
    return ccw(np.array([c + e * s_ + n * t for s_, t in [(-L / 2, n0), (L / 2, n0), (L / 2, n1), (-L / 2, n1)]]))

def jack_box(L, n0, n1):
    c, n, e = JACK_C, JACK_N, JACK_E
    return ccw(np.array([c + e * s_ + n * t for s_, t in [(-L / 2, n0), (L / 2, n0), (L / 2, n1), (-L / 2, n1)]]))

def _ctrl_local(man):
    """map a part built in plate-local coords (x along the rim, y out of the body, z up) onto the rim face"""
    ang = math.degrees(math.atan2(CTRL_E[1], CTRL_E[0]))
    return man.mirror([0, 1, 0]).rotate([0, 0, ang]).translate([*CTRL_C, Z_SPLIT])   # local +y -> outward normal

def _rrect_xz(L, H, r, y0, y1, x0=0.0, z0=0.0):
    """rounded rectangle in the face plane (x along, z up), extruded along local y"""
    pts = [(x0 + sx * (L / 2 - r) + r * math.cos(t), y, z0 + sz * (H / 2 - r) + r * math.sin(t))
           for sx in (-1, 1) for sz in (-1, 1) for t in np.linspace(0, 2 * math.pi, 17)[:-1] for y in (y0, y1)]
    return hull(pts)

def _ycyl(r, y0, y1, x, z, n=32):
    return m3.Manifold.cylinder(y1 - y0, r, r, n).rotate([-90, 0, 0]).translate([x, y0, z])

def ctrl_cutters(depth=40.0):
    """openings through the rim wall on the control face, in body coordinates: the round window behind the plate
    and the BB-04's body opening. Both are template routes from outside."""
    win = _ycyl(CTRL_WIN_D / 2, -depth, 6.0, CTRL_PLATE_S, 0.0, 64)
    L, H = BB['cut']
    bb = _rrect_xz(L, H, 3.2, -depth, 6.0, BB['s'])       # 1/4 in router bit corners
    return _ctrl_local(win + bb)

def build_control_plate():
    L, H, T = CTRL_PLATE; x0 = CTRL_PLATE_S
    plate = _rrect_xz(L, H, 4.0, 0.0, T - 0.6, x0) + _rrect_xz(L - 1.2, H - 1.2, 3.4, T - 0.6, T, x0)
    for _, x, d in CTRL_HOLES:
        plate = plate - _ycyl(d / 2, -5, 10, x, 0)
    for se, sz in CTRL_SCREWS:
        plate = plate - _ycyl(1.4, -5, 10, se, sz, 16)
    return _ctrl_local(plate)

def build_battery_box():
    """Gotoh BB-04: flange screwed to the rim face, hinged lid in the flange, box body through the wall, 9 V inside"""
    s0 = BB['s']; FL, FH, FT = BB['flange']; bl, bw, bd = BB['body']
    flange = _rrect_xz(FL, FH, 3.0, 0.0, FT, s0) - _rrect_xz(FL - 6.0, FH - 5.0, 2.0, -1, FT + 1, s0)
    body = _rrect_xz(bl, bw, 3.5, -bd, 0.01, s0) - _rrect_xz(bl - 2.4, bw - 2.4, 2.3, -bd + 1.2, 1.0, s0)
    sx_, sz_ = BB['screws']
    heads = union([_ycyl(2.3, FT, FT + 0.9, s0 + a * sx_ / 2, b * sz_ / 2, 16) for a in (-1, 1) for b in (-1, 1)])
    box = flange + body + heads
    lid = _rrect_xz(FL - 6.6, FH - 5.6, 1.7, 0.0, FT + 0.4, s0)
    lid = lid - _rrect_xz(10.0, 4.0, 1.9, FT - 0.4, FT + 1.0, s0 + (FL - 6.6) / 2 - 8.0, 0.0)     # thumb notch for the catch
    lid = lid + _ycyl(1.2, FT + 0.4, FT + 0.7, s0 - (FL - 6.6) / 2 + 5.0, 0.0, 12)               # hinge-end rib
    Lb, Tb, Hb = BATTERY
    bat = m3.Manifold.cube([Lb, Hb, Tb], True).translate([s0, -(Hb / 2 + 2.0), 0.0])           # lies in the box, terminals at one end
    return _ctrl_local(lid), _ctrl_local(box), _ctrl_local(bat)

def build_control_hardware():
    """stacked knob, screw heads, and the pot body behind the plate (for clearance checks)"""
    T = CTRL_PLATE[2]; out = []
    ring = lambda r0, r1, y0, y1: hull([(r0 * math.cos(t), y0, r0 * math.sin(t)) for t in np.linspace(0, 2 * math.pi, 7)[:-1]] +
                                       [(r1 * math.cos(t), y1, r1 * math.sin(t)) for t in np.linspace(0, 2 * math.pi, 7)[:-1]])
    for name, x, d in CTRL_HOLES:      # stacked (concentric) pot: outer knob for volume, inner knob for tone
        out.append(ring(11.0, 10.0, T + 1.0, T + 9.0).translate([x, 0, 0]))
        out.append(ring(7.0, 6.0, T + 9.0, T + 17.0).translate([x, 0, 0]))
        out.append(_ycyl(12.0, -14.0, 0.0, x, 0.0))      # stacked pot body
    for se, sz in CTRL_SCREWS:
        out.append(_ycyl(2.6, T, T + 1.0, se, sz, 16))
    return _ctrl_local(union(out))

# ----------------------------------------------------------------------------- OUTPUT JACK (treble tail facet)
def jack_frame():
    a, b = SHARP[names.index('tailR')], SHARP[names.index('tTail')]
    e = (b - a) / np.linalg.norm(b - a); n = np.array([e[1], -e[0]])
    return (a + b) / 2, n, e
JACK_C, JACK_N, JACK_E = jack_frame()
JACK_HOLE = 12.2          # Switchcraft 152B stereo long-thread jack: 12 mm bushing, panels to 1.25 in
JACK_PAD_L = 34.0

def _jack_local(man):
    ang = math.degrees(math.atan2(JACK_E[1], JACK_E[0]))
    return man.mirror([0, 1, 0]).rotate([0, 0, ang]).translate([*JACK_C, Z_SPLIT])

def jack_cutter():
    return _jack_local(_ycyl(JACK_HOLE / 2, -40.0, 6.0, 0.0, 0.0))

def build_jack():
    """152B in the rim wall: flanged nut and washer outside, threaded barrel through the wall"""
    hexa = hull([(7.5 * math.cos(t), y, 7.5 * math.sin(t)) for t in np.linspace(0, 2 * math.pi, 7)[:-1] for y in (1.0, 4.5)])
    wash = _ycyl(9.0, 0.0, 1.0, 0, 0)
    barrel = _ycyl(6.0, -24.0, 5.5, 0, 0)
    sock = _ycyl(3.3, -10.0, 8.0, 0, 0)
    return _jack_local((wash + hexa + barrel) - sock)

BATTERY = (48.5, 17.5, 26.5)      # 9 V: long, thick, tall
def fb_bottom_z(y):
    u = U_FB_END + (n2b(U_FB_END, 0, -P['fb_t'])[1] - y) * -1 / math.cos(TH)
    return n2b(u, 0, -P['fb_t'])[2]

def build_pickup_bracket(pk_info):
    """1.5 mm strap from the underside of the fingerboard extension down to the pickup's mounting tab"""
    zt = pk_info['z_bottom']; y0 = pk_info['y']; sl = string_slope()
    x0, x1 = 22.0, 40.0; yb = 372.0
    zfb = fb_bottom_z(yb)
    top = hull([(x, y, fb_bottom_z(y) - dz) for x in (x0 - 5, x0 + 5) for y in (yb - 8, yb + 8) for dz in (0.0, 1.5)])
    drop = hull([(x, y, z) for x in (x0 - 5, x0 + 5) for y in (yb - 8, yb - 6.5) for z in (zt - 1.5 + sl * (yb - 8 - y0), zfb)])
    run = hull([(x, y, zt - 1.5 + sl * (y - y0) + dz) for x in (x0 - 5, x0 + 5) for y in (y0, yb - 6.5) for dz in (0.0, 1.5)])
    arm = hull([(x, y, zt - 1.5 + sl * (y - y0) + dz) for x in (x0 - 5, x1) for y in (y0 - 5, y0 + 5) for dz in (0.0, 1.5)])
    return union([top, drop, run, arm])

WIRE_HOLE = (12.0, 384.0, 6.0)          # pickup lead into the body, hidden under the fingerboard extension

# ---- alternative pickup: Curtis Novak DeArmond-style gold foil in the "mustache" casing, on a bent bracket.
# PLACEHOLDER sizes from typical DeArmond mustache pickups -- measure the one you buy and update these.
GOLDFOIL = dict(body_L=78.0, body_W=34.0, body_H=9.0,      # centre housing: across strings, along strings, tall
                span=124.0, ear_W=14.0, ear_t=1.6,         # mustache casing tip to tip, ear width, flange thickness
                screw_x=56.0, screw_d=3.2,                 # mounting screw positions (+-) and hole size
                spring_gap=3.0)                            # spring space between the casing ears and the shelf, set by the screws
GF_BRACKET = dict(t=1.5, w=20.0, tab_L=18.0, shelf_W=14.0, screw_d=3.0, fb_screws=(-6.0, 6.0))   # brass or aluminium

def _goldfoil_place():
    g = GOLDFOIL
    fb_end_y = n2b(U_FB_END, 0, 0)[1]
    yc = fb_end_y - g['body_W'] / 2 - 2.0
    ztop = string_z_at(yc) - P['pickup']['string_gap']
    return yc, ztop

def mustache_outline(n=40):
    """plan outline of the mustache casing: full width over the body, sweeping to rounded ears"""
    g = GOLDFOIL; hw = g['span'] / 2; bw = g['body_W'] / 2 + 3.0; ew = g['ear_W'] / 2
    xs = np.linspace(-hw + ew, hw - ew, n)
    def half(x):
        t = min(1.0, max(0.0, (abs(x) - g['body_L'] / 2 + 6) / (hw - ew - g['body_L'] / 2 + 6)))
        return bw + (ew - bw) * (3 * t * t - 2 * t * t * t)
    top = [(x, half(x)) for x in xs]; bot = [(x, -half(x)) for x in xs[::-1]]
    capR = [(hw - ew + ew * math.cos(a), ew * math.sin(a)) for a in np.linspace(-math.pi / 2, math.pi / 2, 9)]
    capL = [(-hw + ew - ew * math.cos(a), -ew * math.sin(a)) for a in np.linspace(-math.pi / 2, math.pi / 2, 9)]
    return np.array(capR + top[::-1][1:-1] + capL + bot[1:-1]) if False else ccw(np.array(top + capL[::-1] + bot + capR[::-1]))

def build_goldfoil():
    """gold foil pickup in a mustache casing; ears screw down to the bracket shelf over springs"""
    g = GOLDFOIL; yc, ztop = _goldfoil_place(); sl = string_slope()
    zb = ztop - g['body_H']
    body = hull([(x, y, z) for x, y in rrect(g['body_L'], g['body_W'], 4.0) for z in (0.0, g['body_H'] - 1.0)] +
                [(x * (g['body_L'] - 3) / g['body_L'], y * (g['body_W'] - 3) / g['body_W'], g['body_H']) for x, y in rrect(g['body_L'], g['body_W'], 4.0)])
    flange = prism(mustache_outline(), 0.0, g['ear_t'])
    for sx in (-1, 1):
        flange = flange - m3.Manifold.cylinder(10, g['screw_d'] / 2, g['screw_d'] / 2, 20).translate([sx * g['screw_x'], 0, -5])
    foil = m3.Manifold.cube([g['body_L'] - 14, g['body_W'] - 12, 0.4], True).translate([0, 0, g['body_H'] + 0.2])
    ang = math.degrees(math.atan(sl))
    place = lambda man: man.rotate([ang, 0, 0]).translate([0, yc, zb])
    info = dict(model='Novak DeArmond-style gold foil, mustache casing (placeholder size)', y=float(yc), z_top=float(ztop), z_bottom=float(zb))
    return place(body + flange), place(foil), info

def build_goldfoil_bracket(info):
    """one bent strap: tab screwed under the fingerboard extension, a drop just past the fingerboard end
    (under the board, clear of the pickup), and a shelf running across under the pickup's ears"""
    g, b = GOLDFOIL, GF_BRACKET; t = b['t']; sl = string_slope()
    yc = info['y']; zb = info['z_bottom']
    fb_end_y = n2b(U_FB_END, 0, 0)[1]
    y_drop = fb_end_y + 1.0                                  # just under the board end, past the casing
    z_shelf_top = zb - g['spring_gap']                       # shelf top, under the ears
    zs = lambda y: z_shelf_top + sl * (y - yc)
    w2 = b['w'] / 2
    tab = hull([(x, y, fb_bottom_z(y) - dz) for x in (-w2, w2) for y in (y_drop, y_drop + b['tab_L']) for dz in (0.0, t)])
    drop = hull([(x, y, z) for x in (-w2, w2) for y in (y_drop, y_drop + t) for z in (zs(y_drop) - t, fb_bottom_z(y_drop))])
    arm = hull([(x, y, zs(y) - dz) for x in (-w2, w2) for y in (yc - b['shelf_W'] / 2, y_drop + t) for dz in (0.0, t)])
    span = g['span'] / 2 + 2.0
    shelf = hull([(x, y, zs(y) - dz) for x in (-span, span) for y in (yc - b['shelf_W'] / 2, yc + b['shelf_W'] / 2) for dz in (0.0, t)])
    br = union([tab, drop, arm, shelf])
    for sx in (-1, 1):
        br = br - m3.Manifold.cylinder(60, b['screw_d'] / 2, b['screw_d'] / 2, 16).translate([sx * g['screw_x'], yc, zs(yc) - 30])
    for x in b['fb_screws']:
        y = y_drop + b['tab_L'] / 2
        br = br - m3.Manifold.cylinder(30, 1.5, 1.5, 16).translate([x, y, fb_bottom_z(y) - 15])
    # clearance from the shelf's underside to the top
    xs = np.linspace(-span, span, 31); ys = np.linspace(yc - b['shelf_W'] / 2, y_drop + t, 7)
    clr = min(zs(y) - t - top_height_at(x, y) for x in xs for y in ys if abs(x) <= w2 or abs(y - yc) <= b['shelf_W'] / 2)
    info['bracket_clearance_to_top'] = float(clr)
    info['flat_pattern'] = dict(tab=b['tab_L'], drop=float(fb_bottom_z(y_drop) - zs(y_drop)), arm=float(y_drop - (yc - b['shelf_W'] / 2)),
                                shelf_span=2 * span, width=b['w'], shelf_W=b['shelf_W'], t=t)
    return br

MORTISE = None

def build_strings(anchor):
    xs_b = np.linspace(-26.2, 26.2, 6); xs_n = np.linspace(-17.5, 17.5, 6)
    gauges = [0.56, 0.45, 0.36, 0.30, 0.22, 0.165]   # radius-ish (13s flatwound, scaled)
    out, lines = [], []
    for i in range(6):
        a = np.array([xs_b[i] * 0.86, anchor[1] + 4, anchor[2] + 2.6])
        s = np.array([xs_b[i], SADDLE[1], SADDLE[2]])
        nut = np.array([xs_n[i], NUT_PT[1], NUT_PT[2]])
        side = -1 if i < 3 else 1
        k = i if i < 3 else 5 - i
        post = h2b(side * TUNER_POS[k][0], TUNER_POS[k][1], 1.6 + 10)
        pts = [a, s, nut, post]
        lines.append([p.tolist() for p in pts])
        for p, q in zip(pts[:-1], pts[1:]):
            out.append(hull([p + d for d in sphere_pts(gauges[i])] + [q + d for d in sphere_pts(gauges[i])]))
    return union(out), lines

def sphere_pts(r):
    return [np.array([r * math.cos(a), r * math.sin(a) * math.cos(b), r * math.sin(a) * math.sin(b)])
            for a in np.linspace(0, math.pi, 4) for b in np.linspace(0, 2 * math.pi, 6)[:-1]]

def build_strap_buttons():
    b1 = m3.Manifold.cylinder(12, 5.0, 5.0, 24).rotate([90, 0, 0]).translate([0, -1, Z_BACK_IN + 8])
    hn = SHARP[names.index('bHorn')]
    return b1

# ----------------------------------------------------------------------------- export
def to_trimesh(man):
    import trimesh
    mh = man.simplify(0.02).to_mesh(); v = np.array(mh.vert_properties)[:, :3]; f = np.array(mh.tri_verts)
    t = trimesh.Trimesh(v, f, process=True)
    if t.body_count > 1:            # drop zero-volume slivers that coplanar CSG faces leave behind
        parts = [p for p in t.split(only_watertight=False) if abs(p.volume) > 5.0]
        if parts: t = trimesh.util.concatenate(parts) if len(parts) > 1 else parts[0]
    if not t.is_watertight:          # needle triangles from the crease blends: weld and collapse them
        import meshfix; t = meshfix.clean(t)
    parts = t.split(only_watertight=False)
    if len(parts) > 1:               # stray specks (a few faces each) that show as dots in the viewer
        big = [p for p in parts if len(p.faces) > 20 or abs(p.volume) > 5.0]
        if big: t = trimesh.util.concatenate(big) if len(big) > 1 else big[0]
    return t

def main():
    global MORTISE
    import trimesh
    MORTISE = mortise()
    parts = {}
    print('top...');  parts['top'] = build_top()
    band = m3.Manifold.extrude(OUTL_CS - OUTL_CS.offset(-BINDING_W, m3.JoinType.Round, circular_segments=96), 2 * BIG).translate([0, 0, -BIG])
    parts['top_binding'] = parts['top'] ^ band          # left undyed: faux binding
    parts['top'] = parts['top'] - band
    print('back...'); parts['back'] = build_back()
    print('rim...');  parts['rim'] = build_rim()
    print('neck...'); parts['neck'] = build_neck() - tuner_holes()
    parts['fingerboard'] = build_fingerboard()
    parts['faceplate'] = build_faceplate()
    parts['truss_cover'] = build_truss_cover()
    parts['frets'] = build_frets()
    parts['nut'] = build_nut()
    parts['side_dots'] = build_side_dots()
    parts['tuners'] = build_tuners()
    br, bridge_h = build_bridge(); parts['bridge'] = br
    tp, anchor = build_tailpiece(); parts['tailpiece'] = tp
    pk_cover, pk_plate, pk_poles, pk_info = build_pickup()
    parts['pickup'] = pk_cover
    parts['pickup_bracket'] = build_pickup_bracket(pk_info)
    gf_body, gf_foil, gf_info = build_goldfoil()
    parts['goldfoil'] = gf_body; parts['goldfoil_foil'] = gf_foil
    parts['goldfoil_bracket'] = build_goldfoil_bracket(gf_info)
    parts['control_plate'] = build_control_plate()
    parts['controls'] = build_control_hardware()
    parts['battery_door'], parts['battery_sleeve'], parts['battery'] = build_battery_box()
    parts['jack'] = build_jack()
    parts['strap_buttons'] = build_strap_buttons()
    st, string_lines = build_strings(anchor); parts['strings'] = st

    for k, v in parts.items():
        if v.is_empty(): print('EMPTY', k)

    MAT = dict(goldfoil='nickel', goldfoil_foil='gold', goldfoil_bracket='brass', top_binding='spruce', truss_cover='nickel', badge='nickel', badge_qr='black', top='spruce_black', back='maple', rim='mahogany', neck='mahogany', fingerboard='richlite', faceplate='maple',
               frets='nickel', nut='bone', side_dots='bone', tuners='nickel', bridge='ebony', tailpiece='nickel',
               pickup='black', pickup_plate='brass', pickup_poles='nickel', pickup_bracket='nickel', control_plate='black', controls='black', battery='steel', battery_door='black', battery_sleeve='black',
               jack='nickel', strap_buttons='nickel', strings='steel')
    GROUP = dict(goldfoil='Pickup', goldfoil_foil='Pickup', goldfoil_bracket='Pickup', top_binding='Top', truss_cover='Neck', badge='Neck', badge_qr='Neck', top='Top', back='Back', rim='Rim', neck='Neck', fingerboard='Neck', faceplate='Neck', frets='Neck',
                 nut='Neck', side_dots='Neck', tuners='Hardware', bridge='Hardware', tailpiece='Hardware',
                 pickup='Pickup', pickup_plate='Pickup', pickup_poles='Pickup', pickup_bracket='Pickup', control_plate='Hardware', controls='Hardware', battery='Hardware', battery_door='Hardware', battery_sleeve='Hardware',
                 jack='Hardware', strap_buttons='Hardware', strings='Strings')
    EXPLODE = dict(Top=[0, 0, 70], Back=[0, 0, -70], Rim=[0, 0, 0], Neck=[0, 120, 45], Hardware=[0, 0, 120],
                   Pickup=[0, 0, 150], Strings=[0, 0, 200])
    DENS = dict(gold=8400, brass=8500, spruce=462, spruce_black=462,   # top: Port Orford cedar (SG 0.43 vs Sitka 0.40)
                 spruce_back=430, maple=550, mahogany=590, ebony=1150, richlite=1400, bone=1900, nickel=8400, steel=7850, black=1400)

    scene = trimesh.Scene()
    vols = {}
    for k, man in parts.items():
        tm = to_trimesh(man)
        vols[k] = man.volume()
        tm.metadata['name'] = k
        scene.add_geometry(tm, node_name=k, geom_name=k,
                           metadata=dict(material=MAT[k], group=GROUP[k], explode=EXPLODE[GROUP[k]]))
        if k == 'top':
            to_trimesh(build_top()).export(os.path.join(OUT, 'top.stl'))   # one solid part (dye mask is cosmetic)
            _bc = globals()['build_collars']; globals()['build_collars'] = lambda bi: m3.Manifold()
            to_trimesh(build_top()).export(os.path.join(OUT, 'top_fea.stl'))   # for fea.py: lips left off (conservative)
            globals()['build_collars'] = _bc
        elif k in ('back', 'rim', 'neck', 'fingerboard', 'bridge'):
            tm.export(os.path.join(OUT, f'{k}.stl'))
    scene.export(os.path.join(OUT, 'juris_c75a.glb'))
    # extras per node are not always carried by trimesh; also write a sidecar json
    meta = {k: dict(material=MAT[k], group=GROUP[k], explode=EXPLODE[GROUP[k]]) for k in parts}

    # masses (wood exact from volume; hardware nominal where solids are simplified)
    mass = {k: vols[k] * 1e-9 * DENS[MAT[k]] * 1000 for k in parts}
    nominal = dict(goldfoil=110, goldfoil_foil=0, goldfoil_bracket=12, tuners=6 * 33, tailpiece=95, pickup=85, battery=45, frets=22, strings=10, strap_buttons=8, controls=50, control_plate=6, battery_door=4, battery_sleeve=14, jack=12, pickup_bracket=10, pickup_plate=0, pickup_poles=0, badge=8, badge_qr=0, truss_cover=3)
    for k, g in nominal.items(): mass[k] = g
    total_g = sum(v for k, v in mass.items() if not k.startswith('goldfoil'))     # gold foil is the alternative pickup

    # air volume & Helmholtz estimate
    outer_body = (BODY_PRISM ^ terrain_below(TOP_TRIS, TOP_HEX, Z_EDGE, Z_PLAT, Z_PLAT, Z_EDGE)
                  ^ terrain_above(BACK_TRIS, BACK_HEX, 0.0, -P['back_arch'], -P['back_arch'], 0.0))
    air = outer_body.volume() - vols['top'] - vols['top_binding'] - vols['back'] - vols['rim'] - (parts['neck'] ^ outer_body).volume()
    A = sum(abs(signed_area(s)) for s in SLOTS)
    t_hole = top_z_outer(0) - top_z_inner(0)
    # slots sit on a fall ~ halfway: thickness ~ 3.9
    t_hole = P['top_fall_t_at_plateau'] + 0.1 + P['collar_h']      # plate at the slots + collar lip
    a_eq = math.sqrt(A / 2 / math.pi)
    L_eff = t_hole + 1.7 * a_eq
    c = 343e3
    f_h = c / (2 * math.pi) * math.sqrt(A / (air * L_eff))

    fb_end = n2b(U_FB_END, 0, 0)
    fb_end_under = n2b(U_FB_END, 0, -P['fb_t'])
    fret12_u = FRET_U[11]
    report = dict(
        units='mm',
        body_length=float(SHARP[:, 1].max()), lower_bout=float(OUTL[:, 0].max() - OUTL[:, 0].min()),
        body_length_in=float(OUTL[:, 1].max() / IN), lower_bout_in=float((OUTL[:, 0].max() - OUTL[:, 0].min()) / IN),
        depth_at_edge=Z_EDGE, depth_max=Z_PLAT + P['back_arch'],
        neck_joint_fret=float(-12 * math.log2(1 - U_NUT / P['scale'])),
        saddle_z=float(SADDLE[2]), plateau_z=Z_PLAT, bridge_height=float(bridge_h), bridge_height_in=float(bridge_h / IN),
        string_break_angle_deg=float(math.degrees(math.atan2(SADDLE[2] - (anchor[2] + 2.6), SADDLE[1] - anchor[1] - 4))),
        nut_vs_top_edge=float(NUT_PT[2] - Z_EDGE),
        fb_end_y=float(fb_end[1]), fb_extension_clearance=float(fb_end_under[2] - top_height_at(0, fb_end_under[1])),
        pickup=pk_info, goldfoil=gf_info,
        soundhole_area_each=float(A / 2), air_volume_l=float(air * 1e-6), helmholtz_hz_est=float(f_h),
        masses_g={k: round(v, 1) for k, v in mass.items()}, total_mass_g=round(total_g), total_mass_lb=round(total_g / 453.6, 2),
        manifold_ok={k: bool(v.status() == m3.Error.NoError and not v.is_empty()) for k, v in parts.items()},
        params=P,
    )
    json.dump(report, open(os.path.join(OUT, 'report.json'), 'w'), indent=1, default=float)
    json.dump(dict(meta=meta, strings=string_lines), open(os.path.join(OUT, 'meta.json'), 'w'))
    write_dxf()
    print(json.dumps({k: v for k, v in report.items() if k not in ('params', 'masses_g', 'manifold_ok')}, indent=1, default=float))
    print('mass', report['masses_g'])
    print('manifold', report['manifold_ok'])

def write_dxf():
    import ezdxf
    doc = ezdxf.new('R2010'); doc.units = ezdxf.units.MM; ms = doc.modelspace()
    for lay, col in [('OUTLINE', 7), ('RIM_INNER', 1), ('LINING_INNER', 2), ('BLOCKS', 3), ('TOP_PLATEAU', 4),
                     ('TOP_FOLDS', 5), ('TOP_BEAM', 6), ('SOUNDHOLES', 1), ('SOUNDHOLE_LIPS', 6), ('BACK_PLATEAU', 4), ('BACK_FOLDS', 5), ('TOP_BARS', 6), ('TOP_DYE_MASK', 1), ('CENTERLINE', 8), ('BRIDGE_LINE', 8)]:
        doc.layers.add(lay, color=col)
    def poly(pts, layer, closed=True): ms.add_lwpolyline([tuple(map(float, p)) for p in pts], close=closed, dxfattribs={'layer': layer})
    poly(OUTL, 'OUTLINE')
    for cs, lay in [(OUTL_CS.offset(-P['rim_wall'], m3.JoinType.Round), 'RIM_INNER'),
                    (OUTL_CS.offset(-(P['rim_wall'] + P['lining_w']), m3.JoinType.Round), 'LINING_INNER')]:
        for c in cs.to_polygons(): poly(c, lay)
    poly(TOP_HEX, 'TOP_PLATEAU'); poly(offset_polygon(ccw(TOP_HEX), P['beam_w']), 'TOP_BEAM')
    poly(BACK_HEX, 'BACK_FOLDS'); poly(BACK_FOOT, 'BACK_FOLDS')
    poly(ARCH, 'TOP_FOLDS')
    for tris, hexa, lay, outer in [(TOP_TRIS, TOP_HEX, 'TOP_FOLDS', ARCH), (BACK_TRIS, BACK_HEX, 'BACK_FOLDS', BACK_FOOT)]:
        seen = set()
        for t in tris:
            for a, b in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])]:
                if a[0] == 'o' and b[0] == 'o': continue
                if a[0] == 'p' and b[0] == 'p': continue
                key = tuple(sorted([a, b]))
                if key in seen: continue
                seen.add(key)
                pa = outer[a[1]] if a[0] == 'o' else hexa[a[1]]; pb = outer[b[1]] if b[0] == 'o' else hexa[b[1]]
                ms.add_line(tuple(pa), tuple(pb), dxfattribs={'layer': lay})
    for s in SLOTS: poly(s, 'SOUNDHOLES')
    for c in collar_cs(): poly(np.array(c.exterior.coords)[:-1], 'SOUNDHOLE_LIPS')
    for c in OUTL_CS.offset(-BINDING_W, m3.JoinType.Round, circular_segments=96).to_polygons(): poly(c, 'TOP_DYE_MASK')
    for ra, rb in bar_segments() + [cross_bar_segment()]: ms.add_line(tuple(ra), tuple(rb), dxfattribs={'layer': 'TOP_BARS'})
    ms.add_line((0, -10), (0, 480), dxfattribs={'layer': 'CENTERLINE'})
    ms.add_line((-60, P['bridge_y']), (60, P['bridge_y']), dxfattribs={'layer': 'BRIDGE_LINE'})
    doc.saveas(os.path.join(OUT, 'juris_c75a_templates.dxf'))

if __name__ == '__main__':
    main()
