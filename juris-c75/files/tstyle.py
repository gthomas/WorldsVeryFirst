"""Juris C-75T — a T-style solid body on the C-75 outline, with the C-75's bevelled edges.
Bolt-on 25.5" neck (stock T-style neck, Fender pocket), Tele bridge pickup in a vintage-style bridge plate,
neck single coil screwed to the body (no pickguard), control plate with 3-way switch, volume and tone.
Body frame as the C-75A: x+ treble, y+ toward the neck (y=0 tail), z+ up (z=0 back). Units mm.
Run: python tstyle.py  -> out_t/ (GLB, STL, DXF, report.json)"""
import json, math, os
import numpy as np
import manifold3d as m3
import build as B                       # shared CSG helpers (prism, hull, union, loft, ...)
import outline as O
import head_c75 as H

IN = 25.4
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out_t'); os.makedirs(OUT, exist_ok=True)

P = dict(
    T=44.45,                 # body thickness, 1-3/4"
    bevel_drop=11.0,         # how far each bevel falls at the outline, top and back (C-75 side elevation: ~0.45")
    scale=647.7,             # 25.5"
    comp=2.5,                # average saddle compensation
    frets=21, fret_body=16,  # body edge at the 16th fret (Tele)
    pocket_w=56.1, pocket_L=76.2, pocket_d=15.9, pocket_r=3.2,      # 2.2" heel + 0.2 mm, 3", 5/8"
    heel_t=25.4,             # neck thickness at the heel including the fingerboard
    fb_t=6.35, fb_radius=304.8, nut_w=42.0, heel_w=55.9,           # 12" radius, 1.65" nut
    neck_d1=21.0, neck_d12=22.5,
    head_angle=13.0,         # tilt-back headstock with a volute, as on the C-75 (so the neck is custom, not a stock paddle neck)
    action12=2.0, fret_h=1.2,
    string_spacing_bridge=54.0,  # 2-1/8"
    neck_half=35.0,          # body half-width at the neck end: pocket walls ~7 mm
    tail_half=16.0,
)

# ----------------------------------------------------------------------------- OUTLINE
O.TAIL_HALF = P['tail_half']; O.NECK_HALF = P['neck_half']
names, SHARP = O.sharp_outline()
RADII = dict(O.RADII); RADII.update(neckT=0, neckB=0)
OUTL = B.ccw(O.fillet(names, SHARP, RADII))
Y_END = SHARP[names.index('neckT')][1]
T = P['T']

# ----------------------------------------------------------------------------- BEVELS
# Inner bevel lines traced from Juris_C_75_Body_2D.pdf (body mm), as straight runs. The drawing is asymmetric
# (the bass bevel is wider through the lower bout and wraps onto the tail), and that is kept.
BASS_CHAIN = [(-55.4, 15.2), (-140.2, 137.3), (-140.9, 158.8), (-85.1, 274.2), (-84.0, 308.8), (-111.4, 387.9),
              (-82.5, 433.4), (-77.6, 433.6), (-42.0, 421.5)]
TREB_CHAIN = [(97.0, 25.9), (144.4, 119.3), (143.2, 146.2), (93.5, 256.0), (87.8, 286.5), (113.5, 359.4),
              (105.3, 367.7), (66.0, 359.0), (57.6, 358.9), (40.0, 383.0)]
# straight runs between those points (the short ones in between are the drawing's fillets): (start, end) index pairs
BASS_RUNS = [(0, 1), (2, 3), (4, 5), (5, 6), (7, 8)]
TREB_RUNS = [(0, 1), (2, 3), (4, 5), (6, 7), (8, 9)]

def _isect(p, d, q, e):
    A = np.array([d, -e]).T; t = np.linalg.solve(A, np.asarray(q) - np.asarray(p)); return np.asarray(p) + t[0] * np.asarray(d)

def bevel_lines(chain, runs):
    """merge each run's line with its neighbours: corners are where consecutive run lines meet"""
    C = np.array(chain, float)
    L = [(C[a], C[b] - C[a]) for a, b in runs]
    verts = [C[runs[0][0]]]
    for i in range(len(L) - 1):
        verts.append(_isect(L[i][0], L[i][1], L[i + 1][0], L[i + 1][1]))
    verts.append(C[runs[-1][1]])
    return np.array(verts)

def _to_edge(p, d, outl, back=0.3):
    """from p along d to the first crossing of the outline (stops `back` mm short, inside)"""
    from shapely.geometry import LineString, Point
    d = np.asarray(d, float) / np.linalg.norm(d)
    g = LineString([p, p + d * 400]).intersection(outl.exterior)
    pts = [g] if g.geom_type == 'Point' else list(getattr(g, 'geoms', []))
    pts = [np.array(q.coords[0]) for q in pts if q.geom_type == 'Point']
    if not pts: return np.asarray(p, float)
    q = min(pts, key=lambda q: np.linalg.norm(q - p)); return q - d * back

def bevel_chains(outl_pts, sharp, nm):
    """flat bevel facets for the lower body and upper-bout sides (the traced runs up to the upper-bout corners),
    each run carried out to the body edge so it tapers to nothing. The horns and the cutaway are done by
    `edge_band` instead, as on the built guitar."""
    from shapely.geometry import Polygon
    outl = Polygon(outl_pts)
    bv = bevel_lines(BASS_CHAIN, BASS_RUNS[:3]); tv = bevel_lines(TREB_CHAIN, TREB_RUNS[:3])
    bv[0] = _to_edge(bv[1], bv[0] - bv[1], outl); tv[0] = _to_edge(tv[1], tv[0] - tv[1], outl)
    bv[-1] = _to_edge(bv[-2], bv[-1] - bv[-2], outl); tv[-1] = _to_edge(tv[-2], tv[-1] - tv[-2], outl)
    return bv, tv

def bevel_band(outl_pts, sharp, nm, z_top, down, drop=None, step=1.5):
    """The C-75 bevel as one band that follows the outline all the way round. Its width at each point is the
    distance to the traced inner bevel line (wide on the bass side, narrower on the treble, wrapping onto the
    tail as drawn); round the horns, the cutaway and into the neck it follows the edge and narrows, as on the
    built guitar. The bevel falls `drop` at the edge (less where it is narrower than 15 mm, so it fades out)."""
    from shapely.geometry import Polygon, Point, LineString
    D = P['bevel_drop'] if drop is None else drop
    O_ = np.asarray(outl_pts, float); ring = Polygon(O_).exterior; L = ring.length
    sp = lambda name: ring.project(Point(sharp[nm.index(name)]))
    outl = Polygon(O_)
    bch = bevel_lines(BASS_CHAIN, BASS_RUNS); tch = bevel_lines(TREB_CHAIN, TREB_RUNS)[:4]
    bch[0] = _to_edge(bch[1], bch[0] - bch[1], outl); tch[0] = _to_edge(tch[1], tch[0] - tch[1], outl)
    bL, tL = LineString(bch), LineString(tch)
    s_tT, s_tH, s_tC, s_nT = ring.project(Point(tch[0])), sp('tHorn'), sp('tCut'), sp('neckT')
    s_nB, s_bC, s_bH, s_bT = sp('neckB'), sp('bCut'), sp('bHorn'), ring.project(Point(bch[0]))
    pt = lambda sv: np.array(ring.interpolate(sv % L).coords[0])
    w_tH = tL.distance(Point(pt(s_tH))); w_bH = bL.distance(Point(pt(s_bH)))
    def width(sv):
        p = Point(pt(sv))
        if s_tT <= sv <= s_tH: return tL.distance(p)
        if s_tH < sv <= s_nT: return float(np.interp(sv, [s_tH, s_tC, s_nT], [w_tH, 10.0, 3.5]))
        if s_nB <= sv < s_bH: return float(np.interp(sv, [s_nB, s_bC, s_bH], [3.5, 7.0, w_bH]))
        if s_bH <= sv <= s_bT: return bL.distance(p)
        return 0.0
    st = np.arange(0.0, L, step)
    pts = []; _jit = np.random.default_rng(7 if down > 0 else 11)
    for sv in st:
        p = pt(sv); a = pt(sv - 0.8); b = pt(sv + 0.8)
        t = (b - a) / np.linalg.norm(b - a); n_in = np.array([-t[1], t[0]])
        w = width(sv)
        if w < 0.3:
            pts.append(None); continue
        d = min(D, w * D / 15.0)
        A = p + n_in * w + _jit.normal(0, 0.004, 2); Bo = p - n_in * 3.0; zl = z_top - down * (d + 3.0 * d / w)   # tiny jitter: stations that fan to one inner corner must not share an edge
        pts.append([(A[0], A[1], z_top + down * 0.001), (Bo[0], Bo[1], zl), (A[0], A[1], z_top + down * 60), (Bo[0], Bo[1], z_top + down * 60)])
    pieces = [B.hull(pts[i] + pts[(i + 1) % len(pts)]) for i in range(len(pts)) if pts[i] is not None and pts[(i + 1) % len(pts)] is not None]
    return B.union(pieces)

BASS_V, TREB_V = bevel_chains(OUTL, SHARP, names)

def _outward(p, q, body_side_pt):
    d = (q - p) / np.linalg.norm(q - p); n = np.array([-d[1], d[0]])
    return n if np.dot(body_side_pt - p, n) < 0 else -n

def bevel_cutters(verts, z_top, down, outl_pts=None):
    """one plane per run: it contains the inner line at z_top and falls `bevel_drop` by the time it reaches the
    outline. Each plane only acts inside its own sector (bounded by the bisectors at its ends), so concave runs
    don't over-cut their neighbours. `down` = +1 for the top (falls toward the back), -1 for the back."""
    from shapely.geometry import Polygon, LineString
    outl = Polygon(OUTL if outl_pts is None else outl_pts); centre = np.array([0.0, 220.0])
    ns = [_outward(verts[i], verts[i + 1], centre) for i in range(len(verts) - 1)]
    cut = []
    for i in range(len(verts) - 1):
        p, q, n = verts[i], verts[i + 1], ns[i]
        bp = n if i == 0 else (ns[i - 1] + n) / np.linalg.norm(ns[i - 1] + n)
        bq = n if i == len(ns) - 1 else (ns[i + 1] + n) / np.linalg.norm(ns[i + 1] + n)
        # width: distance from the inner line to the outline, measured outward from the run's middle
        ws = []
        for t in (0.25, 0.5, 0.75):
            m_ = p + (q - p) * t
            ray = LineString([m_, m_ + n * 200]).intersection(outl.exterior)
            pts = [ray] if ray.geom_type == 'Point' else list(getattr(ray, 'geoms', []))
            if pts: ws.append(min(np.linalg.norm(np.array(g.coords[0]) - m_) for g in pts))
        w = max(max(ws) if ws else 20.0, 8.0)          # never steeper than the drop over 8 mm
        slope = P['bevel_drop'] / w
        R = 400.0
        sector = [p, q, q + bq * R, p + bp * R]
        zs = lambda pt: z_top - down * slope * max(0.0, float(np.dot(pt - p, n)))
        top = z_top + down * 200
        pts3 = [(x, y, zs(np.array([x, y]))) for x, y in sector] + [(x, y, top) for x, y in sector]
        cut.append(B.hull(pts3))
    return B.union(cut)

def body_blank():
    return B.prism(OUTL, 0.0, T)

# ----------------------------------------------------------------------------- NECK GEOMETRY (straight, no angle)
FB_TOP_Z = T - P['pocket_d'] + P['heel_t']                 # fingerboard surface (centre) at the heel
FRET_X = lambda n: P['scale'] * (1 - 2 ** (-n / 12.0))       # distance from the nut
Y_NUT = Y_END + FRET_X(P['fret_body'])
Y_SADDLE = Y_NUT - (P['scale'] + P['comp'])
Y_HEEL = Y_END - P['pocket_L']
Y_FB_END = Y_NUT - FRET_X(P['frets']) - 6.0
Z_NUT = FB_TOP_Z + P['fret_h'] + 0.5
Z12 = FB_TOP_Z + P['fret_h'] + P['action12']
Z_SADDLE = Z_NUT + (Z12 - Z_NUT) * (Y_NUT - Y_SADDLE) / (FRET_X(12))

def string_z(y): return Z_NUT + (Z12 - Z_NUT) * (Y_NUT - y) / FRET_X(12)

# ----------------------------------------------------------------------------- ROUTES
def rrect_at(cx, cy, L, W, r, ang=0.0):
    pts = B.rrect(L, W, r, n=8)
    c, s = math.cos(math.radians(ang)), math.sin(math.radians(ang))
    return B.ccw(np.array([(cx + x * c - y * s, cy + x * s + y * c) for x, y in pts]))

def neck_pocket():
    w, L, d, r = P['pocket_w'], P['pocket_L'], P['pocket_d'], P['pocket_r']
    rect = np.array([(-w / 2, Y_HEEL), (w / 2, Y_HEEL), (w / 2, Y_END + 20), (-w / 2, Y_END + 20)])
    pk = B.prism(rect, T - d, T + 50)
    # dog-bone reliefs so a square T-style heel seats against the pocket end with a 1/4" cutter
    for sx in (-1, 1):
        pk = pk + m3.Manifold.cylinder(d + 50, r, r, 32).translate([sx * (w / 2 - r * 0.6), Y_HEEL + r * 0.6, T - d])
    return pk

NECK_SCREWS = [(sx * 19.0, Y_END - 12.7 - k * 50.8) for sx in (-1, 1) for k in (0, 1)]   # drill to your neck plate

BRIDGE = dict(W=82.5, L=92.0, back=25.0, ang=0.0)             # vintage 3-saddle plate stand-in: rear edge 25 mm behind the saddles
BRIDGE_PU = dict(cy=Y_SADDLE + 35.0, L=76.0, W=41.0, d=19.0, ang=-8.0)   # route for a Tele bridge pickup (baseplate 2.87 x 1.51 in), slanted: treble toward the bridge
NECK_PU = dict(cy=Y_HEEL - 22.0, L=76.0, W=22.0, d=19.0)       # route for a Tele neck pickup (cover 2.56 x 0.62 in) with its two mounting tabs; screwed to the body on springs
FERRULE_Y = Y_SADDLE - 20.0
CTRL = dict(cx=101.0, cy=176.0, L=139.7, W=31.75, cav_L=118.0, cav_W=28.0, d=36.0)     # 5-1/2 x 1-1/4" plate

def ctrl_angle():
    # parallel to the treble lower-bout bevel line
    a, b = TREB_V[1], TREB_V[2]
    return math.degrees(math.atan2(b[1] - a[1], b[0] - a[0]))

def routes():
    r = []
    r.append(neck_pocket())
    r.append(B.prism(rrect_at(0, BRIDGE_PU['cy'], BRIDGE_PU['L'], BRIDGE_PU['W'], 6.0, BRIDGE_PU['ang']), T - BRIDGE_PU['d'], T + 50))
    r.append(B.prism(rrect_at(0, NECK_PU['cy'], NECK_PU['L'], NECK_PU['W'], 12.0), T - NECK_PU['d'], T + 50))
    ang = ctrl_angle()
    r.append(B.prism(rrect_at(CTRL['cx'], CTRL['cy'], CTRL['cav_L'], CTRL['cav_W'], 6.0, ang), T - CTRL['d'], T + 50))
    # string-through: 1/8" holes and ferrule counterbores on the back
    xs = np.linspace(-P['string_spacing_bridge'] / 2, P['string_spacing_bridge'] / 2, 6)
    for x in xs:
        r.append(m3.Manifold.cylinder(T + 20, 1.6, 1.6, 16).translate([x, FERRULE_Y, -10]))
        r.append(m3.Manifold.cylinder(6.5, 4.1, 4.1, 24).translate([x, FERRULE_Y, -1]))
    # neck screws (through, 3/16" clearance) from the back
    for x, y in NECK_SCREWS:
        r.append(m3.Manifold.cylinder(T + 20, 2.4, 2.4, 20).translate([x, y, -10]))
    return B.union(r)

# ----------------------------------------------------------------------------- BODY
def build_body(bevel=True):
    body = body_blank()
    if bevel:
        body = body - bevel_band(OUTL, SHARP, names, T, +1) - bevel_band(OUTL, SHARP, names, 0.0, -1)
    # keep the pocket walls full height: bevels stop short of the neck end (the chains end before x = +-40)
    return body - routes()

# ----------------------------------------------------------------------------- HARDWARE STAND-INS
def build_bridge():
    """vintage-style 3-saddle plate stand-in: flat plate with raised sides and back, pickup opening, 3 barrel saddles"""
    w, L = BRIDGE['W'], BRIDGE['L']; y0 = Y_SADDLE - BRIDGE['back']; t = 1.2; hw = 12.0
    rect = lambda x0, x1, ya, yb: np.array([(x0, ya), (x1, ya), (x1, yb), (x0, yb)])
    base = B.prism(rect(-w / 2, w / 2, y0, y0 + L), T, T + t)
    sides = B.prism(rect(-w / 2, -w / 2 + t, y0, y0 + 42), T, T + hw) + B.prism(rect(w / 2 - t, w / 2, y0, y0 + 42), T, T + hw)
    back = B.prism(rect(-w / 2, w / 2, y0, y0 + t), T, T + hw)
    hole = B.prism(rrect_at(0, BRIDGE_PU['cy'], BRIDGE_PU['L'] - 4, BRIDGE_PU['W'] - 4, 5.0, BRIDGE_PU['ang']), T - 1, T + 5)
    plate = (base + sides + back) - hole
    saddles = []
    sl = (w - 2 * t - 4) / 3
    for k in range(3):
        x0 = -w / 2 + t + 2 + k * sl
        saddles.append(m3.Manifold.cylinder(sl - 1.5, 3.2, 3.2, 24).rotate([0, 90, 0]).translate([x0 + 0.75, Y_SADDLE, Z_SADDLE - 3.2]))
    return plate + B.union(saddles)

def _rot_about(pts, cx, cy, ang):
    c, s_ = math.cos(math.radians(ang)), math.sin(math.radians(ang))
    return np.array([(cx + x * c - y * s_, cy + x * s_ + y * c) for x, y in pts])

def build_bridge_pickup():
    """Tele bridge pickup: copper-plated steel baseplate (2.87 x 1.51 in, three mounting ears), black flatwork
    bobbin, six alnico poles; slanted in the bridge plate. Top sits 3 mm under the strings."""
    p = BRIDGE_PU; cy = p['cy']; ang = p['ang']
    ztop = string_z(cy) - 3.0; zb = T - 14.0
    plate = rrect_at(0, cy, 72.9, 38.4, 5.0, ang)
    ears = [_rot_about([(0, -38.4 / 2 - 3.0)], 0, cy, ang)[0], _rot_about([(-21.8, 38.4 / 2 + 2.5)], 0, cy, ang)[0], _rot_about([(21.8, 38.4 / 2 + 2.5)], 0, cy, ang)[0]]
    base = B.prism(plate, zb, zb + 1.0) + B.union([m3.Manifold.cylinder(1.0, 4.5, 4.5, 20).translate([e[0], e[1], zb]) for e in ears])
    bobbin = B.prism(rrect_at(0, cy, 69.0, 19.5, 6.0, ang), zb + 1.0, ztop)
    flat = B.prism(rrect_at(0, cy, 71.0, 21.5, 6.5, ang), ztop - 1.6, ztop)
    poles = B.union([m3.Manifold.cylinder(1.2, 2.5, 2.5, 16).translate([*_rot_about([(x, 0)], 0, cy, ang)[0], ztop - 1.0])
                     for x in np.linspace(-26.25, 26.25, 6)])
    return bobbin + flat, poles, base

def build_neck_pickup():
    """Tele neck pickup: chrome cover 2.56 x 0.62 in with round ends, two mounting tabs; top 3 mm under the strings"""
    p = NECK_PU; cy = p['cy']; ztop = string_z(cy) - 3.0
    cover = B.prism(rrect_at(0, cy, 65.1, 15.9, 7.9, 0.0), T - 6.0, ztop - 1.5) + \
            B.hull([(x, cy + y, z) for x, y in B.rrect(62.0, 13.0, 6.4) for z in (ztop - 1.5, ztop)])
    tabs = B.prism(rrect_at(0, cy, 75.0, 11.0, 4.0, 0.0), T - 6.0, T - 5.0)
    return cover + tabs

def build_control_plate():
    ang = ctrl_angle()
    plate = B.prism(rrect_at(CTRL['cx'], CTRL['cy'], CTRL['L'], CTRL['W'], 6.0, ang), T, T + 1.6)
    c, s = math.cos(math.radians(ang)), math.sin(math.radians(ang))
    at = lambda u: (CTRL['cx'] + u * c, CTRL['cy'] + u * s)
    knobs = []
    for u, kind in ((-45.0, 'switch'), (0.0, 'vol'), (38.0, 'tone')):
        x, y = at(u)
        if kind == 'switch':
            knobs.append(m3.Manifold.cylinder(14, 1.6, 1.6, 12).translate([x, y, T + 1.6]) + m3.Manifold.cylinder(10, 5, 5, 24).translate([x, y, T + 15.6]))
        else:
            knobs.append(m3.Manifold.cylinder(17, 9.5, 9.5, 40).translate([x, y, T + 1.6]))
    return plate, B.union(knobs)

# ----------------------------------------------------------------------------- NECK (stock T-style stand-in, Juris head cut from a paddle)
def build_neck():
    """reference model of a stock 25.5" T-style neck: square heel 2.2 x 3", 9.5" board, 21 frets, flat headstock
    with the Juris outline cut from a paddle head"""
    fbt = P['fb_t']; z_fbp = FB_TOP_Z - fbt          # fingerboard gluing plane = neck wood top (flat neck)
    def width(y):
        t = (Y_NUT - y) / (Y_NUT - Y_FB_END)
        return P['nut_w'] + (P['heel_w'] - P['nut_w']) * min(1, max(0, t))
    def depth(y):
        s = Y_NUT - y; d = P['neck_d1'] + (P['neck_d12'] - P['neck_d1']) * (s - FRET_X(1)) / (FRET_X(12) - FRET_X(1))
        return d - fbt
    import neck_smooth as NS
    heel_below = z_fbp - (T - P['pocket_d'])                    # heel thickness under the gluing plane
    neck, face, tilt = NS.build(Y_NUT, Y_END, z_fbp, width, depth, head_angle_deg=P['head_angle'], head_T=14.0, volute=4.0, holes_on=build_neck.holes,
                             heel_extra=max(0.0, heel_below - depth(Y_END)), heel_w=P['heel_w'], heel_len=40.0, heel_p=0.3)
    # square heel inside the pocket only (flush with the body end, nothing sticks out)
    heel = B.prism(np.array([(-P['heel_w'] / 2, Y_HEEL), (P['heel_w'] / 2, Y_HEEL), (P['heel_w'] / 2, Y_END), (-P['heel_w'] / 2, Y_END)]),
                   T - P['pocket_d'], z_fbp)
    neck = neck + heel
    for x, y in NECK_SCREWS:
        neck = neck - m3.Manifold.cylinder(25, 1.6, 1.6, 12).translate([x, y, T - P['pocket_d'] - 1])
    # truss rod channel, 1/4 x 3/8 in, adjusting at the headstock through the access slot
    neck = neck - B.prism(np.array([(-3.175, Y_HEEL + 12), (3.175, Y_HEEL + 12), (3.175, Y_NUT + 22), (-3.175, Y_NUT + 22)]), z_fbp - 9.5, z_fbp + 1)
    build_neck.face = face; build_neck.tilt = tilt
    return neck

build_neck.holes = True

def build_fingerboard():
    fbt = P['fb_t']; R = P['fb_radius']; z0 = FB_TOP_Z - fbt
    def w(y):
        t = (Y_NUT - y) / (Y_NUT - Y_FB_END); return P['nut_w'] + (P['heel_w'] - P['nut_w']) * t
    pts = []
    for y in (Y_NUT, Y_FB_END):
        hw = w(y) / 2
        for x in np.linspace(-hw, hw, 13):
            zt = FB_TOP_Z - (R - math.sqrt(R * R - x * x))
            pts += [(x, y, z0), (x, y, zt)]
    fb = B.hull(pts)
    frets = []
    for n in range(1, P['frets'] + 1):
        y = Y_NUT - FRET_X(n); hw = w(y) / 2 - 0.5
        frets.append(B.hull([(x, y + dy, FB_TOP_Z - (R - math.sqrt(R * R - x * x)) + dz) for x in (-hw, 0, hw) for dy in (-1.2, 1.2) for dz in (-0.2, P['fret_h'] * (1 if abs(dy) < 1 else 0.4))]))
    return fb, B.union(frets)

def head_pt(x, s, n):
    """point on the tilted headstock: x across, s from the nut line, n above the face"""
    zf = FB_TOP_Z - P['fb_t']; a = math.radians(P['head_angle'])
    return np.array([x, Y_NUT + s * math.cos(a) - n * math.sin(-a), zf - s * math.sin(a) + n * math.cos(a)])

def build_tuners():
    out = []
    z_face = FB_TOP_Z - P['fb_t']
    for v, s_, sv in H.TUNER_LIST:
        y = Y_NUT + s_
        out.append(m3.Manifold.cylinder(18, 3.0, 3.0, 20).translate([v, y, z_face]))
        out.append(B.hull([(v + sv * dx, y + dy, z_face - 14 - dz) for dx in (2, 26) for dy in (-7, 7) for dz in (0, 14)]))
    return build_neck.tilt(B.union(out))

def build_strings():
    xs_b = np.linspace(-26.25, 26.25, 6); xs_n = np.linspace(-17.5, 17.5, 6)
    out = []; lines = []
    for i in range(6):
        a = np.array([xs_b[i], FERRULE_Y, T + 1.5]); s = np.array([xs_b[i], Y_SADDLE, Z_SADDLE]); nut = np.array([xs_n[i], Y_NUT, Z_NUT])
        side = -1 if i < 3 else 1; k = i if i < 3 else 5 - i
        post = head_pt(side * H.TUNER_POS[k][0], H.TUNER_POS[k][1], 8.0)
        pts = [a, s, nut, post]; lines.append([p.tolist() for p in pts])
        for p, q in zip(pts[:-1], pts[1:]):
            out.append(B.hull([p + d for d in B.sphere_pts(0.3)] + [q + d for d in B.sphere_pts(0.3)]))
    return B.union(out), lines

def build_jack():
    """Tele cup jack on the treble tail edge (drilled by hand, 7/8")"""
    a, b = SHARP[names.index('tailR')], SHARP[names.index('tTail')]
    m = a + (b - a) * 0.55; d = (b - a) / np.linalg.norm(b - a); n = np.array([d[1], -d[0]])
    ang = math.degrees(math.atan2(n[1], n[0]))
    cup = m3.Manifold.cylinder(3, 11.1, 11.1, 32).rotate([0, 90, 0]).rotate([0, 0, ang]).translate([m[0], m[1], T / 2])
    bore = m3.Manifold.cylinder(30, 11.1, 11.1, 32).rotate([0, 90, 0]).rotate([0, 0, ang]).translate([*(m - n * 27), T / 2])
    return cup, bore, m, n

# ----------------------------------------------------------------------------- MAIN
def main():
    import trimesh
    parts = {}
    body = build_body()
    cup, bore, jm, jn = build_jack()
    parts['body'] = body - bore
    parts['bridge'] = build_bridge()
    bpu, poles, bbase = build_bridge_pickup(); parts['bridge_pickup'] = bpu; parts['pickup_poles'] = poles; parts['bridge_pickup_base'] = bbase
    parts['neck_pickup'] = build_neck_pickup()
    plate, knobs = build_control_plate(); parts['control_plate'] = plate; parts['knobs'] = knobs
    parts['jack'] = cup
    parts['neck'] = build_neck(); parts['faceplate'] = build_neck.face
    fb, frets = build_fingerboard(); parts['fingerboard'] = fb; parts['frets'] = frets
    parts['tuners'] = build_tuners()
    st, lines = build_strings(); parts['strings'] = st
    MAT = dict(bridge_pickup_base='brass', faceplate='maple', body='ash', bridge='nickel', bridge_pickup='black', pickup_poles='nickel', neck_pickup='nickel',  control_plate='nickel',
               knobs='nickel', jack='nickel', neck='maple', fingerboard='ebony', frets='nickel', tuners='nickel', strings='steel')
    GROUP = dict(bridge_pickup_base='Pickups', faceplate='Neck', body='Body', bridge='Hardware', bridge_pickup='Pickups', pickup_poles='Pickups', neck_pickup='Pickups', control_plate='Hardware',
                 knobs='Hardware', jack='Hardware', neck='Neck', fingerboard='Neck', frets='Neck', tuners='Hardware', strings='Strings')
    EXPLODE = dict(Body=(0, 0, 0), Neck=(0, 60, 40), Pickups=(0, 0, 50), Hardware=(0, 0, 80), Strings=(0, 0, 110))
    scene = trimesh.Scene(); man_ok = {}; vols = {}
    for k, m in parts.items():
        tm = B.to_trimesh(m); man_ok[k] = bool(tm.is_watertight); vols[k] = m.volume()
        scene.add_geometry(tm, node_name=k, geom_name=k, metadata=dict(material=MAT[k], group=GROUP[k], explode=EXPLODE[GROUP[k]]))
        if k in ('body', 'neck'): tm.export(os.path.join(OUT, f'{k}.stl'))
    scene.export(os.path.join(OUT, 'juris_c75t.glb'))
    DENS = dict(ash=480, maple=700)
    mass = dict(body=vols['body'] * 480e-6, neck=vols['neck'] * 700e-6 + vols['fingerboard'] * 1150e-6)
    hw = dict(bridge=130, pickups=150, control=110, tuners=198, strings=10, frets=20, jack=10, neckplate=40)
    total = sum(mass.values()) + sum(hw.values())
    report = dict(units='mm', params=P, body_length=float(SHARP[:, 1].max()), lower_bout=float(OUTL[:, 0].max() - OUTL[:, 0].min()),
                  y_end=float(Y_END), y_nut=float(Y_NUT), y_saddle=float(Y_SADDLE), fb_top_above_body=float(FB_TOP_Z - T),
                  string_above_body_at_saddle=float(Z_SADDLE - T), saddle_z=float(Z_SADDLE),
                  masses_g={k: round(v, 1) for k, v in mass.items()}, hardware_g=hw, total_mass_g=round(total), total_mass_lb=round(total / 453.6, 2),
                  manifold=man_ok, bevel=dict(bass=BASS_V.round(1).tolist(), treble=TREB_V.round(1).tolist()),
                  ctrl_angle_deg=ctrl_angle(), jack=dict(at=jm.tolist(), normal=jn.tolist()))
    json.dump(report, open(os.path.join(OUT, 'report.json'), 'w'), indent=1)
    meta = {k: dict(material=MAT[k], group=GROUP[k], explode=list(EXPLODE[GROUP[k]])) for k in parts}
    json.dump(dict(meta=meta, strings=lines), open(os.path.join(OUT, 'meta.json'), 'w'))
    print(json.dumps({k: report[k] for k in ('body_length', 'lower_bout', 'y_saddle', 'fb_top_above_body', 'string_above_body_at_saddle', 'total_mass_lb')}))
    print('manifold', man_ok)
    return parts

if __name__ == '__main__':
    main()
