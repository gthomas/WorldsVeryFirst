"""Juris C-75 "The Plaintiff" — the original chambered set-neck guitar, modelled from its technical diagram and the
2D body drawing (Juris_C_75_Body_2D.pdf): 2" chambered mahogany body with a maple cap and the bevelled edges,
24.75" set neck at 4 degrees, '59 roundback, two surface-mount T-Armond style pickups, a floating bridge on a
flat-bottom base, and the geometric trapeze tailpiece. Display model (chambers are representative).
Run: python c75.py -> out_c/"""
import json, math, os
import numpy as np
import manifold3d as m3
import build as B
import head_c75 as H
import outline as O
import tstyle as TS                     # bevel chains + cutter

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out_c'); os.makedirs(OUT, exist_ok=True)
IN = 25.4
P = dict(T=50.8, cap=6.35, bevel_drop=11.4, scale=628.65, comp=2.5, frets=22, neck_angle=4.0, head_angle=13.0,
         fb_t=6.35, fb_radius=304.8, nut_w=42.86, fb_end_w=56.4, neck_d1=21.8, neck_d12=24.4, fret_h=1.2, action12=2.0,
         bridge_target=24.0, tail_half=16.0, neck_half=28.0)

# ---- outline (original C-75: narrow tail flat, neck-width tongue)
O.TAIL_HALF = P['tail_half']; O.NECK_HALF = P['neck_half']
names, SHARP = O.sharp_outline()
RAD = dict(O.RADII); RAD.update(neckT=0, neckB=0)
OUTL = B.ccw(O.fillet(names, SHARP, RAD))
Y_END = SHARP[names.index('neckT')][1]
T = P['T']

# drawing pixel (pdf.png, 720 x 1080 preview) -> body mm, same calibration as the bevel trace
def _px_map():
    import cv2
    c = np.load('/home/claude/contour.npy'); x0, y0, w, h = cv2.boundingRect(c[:, None, :].astype(np.int32))
    sx, sy = w / 13.5, h / 18.0; cx = x0 + w / 2; ybot = y0 + h
    L = {k: O._line(k) for k in O._E}
    t_tail = O._at_x(L['t5'], O.TAIL_HALF); b_tail = O._at_x(L['b5'], -O.TAIL_HALF); ytail = (t_tail[1] + b_tail[1]) / 2
    ymin_in = min(min(a[1], b[1]) for a, b in O._E.values())
    k = 2400 / 720
    return lambda px, py: np.array([(px * k - cx) / sx * IN, ((ybot - py * k) / sy + ymin_in) * IN - ytail])
PX = _px_map()

NECK_PU = PX(386, 323); BRIDGE_PU = PX(386, 513); BRIDGE_C = PX(386, 587)
Y_SADDLE = float(BRIDGE_C[1])
KNOBS = [(PX(557, 621), 'vol'), (PX(586, 710), 'tone')]     # the built body has these two (the drawing's other two are not drilled)
TOGGLE = PX(246, 225)
TONE_SW = np.array([85.9, 142.9])   # Gretsch-style 3-way tone switch, between the knobs and in toward the centre (as in the photo)
GUARD = np.array([PX(446, 268), PX(493, 300), PX(578, 548), PX(467, 532), PX(446, 470)])
TP_BAR_Y = float(PX(386, 675)[1])

# ---- neck geometry: built straight (fb top at Z_J over the body end), then tilted back about the joint
TH = math.radians(P['neck_angle'])
FX = lambda n: P['scale'] * (1 - 2 ** (-n / 12.0))
Y_NUT_PLAN = Y_SADDLE + P['scale'] + P['comp']                # nut position along the body (plan)
U_NUT = (Y_NUT_PLAN - Y_END) / math.cos(TH)                  # straight-frame distance from the joint to the nut

def rot_pts(pts, ang, pivot):
    """rotate points about an x-parallel axis through pivot (positive = nose up toward +z for +y)"""
    p = np.asarray(pts, float) - pivot; c, s = math.cos(ang), math.sin(ang)
    out = p.copy(); out[..., 1] = p[..., 1] * c - p[..., 2] * s; out[..., 2] = p[..., 1] * s + p[..., 2] * c
    return out + pivot

def rot_man(man, ang, pivot):
    return man.translate([-pivot[0], -pivot[1], -pivot[2]]).rotate([math.degrees(ang), 0, 0]).translate(list(pivot))

def solve_joint():
    """height of the fingerboard over the body at the joint that gives the target bridge height"""
    def bridge_h(zj):
        piv = np.array([0, Y_END, zj])
        nut = rot_pts([0, Y_END + U_NUT, zj + P['fret_h'] + 0.5], -TH, piv)
        f12 = rot_pts([0, Y_END + U_NUT - FX(12), zj + P['fret_h'] + P['action12']], -TH, piv)
        z = nut[2] + (f12[2] - nut[2]) * (Y_SADDLE - nut[1]) / (f12[1] - nut[1])
        return z - T
    lo, hi = T - 20, T + 40
    for _ in range(60):
        mid = (lo + hi) / 2
        if bridge_h(mid) < P['bridge_target']: lo = mid
        else: hi = mid
    return (lo + hi) / 2, bridge_h((lo + hi) / 2)
Z_J, BRIDGE_H = solve_joint()
PIV = np.array([0.0, Y_END, Z_J])
Y_N = Y_END + U_NUT                                           # nut, straight frame

def string_z_body(y):
    nut = rot_pts([0, Y_N, Z_J + P['fret_h'] + 0.5], -TH, PIV); f12 = rot_pts([0, Y_N - FX(12), Z_J + P['fret_h'] + P['action12']], -TH, PIV)
    return nut[2] + (f12[2] - nut[2]) * (y - nut[1]) / (f12[1] - nut[1])

# ---- body
def build_body():
    TS.P['bevel_drop'] = P['bevel_drop']
    body = B.prism(OUTL, 0.0, T)
    body = body - TS.bevel_band(OUTL, SHARP, names, T, +1) - TS.bevel_band(OUTL, SHARP, names, 0.0, -1)
    body = body - TS.bevel_band(OUTL, SHARP, names, T, +1, round_r=1.6) - TS.bevel_band(OUTL, SHARP, names, 0.0, -1, round_r=1.6)
    # chambers in both wings under the cap, leaving a solid centre block under the pickups and bridge
    from shapely.geometry import Polygon, box
    inner = Polygon(np.r_[TS.bevel_chains(OUTL, SHARP, names)[0], TS.bevel_chains(OUTL, SHARP, names)[1][::-1]]).buffer(0)
    ch = []
    for side in (-1, 1):
        wing = inner.intersection(box(45, 30, 200, 380) if side > 0 else box(-200, 30, -45, 400)).buffer(-12, join_style=1)
        for kp, _ in KNOBS + [(TONE_SW, 'sw')]:
            wing = wing.difference(Polygon(B.rrect(70, 70, 10)).buffer(0).__class__([(kp[0] - 25, kp[1] - 25), (kp[0] + 25, kp[1] - 25), (kp[0] + 25, kp[1] + 25), (kp[0] - 25, kp[1] + 25)]))
        geoms = getattr(wing, 'geoms', [wing])
        for g in geoms:
            if g.area > 500: ch.append(B.prism(B.ccw(np.array(g.exterior.coords)[:-1]), 6.0, T - P['cap']))
    body = body - B.union(ch) - pickup_route(NECK_PU) - pickup_route(BRIDGE_PU)
    cap = body ^ B.prism(OUTL, T - P['cap'], T + 1)
    return body - B.prism(OUTL, T - P['cap'], T + 1), cap

# ---- neck (straight frame, then tilted)
def build_neck():
    fbt = P['fb_t']; zgp = Z_J - fbt                          # gluing plane at the joint
    y_fb_end = Y_N - FX(P['frets']) - 6.0
    def width(y):
        t = (Y_N - y) / (Y_N - y_fb_end); return P['nut_w'] + (P['fb_end_w'] - P['nut_w']) * min(1, max(0, t))
    def depth(y):
        s = Y_N - y; return P['neck_d1'] + (P['neck_d12'] - P['neck_d1']) * (s - FX(1)) / (FX(12) - FX(1)) - fbt
    HT = 14.0; HA = math.radians(P['head_angle']); hpiv = np.array([0, Y_N, zgp])
    import neck_smooth as NS
    # heel flares deeper toward the body end so it meets the end face about 6 mm above the back
    heel_extra = max(0.0, (zgp - 8.0) - depth(Y_END))
    neck, face, _ = NS.build(Y_N, Y_END, zgp, width, depth, head_angle_deg=P['head_angle'], head_T=HT, volute=4.0,
                             heel_extra=heel_extra, heel_w=P['fb_end_w'], heel_len=48.0, heel_p=0.45)
    tenon = B.hull([(sx * 26.0, y, z) for sx in (-1, 1) for y in (Y_END - 70, Y_END - 0.5) for z in (zgp, T - 32)])   # glued in, hidden in the body
    neck = neck + tenon
    neck_face = face
    tun = []
    for v, s, sv in H.TUNER_LIST:
        p = np.array([v, Y_N + s, zgp])
        post = rot_man(m3.Manifold.cylinder(16, 3.0, 3.0, 20).translate([v, Y_N + s, zgp]), -HA, hpiv)
        btn = rot_man(B.hull([(v + sv * dx, Y_N + s + dy, zgp - HT - dz) for dx in (2, 26) for dy in (-7, 7) for dz in (0, 13)]), -HA, hpiv)
        tun += [post, btn]
    badge = rot_man(B.prism(B.rrect(24, 30, 2.0) + np.array([0, Y_N + 160]), zgp + 0.5, zgp + 1.3), -HA, hpiv)
    tcov, tscr = build_truss_cover(hpiv, zgp, HA)
    # fingerboard + frets
    R = P['fb_radius']; pts = []
    for y in (Y_N, y_fb_end):
        hw = width(y) / 2
        for x in np.linspace(-hw, hw, 13):
            pts += [(x, y, zgp), (x, y, Z_J - (R - math.sqrt(R * R - x * x)))]
    fb = B.hull(pts)
    frets = []
    for n in range(1, P['frets'] + 1):
        y = Y_N - FX(n); hw = width(y) / 2 - 0.5
        frets.append(B.hull([(x, y + dy, Z_J - (R - math.sqrt(R * R - x * x)) + dz) for x in (-hw, 0, hw) for dy in (-1.2, 1.2) for dz in (-0.2, P['fret_h'] * (1 if abs(dy) < 1 else .4))]))
    nut = B.hull([(x, y, z) for x in (-21.4, 21.4) for y in (Y_N, Y_N + 4) for z in (zgp, Z_J + P['fret_h'] + 0.5)])
    tilt = lambda m: rot_man(m, -TH, PIV)
    tuner_posts = [rot_pts(rot_pts([v, Y_N + s, zgp + 8], -HA, hpiv), -TH, PIV) for v, s, _ in H.TUNER_LIST]
    return dict(neck=tilt(neck), faceplate=tilt(neck_face), fingerboard=tilt(fb), frets=tilt(B.union(frets)), nut=tilt(nut), tuners=tilt(B.union(tun)), badge=tilt(badge), truss_cover=tilt(tcov), truss_screws=tilt(tscr)), tuner_posts

# ---- hardware
HB_ROUTE = dict(L=73.0, W=38.0, ear_L=89.0, ear_W=10.0, d=19.0)       # as routed on the built body

def pickup_route(c):
    y = float(c[1]); r = HB_ROUTE
    return B.prism(B.rrect(r['L'], r['W'], 3.2) + np.array([0, y]), T - r['d'], T + 5) + \
           B.prism(B.rrect(r['ear_L'], r['ear_W'], 3.0) + np.array([0, y]), T - r['d'], T + 5)

def _circ(r, n=72): return np.array([(r * math.cos(t), r * math.sin(t)) for t in np.linspace(0, 2 * math.pi, n, endpoint=False)])
def _chull(rings): return B.hull([(float(x), float(y), z) for pts, z in rings for x, y in pts])
def _sp(r, n=10):
    return [np.array([r * math.sin(a) * math.cos(b), r * math.sin(a) * math.sin(b), r * math.cos(a)]) for a in np.linspace(0, math.pi, n) for b in np.linspace(0, 2 * math.pi, 2 * n, endpoint=False)]

def build_pickup(c):
    """TV Jones T-Armond, wood mount, no ears: black base plate screwed to the top at its corners, chrome cover
    with two staggered rows of six flat pole pieces and the logo bar"""
    y = float(c[1]); zt = string_z_body(y) - 3.5
    plate = B.prism(B.rrect(86, 46, 2.5) + np.array([0, y]), T, T + 2.0)
    screws = B.union([B.hull([(sx * 39.5 + dx, y + sy * 19.5 + dy, T + 2.0 + dz) for dx in (-1.6, 1.6) for dy in (-1.6, 1.6) for dz in (0, 0.7)])
                      for sx in (-1, 1) for sy in (-1, 1)])
    cover = _chull([(B.rrect(70.6, 38.6, 4.5) + np.array([0, y]), T + 2.0), (B.rrect(70.6, 38.6, 4.5) + np.array([0, y]), zt - 1.2),
                    (B.rrect(68.2, 36.2, 3.4) + np.array([0, y]), zt)])
    xs = np.linspace(-26.2, 26.2, 6)
    poles = B.union([B.prism(B.rrect(6.4, 4.6, 1.2) + np.array([x + off, y + yo]), zt - 0.1, zt + 0.45)
                     for yo, off in ((-6.5, -2.0), (6.5, 2.0)) for x in xs])
    logo = B.prism(B.rrect(22, 3.2, 1.0) + np.array([-18, y]), zt - 0.1, zt + 0.25)
    return plate, cover, poles + logo + screws

def build_bridge():
    """rosewood base on the top; brass Compton one-piece compensated bridge (dog-bone ends over the posts) sitting
    on brass thumbwheels"""
    y = Y_SADDLE; L, W, H = 135.0, 14.95, 11.75
    base = B.hull([(x * L / 2, y + yy * W / 2, T + z) for x in (-1, 1) for yy in (-1, 1) for z in (0, 3)] +
                  [(x * 30, y + yy * W / 2, T + H) for x in (-1, 1) for yy in (-1, 1)])
    zw = T + H + 1.0
    wheels = B.union([m3.Manifold.cylinder(3.2, 7.4, 7.4, 40).translate([sx * 37, y, zw]) for sx in (-1, 1)])
    posts = B.union([m3.Manifold.cylinder(BRIDGE_H - H - 2.5, 2.0, 2.0, 16).translate([sx * 37, y, T + H]) for sx in (-1, 1)])
    zb = zw + 3.2; ztop = T + BRIDGE_H
    lobes = B.union([m3.Manifold.cylinder(ztop - 2.2 - zb, 7.2, 7.2, 48).translate([sx * 37, y, zb]) for sx in (-1, 1)])
    bar = B.hull([(sx * 36, y + dy, zb + 0.3) for sx in (-1, 1) for dy in (-4.2, 4.2)] + [(sx * 36, y + dy, ztop - 3.2) for sx in (-1, 1) for dy in (-3.2, 3.2)])
    # compensated crown: set back toward the tail for the wound strings, the B set back from the plain G and E
    xs = np.linspace(-26.2, 26.2, 6); off = np.array([2.4, 1.8, 1.2, 0.2, 1.5, 0.0]) - 1.0
    crown = []
    for i in range(5):
        pts = []
        for x, o in ((xs[i], off[i]), (xs[i + 1], off[i + 1])):
            pts += [(x, y + o - 0.6, ztop), (x, y + o + 0.6, ztop), (x, y - 3.0, ztop - 3.4), (x, y + 3.0, ztop - 3.4)]
        if i == 3: pts = [p for p in pts if p[0] == xs[3]] + [(xs[3] + 0.6, y + off[4] - 0.6, ztop), (xs[3] + 0.6, y + off[4] + 0.6, ztop), (xs[3] + 0.6, y - 3.0, ztop - 3.4), (xs[3] + 0.6, y + 3.0, ztop - 3.4)]
        crown.append(B.hull(pts))
    ends = B.union([B.hull([(sx * 26.2, y + o - 0.6, ztop - 0.6), (sx * 26.2, y + o + 0.6, ztop - 0.6), (sx * 32, y - 2.6, ztop - 3.2), (sx * 32, y + 2.6, ztop - 3.2),
                            (sx * 26.2, y - 3.0, ztop - 3.4), (sx * 26.2, y + 3.0, ztop - 3.4)]) for sx, o in ((-1, off[0]), (1, off[5]))])
    holes = B.union([m3.Manifold.cylinder(30, 2.1, 2.1, 16).translate([sx * 37, y, zb - 1]) for sx in (-1, 1)])
    compton = (lobes + bar + B.union(crown) + ends) - holes
    return base, wheels + posts, compton

def build_tailpiece():
    """Allparts short trapeze: chrome hinge plate over the tail, one bent wire loop, and a solid string bar held to
    the wire ends with hex nuts"""
    yb = TP_BAR_Y; zbar0, zbar1 = T + 7.0, T + 17.0
    bar = B.hull([(sx * (39 - r), yb + sy * (8 - r), z) for sx in (-1, 1) for sy in (-1, 1) for z in (zbar0, zbar1) for r in (0, 1.5)
                  if (r == 0) != (z == zbar1)] + [(sx * 37.5, yb + sy * 6.5, zbar1) for sx in (-1, 1) for sy in (-1, 1)])
    sl = B.union([B.prism(B.rrect(1.6, 7, 0.7) + np.array([x, yb + 5.0]), zbar1 - 3, zbar1 + 1) for x in np.linspace(-22.5, 22.5, 6)])
    bar = bar - sl
    d = TS.P['bevel_drop'] if False else min(P['bevel_drop'], 12.0 * P['bevel_drop'] / 15.0)
    rod = lambda a, b: B.hull([tuple(np.array(a) + q) for q in _sp(2.3, 8)] + [tuple(np.array(b) + q) for q in _sp(2.3, 8)])
    hz = T + 3.2
    wire = B.union([rod((-15, 20, hz), (15, 20, hz))] + [rod((sx * 15, 20, hz), (sx * 33, yb - 2, zbar0 + 2)) for sx in (-1, 1)] +
                   [rod((sx * 33, yb - 2, zbar0 + 2), (sx * 33, yb + 2, zbar0 + 2)) for sx in (-1, 1)])
    nuts = B.union([m3.Manifold.cylinder(4.5, 4.4, 4.4, 6).translate([sx * 33, yb - 2, zbar1 - 0.3]) +
                    m3.Manifold.cylinder(6.5, 1.9, 1.9, 16).translate([sx * 33, yb - 2, zbar1]) for sx in (-1, 1)])
    # hinge plate: flat on the top, down the tail bevel, then down the tail edge
    W1, W0 = 36.0, 44.0
    flat = B.hull([(sx * W1 / 2, yy, T + z) for sx in (-1, 1) for yy in (12.0, 26.0) for z in (0.02, 1.6)])
    knuckle = B.hull([tuple(np.array((x, 22.0, hz)) + q) for x in (-17, 17) for q in _sp(3.0, 8)])
    slope = B.hull([(sx * W0 / 2, 0.4, T - d + z) for sx in (-1, 1) for z in (0.02, 1.6)] + [(sx * W1 / 2, 12.0, T + z) for sx in (-1, 1) for z in (0.02, 1.6)])
    side = B.hull([(sx * W0 / 2, yy, z) for sx in (-1, 1) for yy in (-1.4, 0.2) for z in (T - d - 26, T - d + 1.6)])
    screws = B.union([m3.Manifold.cylinder(1.6, 2.6, 2.4, 16).rotate([90, 0, 0]).translate([sx * 13, -1.3, T - d - 14 + dz]) for sx in (-1, 1) for dz in (0, -8)])
    return bar + wire + nuts + flat + knuckle + slope + side + screws, np.array([0, yb + 5.0, zbar1 - 0.5])

def _knob(p, h=15.0, r=9.6):
    """amber speckled ribbed knob: straight fluted skirt, softened top edge, set-screw boss"""
    body = _chull([(_circ(r) + p, T + 0.2), (_circ(r) + p, T + h - 1.8), (_circ(r - 0.7) + p, T + h - 0.6), (_circ(r - 1.9) + p, T + h)])
    ribs = B.union([B.prism(B.rrect(1.1, 1.6, 0.5) @ np.array([[math.cos(a), math.sin(a)], [-math.sin(a), math.cos(a)]]) + p + r * np.array([math.cos(a), math.sin(a)]), T, T + h - 2.0)
                    for a in np.linspace(0, 2 * math.pi, 30, endpoint=False)])
    return body - ribs

def _toggle(p, throw=1):
    """full-size toggle: chrome skirted dress nut over the hex nut, threaded bushing, lever with a Gretsch-style tip"""
    x, y = p
    skirt = m3.Manifold.cylinder(2.4, 9.0, 5.2, 48).translate([x, y, T])
    hexn = m3.Manifold.cylinder(1.8, 4.8, 4.8, 6).translate([x, y, T + 2.4])
    bush = m3.Manifold.cylinder(5.0, 3.1, 3.1, 24).translate([x, y, T + 2.4])
    a = math.radians(14) * throw; top = np.array([x, y + 17 * math.sin(a), T + 7.0 + 17 * math.cos(a)])
    lever = B.hull([tuple(np.array([x, y, T + 6.5]) + q) for q in _sp(1.1, 6)] + [tuple(top + q) for q in _sp(0.9, 6)])
    tip = B.hull([tuple(top + np.array([0, -1.0 * math.sin(a), -1.0 * math.cos(a)]) + q) for q in _sp(2.3, 10)] +
                 [tuple(top + np.array([0, 5.0 * math.sin(a), 5.0 * math.cos(a)]) + q) for q in _sp(2.6, 10)])
    return skirt + hexn + bush + lever + tip

def build_controls():
    knobs = B.union([_knob(np.array(p)) for p, _ in KNOBS])
    switches = _toggle(TONE_SW, 1) + _toggle(TOGGLE, -1)
    return knobs, switches

def build_jack():
    """Electrosocket on the treble side of the lower bout: round plate with two screws, plug socket in the middle"""
    a, b = SHARP[names.index('tailR')], SHARP[names.index('tTail')]
    m = a + (b - a) * 0.55; dd = (b - a) / np.linalg.norm(b - a); n = np.array([dd[1], -dd[0]])
    ang = math.degrees(math.atan2(n[1], n[0]))
    place = lambda man, out: man.rotate([0, 90, 0]).rotate([0, 0, ang]).translate([m[0] + n[0] * out, m[1] + n[1] * out, T / 2])
    plate = place(m3.Manifold.cylinder(1.5, 13.5, 13.0, 64), 0.0)
    rim = place(m3.Manifold.cylinder(2.6, 6.0, 6.0, 40), 0.0) - place(m3.Manifold.cylinder(4, 4.4, 4.4, 32), -0.5)
    screws = B.union([place(m3.Manifold.cylinder(0.9, 1.7, 1.4, 16).translate([0, 0, 9.5 * s]), 1.5) for s in (-1, 1)])
    bore = place(m3.Manifold.cylinder(40, 11.1, 11.1, 40), -38.0)
    return plate + rim + screws, bore

def build_truss_cover(hpiv, zgp, HA):
    """rosewood bell over the truss rod access slot, two screws"""
    bell = [(-6.0, 2.0), (6.0, 2.0)] + [(13.0 * math.sin(t) * (0.55 + 0.45 * t / (math.pi / 2)), 30.0 + 6.0 * math.cos(t)) for t in np.linspace(math.pi / 2, -math.pi / 2, 15)][::-1]
    bell = [(-6.0, 2.0), (6.0, 2.0), (12.0, 26.0), (13.0, 31.0), (11.5, 35.0), (6.0, 37.5), (0.0, 38.2), (-6.0, 37.5), (-11.5, 35.0), (-13.0, 31.0), (-12.0, 26.0)]
    pp = B.ccw(np.array(bell) + np.array([0, Y_N]))
    cov = B.prism(pp, zgp + 0.6, zgp + 3.0)
    scr = B.union([B.hull([(0 + dx, Y_N + s + dy, zgp + 3.0 + dz) for dx in (-1.5, 1.5) for dy in (-1.5, 1.5) for dz in (0, 0.6)]) for s in (5.0, 33.0)])
    return rot_man(cov, -HA, hpiv), rot_man(scr, -HA, hpiv)

def build_strings(anchor, posts):
    xs_b = np.linspace(-26.2, 26.2, 6); xs_n = np.linspace(-17.5, 17.5, 6)
    nut_z = string_z_body(Y_NUT_PLAN); sad_z = T + BRIDGE_H
    out = []
    for i in range(6):
        a = np.array([xs_b[i] * 0.86, anchor[1] + 4, anchor[2]]); s = np.array([xs_b[i], Y_SADDLE, sad_z])
        nut = rot_pts([xs_n[i], Y_N, Z_J + P['fret_h'] + 0.5], -TH, PIV)
        side = -1 if i < 3 else 1; k = i if i < 3 else 5 - i
        idx = [j for j, (v, s_, sv) in enumerate(H.TUNER_LIST) if sv == side][k] if False else None
        post = [p for (v, s_, sv), p in zip(H.TUNER_LIST, posts) if sv == side and abs(abs(v) - H.TUNER_POS[k][0]) < 1e-6 and abs(s_ - H.TUNER_POS[k][1]) < 1e-6][0]
        pts = [a, s, nut, post]
        for p, q in zip(pts[:-1], pts[1:]):
            out.append(B.hull([p + d for d in B.sphere_pts(0.3)] + [q + d for d in B.sphere_pts(0.3)]))
    return B.union(out)

def main():
    parts = {}
    body, cap = build_body(); parts['body'] = body; parts['cap'] = cap
    nk, posts = build_neck(); parts.update(nk)
    for nm, c in (('neck_pickup', NECK_PU), ('bridge_pickup', BRIDGE_PU)):
        ring, cover, poles = build_pickup(c)
        parts[nm + '_ring'] = ring; parts[nm] = cover; parts[nm + '_poles'] = poles
    br, wh, cp = build_bridge(); parts['bridge'] = br; parts['thumbwheels'] = wh; parts['compton'] = cp
    tp, anchor = build_tailpiece(); parts['tailpiece'] = tp
    parts['knobs'], parts['switches'] = build_controls()
    parts['jack'] = build_jack()[0]
    parts['strings'] = build_strings(anchor, posts)
    MAT = dict(faceplate='maple', body='mahogany_nat', cap='maple', cap_edge='maple_nat', neck='mahogany', fingerboard='ebony_brown', frets='fretgold', nut='bone', tuners='nickel',
               badge='nickel', truss_cover='rosewood', truss_screws='nickel', neck_pickup_ring='black_plastic', neck_pickup='chrome', neck_pickup_poles='chrome',
               bridge_pickup_ring='black_plastic', bridge_pickup='chrome', bridge_pickup_poles='chrome', bridge='rosewood', thumbwheels='brass', compton='brass',
               tailpiece='chrome', knobs='amber', switches='chrome', jack='chrome', strings='steel')
    GROUP = {k: ('Body' if k in ('body', 'cap') else 'Neck' if k in ('neck', 'faceplate', 'fingerboard', 'frets', 'nut', 'badge', 'truss_cover', 'truss_screws') else
                 'Pickups' if 'pickup' in k else 'Strings' if k == 'strings' else 'Hardware') for k in parts}
    EXPLODE = dict(Body=(0, 0, 0), Neck=(0, 60, 40), Pickups=(0, 0, 50), Hardware=(0, 0, 80), Strings=(0, 0, 110))
    import trimesh
    scene = trimesh.Scene(); ok = {}; vols = {}
    meshes = {k: B.to_trimesh(m) for k, m in parts.items()}
    # the stain only takes on the flat of the cap; where the bevels cut through it the maple shows as a light
    # line round the edge (the built guitar's faux binding)
    z0 = T - P['cap']; stripe = parts['cap'] ^ B.prism(OUTL, z0 - 1.0, z0 + 2.6)
    meshes['cap'] = B.to_trimesh(parts['cap'] - B.prism(OUTL, z0 - 1.0, z0 + 2.6)); meshes['cap_edge'] = B.to_trimesh(stripe)
    GROUP['cap_edge'] = 'Body'
    for k, tm in meshes.items():
        ok[k] = bool(tm.is_watertight) if k not in ('cap', 'cap_edge') else True; vols[k] = parts[k].volume() if k in parts else 0
        scene.add_geometry(tm, node_name=k, geom_name=k, metadata=dict(material=MAT[k], group=GROUP[k], explode=EXPLODE[GROUP[k]]))
    scene.export(os.path.join(OUT, 'juris_c75.glb'))
    mass = vols['body'] * 560e-6 + vols['cap'] * 700e-6 + vols['neck'] * 560e-6 + vols['fingerboard'] * 1150e-6 + 6 * 33 + 300 + 95 + 150
    nut_y = rot_pts([0, Y_N, Z_J], -TH, PIV)[1]
    joint_fret = 12 * math.log2(P['scale'] / (P['scale'] - (nut_y - Y_END)))
    report = dict(units='mm', params=P, body_length=float(SHARP[:, 1].max()), lower_bout=float(OUTL[:, 0].max() - OUTL[:, 0].min()),
                  bridge_height=float(BRIDGE_H), fb_above_body_at_joint=float(Z_J - T), neck_joint_fret=float(joint_fret), y_saddle=Y_SADDLE,
                  total_mass_g=round(mass), total_mass_lb=round(mass / 453.6, 2), manifold=ok)
    json.dump(report, open(os.path.join(OUT, 'report.json'), 'w'), indent=1)
    json.dump(dict(meta={k: dict(material=MAT[k], group=GROUP[k], explode=list(EXPLODE[GROUP[k]])) for k in meshes}), open(os.path.join(OUT, 'meta.json'), 'w'))
    print(json.dumps({k: report[k] for k in ('body_length', 'bridge_height', 'fb_above_body_at_joint', 'neck_joint_fret', 'y_saddle', 'total_mass_lb')}))
    print('manifold', ok)

if __name__ == '__main__':
    main()
