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
    TS.OUTL = OUTL; TS.P['bevel_drop'] = P['bevel_drop']
    body = B.prism(OUTL, 0.0, T)
    for verts in (TS.BASS_V, TS.TREB_V):
        body = body - TS.bevel_cutters(verts, T, +1) - TS.bevel_cutters(verts, 0.0, -1)
    # chambers in both wings under the cap, leaving a solid centre block under the pickups and bridge
    from shapely.geometry import Polygon, box
    inner = Polygon(np.r_[TS.BASS_V, TS.TREB_V[::-1]]).buffer(0)
    ch = []
    for side in (-1, 1):
        wing = inner.intersection(box(45, 30, 200, 380) if side > 0 else box(-200, 30, -45, 400)).buffer(-12, join_style=1)
        for kp, _ in KNOBS:
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
    return dict(neck=tilt(neck), faceplate=tilt(neck_face), fingerboard=tilt(fb), frets=tilt(B.union(frets)), nut=tilt(nut), tuners=tilt(B.union(tun)), badge=tilt(badge)), tuner_posts

# ---- hardware
HB_ROUTE = dict(L=73.0, W=38.0, ear_L=89.0, ear_W=10.0, d=19.0)       # as routed on the built body

def pickup_route(c):
    y = float(c[1]); r = HB_ROUTE
    return B.prism(B.rrect(r['L'], r['W'], 3.2) + np.array([0, y]), T - r['d'], T + 5) + \
           B.prism(B.rrect(r['ear_L'], r['ear_W'], 3.0) + np.array([0, y]), T - r['d'], T + 5)

def build_pickup(c):
    """surface-mount T-Armond style: flange on the top over the humbucker route, cover with six poles"""
    y = float(c[1]); zt = string_z_body(y) - 3.5
    flange = B.prism(B.rrect(95, 44, 6.0) + np.array([0, y]), T, T + 1.5)
    cover = B.prism(B.rrect(70, 38, 4.0) + np.array([0, y]), T + 1, zt)
    poles = B.union([m3.Manifold.cylinder(1.0, 2.6, 2.6, 16).translate([x, y, zt]) for x in np.linspace(-26.25, 26.25, 6)])
    return flange, cover, poles

def build_bridge():
    y = Y_SADDLE; L, W, H = 135.0, 14.95, 11.75
    base = B.hull([(x * L / 2, y + yy * W / 2, T + z) for x in (-1, 1) for yy in (-1, 1) for z in (0, 3)] +
                  [(x * 30, y + yy * W / 2, T + H) for x in (-1, 1) for yy in (-1, 1)])
    top = B.hull([(x, y + dy, T + H + 6 + dz) for x in (-44, 44) for dy in (-5, 5) for dz in (0, BRIDGE_H - H - 6)])
    wheels = B.union([m3.Manifold.cylinder(3.0, 7.5, 7.5, 32).translate([sx * 37, y, T + H + 1]) for sx in (-1, 1)])
    return base + top, wheels

def build_tailpiece():
    yb = TP_BAR_Y; zt = T + 4.0
    bar = B.hull([(x, yb + dy, zt + dz) for x in (-37, 37) for dy in (-3, 9) for dz in (-2.5, 2.5)])
    arms = B.union([B.hull([(sx * 34, yb, zt), (sx * 30, yb, zt), (sx * 11, 6.0, T + 2.5), (sx * 15, 6.0, T + 2.5)] +
                           [(sx * 34, yb, zt + 3), (sx * 30, yb, zt + 3), (sx * 11, 6.0, T + 5.5), (sx * 15, 6.0, T + 5.5)]) for sx in (-1, 1)])
    hinge = B.hull([(sx * 15, y, z) for sx in (-1, 1) for y in (-6.0, 8.0) for z in (T + 2.5, T + 6)]) + \
            B.hull([(sx * 13, y, z) for sx in (-1, 1) for y in (-6.0, -1.0) for z in (T - 30, T + 6)])
    mid = ((yb + 6) / 2)
    dia = B.hull([(0, mid - 22, zt), (0, mid + 22, zt), (10, mid, zt), (-10, mid, zt), (0, mid - 22, zt + 2.5), (0, mid + 22, zt + 2.5)])
    return bar + arms + hinge + dia, np.array([0, yb, zt + 2.5])

def build_controls():
    out = []
    for p, kind in KNOBS:
        r = 9.5 if kind == 'vol' else 8.0
        out.append(m3.Manifold.cylinder(17, r, r * 0.92, 40).translate([p[0], p[1], T]))
    out.append(m3.Manifold.cylinder(2.0, 7.0, 7.0, 24).translate([TOGGLE[0], TOGGLE[1], T]) + m3.Manifold.cylinder(18, 1.4, 1.4, 12).translate([TOGGLE[0], TOGGLE[1], T]) +
               m3.Manifold.sphere(3.2, 16).translate([TOGGLE[0], TOGGLE[1], T + 19]))
    return B.union(out)

def build_guard():
    return B.prism(B.ccw(GUARD), T, T + 2.3)

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
    br, wh = build_bridge(); parts['bridge'] = br; parts['thumbwheels'] = wh
    tp, anchor = build_tailpiece(); parts['tailpiece'] = tp
    parts['controls'] = build_controls(); parts['pickguard'] = build_guard()
    parts['strings'] = build_strings(anchor, posts)
    MAT = dict(faceplate='maple', body='mahogany', cap='maple', neck='maple', fingerboard='ebony', frets='nickel', nut='bone', tuners='nickel', badge='nickel',
               neck_pickup_ring='nickel', neck_pickup='nickel', neck_pickup_poles='nickel', bridge_pickup_ring='nickel', bridge_pickup='nickel',
               bridge_pickup_poles='nickel', bridge='ebony', thumbwheels='nickel', tailpiece='nickel', controls='black', pickguard='black', strings='steel')
    GROUP = {k: ('Body' if k in ('body', 'cap') else 'Neck' if k in ('neck', 'faceplate', 'fingerboard', 'frets', 'nut', 'badge') else
                 'Pickups' if 'pickup' in k else 'Strings' if k == 'strings' else 'Hardware') for k in parts}
    EXPLODE = dict(Body=(0, 0, 0), Neck=(0, 60, 40), Pickups=(0, 0, 50), Hardware=(0, 0, 80), Strings=(0, 0, 110))
    import trimesh
    scene = trimesh.Scene(); ok = {}; vols = {}
    for k, m in parts.items():
        tm = B.to_trimesh(m); ok[k] = bool(tm.is_watertight); vols[k] = m.volume()
        scene.add_geometry(tm, node_name=k, geom_name=k, metadata=dict(material=MAT[k], group=GROUP[k], explode=EXPLODE[GROUP[k]]))
    scene.export(os.path.join(OUT, 'juris_c75.glb'))
    mass = vols['body'] * 560e-6 + vols['cap'] * 700e-6 + vols['neck'] * 560e-6 + vols['fingerboard'] * 1150e-6 + 6 * 33 + 300 + 95 + 150
    nut_y = rot_pts([0, Y_N, Z_J], -TH, PIV)[1]
    joint_fret = 12 * math.log2(P['scale'] / (P['scale'] - (nut_y - Y_END)))
    report = dict(units='mm', params=P, body_length=float(SHARP[:, 1].max()), lower_bout=float(OUTL[:, 0].max() - OUTL[:, 0].min()),
                  bridge_height=float(BRIDGE_H), fb_above_body_at_joint=float(Z_J - T), neck_joint_fret=float(joint_fret), y_saddle=Y_SADDLE,
                  total_mass_g=round(mass), total_mass_lb=round(mass / 453.6, 2), manifold=ok)
    json.dump(report, open(os.path.join(OUT, 'report.json'), 'w'), indent=1)
    json.dump(dict(meta={k: dict(material=MAT[k], group=GROUP[k], explode=list(EXPLODE[GROUP[k]])) for k in parts}), open(os.path.join(OUT, 'meta.json'), 'w'))
    print(json.dumps({k: report[k] for k in ('body_length', 'bridge_height', 'fb_above_body_at_joint', 'neck_joint_fret', 'y_saddle', 'total_mass_lb')}))
    print('manifold', ok)

if __name__ == '__main__':
    main()
