"""Smooth neck solid for the C-75 and C-75T, built in a straight neck frame: y along the neck, z up, the
fingerboard gluing plane at z = zgp. The shaft carve flows into the headstock as a series of vertical slices
(round '59-style section -> square headstock section), with the back blended into the head by a soft minimum
and an optional volute, the same way the C-75A's volute is made. The heel end flares wider and deeper as the
neck meets the body instead of ending in a block. Returns manifolds in the straight frame; callers tilt them."""
import math
import numpy as np
import manifold3d as m3
import build as B
import head_c75 as H

_HR = np.array(H._RIGHT[:-1])
def head_half(s): return float(np.interp(s, _HR[:, 1], _HR[:, 0]))

def _smin(a, b, k):
    return -k * math.log(math.exp(-a / k) + math.exp(-b / k))

def _section(y, half, top, bot, p, n=33):
    """convex slice: flat top at `top`, back curve down to `bot`; p -> 0 squares the section off"""
    D = top - bot
    return [(half * math.cos(a), y, top - D * (math.sin(a) ** p)) for a in np.linspace(0, math.pi, n)]

def build(Y_N, Y_body, zgp, width, depth, head_angle_deg=0.0, head_T=14.0, volute=0.0,
          heel_extra=0.0, heel_w=None, heel_len=45.0, heel_p=0.35, s_plate=32.0):
    """Y_N: nut line. Y_body: where the neck meets the body (the slices stop just inside it).
    width(y), depth(y): shaft width and depth below the gluing plane. heel_extra/heel_w: how much deeper/wider the
    neck gets over the last `heel_len` mm before the body; heel_p: how square the section gets there."""
    HA = math.radians(head_angle_deg); tH = math.tan(HA); cH = math.cos(HA)
    face = lambda s: zgp if s <= 0 else zgp - s * tH
    head_back = lambda s: face(s) - head_T / cH
    rings = []
    ys = np.r_[np.linspace(Y_body - 1.0, Y_body + heel_len, 22), np.linspace(Y_body + heel_len + 6, Y_N - 40, 30),
               np.linspace(Y_N - 38, Y_N + s_plate + 3, 46)]
    for y in ys:
        s = y - Y_N
        # shaft
        yy = min(y, Y_N)
        half = width(yy) / 2; d = depth(yy); p = 0.85
        t = float(np.clip((Y_body + heel_len - y) / heel_len, 0, 1)); t = t * t * (3 - 2 * t)      # 1 at the body
        if t > 0:
            hw = (heel_w or width(Y_body)) / 2
            half = half + (hw - half) * t; d = d + heel_extra * t; p = 0.85 + (heel_p - 0.85) * t
        top = face(s); bot = zgp - d
        if s > -30:
            k = 0.05 + 6.0 * math.exp(-(s / 20.0) ** 2)                # soft blend only around the break
            bot = _smin(bot, head_back(s), k) if s < s_plate + 3 else head_back(s) + 0.05
            bot -= volute * math.exp(-((s - 6.0) / 13.0) ** 2)
            if s > 0:
                half = head_half(s) - 0.02
            u = float(np.clip(s / s_plate, 0, 1)); u = u * u * (3 - 2 * u)
            p = p * (1 - u) + 0.04 * u                                    # round neck -> square headstock section
        rings.append(np.array(_section(y, half, top - 0.01, bot, p)))
    body = B.union([B.hull(list(a) + list(b)) for a, b in zip(rings[:-1], rings[1:])])
    # headstock plate from s_plate on: the exact outline, tilted about the nut
    from shapely.geometry import Polygon, box
    poly = Polygon(H.HEAD_OUTLINE).intersection(box(-200, s_plate - 2.0, 200, 400))
    pp = B.ccw(np.array(poly.exterior.coords)[:-1] + np.array([0, Y_N]))
    piv = [0.0, Y_N, zgp]
    tilt = lambda m: m.translate([-piv[0], -piv[1], -piv[2]]).rotate([-head_angle_deg, 0, 0]).translate(piv)
    plate = tilt(B.prism(pp, zgp - head_T, zgp))
    neck = body + plate
    # tuner holes, truss rod slot
    holes = B.union([tilt(m3.Manifold.cylinder(60, H.TUNER_HOLE / 2, H.TUNER_HOLE / 2, 28).translate([v, Y_N + s, zgp - 40]))
                     for v, s, _ in H.TUNER_LIST])
    slot = tilt(B.prism(H.TRUSS_SLOT + np.array([0, Y_N]), zgp - 7, zgp + 2))
    # headstock face skin (finished to match the top): the face from just past the nut to the tip
    fpoly = Polygon(H.HEAD_OUTLINE).intersection(box(-200, 3.0, 200, 400)).buffer(-0.3)
    fp = B.ccw(np.array(fpoly.exterior.coords)[:-1] + np.array([0, Y_N]))
    face_skin = tilt(B.prism(fp, zgp - 0.01, zgp + 0.6)) - holes - slot
    return neck - holes - slot, face_skin, tilt
