"""Headstock of the C-75 from Chuck's model (WyreCustomGuitars, Juris.dxf), measured in the plane of the
headstock face: concave shoulders flaring out from the neck, waisted sides, wide top corners and a gabled tip;
3+3 tuners, a pointed truss rod cover with three screws over the access slot. s = distance from the nut line
along the head (Chuck's outline starts at the back of the 4.5 mm nut), x = across. mm.
The previous outline, traced from a photo, is kept in head_c75_photo.py."""
import numpy as np

# right half from the nut to the tip; the left half mirrors it
_RIGHT = [(21.4, 0.0), (21.4, 4.5), (21.43, 6.0), (21.49, 7.5), (21.64, 9.0), (21.87, 10.5), (22.2, 12.0), (22.62, 13.5), (23.13, 15.0), (23.75, 16.5), (24.48, 18.0), (25.34, 19.5), (26.33, 21.0), (27.48, 22.5), (28.82, 24.0), (30.38, 25.5), (32.25, 27.0), (34.53, 28.5), (37.28, 30.0), (38.21, 31.5), (38.05, 33.0), (37.66, 34.5), (36.41, 39.5), (35.32, 44.5), (34.36, 49.5), (33.53, 54.5), (32.83, 59.5), (32.26, 64.5), (31.82, 69.5), (31.52, 74.5), (31.33, 79.5), (31.28, 84.5), (31.35, 89.5), (31.55, 94.5), (31.88, 99.5), (32.33, 104.5), (32.92, 109.5), (33.63, 114.5), (34.48, 119.5), (35.46, 124.5), (36.58, 129.5), (37.84, 134.5), (39.23, 139.5), (40.68, 144.5), (40.62, 146.5), (38.14, 148.5), (34.89, 150.5), (31.64, 152.5), (28.39, 154.5), (25.14, 156.5), (21.89, 158.5), (18.65, 160.5), (15.4, 162.5), (12.15, 164.5), (8.9, 166.5), (5.65, 168.5), (1.42, 171.1), (0.0, 171.6)]
HEAD_OUTLINE = np.array(_RIGHT + [(-x, s) for x, s in _RIGHT[-2::-1]], float)     # CCW: up the treble side, down the bass side
TUNER_POS = [(21.9, 48.0), (18.6, 83.6), (21.8, 119.2)]                             # (|x|, s), nearest the nut first
TUNER_LIST = [(sv * x, s, sv) for x, s in TUNER_POS for sv in (-1, 1)]
TRUSS_SLOT = np.array([(-4.9, 8.0), (4.9, 8.0), (4.9, 31.0), (-4.9, 31.0)])
TRUSS_COVER = np.array([(-13.37, 4.5), (13.36, 4.5), (9.22, 45.53), (0.0, 54.5), (-9.23, 45.53)])
TRUSS_SCREWS = [(-9.5, 7.5), (9.5, 7.5), (0.0, 49.0)]
TUNER_HOLE = 10.0                                                                  # as drawn (10 mm bushings)
