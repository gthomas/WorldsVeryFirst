"""Headstock of the built C-75 (from the photo of the unfinished guitar, scaled from the 13.5" lower bout):
waisted sides, a gabled top, 3+3 tuner holes in near-straight columns (middle pair pulled in with the waist),
truss rod access slot just above the nut. s = distance from the nut line along the head, x = across. mm."""
import numpy as np

# right half from the nut to the tip; the left half mirrors it
_RIGHT = [(21.4, 0.0), (28.0, 3.0), (36.0, 6.5), (43.0, 12.0), (41.8, 30.0), (39.6, 55.0), (37.8, 80.0), (37.6, 95.0),
          (38.6, 115.0), (41.2, 140.0), (45.0, 163.0), (0.0, 190.0)]
HEAD_OUTLINE = np.array(_RIGHT + [(-x, s) for x, s in _RIGHT[-2::-1]], float)     # CCW: up the treble side, down the bass side
TUNER_POS = [(24.0, 52.4), (20.7, 93.0), (24.4, 133.7)]                             # (|x|, s), nearest the nut first
TUNER_LIST = [(sv * x, s, sv) for x, s in TUNER_POS for sv in (-1, 1)]
TRUSS_SLOT = np.array([(-4.9, 6.0), (4.9, 6.0), (4.9, 29.0), (-4.9, 29.0)])
TUNER_HOLE = 8.73                                                                  # 11/32" for Grover 133N
