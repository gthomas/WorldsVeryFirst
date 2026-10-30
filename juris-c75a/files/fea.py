"""Static FEA of the Juris C-75A top under string load.
Quadratic tets (P2), orthotropic quartersawn spruce, top fixed where it is glued to the rim/linings.
Units: mm, N, MPa."""
import json, math, sys, time
import numpy as np, trimesh, tetgen
from skfem import MeshTet, ElementTetP2, ElementVector, Basis, FacetBasis, BilinearForm, LinearForm, asm, condense
from skfem.helpers import sym_grad
from shapely.geometry import Polygon, Point
import build as B

TENSION_N = float(sys.argv[1]) if len(sys.argv) > 1 else 850.0     # 13-56 flatwound set, approx.
MAXVOL = float(sys.argv[2]) if len(sys.argv) > 2 else 8.0

# ---- load
sad, nut = B.SADDLE, B.NUT_PT
rep = json.load(open('out/report.json'))
a_tail = math.radians(rep['string_break_angle_deg'])
a_neck = math.atan2(nut[2] - sad[2], nut[1] - sad[1])          # negative: strings fall toward the nut
F = TENSION_N * (math.sin(a_tail) + math.sin(-a_neck))

# ---- mesh
tm = trimesh.load('out/top.stl')
tg = tetgen.TetGen(tm.vertices, tm.faces)
_r = tg.tetrahedralize(order=1, mindihedral=12, minratio=1.6, maxvolume=MAXVOL); nodes, elems = _r[0], _r[1]
m = MeshTet(nodes.T.copy(), elems.T.copy())
e = ElementVector(ElementTetP2())
ib = Basis(m, e, intorder=3)
print(f'mesh: {m.nelements} tets, {ib.N} dofs', flush=True)

# ---- quartersawn spruce: x = radial (across the width), y = longitudinal (grain), z = tangential (thickness)
EL, ER, ET = 11000.0, 900.0, 500.0
GLR, GLT, GRT = 700.0, 600.0, 35.0
nLR, nLT, nRT = 0.37, 0.42, 0.47
S = np.zeros((6, 6))
S[0, 0], S[1, 1], S[2, 2] = 1 / ER, 1 / EL, 1 / ET
S[0, 1] = S[1, 0] = -nLR / EL
S[0, 2] = S[2, 0] = -nRT / ER
S[1, 2] = S[2, 1] = -nLT / EL
S[3, 3], S[4, 4], S[5, 5] = 1 / GLT, 1 / GRT, 1 / GLR      # yz, xz, xy
C = np.linalg.inv(S)

def voigt(E):
    return [E[0, 0], E[1, 1], E[2, 2], 2 * E[1, 2], 2 * E[0, 2], 2 * E[0, 1]]

@BilinearForm
def stiffness(u, v, w):
    eu, ev = voigt(sym_grad(u)), voigt(sym_grad(v))
    out = 0
    for i in range(6):
        for j in range(6):
            if C[i, j] != 0:
                out = out + C[i, j] * eu[j] * ev[i]
    return out

t0 = time.time()
K = asm(stiffness, ib)
print('assembled', round(time.time() - t0, 1), 's', flush=True)

# ---- supports: glued footprint on the rim + glue ledges
ledge_in = Polygon(B.OUTL).buffer(-(B.P['rim_wall'] + B.P['lining_w']))
# plus the neck block, which the top also sits on
nb = Polygon([(-44, 392), (44, 392), (44, 470), (-44, 470)]).intersection(Polygon(B.OUTL))
def glued(x):
    return np.array([(abs(z - B.Z_RIM_TOP) < 0.05) and ((not ledge_in.contains(Point(px, py))) or nb.contains(Point(px, py)))
                     for px, py, z in x.T])
fixed = ib.get_dofs(lambda x: glued(x)).all()
print('fixed dofs', len(fixed), flush=True)

# ---- bridge feet on the plateau
feet = [Polygon([(sx * 22, B.P['bridge_y'] - 5.0), (sx * 67.5, B.P['bridge_y'] - 3.0),
                 (sx * 67.5, B.P['bridge_y'] + 3.0), (sx * 22, B.P['bridge_y'] + 5.0)]) for sx in (-1, 1)]
foot_area = sum(f.area for f in feet)
p = F / foot_area
def on_feet(x):
    return np.array([(abs(z - B.Z_PLAT) < 0.05) and any(f.contains(Point(px, py)) for f in feet) for px, py, z in x.T])
facets = m.facets_satisfying(lambda x: on_feet(x), boundaries_only=True)
fb = FacetBasis(m, e, facets=facets, intorder=3)

@LinearForm
def unit(v, w):
    return -v[2]
f1 = asm(unit, fb)
p = F / (-f1.sum())            # loaded facets overhang the foot outline a little; scale so the total is exactly F
f = p * f1
applied = -f.sum()   # check
print(f'downforce {F:.1f} N, feet area {foot_area:.0f} mm2, assembled load {applied:.1f} N', flush=True)

# ---- solve
import pyamg
from scipy.sparse.linalg import cg
Kc, fc, xI, I = condense(K, f, D=fixed)
# rigid-body modes as near-nullspace for the AMG
X_ = ib.doflocs[:, I] if False else None
dl = ib.doflocs
comp = np.zeros(ib.N, int)
for c_ in range(3): comp[ib.nodal_dofs[c_]] = c_
for c_ in range(3):
    if hasattr(ib, 'edge_dofs') and ib.edge_dofs.size: comp[ib.edge_dofs[c_]] = c_
xyz = dl
Bm = np.zeros((ib.N, 6))
for c_ in range(3): Bm[comp == c_, c_] = 1.0
x_, y_, z_ = xyz
Bm[comp == 0, 3] = -y_[comp == 0]; Bm[comp == 1, 3] = x_[comp == 1]
Bm[comp == 1, 4] = -z_[comp == 1]; Bm[comp == 2, 4] = y_[comp == 2]
Bm[comp == 2, 5] = -x_[comp == 2]; Bm[comp == 0, 5] = z_[comp == 0]
import pypardiso
t0 = time.time()
uI = pypardiso.spsolve(Kc.tocsr(), fc)
r = np.linalg.norm(Kc @ uI - fc) / np.linalg.norm(fc)
print('solved (pardiso)', round(time.time() - t0, 1), 's, rel residual', r, flush=True)
u = np.zeros(K.shape[0]); u[I] = uI

# ---- results
U = u[ib.nodal_dofs].T if hasattr(ib, 'nodal_dofs') else None
uz = u[ib.nodal_dofs[2]]; ux = u[ib.nodal_dofs[0]]; uy = u[ib.nodal_dofs[1]]
X = m.p
bridge_nodes = on_feet(X)
out = dict(tension_N=TENSION_N, downforce_N=F, tets=int(m.nelements), dofs=int(ib.N),
           max_down_mm=float(-uz.min()), bridge_down_mm=float(-uz[bridge_nodes].mean()),
           where_max=X[:, np.argmin(uz)].round(1).tolist())

# stresses at quadrature points
du = ib.interpolate(u)
eps = sym_grad(du)
ev = np.array(voigt(eps))                    # 6 x nelem x nqp
sig = np.einsum('ij,jkl->ikl', C, ev)
def pct(a, q): return float(np.percentile(np.abs(a), q))
names = ['xx_cross', 'yy_grain', 'zz_thick', 'yz_LT', 'xz_RT', 'xy_LR']
out['stress_MPa_p99'] = {n: pct(sig[i], 99.0) for i, n in enumerate(names)}
out['stress_MPa_max'] = {n: float(np.abs(sig[i]).max()) for i, n in enumerate(names)}
# allowable (short-term, clear spruce, approx.)
allow = dict(xx_cross=4.0, yy_grain=35.0, zz_thick=4.0, yz_LT=6.0, xz_RT=1.5, xy_LR=6.0)
out['safety_factor_p99'] = {n: allow[n] / max(out['stress_MPa_p99'][n], 1e-9) for n in names}
json.dump(out, open('out/fea.json', 'w'), indent=1)
print(json.dumps(out, indent=1))

# deflection field for plotting: top-surface nodes
top_mask = np.array([z > B.Z_RIM_TOP + 2.5 or (abs(z - B.Z_EDGE) < 0.05) for z in X[2]])
np.savez('out/fea_field.npz', x=X[0], y=X[1], z=X[2], uz=uz)
