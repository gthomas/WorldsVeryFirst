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
import os
MESH = os.environ.get('FEA_MESH', 'out/top_fea.stl')
TAG = os.environ.get('FEA_TAG', '')
OUTJ = f'out/fea_{TAG}.json' if TAG else 'out/fea.json'
NMODES = int(os.environ.get('FEA_MODES', '4'))
RHO = 430e-12          # spruce, tonne/mm^3 (N, mm, s units)
BRIDGE_G = 37.0        # base + top, carried on the feet

# ---- load
sad, nut = B.SADDLE, B.NUT_PT
rep = json.load(open('out/report.json'))
a_tail = math.radians(rep['string_break_angle_deg'])
a_neck = math.atan2(nut[2] - sad[2], nut[1] - sad[1])          # negative: strings fall toward the nut
F = TENSION_N * (math.sin(a_tail) + math.sin(-a_neck))

# ---- mesh
tm = trimesh.load(MESH)   # top without the soundhole lips: they only add stiffness, and tetgen can't mesh their joint
tg = tetgen.TetGen(tm.vertices, tm.faces)
_r = tg.tetrahedralize(order=1, mindihedral=12, minratio=float(os.environ.get('FEA_MINRATIO', '2.0')), maxvolume=MAXVOL); nodes, elems = _r[0], _r[1]
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
# Load every plateau-top facet by the exact area it shares with the feet, so the load does not
# depend on how tetgen happened to triangulate the plateau. For a P2 triangle under uniform
# pressure the corner nodes carry nothing and each mid-edge node carries a third.
from shapely.geometry import Polygon as _Pg
from shapely.ops import unary_union
feet_u = unary_union(feet)
bf = m.boundary_facets()
fz = m.p[2][m.facets[:, bf]]
plat = bf[np.all(np.abs(fz - B.Z_PLAT) < 0.05, axis=0)]
f = np.zeros(ib.N); loaded_area = 0.0
for fi in plat:
    tri = m.p[:2, m.facets[:, fi]].T
    a = _Pg(tri).intersection(feet_u).area
    if a <= 0: continue
    loaded_area += a
    for ed in m.f2e[:, fi]:
        f[ib.edge_dofs[2, ed]] -= p * a / 3.0
applied = -f.sum()
print(f'downforce {F:.1f} N, feet area {foot_area:.0f} mm2, loaded {loaded_area:.0f} mm2, assembled load {applied:.1f} N', flush=True)

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
ONLY_MODES = os.environ.get("FEA_ONLY_MODES") == "1"
del K; import gc; gc.collect()
if ONLY_MODES:
    out = json.load(open(OUTJ)) if os.path.exists(OUTJ) else {}
else:
    import pypardiso
    t0 = time.time()
    uI = pypardiso.spsolve(Kc.tocsr(), fc)
    r = np.linalg.norm(Kc @ uI - fc) / np.linalg.norm(fc)
    print('solved (pardiso)', round(time.time() - t0, 1), 's, rel residual', r, flush=True)
    u = np.zeros(ib.N); u[I] = uI

    # ---- results
    U = u[ib.nodal_dofs].T if hasattr(ib, 'nodal_dofs') else None
    uz = u[ib.nodal_dofs[2]]; ux = u[ib.nodal_dofs[0]]; uy = u[ib.nodal_dofs[1]]
    X = m.p
    bridge_nodes = on_feet(X)
    out = dict(tension_N=TENSION_N, downforce_N=F, tets=int(m.nelements), dofs=int(ib.N),
               max_down_mm=float(-uz.min()), bridge_down_mm=float((f @ u) / max(-f.sum(), 1e-9)),
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
    # ---- stress at the fold lines (creases between facets) vs. away from them
    from shapely.geometry import LineString, MultiLineString
    segs = [(B.TOP_HEX[i], B.TOP_HEX[(i + 1) % len(B.TOP_HEX)]) for i in range(len(B.TOP_HEX))]
    segs += [(B.ARCH[i], B.ARCH[(i + 1) % len(B.ARCH)]) for i in range(len(B.ARCH))]
    for t_ in B.TOP_TRIS:
        for a_, b_ in [(t_[0], t_[1]), (t_[1], t_[2]), (t_[2], t_[0])]:
            if a_[0] != b_[0]:
                pa = B.ARCH[a_[1]] if a_[0] == 'o' else B.TOP_HEX[a_[1]]; pb = B.ARCH[b_[1]] if b_[0] == 'o' else B.TOP_HEX[b_[1]]
                segs.append((pa, pb))
    folds = MultiLineString([[tuple(p), tuple(q)] for p, q in segs])
    util = np.max(np.abs(sig) / np.array([allow[n] for n in names])[:, None, None], axis=0)     # nelem x nqp
    gx = ib.mapping.F(ib.quadrature[0])                                                          # 3 x nelem x nqp
    shell_q = gx[2] > B.Z_RIM_TOP + 0.5                     # leave out the glued skirt seat
    xs, ys, us = gx[0][shell_q], gx[1][shell_q], util[shell_q]
    import shapely
    dist = shapely.distance(folds, shapely.points(np.c_[xs, ys]))
    near = dist < 4.0
    hi = us >= np.percentile(us, 99.9)
    out['folds'] = dict(util_p99_near=float(np.percentile(us[near], 99)), util_p99_away=float(np.percentile(us[~near], 99)),
                        util_max_near=float(us[near].max()), util_max_away=float(us[~near].max()),
                        share_of_area_near=float(near.mean()), share_of_top01pct_near=float(near[hi].mean()))
    print('folds', {k: round(v, 3) for k, v in out['folds'].items()}, flush=True)
    # ---- natural frequencies of the top as glued to a rigid rim, carrying the bridge
if NMODES and (ONLY_MODES or os.environ.get('FEA_MODES_INLINE') == '1'):
    from skfem import BilinearForm as _BF
    from skfem.helpers import dot as _dot
    from scipy.sparse.linalg import eigsh, LinearOperator
    import scipy.sparse as sp
    @_BF
    def _mass(u, v, w): return RHO * _dot(u, v)
    M = asm(_mass, ib)
    # bridge mass, spread over the feet in proportion to contact area (all three directions)
    mb = np.zeros(ib.N)
    for fi in plat:
        tri = m.p[:2, m.facets[:, fi]].T
        a_ = _Pg(tri).intersection(feet_u).area
        if a_ <= 0: continue
        for ed in m.f2e[:, fi]:
            for c_ in range(3): mb[ib.edge_dofs[c_, ed]] += BRIDGE_G * 1e-6 * (a_ / foot_area) / 3.0
    M = M + sp.diags(mb)
    import gc
    try: pypardiso.ps.free_memory(everything=True)     # drop the static factorization (memory)
    except Exception: pass
    Mc = M[I][:, I].tocsr(); Kcc = Kc.tocsr()
    del M; import gc; gc.collect()
    sigma = 0.0     # K is positive definite once clamped, so shift-invert about zero finds the lowest modes
    A = sp.triu(Kcc, format='csr'); del Kc; gc.collect()
    solver = pypardiso.PyPardisoSolver(mtype=2); solver.factorize(A)      # symmetric positive definite: half the memory
    OP = LinearOperator(Kcc.shape, matvec=lambda x: solver.solve(A, x), dtype=float)
    t0 = time.time()
    vals, vecs = eigsh(Kcc, k=NMODES, M=Mc, sigma=sigma, which='LM', OPinv=OP)
    order = np.argsort(vals); vals = vals[order]; vecs = vecs[:, order]
    fr = np.sqrt(np.maximum(vals, 0)) / (2 * np.pi)
    modes = []
    for k in range(NMODES):
        ph = np.zeros(ib.N); ph[I] = vecs[:, k]
        z_ = ph[ib.nodal_dofs[2]]; s_ = np.sign(z_[np.argmax(np.abs(z_))])
        # shape label: how much of the plate moves in phase (monopole) vs. split (dipole)
        top_n = np.abs(m.p[2] - B.Z_PLAT) < 30
        zz = z_[top_n] * s_
        lr = np.sign(m.p[0][top_n]); ud = np.sign(m.p[1][top_n] - B.P['bridge_y'])
        sc = dict(mono=abs(zz.sum()), lr=abs((zz * lr).sum()), ud=abs((zz * ud).sum()))
        modes.append(dict(hz=float(fr[k]), shape=max(sc, key=sc.get)))
    out['modes'] = modes
    print('modes', [(round(md['hz']), md['shape']) for md in modes], round(time.time() - t0, 1), 's', flush=True)
json.dump(out, open(OUTJ, 'w'), indent=1)
print(json.dumps(out, indent=1))

# deflection field for plotting: top-surface nodes
if not ONLY_MODES:
    np.savez(OUTJ.replace('.json', '_field.npz'), x=X[0], y=X[1], z=X[2], uz=uz)
