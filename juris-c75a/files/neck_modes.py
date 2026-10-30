"""Natural frequencies of the neck + fingerboard + headstock face, with the six tuners as added mass.
The neck is held rigid where it sits in the body (tenon, tongue and heel face), so these are upper bounds:
the real body flexes and pulls them down somewhat. Quadratic tets, orthotropic woods, grain along the neck.
Output: out/neck_modes.json, with the nearest note to each mode and where that note falls on the neck.
Units: mm, N, tonne, s."""
import json, math, os, time
import numpy as np, trimesh, tetgen, scipy.sparse as sp
from skfem import MeshTet, ElementTetP2, ElementVector, Basis, BilinearForm, asm, condense
from skfem.helpers import sym_grad, dot
import build as B

NMODES = int(os.environ.get('NECK_MODES', '8'))
TUNER_G = 33.0                         # Grover 133N, each, with bushing and screw

# ---- geometry: one solid, element material picked by which part its centroid falls in
B.fret_slots = lambda *a, **k: B.m3.Manifold()          # frets fill their slots
neck = B.build_neck() + B.truss_channel()                # rod fills its channel; tuner posts fill their holes
fb = B.build_fingerboard()
face = B.build_faceplate()
whole = neck + fb + face
tm = B.to_trimesh(whole)
tm.export('out/neck_fea.stl')
# the lofted volute and heel defeat the tet mesher, so build the solid from 2 mm cubes (trilinear hexes)
PITCH = float(os.environ.get('NECK_PITCH', '2.0'))
vox = tm.voxelized(PITCH).fill()
idx = np.argwhere(vox.matrix)
origin = vox.transform[:3, 3] - PITCH / 2            # corner of voxel (0,0,0)
corners = np.array([(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)])
allc = (idx[:, None, :] + corners[None]).reshape(-1, 3)
uniq, inv = np.unique(allc, axis=0, return_inverse=True)
from skfem import MeshHex, ElementHex1
m = MeshHex((uniq * PITCH + origin).T.copy(), inv.reshape(-1, 8).T.copy())
e = ElementVector(ElementHex1())
ib = Basis(m, e, intorder=2)
print(f'mesh: {m.nelements} hexes at {PITCH} mm, {ib.N} dofs', flush=True)
cent = m.p[:, m.t].mean(axis=1).T
fb_tm, face_tm = B.to_trimesh(fb), B.to_trimesh(face)
is_eb = fb_tm.contains(cent) | face_tm.contains(cent)

# neck frame axes in body coordinates: L along the neck, R across (x), T through the neck depth
L = B.n2b(1, 0, 0) - B.n2b(0, 0, 0); L /= np.linalg.norm(L)
Rx = np.array([1.0, 0, 0]); Tn = np.cross(L, Rx); Tn /= np.linalg.norm(Tn)
Q = np.array([Rx, L, Tn])               # rows: local axes (R, L, T) in body coords

def ortho_C(EL, ER, ET, GLR, GLT, GRT, nLR, nLT, nRT):
    """stiffness in the local (R, L, T) frame, Voigt order xx yy zz yz xz xy with x=R, y=L, z=T"""
    S = np.zeros((6, 6))
    S[0, 0], S[1, 1], S[2, 2] = 1 / ER, 1 / EL, 1 / ET
    S[0, 1] = S[1, 0] = -nLR / EL; S[0, 2] = S[2, 0] = -nRT / ER; S[1, 2] = S[2, 1] = -nLT / EL
    S[3, 3], S[4, 4], S[5, 5] = 1 / GLT, 1 / GRT, 1 / GLR
    return np.linalg.inv(S)

def rotate_C(C, Q):
    """rotate a Voigt stiffness from local to body axes (Bond transformation)"""
    a = Q.T                               # body = a @ local
    M = np.zeros((6, 6))
    pairs = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
    for I, (i, j) in enumerate(pairs):
        for J, (k, l) in enumerate(pairs):
            if J < 3: M[I, J] = a[i, k] * a[j, l]
            else: M[I, J] = a[i, k] * a[j, l] + a[i, l] * a[j, k]
    return M @ C @ M.T

# Honduran mahogany and ebony, typical published values (MPa); density tonne/mm^3
C_mah = rotate_C(ortho_C(10100, 1080, 660, 880, 630, 210, 0.31, 0.39, 0.60), Q)
C_ebo = rotate_C(ortho_C(17000, 2500, 1600, 1500, 1300, 500, 0.35, 0.40, 0.55), Q)
RHO_MAH, RHO_EBO = 590e-12, 1150e-12

def voigt(E):
    return [E[0, 0], E[1, 1], E[2, 2], 2 * E[1, 2], 2 * E[0, 2], 2 * E[0, 1]]

eb_field = ib.zeros()          # placeholder to get shapes
eb_el = is_eb.astype(float)

@BilinearForm
def stiffness(u, v, w):
    eu, ev = voigt(sym_grad(u)), voigt(sym_grad(v))
    k = w['eb']
    out = 0
    for i in range(6):
        for j in range(6):
            c = C_mah[i, j] * (1 - k) + C_ebo[i, j] * k
            if np.any(c != 0):
                out = out + c * eu[j] * ev[i]
    return out

@BilinearForm
def mass(u, v, w):
    k = w['eb']
    return (RHO_MAH * (1 - k) + RHO_EBO * k) * dot(u, v)

ebq = np.repeat(eb_el[:, None], ib.X.shape[-1], axis=1)
t0 = time.time()
K = asm(stiffness, ib, eb=ebq)
M = asm(mass, ib, eb=ebq)
print('assembled', round(time.time() - t0, 1), 's', flush=True)

# tuners: lumped mass on the dofs round each post hole
X = ib.doflocs
mt = np.zeros(ib.N)
for v, s, _ in B.TUNER_LIST:
    a = B.h2b(v, s, 0.0); b = B.h2b(v, s, -10.0); ax = (b - a) / np.linalg.norm(b - a)
    d = X.T - a; r = np.linalg.norm(d - np.outer(d @ ax, ax), axis=1)
    sel = np.where(r < B.P['tuner_hole'] / 2 + 2.0)[0]
    sel_nodes = np.unique(sel // 3 if False else sel)
    mt[sel] += TUNER_G * 1e-6 / max(len(sel) / 3, 1)
M = M + sp.diags(mt)
print('tuner mass on', int((mt > 0).sum()), 'dofs, total', round(mt.sum() / 3 * 1e6, 1), 'g', flush=True)

# ---- held rigid where it sits in the body: everything on the body side of the neck-end plane
fixed = np.where(ib.doflocs[1] < B.Y_NECK_END + 0.05)[0]
I = np.setdiff1d(np.arange(ib.N), fixed)
Kc = K[I][:, I].tocsr(); Mc = M[I][:, I].tocsr()
del K, M
import pypardiso
from scipy.sparse.linalg import eigsh, LinearOperator
A = sp.triu(Kc, format='csr')
solver = pypardiso.PyPardisoSolver(mtype=2); solver.factorize(A)
OP = LinearOperator(Kc.shape, matvec=lambda x: solver.solve(A, x), dtype=float)
t0 = time.time()
vals, vecs = eigsh(Kc, k=NMODES, M=Mc, sigma=0.0, which='LM', OPinv=OP)
o = np.argsort(vals); vals, vecs = vals[o], vecs[:, o]
fr = np.sqrt(np.maximum(vals, 0)) / (2 * np.pi)
print('solved', round(time.time() - t0, 1), 's', flush=True)

# ---- classify each mode by how the headstock tip moves, and find nearby notes on the neck
NAMES = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']
OPEN = {'E2': 40, 'A2': 45, 'D3': 50, 'G3': 55, 'B3': 59, 'E4': 64}       # MIDI numbers
def note(f):
    mid = 69 + 12 * math.log2(f / 440.0); n = round(mid)
    return f'{NAMES[n % 12]}{n // 12 - 1}', n, (mid - n) * 100
tip = B.h2b(0, 190, -B.HEAD_T / 2)
tip_dofs = np.argsort(np.linalg.norm(X.T - tip, axis=1))[:60]
modes = []
for k in range(NMODES):
    ph = np.zeros(ib.N); ph[I] = vecs[:, k]
    comp = np.array([np.abs(ph[tip_dofs[(tip_dofs % 3) == c]]).mean() for c in range(3)]) if False else None
    # displacement components at headstock-tip dofs (vector dofs are interleaved x,y,z per location)
    ux = ph[ib.nodal_dofs[0]]; uy = ph[ib.nodal_dofs[1]]; uz = ph[ib.nodal_dofs[2]]
    tipn = np.argsort(np.linalg.norm(m.p.T - tip, axis=1))[:20]
    side = np.abs(ux[tipn]).mean(); vert = np.abs(np.c_[uy[tipn], uz[tipn]] @ Tn[1:]).mean()
    # torsion: opposite vertical motion on the two edges of the headstock
    l_ = B.h2b(-38, 150, -7); r_ = B.h2b(38, 150, -7)
    ln = np.argsort(np.linalg.norm(m.p.T - l_, axis=1))[:10]; rn = np.argsort(np.linalg.norm(m.p.T - r_, axis=1))[:10]
    tw = abs((np.c_[uy[ln], uz[ln]] @ Tn[1:]).mean() - (np.c_[uy[rn], uz[rn]] @ Tn[1:]).mean()) / 2
    kind = max([('up-and-down bend', vert), ('side-to-side bend', side), ('twist', tw)], key=lambda t: t[1])[0]
    nm, n, cents = note(fr[k])
    where = [f'{s} fret {n - o}' for s, o in OPEN.items() if 0 <= n - o <= 15]
    modes.append(dict(hz=float(fr[k]), kind=kind, note=nm, cents=float(cents), where=where))
    print(f'{fr[k]:7.1f} Hz  {kind:18s}  {nm:4s} {cents:+5.0f}c  {", ".join(where)}', flush=True)
json.dump(dict(tuner_g=TUNER_G, modes=modes, note='neck held rigid at the body; real values somewhat lower'),
          open('out/neck_modes.json', 'w'), indent=1)
