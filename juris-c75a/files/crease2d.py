"""Plane-strain section across one fold of the top: a 3 mm plate bent through theta, sharp vs. a radius on the
re-entrant (inside) face. Clamp one arm, put a unit bending moment (per mm of fold length) on the other, and compare
the peak stress in the fold zone with the nominal plate bending stress 6M/t^2. Section plane = radial x thickness
(cross-grain), quartersawn: E_R, E_T, G_RT."""
import sys, math, json, numpy as np, gmsh
from skfem import MeshTri, ElementTriP2, ElementVector, Basis, FacetBasis, BilinearForm, LinearForm, asm, condense, solve
from skfem.helpers import sym_grad, dot

def section(theta_deg, R, t=3.0, L=30.0, h_fine=0.08, h=0.6):
    th = math.radians(theta_deg)
    # arm directions: left arm goes down-left, right arm down-right (ridge on top: the re-entrant side is underneath)
    dL = np.array([-math.cos(th / 2), -math.sin(th / 2)]); dR = np.array([math.cos(th / 2), -math.sin(th / 2)])
    nL = np.array([-dL[1], dL[0]]) * -1; nR = np.array([dR[1], -dR[0]]) * -1      # inward (downward-ish) normals
    nL = nL if nL[1] < 0 else -nL; nR = nR if nR[1] < 0 else -nR
    top = [dL * L, np.zeros(2), dR * L]                               # outer (convex) face, crisp
    # inner face: offset lines by t; their corner is the re-entrant one
    bl0, bl1 = dL * L + nL * t, nL * t; br0, br1 = nR * t, dR * L + nR * t
    # intersection of the two inner lines
    A = np.array([dL, -dR]).T; s_ = np.linalg.solve(A, br0 - bl1)
    corner = bl1 + dL * s_[0]
    gmsh.initialize(); gmsh.option.setNumber('General.Terminal', 0); gmsh.model.add('s')
    occ = gmsh.model.occ
    if R <= 0:
        pts = [top[0], top[1], top[2], dR * L + nR * t, corner, dL * L + nL * t]
    else:
        # fillet: circle of radius R tangent to both inner lines, on the material side (above the corner)
        half = (math.pi - th) / 2                                     # half the inner angle seen from the void... 
        bis = np.array([0.0, 1.0])                                    # bisector points up, into the material
        dist = R / math.sin(math.radians(90) - th / 2) if False else R / math.cos(th / 2)
        # the inner lines make angle (pi - th) at the corner, opening downward; circle sits below the corner (void side)
        ctr = corner - bis * (R / math.sin((math.pi - th) / 2))
        tL = corner + dL * (R / math.tan((math.pi - th) / 2)); tR = corner + dR * (R / math.tan((math.pi - th) / 2))
        arc = [ctr + R * np.array([math.cos(a), math.sin(a)]) for a in np.linspace(math.atan2(*(tR - ctr)[::-1]), math.atan2(*(tL - ctr)[::-1]), 60)]
        pts = [top[0], top[1], top[2], dR * L + nR * t] + arc + [dL * L + nL * t]
    P = [occ.addPoint(p[0], p[1], 0) for p in pts]
    Ls = [occ.addLine(P[i], P[(i + 1) % len(P)]) for i in range(len(P))]
    occ.addPlaneSurface([occ.addCurveLoop(Ls)]); occ.synchronize()
    f = gmsh.model.mesh.field
    d = f.add('Distance'); f.setNumbers(d, 'PointsList', [P[1]]); f.setNumber(d, 'Sampling', 200)
    th_ = f.add('Threshold'); f.setNumber(th_, 'InField', d); f.setNumber(th_, 'SizeMin', h_fine); f.setNumber(th_, 'SizeMax', h)
    f.setNumber(th_, 'DistMin', max(2.0, 0.6 * R * math.sin(th / 2) + 1.0)); f.setNumber(th_, 'DistMax', 12.0); f.setAsBackgroundMesh(th_)
    gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary', 0)
    gmsh.model.mesh.generate(2)
    nt, nc, _ = gmsh.model.mesh.getNodes(); X = nc.reshape(-1, 3)[:, :2]; idx = {g: i for i, g in enumerate(nt)}
    et, tg, cn = gmsh.model.mesh.getElements(2); T = np.array([idx[c] for c in cn[0]]).reshape(-1, 3)
    gmsh.finalize()
    u, inv = np.unique(T, return_inverse=True); X = X[u]; T = inv.reshape(-1, 3)
    return MeshTri(X.T.copy(), T.T.copy()), dL, dR, nL, nR, L, t, corner

def run(theta, R, Er=900.0 * 1.083, Et=500.0 * 1.083, Grt=35.0 * 1.083, nu=0.47, sign=1):
    m, dL, dR, nL, nR, L, t, corner = section(theta, R)
    # plane strain, orthotropic in the section (x = radial/across, y = thickness/tangential)
    S = np.array([[1 / Er, -nu / Er, 0], [-nu / Er, 1 / Et, 0], [0, 0, 1 / Grt]]); C = np.linalg.inv(S)
    e = ElementVector(ElementTriP2()); ib = Basis(m, e)
    @BilinearForm
    def a(u, v, w):
        eu, ev = sym_grad(u), sym_grad(v)
        E = [eu[0, 0], eu[1, 1], 2 * eu[0, 1]]; F = [ev[0, 0], ev[1, 1], 2 * ev[0, 1]]
        return sum(C[i, j] * E[j] * F[i] for i in range(3) for j in range(3))
    K = asm(a, ib)
    x = m.p
    sL = x.T @ dL; sR = x.T @ dR
    clamp = ib.get_dofs(lambda X: (X.T @ dL) > L - 1e-3).all()
    endR = m.facets_satisfying(lambda X: (X.T @ dR) > L - 1e-3)
    fb = FacetBasis(m, e, facets=endR)
    mid = dR * L + nR * t / 2
    I_ = t ** 3 / 12.0
    @LinearForm
    def load(v, w):        # linear axial traction = unit moment per mm
        X = w.x; zeta = (X[0] - mid[0]) * nR[0] + (X[1] - mid[1]) * nR[1]
        sig = sign * 1.0 * zeta / I_
        return sig * (v[0] * dR[0] + v[1] * dR[1])
    f = asm(load, fb)
    u = solve(*condense(K, f, D=clamp))
    du = ib.interpolate(u); eps = sym_grad(du)
    ev = np.array([eps[0, 0], eps[1, 1], 2 * eps[0, 1]]); sig = np.einsum('ij,jkl->ikl', C, ev)
    gx = ib.mapping.F(ib.quadrature[0])
    r = np.hypot(gx[0] - corner[0] * 0, gx[1] - 0)            # distance from the fold line (outer apex at origin)
    zone = np.hypot(gx[0], gx[1]) < 12.0
    s_tan = np.abs(sig[0])     # cross-grain normal stress (along the section)
    nominal = 6.0 / t ** 2
    far = (np.hypot(gx[0], gx[1]) > 18) & (np.hypot(gx[0], gx[1]) < 24)
    return dict(theta=theta, R=R, sign=sign, peak_ratio=float(s_tan[zone].max() / nominal),
                p999_ratio=float(np.percentile(s_tan[zone], 99.9) / nominal), far_ratio=float(s_tan[far].max() / nominal),
                shear_peak=float(np.abs(sig[2][zone]).max() / nominal), nodes=int(m.p.shape[1]))

if __name__ == '__main__':
    out = []
    for theta in (8.0, 15.0, 22.0):
        for R in (0.0, 10.0, 25.0):
            for sg in (1, -1):
                r = run(theta, R, sign=sg); out.append(r); print(r, flush=True)
    json.dump(out, open('out/crease2d.json', 'w'), indent=1)
