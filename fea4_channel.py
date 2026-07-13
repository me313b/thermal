"""fea4_channel.py - duct coefficients for the plate-channel architecture.

Fully developed laminar analysis of the true lens-shaped subchannel (cell
crown vs flat plate over one pitch), by the influence-coefficient method
for doubly connected ducts (Shah & London):

    [dT_cell-bulk ]   [Rcc  Rcp] [q'_cell ]
    [dT_plate-bulk] = [Rpc  Rpp] [q'_plate]      K per (W/m), per pitch

plus fRe for the hydraulics. Velocity: lap w = -1 with no-slip on cell
and plate, symmetry at the pitch cuts and the row centreline. Energy:
k lap T = (w/wbar) q'_tot/A with uniform-flux (H2) Neumann walls,
flow-weighted mean pinned to zero; reported wall temperatures are
extrapolated the half cell to the true wall. The pass-through v1 study
(stagnant slot) gave lambda ~ 0.13 and is superseded: oil is the primary
carrier, plates the secondary sink.

Validation gates (parallel-plate limit): fRe = 96, one-side-heated
H2 Nu = 5.385.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def _assemble(oil, dx, dy, masked="dirichlet", top="none",
              bottom="none"):
    """Laplacian on the oil mask. masked: BC where a masked (cell)
    neighbour is hit. top/bottom: BC at the j = m-1 / j = 0 domain edge.
    'dirichlet' -> value 0 via ghost; 'none' -> symmetry/Neumann."""
    n, m = oil.shape
    idx = -np.ones((n, m), int)
    idx[oil] = np.arange(oil.sum())
    N = int(oil.sum())
    rows, cols, vals = [], [], []
    for i in range(n):
        for j in range(m):
            if not oil[i, j]:
                continue
            p = idx[i, j]
            diag = 0.0
            for di, dj, h in ((1, 0, dx), (-1, 0, dx),
                              (0, 1, dy), (0, -1, dy)):
                ii, jj = i + di, j + dj
                w = 1.0 / h ** 2
                if 0 <= ii < n and 0 <= jj < m:
                    if oil[ii, jj]:
                        rows.append(p); cols.append(idx[ii, jj])
                        vals.append(w); diag -= w
                    elif masked == "dirichlet":
                        diag -= 2 * w
                else:
                    edge = (top if jj >= m else
                            bottom if jj < 0 else "none")
                    if edge == "dirichlet":
                        diag -= 2 * w
            rows.append(p); cols.append(p); vals.append(diag)
    return sp.csr_matrix((vals, (rows, cols)), shape=(N, N)), idx


def _cell_faces(oil, dx, dy):
    n, m = oil.shape
    out = []
    for i in range(n):
        for j in range(m):
            if not oil[i, j]:
                continue
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                ii, jj = i + di, j + dj
                if 0 <= ii < n and 0 <= jj < m and not oil[ii, jj]:
                    out.append((i, j, dy if di else dx,
                                dx if di else dy))
                    # (face length, normal spacing)
    return out


def duct_coefficients(D=0.021, pitch_x=0.0215, s=0.002, k_oil=0.13,
                      n=160, rect=False):
    Lx = pitch_x
    Ly = s if rect else D / 2 + s
    nx = n
    ny = max(int(round(n * Ly / Lx)), 24)
    dx, dy = Lx / nx, Ly / ny
    x = (np.arange(nx) + 0.5) * dx
    y = (np.arange(ny) + 0.5) * dy
    X, Y = np.meshgrid(x, y, indexing="ij")
    oil = (np.ones_like(X, bool) if rect
           else np.hypot(X - pitch_x / 2, Y) > D / 2)
    N = int(oil.sum())
    Acs = N * dx * dy

    # ---- velocity ----
    Av, idx = _assemble(oil, dx, dy, masked="dirichlet", top="dirichlet",
                        bottom="dirichlet" if rect else "none")
    w = spla.spsolve(Av, -np.ones(N))
    wf = np.zeros_like(X); wf[oil] = w
    wbar = w.mean()
    cf = ([(i, 0, dx, dy) for i in range(nx)] if rect
          else _cell_faces(oil, dx, dy))
    pf = [(i, ny - 1, dx, dy) for i in range(nx) if oil[i, ny - 1]]
    Pc = sum(f[2] for f in cf)
    Pp = sum(f[2] for f in pf)
    Dh = 4 * Acs / (Pc + Pp)
    fRe = 2 * Dh ** 2 / wbar

    # ---- energy: pure Neumann operator ----
    Ae, _ = _assemble(oil, dx, dy, masked="neumann")
    wvec = wf[oil]
    Abig = sp.vstack([
        sp.hstack([Ae, sp.csr_matrix(wvec.reshape(-1, 1))]),
        sp.hstack([sp.csr_matrix(wvec.reshape(1, -1)),
                   sp.csr_matrix((1, 1))])]).tocsr()

    def solveT(heated, other):
        per = sum(f[2] for f in heated)
        b = (wvec / wbar) * (1.0 / Acs) / k_oil
        for (i, j, L, hn) in heated:
            b[idx[i, j]] -= (1.0 / per) / k_oil * L / (dx * dy)
        T = spla.spsolve(Abig, np.append(b, 0.0))[:-1]
        Tf = np.zeros_like(X); Tf[oil] = T
        def wT(faces, q_per):
            tot = sum(f[2] for f in faces)
            v = sum((Tf[i, j] + q_per * (hn / 2) / k_oil) * L
                    for (i, j, L, hn) in faces)
            return v / tot
        return wT(heated, 1.0 / per), wT(other, 0.0)

    Tc_A, Tp_A = solveT(cf, pf)
    Tp_B, Tc_B = solveT(pf, cf)
    return dict(fRe=fRe, Dh=Dh, A=Acs, Pc=Pc, Pp=Pp, wbar=wbar,
                Rcc=Tc_A, Rpc=Tp_A, Rpp=Tp_B, Rcp=Tc_B)


if __name__ == "__main__":
    v = duct_coefficients(rect=True, s=0.002, pitch_x=0.02, n=200)
    Nu = (2 * 0.002) / (v["Rpp"] * 0.13 * v["Pp"])
    print(f"parallel-plate limit: fRe = {v['fRe']:.1f} (target 96.0), "
          f"one-side H2 Nu = {Nu:.3f} (target 5.385)")
    print(f"\n{'s [mm]':>7} {'fRe':>7} {'Dh[mm]':>7} {'h_c':>7} "
          f"{'h_p':>7} {'Rcp/Rcc':>8}")
    for s in (0.0015, 0.002, 0.003):
        c = duct_coefficients(s=s, n=160)
        h_c = 1.0 / (c["Rcc"] * c["Pc"])
        h_p = 1.0 / (c["Rpp"] * c["Pp"])
        print(f"{s*1000:7.1f} {c['fRe']:7.1f} {c['Dh']*1000:7.2f} "
              f"{h_c:7.0f} {h_p:7.0f} {c['Rcp']/c['Rcc']:8.3f}")
    for n in (120, 160, 220):
        c = duct_coefficients(n=n)
        print(f"grid {n}: Rcc = {c['Rcc']:.4f}  fRe = {c['fRe']:.1f}")


# ------------------------------------------------------------------ #
#  Graetz kernel: developing transport with two isothermal walls      #
# ------------------------------------------------------------------ #
def graetz_kernel(D=0.021, pitch_x=0.0215, s=0.002, k_oil=0.13,
                  n=110, nz=140, zstar_max=0.5, rect=False):
    """Developing transport with two isothermal walls, marched in
    z* = z / (Dh Re Pr) (i.e. dz* step h with operator scaled on Dh²).
    Two unit cases by superposition: A (T_c=1, T_p=0, Tb0=0) and
    B (T_c=0, T_p=1, Tb0=0). At each z*, [q'_c; q'_p] and the drivers
    [T_c - Tb; T_p - Tb] of both cases give the exact 2x2 kernel a(z*):
        q' = a(z*) . dT      [W/m.K per pitch]
    FD limit checks: rect symmetric-T walls -> Nu = 7.541 each;
    lens FD a_cp reproduces the conduction short-circuit."""
    Lx = pitch_x
    Ly = s if rect else D / 2 + s
    nx = n
    ny = max(int(round(n * Ly / Lx)), 24)
    dx, dy = Lx / nx, Ly / ny
    x = (np.arange(nx) + 0.5) * dx
    y = (np.arange(ny) + 0.5) * dy
    X, Y = np.meshgrid(x, y, indexing="ij")
    oil = (np.ones_like(X, bool) if rect
           else np.hypot(X - pitch_x / 2, Y) > D / 2)
    N = int(oil.sum())
    Acs = N * dx * dy
    Av, idx = _assemble(oil, dx, dy, masked="dirichlet",
                        top="dirichlet",
                        bottom="dirichlet" if rect else "none")
    w = spla.spsolve(Av, -np.ones(N))
    wbar = w.mean()
    wn = w / wbar
    cfa = ([(i, 0, dx, dy) for i in range(nx)] if rect
           else _cell_faces(oil, dx, dy))
    pfa = [(i, ny - 1, dx, dy) for i in range(nx) if oil[i, ny - 1]]
    Pc = sum(f[2] for f in cfa); Pp = sum(f[2] for f in pfa)
    Dh = 4 * Acs / (Pc + Pp)

    def operator(Tc, Tp):
        A, _ = _assemble(oil, dx, dy, masked="neumann")
        A = A.tolil(); b = np.zeros(N)
        for (i, j, L, hn) in cfa:
            p = idx[i, j]
            A[p, p] -= 2.0 / hn ** 2
            b[p] += 2.0 / hn ** 2 * Tc
        for (i, j, L, hn) in pfa:
            p = idx[i, j]
            A[p, p] -= 2.0 / hn ** 2
            b[p] += 2.0 / hn ** 2 * Tp
        return A.tocsr(), b            # lap T sign: A T + b = source

    zs = np.geomspace(3e-4, zstar_max, nz)
    dzs = np.diff(np.concatenate([[0.0], zs]))
    Q = np.zeros((2, nz, 2))           # [case, z, (qc,qp)]
    DT = np.zeros((2, nz, 2))          # drivers
    for case, (Tc, Tp) in enumerate([(1.0, 0.0), (0.0, 1.0)]):
        A, b = operator(Tc, Tp)
        T = np.zeros(N)
        lu, prev = None, None
        for kz, h in enumerate(dzs):
            step = h * Dh ** 2
            if prev != h:
                M = (sp.diags(wn) / step - A).tocsc()
                lu = spla.splu(M)
                prev = h
            T = lu.solve(wn * T / step + b)
            Tb = float((wn * T).mean())
            qc = sum(2 * k_oil / hn * (Tc - T[idx[i, j]]) * L
                     for (i, j, L, hn) in cfa)
            qp = sum(2 * k_oil / hn * (Tp - T[idx[i, j]]) * L
                     for (i, j, L, hn) in pfa)
            Q[case, kz] = (qc, qp)
            DT[case, kz] = (Tc - Tb, Tp - Tb)
    a = np.zeros((nz, 2, 2))
    for kz in range(nz):
        Dmat = np.array([[DT[0, kz, 0], DT[0, kz, 1]],
                         [DT[1, kz, 0], DT[1, kz, 1]]]).T
        Qmat = np.array([[Q[0, kz, 0], Q[1, kz, 0]],
                         [Q[0, kz, 1], Q[1, kz, 1]]])
        a[kz] = Qmat @ np.linalg.inv(Dmat)
    return dict(zs=zs, a=a, Dh=Dh, A=Acs, Pc=Pc, Pp=Pp, wbar=wbar,
                fRe=2 * Dh ** 2 / wbar)


if __name__ == "__main__" and True:
    # gate 3: rect, symmetric FD limit -> Nu 7.541 on each wall
    g = graetz_kernel(rect=True, s=0.002, pitch_x=0.02, n=140, nz=120)
    aFD = g["a"][-1]
    Nu_sym = (aFD[0, 0] + aFD[0, 1]) * (2 * 0.002) / (0.13 * g["Pc"])
    print(f"gate 3 rect FD both-walls-T Nu = "
          f"{(aFD[0,0]+aFD[0,1])*(2*0.002)/(0.13*g['Pc']):.3f} "
          f"(target 7.541)")
    L = graetz_kernel(n=110, nz=120)
    a0, aF = L["a"][8], L["a"][-1]
    print("lens kernel  z*=%.3g:  a=[[%.2f %.2f];[%.2f %.2f]] W/m.K"
          % (L["zs"][8], a0[0,0], a0[0,1], a0[1,0], a0[1,1]))
    print("lens kernel  FD     :  a=[[%.2f %.2f];[%.2f %.2f]]"
          % (aF[0,0], aF[0,1], aF[1,0], aF[1,1]))
    print("FD cell->plate short-circuit -a_cp = %.2f W/m.K "
          "(v1 conduction gave ~%.2f)"
          % (-aF[0,1], 0.134*158*0.033))
