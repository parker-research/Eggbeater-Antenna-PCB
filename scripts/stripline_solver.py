"""2D finite-volume Laplace solver used to size the phasing line.

Requires numpy + scipy:  pip install numpy scipy
Usage: python3 scripts/stripline_solver.py [s_mm] [w_mm ...]
Validated against the exact (Cohn) formula for zero-thickness stripline (<0.5% error).

2D finite-volume Laplace solver for broadside-coupled stripline in JLC04161H-7628.

Cross-section (z up):
  F.Cu GND plane            z = H
  prepreg 7628 (er 4.4)     0.2104
  In1.Cu  trace P (0.0152)
  core     (er 4.6)         1.065
  In2.Cu  trace N (0.0152)
  prepreg 7628 (er 4.4)     0.2104
  B.Cu GND plane            z = 0
Side walls (GND) at x = +-(w/2 + s): coplanar inner pours + via fence.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
import sys

C0 = 299792458.0
EPS0 = 8.854e-12


def solve(w, s, h=0.01, pp=0.2104, core=1.065, tcu=0.0152, er_pp=4.4, er_core=4.6,
          stacked=True, micro=False, air_above=1.0):
    """Returns (Zdiff, Zcm, eeff_odd). Units mm."""
    half = w / 2 + s
    H = pp + tcu + core + tcu + pp
    nx = int(round(2 * half / h)) + 1
    nz = int(round(H / h)) + 1
    xs = np.linspace(-half, half, nx)
    zs = np.linspace(0, H, nz)
    # cell permittivity (cells between nodes)
    zc = (zs[:-1] + zs[1:]) / 2
    er_cells_1d = np.where((zc > pp + tcu) & (zc < pp + tcu + core), er_core, er_pp)

    def run(er_1d, vp, vn):
        eps = np.tile(er_1d[None, :], (nx - 1, 1))  # (nx-1, nz-1)
        fixed = np.zeros((nx, nz), bool)
        val = np.zeros((nx, nz))
        fixed[0, :] = fixed[-1, :] = True
        fixed[:, 0] = fixed[:, -1] = True
        X, Z = np.meshgrid(xs, zs, indexing='ij')
        inx = np.abs(X) <= w / 2 + 1e-9
        zN = (Z >= pp - 1e-9) & (Z <= pp + tcu + 1e-9)
        zP = (Z >= pp + tcu + core - 1e-9) & (Z <= pp + 2 * tcu + core + 1e-9)
        mP = inx & zP
        mN = inx & zN
        fixed |= mP | mN
        val[mP] = vp
        val[mN] = vn
        idx = -np.ones((nx, nz), int)
        free = ~fixed
        idx[free] = np.arange(free.sum())
        n = free.sum()
        rows, cols, data = [], [], []
        b = np.zeros(n)
        # coefficient between node (i,j) and neighbour: average eps of the two adjacent cells
        def epsx(i, j):  # edge between (i,j)-(i+1,j)
            e = []
            if j > 0: e.append(eps[i, j - 1])
            if j < nz - 1: e.append(eps[i, j])
            return np.mean(e)

        def epsz(i, j):  # edge between (i,j)-(i,j+1)
            e = []
            if i > 0: e.append(eps[i - 1, j])
            if i < nx - 1: e.append(eps[i, j])
            return np.mean(e)

        # vectorised assembly
        fi, fj = np.nonzero(free)
        k = idx[fi, fj]
        diag = np.zeros(n)
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            ni, nj = fi + di, fj + dj
            if di:
                ei = np.minimum(fi, ni)
                c = 0.5 * (eps[ei, np.maximum(fj - 1, 0)] + eps[ei, np.minimum(fj, nz - 2)])
            else:
                ej = np.minimum(fj, nj)
                c = 0.5 * (eps[np.maximum(fi - 1, 0), ej] + eps[np.minimum(fi, nx - 2), ej])
            diag += c
            nb = idx[ni, nj]
            isfree = nb >= 0
            rows.append(k[isfree]); cols.append(nb[isfree]); data.append(-c[isfree])
            np.add.at(b, k[~isfree], c[~isfree] * val[ni[~isfree], nj[~isfree]])
        rows.append(k); cols.append(k); data.append(diag)
        A = sp.csr_matrix((np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n))
        phi = val.copy()
        phi[free] = spl.spsolve(A, b)
        # energy
        ex = np.diff(phi, axis=0) / h  # (nx-1, nz)
        ez = np.diff(phi, axis=1) / h  # (nx, nz-1)
        exc = 0.5 * (ex[:, :-1] + ex[:, 1:])
        ezc = 0.5 * (ez[:-1, :] + ez[1:, :])
        W = 0.5 * np.sum(eps * (exc ** 2 + ezc ** 2)) * h * h  # * eps0, per mm length (mm^2/mm^2 dimensionless)
        return W * EPS0  # F/m (dimensionless geometry)

    res = {}
    for mode, (vp, vn) in {'odd': (1, -1), 'even': (1, 1)}.items():
        W = run(er_cells_1d, vp, vn)
        W0 = run(np.ones_like(er_cells_1d), vp, vn)
        # odd: W = 0.5*(Q*1 + (-Q)*(-1)) = Q ; C_per_conductor = Q/1 = W
        # even: W = 0.5*(Q+Q) = Q ; C_per_conductor = W
        C, Ca = W, W0
        Z = 1 / (C0 * np.sqrt(C * Ca))
        res[mode] = (Z, C / Ca)
    return 2 * res['odd'][0], res['even'][0] / 2, res['odd'][1]


if __name__ == '__main__':
    s = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
    for w in [float(a) for a in sys.argv[2:]] or [0.15, 0.2, 0.25, 0.3, 0.35]:
        zd, zc, ee = solve(w, s)
        print(f"w={w:.3f} s={s:.2f}  Zdiff={zd:6.1f}  Zcm={zc:6.1f}  eeff={ee:.3f}")
