"""
=================================================================
 ehtools2d — 2D kernels (Phase 2: coefficient assembly)
=================================================================

Faithful port of the C ``calculate_coefficients`` (`legacy_eh/2dTfree/tools.h`)
for a **solid filler block** (``tipo 2``: nodes ``NT`` 1-8 and 13). The air-cavity
nodes (0, 9-12) arrive in Phase 6.

Implicit finite-volume scheme, 5 points, conductivity by **harmonic mean**
``kh(a,b)=2ab/(a+b)``. For each node ``P=(i,j)``:

    aP·T_P = aE·T_E + aW·T_W + aN·T_N + aS·T_S + apo·T_P° + S_b

with the y neighbours (N=j-1 outside, S=j+1 inside) **deferred** to the ``d`` term
(line-by-line Gauss-Seidel; the implicit TDMA acts only in x). The solver builds
``a=aP``, ``b=aE``, ``c=aW`` and ``d = aN·T_N + aS·T_S + apo·To + S_b``.

Mesh convention identical to `eh2d`: ``(nx,ny)`` arrays indexed ``[i,j]``,
``j=0`` outside (``Tsa,ho``), ``j=ny-1`` inside (``Tint,hi``); sides
``i=0``/``i=nx-1`` adiabatic.
"""

import numpy as np
from numba import njit, prange

# Node types valid for a solid filler block.
_RELLENA_TYPES = {1, 2, 3, 4, 5, 6, 7, 8, 13}


def _harmonic_faces(k, dx, dy):
    """
    Face conductances by harmonic mean, already scaled by the geometry.
    Returns ``(aN, aS, aE, aW)`` of shape ``(nx,ny)``; each one is 0 on the
    boundary where the neighbour does not exist (matching the coefficients the
    C sets to zero on those boundaries).
    """
    aN = np.zeros_like(k)
    aS = np.zeros_like(k)
    aE = np.zeros_like(k)
    aW = np.zeros_like(k)

    def kh(a, b):
        return 2.0 * a * b / (a + b)

    # North (j-1, towards outside) and South (j+1, towards inside): factor dx/dy.
    aN[:, 1:] = kh(k[:, 1:], k[:, :-1]) * dx / dy
    aS[:, :-1] = kh(k[:, :-1], k[:, 1:]) * dx / dy
    # East (i+1) and West (i-1): factor dy/dx.
    aE[:-1, :] = kh(k[:-1, :], k[1:, :]) * dy / dx
    aW[1:, :] = kh(k[1:, :], k[:-1, :]) * dy / dx
    return aN, aS, aE, aW


def calculate_coefficients_2d(NT, k, rhoc, To, T, Tsa, Tint, ho, hi, dt, dx, dy):
    """
    Assembly of ``a,b,c,d`` for a solid filler block (``tipo 2``).

    Args:
        NT (int array (nx,ny)): node types (must be ⊆ {1-8,13}).
        k, rhoc (float array (nx,ny)): conductivity and heat capacity per node.
        To (float array (nx,ny)): field of the previous time step (``T°``).
        T  (float array (nx,ny)): current field (for the deferred N/S neighbours).
        Tsa, Tint (float): sun-air (outside) and indoor-air temperature.
        ho, hi (float): outdoor and indoor convective coefficients.
        dt, dx, dy (float): time and mesh step.

    Returns:
        (a, b, c, d) arrays (nx,ny). ``a=aP``, ``b=aE``, ``c=aW``,
        ``d = aN·T_N + aS·T_S + apo·To + S_b``.
    """
    NT = np.asarray(NT)
    bad = set(np.unique(NT).tolist()) - _RELLENA_TYPES
    if bad:
        raise NotImplementedError(
            f"calculate_coefficients_2d only supports a solid filler block "
            f"(NT ⊆ {sorted(_RELLENA_TYPES)}); got {sorted(bad)}.")

    k = np.asarray(k, dtype=np.float64)
    rhoc = np.asarray(rhoc, dtype=np.float64)
    To = np.asarray(To, dtype=np.float64)
    T = np.asarray(T, dtype=np.float64)
    nx, ny = k.shape

    apo = rhoc * dx * dy / dt
    aN, aS, aE, aW = _harmonic_faces(k, dx, dy)

    # Deferred y neighbours (on the boundaries the coef. is already 0 → null term).
    TN = np.zeros_like(T)
    TS = np.zeros_like(T)
    TN[:, 1:] = T[:, :-1]   # T_N = T[i, j-1]
    TS[:, :-1] = T[:, 1:]   # T_S = T[i, j+1]

    # aP = apo + (active conductions) ; the inactive ones are already 0.
    a = apo + aN + aS + aE + aW
    d = apo * To + aN * TN + aS * TS

    # Convective boundaries: outside at j=0 (Tsa, ho), inside at j=ny-1 (Tint, hi).
    a[:, 0] += ho * dx
    d[:, 0] += ho * dx * Tsa
    a[:, ny - 1] += hi * dx
    d[:, ny - 1] += hi * dx * Tint

    b = aE.copy()   # = aE (0 at i=nx-1, adiabatic right boundaries/types 2,4,7)
    c = aW.copy()   # = aW (0 at i=0,    adiabatic left boundaries/types 1,3,6)

    return a, b, c, d


def _tdma_rows(a, b, c, d):
    """
    Implicit x TDMA (Thomas) for **all** rows ``j`` at once.
    ``a·T_i = b·T_{i+1} + c·T_{i-1} + d`` (b=aE, c=aW). Vectorised over j;
    the i sweep replicates exactly the C order of operations.
    """
    nx, ny = a.shape
    P = np.zeros((nx, ny))
    Q = np.zeros((nx, ny))
    Tn = np.zeros((nx, ny))

    P[0] = b[0] / a[0]
    Q[0] = d[0] / a[0]
    for i in range(1, nx):
        denom = a[i] - c[i] * P[i - 1]
        P[i] = b[i] / denom
        Q[i] = (d[i] + c[i] * Q[i - 1]) / denom

    Tn[nx - 1] = Q[nx - 1]
    for i in range(nx - 2, -1, -1):
        Tn[i] = P[i] * Tn[i + 1] + Q[i]
    return Tn


def solve_step_2d(NT, k, rhoc, To, Tsa, Tint, ho, hi, dt, dx, dy,
                  La, X, rhoair, cair, tol=1e-10, max_iter=100000):
    """
    Port of ``solve_PQ`` for a solid filler block: **one** time step.

    Inner loop (line-by-line Gauss-Seidel): recomputes coefficients with the most
    recent `T`, solves each row with x TDMA (y neighbours deferred), and repeats
    until ``|error| <= tol`` with ``error = Σ (T-Tn)/T /(nx·ny)`` (signed, like the
    C). Then updates the indoor air `Tint` (lumped node).

    Args:
        To: previous-step field; also the initial condition of the step
            (in the C, at the start of the step ``To == T``).
        Tint: indoor-air temperature at the start of the step.
        La, X, rhoair, cair: indoor-air node parameters.

    Returns:
        dict with ``T`` (solved field), ``Tint`` (updated), ``iters``,
        ``Qin``, ``error``.
    """
    To = np.asarray(To, dtype=np.float64)
    nx, ny = To.shape
    Ti = float(Tint)            # "old" Tint, constant during the inner loop
    T = To.copy()               # initial condition = previous field

    iters = 0
    error = 0.0
    while True:
        iters += 1
        a, b, c, d = calculate_coefficients_2d(
            NT, k, rhoc, To, T, Tsa, Ti, ho, hi, dt, dx, dy)
        Tn = _tdma_rows(a, b, c, d)
        error = float(np.sum((T - Tn) / T) / nx / ny)
        T = Tn
        if abs(error) <= tol or iters >= max_iter:
            break

    # Indoor convective flux and indoor-air update (lumped node).
    Tsurf = T[:, ny - 1]
    Qh = float(np.sum(hi * dt * dx * (Tsurf - Ti)))
    Qin = float(np.sum(np.where(Tsurf > Ti, hi * dx * (Tsurf - Ti), 0.0)))
    Cair = rhoair * cair * La * X
    Tint_new = (Qh + (Cair / dt) * Ti) * dt / Cair

    return {"T": T, "Tint": Tint_new, "iters": iters, "Qin": Qin, "error": error}


# =================================================================
#  Production engine (numba JIT) — Phase 5
# =================================================================
#
# Full-day solver with day-to-day convergence, compiled kernels.
# DELIBERATE differences from the faithful port (Phase 4), documented and
# backed by the Phase 4 regression:
#   1. Indoor-air update with a SINGLE dt (physically correct):
#        Tint += dt·Σ_i hi·dx·(Tsurf_i − Ti)/(ρc·La·X)
#      (the C had a latent dt² that only matches at dt=1). This makes the 2D reduce
#      EXACTLY to the package's 1D at any dt.
#   2. Surface averages to /(nx-1) (half nodes), for Tso and Tsi.
#   3. ho, hi are taken as-is (without the C wall override); the choice of
#      hi is made by the API layer as appropriate.


@njit(cache=True)
def _step_inner(k, rhoc, To, T, Tsa, Tint, ho, hi, dt, dx, dy,
                a, b, c, d, P, Q, Tn, Tnew, tol):
    """Inner loop (line-by-line Gauss-Seidel) of one step; updates ``T`` in place.

    Work buffers (preallocated): ``a,b,c,d`` (nx,ny), ``P,Q,Tn`` (nx),
    ``Tnew`` (nx,ny). Assumes NT ⊆ {1-8,13} (adiabatic sides, convective
    boundaries at j=0/j=ny-1) — valid for a solid filler block.
    """
    nx, ny = k.shape
    iters = 0
    while True:
        iters += 1
        # ---- assemble a,b,c,d ----
        for j in range(ny):
            for i in range(nx):
                kij = k[i, j]
                apo = rhoc[i, j] * dx * dy / dt
                aN = aS = aE = aW = 0.0
                if j > 0:
                    kk = k[i, j - 1]
                    aN = 2.0 * kk * kij / (kk + kij) * dx / dy
                if j < ny - 1:
                    kk = k[i, j + 1]
                    aS = 2.0 * kk * kij / (kk + kij) * dx / dy
                if i < nx - 1:
                    kk = k[i + 1, j]
                    aE = 2.0 * kk * kij / (kk + kij) * dy / dx
                if i > 0:
                    kk = k[i - 1, j]
                    aW = 2.0 * kk * kij / (kk + kij) * dy / dx
                aP = apo + aN + aS + aE + aW
                dd = apo * To[i, j]
                if j > 0:
                    dd += aN * T[i, j - 1]
                if j < ny - 1:
                    dd += aS * T[i, j + 1]
                if j == 0:
                    aP += ho * dx
                    dd += ho * dx * Tsa
                if j == ny - 1:
                    aP += hi * dx
                    dd += hi * dx * Tint
                a[i, j] = aP
                b[i, j] = aE
                c[i, j] = aW
                d[i, j] = dd
        # ---- x TDMA for each row j -> Tnew ----
        for j in range(ny):
            P[0] = b[0, j] / a[0, j]
            Q[0] = d[0, j] / a[0, j]
            for i in range(1, nx):
                denom = a[i, j] - c[i, j] * P[i - 1]
                P[i] = b[i, j] / denom
                Q[i] = (d[i, j] + c[i, j] * Q[i - 1]) / denom
            Tnew[nx - 1, j] = Q[nx - 1]
            for i in range(nx - 2, -1, -1):
                Tnew[i, j] = P[i] * Tnew[i + 1, j] + Q[i]
        # ---- signed error and commit T <- Tnew ----
        error = 0.0
        for j in range(ny):
            for i in range(nx):
                error += (T[i, j] - Tnew[i, j]) / T[i, j]
                T[i, j] = Tnew[i, j]
        error = error / nx / ny
        if abs(error) <= tol:
            break
    return iters


@njit(cache=True)
def solve_day_2d(NT, k, rhoc, Tsa_arr, ho, hi, dt, dx, dy, La, X,
                 rhoair, cair, T0, tol_inner=1e-10, tol_day=5e-4, max_days=60):
    """
    Production engine: runs the full day, repeating it until periodic steady
    state (``mean|T_day−T_prev_day| < tol_day``).

    Returns:
        (Ti_series, Tso_series, Tsi_series, T_field, days, Qin, Qout)
        Series of shape ``(nsteps,)`` (one value per ``Tsa_arr`` step).
    """
    nx, ny = k.shape
    nsteps = Tsa_arr.shape[0]

    T = np.empty((nx, ny))
    To = np.empty((nx, ny))
    Told = np.empty((nx, ny))
    for j in range(ny):
        for i in range(nx):
            T[i, j] = T0
    Tint = T0

    a = np.empty((nx, ny)); b = np.empty((nx, ny))
    c = np.empty((nx, ny)); d = np.empty((nx, ny))
    Tnew = np.empty((nx, ny))
    P = np.empty(nx); Q = np.empty(nx); Tn = np.empty(nx)

    Ti_series = np.empty(nsteps)
    Tso_series = np.empty(nsteps)
    Tsi_series = np.empty(nsteps)

    Cair = rhoair * cair * La * X
    days = 0
    Qin = Qout = 0.0
    C = 1.0e9
    while C > tol_day and days < max_days:
        for j in range(ny):
            for i in range(nx):
                Told[i, j] = T[i, j]
        Qin = Qout = 0.0
        for s in range(nsteps):
            # Tso: outer surface before solving, /(nx-1)
            tso = 0.0
            for i in range(nx):
                tso += T[i, 0]
            Tso_series[s] = tso / (nx - 1)
            # To <- T
            for j in range(ny):
                for i in range(nx):
                    To[i, j] = T[i, j]
            Ti_old = Tint
            _step_inner(k, rhoc, To, T, Tsa_arr[s], Tint, ho, hi, dt, dx, dy,
                        a, b, c, d, P, Q, Tn, Tnew, tol_inner)
            # indoor-air update (ONE dt, physically correct)
            flux = 0.0
            for i in range(nx):
                flux += hi * dx * (T[i, ny - 1] - Ti_old)
            Tint = Ti_old + dt * flux / Cair
            Ti_series[s] = Tint
            # Tsi: inner surface after solving, /(nx-1)
            tsi = 0.0
            for i in range(nx):
                tsi += T[i, ny - 1]
            Tsi_series[s] = tsi / (nx - 1)
            # energy (per unit indoor area)
            e = flux * dt / X
            if e > 0.0:
                Qin += e
            else:
                Qout -= e
        C = 0.0
        for j in range(ny):
            for i in range(nx):
                C += abs(Told[i, j] - T[i, j])
        C = C / nx / ny
        days += 1

    return Ti_series, Tso_series, Tsi_series, T, days, Qin, Qout


# =================================================================
#  Filler block with air cavity (tipo 1) — Phase 6
# =================================================================
#
# Cavity physics: the cavity walls (NT 9-12) convect to the cavity air
# (coef. hh by Nusselt) and radiate among themselves (Stefan-Boltzmann + view
# factors); the cavity air is a lumped node Thueco (NT 0 fixes T=Thueco). Faithful
# port of the tipo==1 branch of solve_PQ + cases 0,9-12 of calculate_coefficients.
# Wall only (beta=90): hh = 0.4005·|ΔT|^0.3033 / e22^0.0901.

_SIGMA = 5.6704e-8


@njit(cache=True)
def _step_hueca(k, rhoc, To, T, Tsa, Tint, Th, ho, hi, dt, dx, dy,
                i1, j1, i2, j2, e22, E,
                Fud, Ful, Fur, Fru, Frd, Frl, Fdl, Fdr, Fdu, Flu, Flr, Fld,
                a, b, c, d, P, Q, Tnew, tol):
    """Inner loop of one step for a filler block with air cavity (wall)."""
    nx, ny = k.shape
    iters = 0
    hh = 1.0
    while True:
        iters += 1
        # --- mean temperatures of the 4 walls (current T) ---
        tup = tdn = tlf = trt = 0.0
        for i in range(i1, i2):
            tup += T[i, j1 - 1]
            tdn += T[i, j2]
        nud = i2 - i1
        for j in range(j1, j2):
            tlf += T[i1 - 1, j]
            trt += T[i2, j]
        nlr = j2 - j1
        tup /= nud; tdn /= nud; tlf /= nlr; trt /= nlr
        # --- Nusselt (wall) ---
        hh = 0.4005 * (abs(tup - tdn) ** 0.3033) / (e22 ** 0.0901)
        # --- radiation between walls (view factors) ---
        Tu = tup + 273.15; Td = tdn + 273.15; Tl = tlf + 273.15; Tr = trt + 273.15
        Tu4 = Tu * Tu * Tu * Tu; Td4 = Td * Td * Td * Td
        Tl4 = Tl * Tl * Tl * Tl; Tr4 = Tr * Tr * Tr * Tr
        sx = dx * E * _SIGMA
        sy = dy * E * _SIGMA
        Qud = sx * (Tu4 - Td4) * Fud
        Qul = sx * (Tu4 - Tl4) * Ful
        Qur = sx * (Tu4 - Tr4) * Fur
        Qru = sy * (Tr4 - Tu4) * Fru
        Qrd = sy * (Tr4 - Td4) * Frd
        Qrl = sy * (Tr4 - Tl4) * Frl
        Qdu = sx * (Td4 - Tu4) * Fdu
        Qdr = sx * (Td4 - Tr4) * Fdr
        Qdl = sx * (Td4 - Tl4) * Fdl
        Qlu = sy * (Tl4 - Tu4) * Flu
        Qlr = sy * (Tl4 - Tr4) * Flr
        Qld = sy * (Tl4 - Td4) * Fld
        # --- assemble a,b,c,d ---
        for j in range(ny):
            for i in range(nx):
                nt_ij = 0
                # node type deduced from the cavity coordinates
                kij = k[i, j]
                apo = rhoc[i, j] * dx * dy / dt
                aN = aS = aE = aW = 0.0
                if j > 0:
                    kk = k[i, j - 1]; aN = 2.0 * kk * kij / (kk + kij) * dx / dy
                if j < ny - 1:
                    kk = k[i, j + 1]; aS = 2.0 * kk * kij / (kk + kij) * dx / dy
                if i < nx - 1:
                    kk = k[i + 1, j]; aE = 2.0 * kk * kij / (kk + kij) * dy / dx
                if i > 0:
                    kk = k[i - 1, j]; aW = 2.0 * kk * kij / (kk + kij) * dy / dx
                in_hole_cols = (i1 <= i) and (i < i2)
                in_hole_rows = (j1 <= j) and (j < j2)
                if in_hole_cols and in_hole_rows:
                    # cavity air: T = Th
                    a[i, j] = 1.0; b[i, j] = 0.0; c[i, j] = 0.0; d[i, j] = Th
                elif in_hole_cols and j == j1 - 1:
                    # top wall (NT 9): cavity to the south
                    aP = apo + aN + hh * dx + aE + aW
                    dd = aN * T[i, j - 1] + hh * dx * Th + apo * To[i, j] - Qur - Qud - Qul
                    a[i, j] = aP; b[i, j] = aE; c[i, j] = aW; d[i, j] = dd
                elif in_hole_cols and j == j2:
                    # bottom wall (NT 10): cavity to the north
                    aP = apo + hh * dx + aS + aE + aW
                    dd = aS * T[i, j + 1] + hh * dx * Th + apo * To[i, j] - Qdl - Qdu - Qdr
                    a[i, j] = aP; b[i, j] = aE; c[i, j] = aW; d[i, j] = dd
                elif in_hole_rows and i == i1 - 1:
                    # left wall (NT 11): cavity to the east
                    aP = apo + aN + aS + hh * dy + aW
                    dd = aN * T[i, j - 1] + aS * T[i, j + 1] + apo * To[i, j] + hh * dy * Th - Qlu - Qlr - Qld
                    a[i, j] = aP; b[i, j] = 0.0; c[i, j] = aW; d[i, j] = dd
                elif in_hole_rows and i == i2:
                    # right wall (NT 12): cavity to the west
                    aP = apo + aN + aS + aE + hh * dy
                    dd = aN * T[i, j - 1] + aS * T[i, j + 1] + apo * To[i, j] + hh * dy * Th - Qrd - Qrl - Qru
                    a[i, j] = aP; b[i, j] = aE; c[i, j] = 0.0; d[i, j] = dd
                else:
                    # standard node (1-8,13)
                    aP = apo + aN + aS + aE + aW
                    dd = apo * To[i, j]
                    if j > 0:
                        dd += aN * T[i, j - 1]
                    if j < ny - 1:
                        dd += aS * T[i, j + 1]
                    if j == 0:
                        aP += ho * dx; dd += ho * dx * Tsa
                    if j == ny - 1:
                        aP += hi * dx; dd += hi * dx * Tint
                    a[i, j] = aP; b[i, j] = aE; c[i, j] = aW; d[i, j] = dd
        # --- x TDMA per row ---
        for j in range(ny):
            P[0] = b[0, j] / a[0, j]
            Q[0] = d[0, j] / a[0, j]
            for i in range(1, nx):
                denom = a[i, j] - c[i, j] * P[i - 1]
                P[i] = b[i, j] / denom
                Q[i] = (d[i, j] + c[i, j] * Q[i - 1]) / denom
            Tnew[nx - 1, j] = Q[nx - 1]
            for i in range(nx - 2, -1, -1):
                Tnew[i, j] = P[i] * Tnew[i + 1, j] + Q[i]
        # --- error and commit ---
        error = 0.0
        for j in range(ny):
            for i in range(nx):
                error += (T[i, j] - Tnew[i, j]) / T[i, j]
                T[i, j] = Tnew[i, j]
        error = error / nx / ny
        if abs(error) <= tol:
            break
    return iters, hh


def _view_factors(a21, e22):
    """Cavity view factors (h=e22 height, l=a21 width), like the C."""
    h, l = e22, a21
    Fur = 0.5 * (1.0 + h / l - (1.0 + (h * h) / (l * l)) ** 0.5)
    Ful = Fur
    Fud = 1.0 - 2.0 * Fur
    Fru = 0.5 * (1.0 + l / h - (1.0 + (l * l) / (h * h)) ** 0.5)
    Frd = Fru
    Frl = 1.0 - 2.0 * Fru
    Fdl, Fdr, Fdu = Ful, Fur, Fud
    Flu, Flr, Fld = Fru, Frl, Frd
    return Fud, Ful, Fur, Fru, Frd, Frl, Fdl, Fdr, Fdu, Flu, Flr, Fld


def solve_step_hueca(NT, k, rhoc, To, Tsa, Tint, Thueco, ho, hi, dt, dx, dy,
                     La, X, rhoair, cair, i1, j1, i2, j2, a21, e22, E, beta,
                     tol=1e-10):
    """
    One time step for a filler block with air cavity (wall, beta=90).
    Returns a dict with ``T, Tint, Thueco, hh, iters``.
    """
    if beta != 90.0:
        raise NotImplementedError("solve_step_hueca: wall only (beta=90) for now.")
    k = np.asarray(k, dtype=np.float64)
    rhoc = np.asarray(rhoc, dtype=np.float64)
    To = np.asarray(To, dtype=np.float64)
    nx, ny = k.shape
    Ti = float(Tint); Th = float(Thueco)
    T = To.copy()
    a = np.empty((nx, ny)); b = np.empty((nx, ny))
    cc = np.empty((nx, ny)); d = np.empty((nx, ny)); Tnew = np.empty((nx, ny))
    P = np.empty(nx); Q = np.empty(nx)
    vf = _view_factors(a21, e22)
    iters, hh = _step_hueca(k, rhoc, To, T, Tsa, Ti, Th, ho, hi, dt, dx, dy,
                            i1, j1, i2, j2, e22, E, *vf,
                            a, b, cc, d, P, Q, Tnew, tol)
    # cavity air (lumped, single dt)
    Qh_hole = 0.0
    for i in range(i1, i2):
        Qh_hole += hh * dx * (T[i, j1 - 1] - Th)
        Qh_hole += hh * dx * (T[i, j2] - Th)
    for j in range(j1, j2):
        Qh_hole += hh * dy * (T[i1 - 1, j] - Th)
        Qh_hole += hh * dy * (T[i2, j] - Th)
    Ch = rhoair * cair * a21 * e22
    Th_new = (Qh_hole + (Ch / dt) * Th) * dt / Ch
    # indoor air (identical to the faithful tipo 2 port: dt²)
    Tsurf = T[:, ny - 1]
    Qh = float(np.sum(hi * dt * dx * (Tsurf - Ti)))
    Cair = rhoair * cair * La * X
    Tint_new = (Qh + (Cair / dt) * Ti) * dt / Cair
    return {"T": T, "Tint": Tint_new, "Thueco": Th_new, "hh": hh, "iters": iters}


@njit(cache=True)
def solve_day_hueca(NT, k, rhoc, Tsa_arr, ho, hi, dt, dx, dy, La, X,
                    rhoair, cair, T0, i1, j1, i2, j2, a21, e22, E,
                    Fud, Ful, Fur, Fru, Frd, Frl, Fdl, Fdr, Fdu, Flu, Flr, Fld,
                    tol_inner, tol_day, max_days):
    """
    Full day with day-to-day convergence for a filler block with air cavity.
    Faithful replica of the C (Tso /nx, Tsi /(nx-1), Tint with dt², Thueco lumped).
    The cavity air ``Th`` and the indoor air ``Tint`` march step by step.

    Returns:
        (Ti_series, Tso_series, Tsi_series, Th_series, T_field, days)
    """
    nx, ny = k.shape
    nsteps = Tsa_arr.shape[0]
    T = np.empty((nx, ny)); To = np.empty((nx, ny)); Told = np.empty((nx, ny))
    for j in range(ny):
        for i in range(nx):
            T[i, j] = T0
    Tint = T0
    Th = T0
    a = np.empty((nx, ny)); b = np.empty((nx, ny)); c = np.empty((nx, ny))
    d = np.empty((nx, ny)); Tnew = np.empty((nx, ny))
    P = np.empty(nx); Q = np.empty(nx)
    Ti_s = np.empty(nsteps); Tso_s = np.empty(nsteps)
    Tsi_s = np.empty(nsteps); Th_s = np.empty(nsteps)
    Cair = rhoair * cair * La * X
    Ch = rhoair * cair * a21 * e22
    days = 0
    C = 1.0e9
    while C > tol_day and days < max_days:
        for j in range(ny):
            for i in range(nx):
                Told[i, j] = T[i, j]
        for s in range(nsteps):
            tso = 0.0
            for i in range(nx):
                tso += T[i, 0]
            Tso_s[s] = tso / nx          # C's Tsout: /nx (faithful replica)
            for j in range(ny):
                for i in range(nx):
                    To[i, j] = T[i, j]
            Ti_old = Tint
            Th_old = Th
            _, hh = _step_hueca(k, rhoc, To, T, Tsa_arr[s], Ti_old, Th_old,
                                ho, hi, dt, dx, dy, i1, j1, i2, j2, e22, E,
                                Fud, Ful, Fur, Fru, Frd, Frl, Fdl, Fdr, Fdu, Flu, Flr, Fld,
                                a, b, c, d, P, Q, Tnew, tol_inner)
            # cavity air
            qh = 0.0
            for i in range(i1, i2):
                qh += hh * dx * (T[i, j1 - 1] - Th_old)
                qh += hh * dx * (T[i, j2] - Th_old)
            for j in range(j1, j2):
                qh += hh * dy * (T[i1 - 1, j] - Th_old)
                qh += hh * dy * (T[i2, j] - Th_old)
            Th = (qh + (Ch / dt) * Th_old) * dt / Ch
            # indoor air (dt², faithful)
            qi = 0.0
            for i in range(nx):
                qi += hi * dt * dx * (T[i, ny - 1] - Ti_old)
            Tint = (qi + (Cair / dt) * Ti_old) * dt / Cair
            Ti_s[s] = Tint
            Th_s[s] = Th
            tsi = 0.0
            for i in range(nx):
                tsi += T[i, ny - 1]
            Tsi_s[s] = tsi / (nx - 1)    # C's max_min: /(nx-1)
        C = 0.0
        for j in range(ny):
            for i in range(nx):
                C += abs(Told[i, j] - T[i, j])
        C = C / nx / ny
        days += 1
    return Ti_s, Tso_s, Tsi_s, Th_s, T, days


@njit(cache=True)
def solve_day_hueca_prod(NT, k, rhoc, Tsa_arr, ho, hi, dt, dx, dy, La, X,
                         rhoair, cair, T0, i1, j1, i2, j2, a21, e22, E,
                         Fud, Ful, Fur, Fru, Frd, Frl, Fdl, Fdr, Fdu, Flu, Flr, Fld,
                         tol_inner, tol_day, max_days):
    """
    **Production** version of the day with air cavity (hollow block / filler
    block with air). Same as :func:`solve_day_hueca` but with the Phase 5
    deliberate corrections: indoor air with a **single `dt`** and surfaces at
    `/(nx-1)`. The cavity physics (radiation + Nusselt + `Thueco`) is identical.
    Also returns `Qin, Qout` (energy per unit indoor area, last day).

    Returns:
        (Ti_series, Tso_series, Tsi_series, Th_series, T_field, days, Qin, Qout)
    """
    nx, ny = k.shape
    nsteps = Tsa_arr.shape[0]
    T = np.empty((nx, ny)); To = np.empty((nx, ny)); Told = np.empty((nx, ny))
    for j in range(ny):
        for i in range(nx):
            T[i, j] = T0
    Tint = T0
    Th = T0
    a = np.empty((nx, ny)); b = np.empty((nx, ny)); c = np.empty((nx, ny))
    d = np.empty((nx, ny)); Tnew = np.empty((nx, ny))
    P = np.empty(nx); Q = np.empty(nx)
    Ti_s = np.empty(nsteps); Tso_s = np.empty(nsteps)
    Tsi_s = np.empty(nsteps); Th_s = np.empty(nsteps)
    Cair = rhoair * cair * La * X
    Ch = rhoair * cair * a21 * e22
    days = 0
    C = 1.0e9
    Qin = Qout = 0.0
    while C > tol_day and days < max_days:
        for j in range(ny):
            for i in range(nx):
                Told[i, j] = T[i, j]
        Qin = Qout = 0.0
        for s in range(nsteps):
            tso = 0.0
            for i in range(nx):
                tso += T[i, 0]
            Tso_s[s] = tso / (nx - 1)        # production: /(nx-1)
            for j in range(ny):
                for i in range(nx):
                    To[i, j] = T[i, j]
            Ti_old = Tint
            Th_old = Th
            _, hh = _step_hueca(k, rhoc, To, T, Tsa_arr[s], Ti_old, Th_old,
                                ho, hi, dt, dx, dy, i1, j1, i2, j2, e22, E,
                                Fud, Ful, Fur, Fru, Frd, Frl, Fdl, Fdr, Fdu, Flu, Flr, Fld,
                                a, b, c, d, P, Q, Tnew, tol_inner)
            # cavity air (single dt, same as the faithful one — already correct)
            qh = 0.0
            for i in range(i1, i2):
                qh += hh * dx * (T[i, j1 - 1] - Th_old)
                qh += hh * dx * (T[i, j2] - Th_old)
            for j in range(j1, j2):
                qh += hh * dy * (T[i1 - 1, j] - Th_old)
                qh += hh * dy * (T[i2, j] - Th_old)
            Th = (qh + (Ch / dt) * Th_old) * dt / Ch
            # indoor air (single dt, physically correct)
            flux = 0.0
            for i in range(nx):
                flux += hi * dx * (T[i, ny - 1] - Ti_old)
            Tint = Ti_old + dt * flux / Cair
            Ti_s[s] = Tint
            Th_s[s] = Th
            tsi = 0.0
            for i in range(nx):
                tsi += T[i, ny - 1]
            Tsi_s[s] = tsi / (nx - 1)
            e = flux * dt / X
            if e > 0.0:
                Qin += e
            else:
                Qout -= e
        C = 0.0
        for j in range(ny):
            for i in range(nx):
                C += abs(Told[i, j] - T[i, j])
        C = C / nx / ny
        days += 1
    return Ti_s, Tso_s, Tsi_s, Th_s, T, days, Qin, Qout


# =================================================================
#  Parallel variant (numba prange) — Phase 7
# =================================================================
#
# The method is already line-by-line Jacobi (rows use a single snapshot of T per
# inner iteration and are solved into a Tnew buffer before T<-Tnew), so the rows
# are INDEPENDENT: parallelising the j sweep with prange does not change the
# algorithm or the convergence, it only spreads rows across threads. Portable: numba
# prange with its threading layer (workqueue by default, no external libs).


@njit(parallel=True, cache=True)
def _step_inner_par(k, rhoc, To, T, Tsa, Tint, ho, hi, dt, dx, dy,
                    a, b, c, d, Tnew, tol):
    """Same as ``_step_inner`` but with the row sweep in ``prange``.

    Each TDMA row j uses local P,Q buffers (one array per prange iteration) to
    avoid races. The error is a sum reduction over j.
    """
    nx, ny = k.shape
    iters = 0
    while True:
        iters += 1
        # ---- assemble (independent rows) ----
        for j in prange(ny):
            for i in range(nx):
                kij = k[i, j]
                apo = rhoc[i, j] * dx * dy / dt
                aN = aS = aE = aW = 0.0
                if j > 0:
                    kk = k[i, j - 1]; aN = 2.0 * kk * kij / (kk + kij) * dx / dy
                if j < ny - 1:
                    kk = k[i, j + 1]; aS = 2.0 * kk * kij / (kk + kij) * dx / dy
                if i < nx - 1:
                    kk = k[i + 1, j]; aE = 2.0 * kk * kij / (kk + kij) * dy / dx
                if i > 0:
                    kk = k[i - 1, j]; aW = 2.0 * kk * kij / (kk + kij) * dy / dx
                aP = apo + aN + aS + aE + aW
                dd = apo * To[i, j]
                if j > 0:
                    dd += aN * T[i, j - 1]
                if j < ny - 1:
                    dd += aS * T[i, j + 1]
                if j == 0:
                    aP += ho * dx; dd += ho * dx * Tsa
                if j == ny - 1:
                    aP += hi * dx; dd += hi * dx * Tint
                a[i, j] = aP; b[i, j] = aE; c[i, j] = aW; d[i, j] = dd
        # ---- TDMA per row (thread-local P,Q) ----
        for j in prange(ny):
            P = np.empty(nx)
            Q = np.empty(nx)
            P[0] = b[0, j] / a[0, j]
            Q[0] = d[0, j] / a[0, j]
            for i in range(1, nx):
                denom = a[i, j] - c[i, j] * P[i - 1]
                P[i] = b[i, j] / denom
                Q[i] = (d[i, j] + c[i, j] * Q[i - 1]) / denom
            Tnew[nx - 1, j] = Q[nx - 1]
            for i in range(nx - 2, -1, -1):
                Tnew[i, j] = P[i] * Tnew[i + 1, j] + Q[i]
        # ---- error (reduction) and commit ----
        error = 0.0
        for j in prange(ny):
            for i in range(nx):
                error += (T[i, j] - Tnew[i, j]) / T[i, j]
                T[i, j] = Tnew[i, j]
        error = error / nx / ny
        if abs(error) <= tol:
            break
    return iters


@njit(cache=True)
def solve_day_2d_par(NT, k, rhoc, Tsa_arr, ho, hi, dt, dx, dy, La, X,
                     rhoair, cair, T0, tol_inner=1e-10, tol_day=5e-4, max_days=60):
    """Parallel version of :func:`solve_day_2d`: the parallelism (``prange`` over
    rows) lives in ``_step_inner_par``; the day/step loop is sequential."""
    nx, ny = k.shape
    nsteps = Tsa_arr.shape[0]
    T = np.empty((nx, ny)); To = np.empty((nx, ny)); Told = np.empty((nx, ny))
    for j in range(ny):
        for i in range(nx):
            T[i, j] = T0
    Tint = T0
    a = np.empty((nx, ny)); b = np.empty((nx, ny))
    c = np.empty((nx, ny)); d = np.empty((nx, ny)); Tnew = np.empty((nx, ny))
    Ti_series = np.empty(nsteps); Tso_series = np.empty(nsteps); Tsi_series = np.empty(nsteps)
    Cair = rhoair * cair * La * X
    days = 0; Qin = Qout = 0.0; C = 1.0e9
    while C > tol_day and days < max_days:
        for j in range(ny):
            for i in range(nx):
                Told[i, j] = T[i, j]
        Qin = Qout = 0.0
        for s in range(nsteps):
            tso = 0.0
            for i in range(nx):
                tso += T[i, 0]
            Tso_series[s] = tso / (nx - 1)
            for j in range(ny):
                for i in range(nx):
                    To[i, j] = T[i, j]
            Ti_old = Tint
            _step_inner_par(k, rhoc, To, T, Tsa_arr[s], Tint, ho, hi, dt, dx, dy,
                            a, b, c, d, Tnew, tol_inner)
            flux = 0.0
            for i in range(nx):
                flux += hi * dx * (T[i, ny - 1] - Ti_old)
            Tint = Ti_old + dt * flux / Cair
            Ti_series[s] = Tint
            tsi = 0.0
            for i in range(nx):
                tsi += T[i, ny - 1]
            Tsi_series[s] = tsi / (nx - 1)
            e = flux * dt / X
            if e > 0.0:
                Qin += e
            else:
                Qout -= e
        C = 0.0
        for j in range(ny):
            for i in range(nx):
                C += abs(Told[i, j] - T[i, j])
        C = C / nx / ny
        days += 1
    return Ti_series, Tso_series, Tsi_series, T, days, Qin, Qout


# =================================================================
#  Joist and filler block — ROOF, N cavities, 3 solids (Phase 8b)
# =================================================================
#
# Generalises the cavity physics (Phase 6) to a roof slab with:
#   - N equal air cavities (each one a lumped node Thueco[c]);
#   - three solid materials (topping, L-shaped joist, filler block) encoded
#     entirely in the per-node k/rhoc fields → the conduction assembly does not
#     change, it only reads k/rhoc;
#   - **roof** Nusselt (Rayleigh, beta=0) in addition to the wall one (beta=90).
# Each node's type is read from NT (1-8,13 standard; 0 air, 9-12 walls) and the
# cavity each air/wall node belongs to, from `cav_of`. Production conventions
# (single dt in Thueco/Tint, surfaces /(nx-1)).
#
# Air properties for the roof Nusselt (identical to the C):
#   gr=9.81, Beta=1/300, nu=1.11e-5, kair=0.0262, alphaair=kair/(rhoair·cair).

_GR = 9.81
_BETA_EXP = 1.0 / 300.0
_NU_AIR = 1.11e-5
_K_AIR = 0.0262


@njit(cache=True)
def _slab_hh(tup, tdn, e22, beta, kair, gr, beta_exp, nu, alphaair):
    """Cavity convective coefficient ``hh``. ``beta=90`` wall, ``beta=0`` roof
    (Rayleigh). ``tup`` = top wall (outside), ``tdn`` = bottom (inside)."""
    if beta == 90.0:
        return 0.4005 * (abs(tup - tdn) ** 0.3033) / (e22 ** 0.0901)
    # roof (beta=0): Rayleigh-Bénard; stable (tdn<=tup) → conduction only.
    if tdn <= tup:
        return kair / e22
    Ra = gr * beta_exp * (tdn - tup) * (e22 ** 3) / nu / alphaair
    dot11 = 1.0 - 1708.0 / Ra
    dot22 = (Ra / 5830.0) ** (1.0 / 3.0) - 1.0
    if dot11 < 0.0:
        dot11 = 0.0
    if dot22 < 0.0:
        dot22 = 0.0
    return kair / e22 * (1.0 + 1.44 * dot11 + dot22)


@njit(cache=True)
def _step_slab(k, rhoc, To, T, Tsa, Tint, Th, ho, hi, dt, dx, dy,
               NT, cav_of, cav_i1, cav_i2, cj1, cj2, n_cav, e22, E, beta,
               kair, gr, beta_exp, nu, alphaair,
               Fud, Ful, Fur, Fru, Frd, Frl, Fdl, Fdr, Fdu, Flu, Flr, Fld,
               a, b, c, d, P, Q, Tnew, hh, Qtop, Qbot, Qleft, Qright, tol):
    """Inner loop (line-by-line Gauss-Seidel) of one step for the N-cavity roof
    slab. ``Th`` (n_cav) constant during the loop; returns ``iters`` and fills
    ``hh`` (n_cav). Updates ``T`` in place."""
    nx, ny = k.shape
    sx = dx * E * _SIGMA
    sy = dy * E * _SIGMA
    iters = 0
    while True:
        iters += 1
        # --- per cavity: mean wall temperatures, hh and net radiation ---
        for cidx in range(n_cav):
            ci1 = cav_i1[cidx]; ci2 = cav_i2[cidx]
            tup = 0.0; tdn = 0.0; tlf = 0.0; trt = 0.0
            for i in range(ci1, ci2):
                tup += T[i, cj1 - 1]
                tdn += T[i, cj2]
            nud = ci2 - ci1
            for j in range(cj1, cj2):
                tlf += T[ci1 - 1, j]
                trt += T[ci2, j]
            nlr = cj2 - cj1
            tup /= nud; tdn /= nud; tlf /= nlr; trt /= nlr
            hh[cidx] = _slab_hh(tup, tdn, e22, beta, kair, gr, beta_exp, nu, alphaair)
            Tu = tup + 273.15; Td = tdn + 273.15
            Tl = tlf + 273.15; Tr = trt + 273.15
            Tu4 = Tu * Tu * Tu * Tu; Td4 = Td * Td * Td * Td
            Tl4 = Tl * Tl * Tl * Tl; Tr4 = Tr * Tr * Tr * Tr
            Qud = sx * (Tu4 - Td4) * Fud; Qul = sx * (Tu4 - Tl4) * Ful
            Qur = sx * (Tu4 - Tr4) * Fur
            Qru = sy * (Tr4 - Tu4) * Fru; Qrd = sy * (Tr4 - Td4) * Frd
            Qrl = sy * (Tr4 - Tl4) * Frl
            Qdu = sx * (Td4 - Tu4) * Fdu; Qdr = sx * (Td4 - Tr4) * Fdr
            Qdl = sx * (Td4 - Tl4) * Fdl
            Qlu = sy * (Tl4 - Tu4) * Flu; Qlr = sy * (Tl4 - Tr4) * Flr
            Qld = sy * (Tl4 - Td4) * Fld
            Qtop[cidx] = Qur + Qud + Qul
            Qbot[cidx] = Qdl + Qdu + Qdr
            Qleft[cidx] = Qlu + Qlr + Qld
            Qright[cidx] = Qrd + Qrl + Qru
        # --- assemble a,b,c,d (driven by NT + cav_of) ---
        for j in range(ny):
            for i in range(nx):
                nt = NT[i, j]
                kij = k[i, j]
                apo = rhoc[i, j] * dx * dy / dt
                aN = aS = aE = aW = 0.0
                if j > 0:
                    kk = k[i, j - 1]; aN = 2.0 * kk * kij / (kk + kij) * dx / dy
                if j < ny - 1:
                    kk = k[i, j + 1]; aS = 2.0 * kk * kij / (kk + kij) * dx / dy
                if i < nx - 1:
                    kk = k[i + 1, j]; aE = 2.0 * kk * kij / (kk + kij) * dy / dx
                if i > 0:
                    kk = k[i - 1, j]; aW = 2.0 * kk * kij / (kk + kij) * dy / dx
                if nt == 0:
                    cc = cav_of[i, j]
                    a[i, j] = 1.0; b[i, j] = 0.0; c[i, j] = 0.0; d[i, j] = Th[cc]
                elif nt == 9:                       # top wall (cavity to the south)
                    cc = cav_of[i, j]; hc = hh[cc]
                    a[i, j] = apo + aN + hc * dx + aE + aW
                    d[i, j] = aN * T[i, j - 1] + hc * dx * Th[cc] + apo * To[i, j] - Qtop[cc]
                    b[i, j] = aE; c[i, j] = aW
                elif nt == 10:                      # bottom wall (cavity to the north)
                    cc = cav_of[i, j]; hc = hh[cc]
                    a[i, j] = apo + hc * dx + aS + aE + aW
                    d[i, j] = aS * T[i, j + 1] + hc * dx * Th[cc] + apo * To[i, j] - Qbot[cc]
                    b[i, j] = aE; c[i, j] = aW
                elif nt == 11:                      # left wall (cavity to the east)
                    cc = cav_of[i, j]; hc = hh[cc]
                    a[i, j] = apo + aN + aS + hc * dy + aW
                    d[i, j] = aN * T[i, j - 1] + aS * T[i, j + 1] + apo * To[i, j] + hc * dy * Th[cc] - Qleft[cc]
                    b[i, j] = 0.0; c[i, j] = aW
                elif nt == 12:                      # right wall (cavity to the west)
                    cc = cav_of[i, j]; hc = hh[cc]
                    a[i, j] = apo + aN + aS + aE + hc * dy
                    d[i, j] = aN * T[i, j - 1] + aS * T[i, j + 1] + apo * To[i, j] + hc * dy * Th[cc] - Qright[cc]
                    b[i, j] = aE; c[i, j] = 0.0
                else:                               # standard (1-8, 13)
                    aP = apo + aN + aS + aE + aW
                    dd = apo * To[i, j]
                    if j > 0:
                        dd += aN * T[i, j - 1]
                    if j < ny - 1:
                        dd += aS * T[i, j + 1]
                    if j == 0:
                        aP += ho * dx; dd += ho * dx * Tsa
                    if j == ny - 1:
                        aP += hi * dx; dd += hi * dx * Tint
                    a[i, j] = aP; b[i, j] = aE; c[i, j] = aW; d[i, j] = dd
        # --- x TDMA per row ---
        for j in range(ny):
            P[0] = b[0, j] / a[0, j]
            Q[0] = d[0, j] / a[0, j]
            for i in range(1, nx):
                denom = a[i, j] - c[i, j] * P[i - 1]
                P[i] = b[i, j] / denom
                Q[i] = (d[i, j] + c[i, j] * Q[i - 1]) / denom
            Tnew[nx - 1, j] = Q[nx - 1]
            for i in range(nx - 2, -1, -1):
                Tnew[i, j] = P[i] * Tnew[i + 1, j] + Q[i]
        # --- error and commit ---
        error = 0.0
        for j in range(ny):
            for i in range(nx):
                error += (T[i, j] - Tnew[i, j]) / T[i, j]
                T[i, j] = Tnew[i, j]
        error = error / nx / ny
        if abs(error) <= tol:
            break
    return iters


@njit(cache=True)
def solve_day_slab_prod(NT, k, rhoc, Tsa_arr, ho, hi, dt, dx, dy, La, X,
                        rhoair, cair, T0, cav_of, cav_i1, cav_i2, cj1, cj2,
                        cavity_width, e22, E, beta,
                        Fud, Ful, Fur, Fru, Frd, Frl, Fdl, Fdr, Fdu, Flu, Flr, Fld,
                        tol_inner, tol_day, max_days):
    """Full day (day-to-day convergence) of the N-air-cavity roof slab. Each cavity
    has its lumped node ``Th[c]``; production conventions (single dt, surfaces
    /(nx-1)). Wall/roof Nusselt according to ``beta``.

    Returns:
        (Ti_series, Tso_series, Tsi_series, Th_mean_series, T_field, days, Qin, Qout)
    """
    nx, ny = k.shape
    nsteps = Tsa_arr.shape[0]
    n_cav = cav_i1.shape[0]
    T = np.empty((nx, ny)); To = np.empty((nx, ny)); Told = np.empty((nx, ny))
    for j in range(ny):
        for i in range(nx):
            T[i, j] = T0
    Tint = T0
    Th = np.empty(n_cav)
    for cidx in range(n_cav):
        Th[cidx] = T0
    a = np.empty((nx, ny)); b = np.empty((nx, ny)); c = np.empty((nx, ny))
    d = np.empty((nx, ny)); Tnew = np.empty((nx, ny))
    P = np.empty(nx); Q = np.empty(nx)
    hh = np.empty(n_cav)
    Qtop = np.empty(n_cav); Qbot = np.empty(n_cav)
    Qleft = np.empty(n_cav); Qright = np.empty(n_cav)
    Ti_s = np.empty(nsteps); Tso_s = np.empty(nsteps)
    Tsi_s = np.empty(nsteps); Th_s = np.empty(nsteps)
    Cair = rhoair * cair * La * X
    Ch = rhoair * cair * cavity_width * e22
    alphaair = _K_AIR / rhoair / cair
    days = 0
    C = 1.0e9
    Qin = Qout = 0.0
    while C > tol_day and days < max_days:
        for j in range(ny):
            for i in range(nx):
                Told[i, j] = T[i, j]
        Qin = Qout = 0.0
        for s in range(nsteps):
            tso = 0.0
            for i in range(nx):
                tso += T[i, 0]
            Tso_s[s] = tso / (nx - 1)
            for j in range(ny):
                for i in range(nx):
                    To[i, j] = T[i, j]
            Ti_old = Tint
            _step_slab(k, rhoc, To, T, Tsa_arr[s], Ti_old, Th, ho, hi, dt, dx, dy,
                       NT, cav_of, cav_i1, cav_i2, cj1, cj2, n_cav, e22, E, beta,
                       _K_AIR, _GR, _BETA_EXP, _NU_AIR, alphaair,
                       Fud, Ful, Fur, Fru, Frd, Frl, Fdl, Fdr, Fdu, Flu, Flr, Fld,
                       a, b, c, d, P, Q, Tnew, hh, Qtop, Qbot, Qleft, Qright, tol_inner)
            # air of each cavity (single dt)
            thsum = 0.0
            for cidx in range(n_cav):
                ci1 = cav_i1[cidx]; ci2 = cav_i2[cidx]; hc = hh[cidx]
                qh = 0.0
                for i in range(ci1, ci2):
                    qh += hc * dx * (T[i, cj1 - 1] - Th[cidx])
                    qh += hc * dx * (T[i, cj2] - Th[cidx])
                for j in range(cj1, cj2):
                    qh += hc * dy * (T[ci1 - 1, j] - Th[cidx])
                    qh += hc * dy * (T[ci2, j] - Th[cidx])
                Th[cidx] = Th[cidx] + dt * qh / Ch
                thsum += Th[cidx]
            Th_s[s] = thsum / n_cav
            # indoor air (single dt)
            flux = 0.0
            for i in range(nx):
                flux += hi * dx * (T[i, ny - 1] - Ti_old)
            Tint = Ti_old + dt * flux / Cair
            Ti_s[s] = Tint
            tsi = 0.0
            for i in range(nx):
                tsi += T[i, ny - 1]
            Tsi_s[s] = tsi / (nx - 1)
            e = flux * dt / X
            if e > 0.0:
                Qin += e
            else:
                Qout -= e
        C = 0.0
        for j in range(ny):
            for i in range(nx):
                C += abs(Told[i, j] - T[i, j])
        C = C / nx / ny
        days += 1
    return Ti_s, Tso_s, Tsi_s, Th_s, T, days, Qin, Qout


def solve_step_slab(NT, k, rhoc, To, Tsa, Tint, Th0, ho, hi, dt, dx, dy,
                    cav_of, cav_i1, cav_i2, cj1, cj2, cavity_width, e22, E, beta,
                    rhoair, cair, La, X, tol=1e-10):
    """One step of the roof slab (for unit testing). Returns a dict with
    ``T, Tint, Thueco (array), hh (array), iters``."""
    k = np.asarray(k, dtype=np.float64)
    rhoc = np.asarray(rhoc, dtype=np.float64)
    To = np.asarray(To, dtype=np.float64)
    NT = np.ascontiguousarray(NT, dtype=np.int64)
    cav_of = np.ascontiguousarray(cav_of, dtype=np.int64)
    cav_i1 = np.ascontiguousarray(cav_i1, dtype=np.int64)
    cav_i2 = np.ascontiguousarray(cav_i2, dtype=np.int64)
    nx, ny = k.shape
    n_cav = cav_i1.shape[0]
    Th = np.full(n_cav, float(Th0))
    a = np.empty((nx, ny)); b = np.empty((nx, ny)); cc = np.empty((nx, ny))
    d = np.empty((nx, ny)); Tnew = np.empty((nx, ny))
    P = np.empty(nx); Q = np.empty(nx)
    hh = np.empty(n_cav)
    Qtop = np.empty(n_cav); Qbot = np.empty(n_cav)
    Qleft = np.empty(n_cav); Qright = np.empty(n_cav)
    vf = _view_factors(cavity_width, e22)
    alphaair = _K_AIR / rhoair / cair
    T = To.copy()
    iters = _step_slab(k, rhoc, To, T, float(Tsa), float(Tint), Th, ho, hi, dt, dx, dy,
                       NT, cav_of, cav_i1, cav_i2, cj1, cj2, n_cav, e22, E, float(beta),
                       _K_AIR, _GR, _BETA_EXP, _NU_AIR, alphaair, *vf,
                       a, b, cc, d, P, Q, Tnew, hh, Qtop, Qbot, Qleft, Qright, tol)
    Ch = rhoair * cair * cavity_width * e22
    for cidx in range(n_cav):
        ci1 = cav_i1[cidx]; ci2 = cav_i2[cidx]; hc = hh[cidx]
        qh = 0.0
        for i in range(ci1, ci2):
            qh += hc * dx * (T[i, cj1 - 1] - Th[cidx])
            qh += hc * dx * (T[i, cj2] - Th[cidx])
        for j in range(cj1, cj2):
            qh += hc * dy * (T[ci1 - 1, j] - Th[cidx])
            qh += hc * dy * (T[ci2, j] - Th[cidx])
        Th[cidx] = Th[cidx] + dt * qh / Ch
    Tsurf = T[:, ny - 1]
    Cair = rhoair * cair * La * X
    flux = float(np.sum(hi * dx * (Tsurf - float(Tint))))
    Tint_new = float(Tint) + dt * flux / Cair
    return {"T": T, "Tint": Tint_new, "Thueco": Th, "hh": hh, "iters": iters}


# =================================================================
#  Parallel variants (numba prange) — optional engine
# =================================================================
#
# The line-by-line sweep is Jacobi (rows use a single snapshot of T per inner
# iteration and are solved into Tnew before T<-Tnew), so the rows are independent:
# parallelising the assembly/TDMA/error over j with prange runs the SAME algorithm
# (same number of iterations, bit-for-bit result, see Phase 7). The per-cavity
# precompute (walls, hh, radiation) is small and left sequential.
# numba auto-detects the cores at runtime (internal threading layer, portable).


@njit(parallel=True, cache=True)
def _step_hueca_par(k, rhoc, To, T, Tsa, Tint, Th, ho, hi, dt, dx, dy,
                    i1, j1, i2, j2, e22, E,
                    Fud, Ful, Fur, Fru, Frd, Frl, Fdl, Fdr, Fdu, Flu, Flr, Fld,
                    a, b, c, d, Tnew, tol):
    """Igual que :func:`_step_hueca` pero con ensamble/TDMA/error en ``prange``."""
    nx, ny = k.shape
    iters = 0
    hh = 1.0
    while True:
        iters += 1
        tup = tdn = tlf = trt = 0.0
        for i in range(i1, i2):
            tup += T[i, j1 - 1]
            tdn += T[i, j2]
        nud = i2 - i1
        for j in range(j1, j2):
            tlf += T[i1 - 1, j]
            trt += T[i2, j]
        nlr = j2 - j1
        tup /= nud; tdn /= nud; tlf /= nlr; trt /= nlr
        hh = 0.4005 * (abs(tup - tdn) ** 0.3033) / (e22 ** 0.0901)
        Tu = tup + 273.15; Td = tdn + 273.15; Tl = tlf + 273.15; Tr = trt + 273.15
        Tu4 = Tu * Tu * Tu * Tu; Td4 = Td * Td * Td * Td
        Tl4 = Tl * Tl * Tl * Tl; Tr4 = Tr * Tr * Tr * Tr
        sx = dx * E * _SIGMA; sy = dy * E * _SIGMA
        Qud = sx * (Tu4 - Td4) * Fud; Qul = sx * (Tu4 - Tl4) * Ful
        Qur = sx * (Tu4 - Tr4) * Fur
        Qru = sy * (Tr4 - Tu4) * Fru; Qrd = sy * (Tr4 - Td4) * Frd
        Qrl = sy * (Tr4 - Tl4) * Frl
        Qdu = sx * (Td4 - Tu4) * Fdu; Qdr = sx * (Td4 - Tr4) * Fdr
        Qdl = sx * (Td4 - Tl4) * Fdl
        Qlu = sy * (Tl4 - Tu4) * Flu; Qlr = sy * (Tl4 - Tr4) * Flr
        Qld = sy * (Tl4 - Td4) * Fld
        for j in prange(ny):
            for i in range(nx):
                kij = k[i, j]
                apo = rhoc[i, j] * dx * dy / dt
                aN = aS = aE = aW = 0.0
                if j > 0:
                    kk = k[i, j - 1]; aN = 2.0 * kk * kij / (kk + kij) * dx / dy
                if j < ny - 1:
                    kk = k[i, j + 1]; aS = 2.0 * kk * kij / (kk + kij) * dx / dy
                if i < nx - 1:
                    kk = k[i + 1, j]; aE = 2.0 * kk * kij / (kk + kij) * dy / dx
                if i > 0:
                    kk = k[i - 1, j]; aW = 2.0 * kk * kij / (kk + kij) * dy / dx
                in_hole_cols = (i1 <= i) and (i < i2)
                in_hole_rows = (j1 <= j) and (j < j2)
                if in_hole_cols and in_hole_rows:
                    a[i, j] = 1.0; b[i, j] = 0.0; c[i, j] = 0.0; d[i, j] = Th
                elif in_hole_cols and j == j1 - 1:
                    a[i, j] = apo + aN + hh * dx + aE + aW
                    d[i, j] = aN * T[i, j - 1] + hh * dx * Th + apo * To[i, j] - Qur - Qud - Qul
                    b[i, j] = aE; c[i, j] = aW
                elif in_hole_cols and j == j2:
                    a[i, j] = apo + hh * dx + aS + aE + aW
                    d[i, j] = aS * T[i, j + 1] + hh * dx * Th + apo * To[i, j] - Qdl - Qdu - Qdr
                    b[i, j] = aE; c[i, j] = aW
                elif in_hole_rows and i == i1 - 1:
                    a[i, j] = apo + aN + aS + hh * dy + aW
                    d[i, j] = aN * T[i, j - 1] + aS * T[i, j + 1] + apo * To[i, j] + hh * dy * Th - Qlu - Qlr - Qld
                    b[i, j] = 0.0; c[i, j] = aW
                elif in_hole_rows and i == i2:
                    a[i, j] = apo + aN + aS + aE + hh * dy
                    d[i, j] = aN * T[i, j - 1] + aS * T[i, j + 1] + apo * To[i, j] + hh * dy * Th - Qrd - Qrl - Qru
                    b[i, j] = aE; c[i, j] = 0.0
                else:
                    aP = apo + aN + aS + aE + aW
                    dd = apo * To[i, j]
                    if j > 0:
                        dd += aN * T[i, j - 1]
                    if j < ny - 1:
                        dd += aS * T[i, j + 1]
                    if j == 0:
                        aP += ho * dx; dd += ho * dx * Tsa
                    if j == ny - 1:
                        aP += hi * dx; dd += hi * dx * Tint
                    a[i, j] = aP; b[i, j] = aE; c[i, j] = aW; d[i, j] = dd
        for j in prange(ny):
            P = np.empty(nx); Q = np.empty(nx)
            P[0] = b[0, j] / a[0, j]; Q[0] = d[0, j] / a[0, j]
            for i in range(1, nx):
                denom = a[i, j] - c[i, j] * P[i - 1]
                P[i] = b[i, j] / denom
                Q[i] = (d[i, j] + c[i, j] * Q[i - 1]) / denom
            Tnew[nx - 1, j] = Q[nx - 1]
            for i in range(nx - 2, -1, -1):
                Tnew[i, j] = P[i] * Tnew[i + 1, j] + Q[i]
        error = 0.0
        for j in prange(ny):
            for i in range(nx):
                error += (T[i, j] - Tnew[i, j]) / T[i, j]
                T[i, j] = Tnew[i, j]
        error = error / nx / ny
        if abs(error) <= tol:
            break
    return iters, hh


@njit(cache=True)
def solve_day_hueca_prod_par(NT, k, rhoc, Tsa_arr, ho, hi, dt, dx, dy, La, X,
                             rhoair, cair, T0, i1, j1, i2, j2, a21, e22, E,
                             Fud, Ful, Fur, Fru, Frd, Frl, Fdl, Fdr, Fdu, Flu, Flr, Fld,
                             tol_inner, tol_day, max_days):
    """Parallel version of :func:`solve_day_hueca_prod` (uses ``_step_hueca_par``)."""
    nx, ny = k.shape
    nsteps = Tsa_arr.shape[0]
    T = np.empty((nx, ny)); To = np.empty((nx, ny)); Told = np.empty((nx, ny))
    for j in range(ny):
        for i in range(nx):
            T[i, j] = T0
    Tint = T0; Th = T0
    a = np.empty((nx, ny)); b = np.empty((nx, ny)); c = np.empty((nx, ny))
    d = np.empty((nx, ny)); Tnew = np.empty((nx, ny))
    Ti_s = np.empty(nsteps); Tso_s = np.empty(nsteps)
    Tsi_s = np.empty(nsteps); Th_s = np.empty(nsteps)
    Cair = rhoair * cair * La * X
    Ch = rhoair * cair * a21 * e22
    days = 0; C = 1.0e9; Qin = Qout = 0.0
    while C > tol_day and days < max_days:
        for j in range(ny):
            for i in range(nx):
                Told[i, j] = T[i, j]
        Qin = Qout = 0.0
        for s in range(nsteps):
            tso = 0.0
            for i in range(nx):
                tso += T[i, 0]
            Tso_s[s] = tso / (nx - 1)
            for j in range(ny):
                for i in range(nx):
                    To[i, j] = T[i, j]
            Ti_old = Tint; Th_old = Th
            _, hh = _step_hueca_par(k, rhoc, To, T, Tsa_arr[s], Ti_old, Th_old,
                                    ho, hi, dt, dx, dy, i1, j1, i2, j2, e22, E,
                                    Fud, Ful, Fur, Fru, Frd, Frl, Fdl, Fdr, Fdu, Flu, Flr, Fld,
                                    a, b, c, d, Tnew, tol_inner)
            qh = 0.0
            for i in range(i1, i2):
                qh += hh * dx * (T[i, j1 - 1] - Th_old)
                qh += hh * dx * (T[i, j2] - Th_old)
            for j in range(j1, j2):
                qh += hh * dy * (T[i1 - 1, j] - Th_old)
                qh += hh * dy * (T[i2, j] - Th_old)
            Th = Th_old + dt * qh / Ch
            flux = 0.0
            for i in range(nx):
                flux += hi * dx * (T[i, ny - 1] - Ti_old)
            Tint = Ti_old + dt * flux / Cair
            Ti_s[s] = Tint; Th_s[s] = Th
            tsi = 0.0
            for i in range(nx):
                tsi += T[i, ny - 1]
            Tsi_s[s] = tsi / (nx - 1)
            e = flux * dt / X
            if e > 0.0:
                Qin += e
            else:
                Qout -= e
        C = 0.0
        for j in range(ny):
            for i in range(nx):
                C += abs(Told[i, j] - T[i, j])
        C = C / nx / ny
        days += 1
    return Ti_s, Tso_s, Tsi_s, Th_s, T, days, Qin, Qout


@njit(parallel=True, cache=True)
def _step_slab_par(k, rhoc, To, T, Tsa, Tint, Th, ho, hi, dt, dx, dy,
                   NT, cav_of, cav_i1, cav_i2, cj1, cj2, n_cav, e22, E, beta,
                   kair, gr, beta_exp, nu, alphaair,
                   Fud, Ful, Fur, Fru, Frd, Frl, Fdl, Fdr, Fdu, Flu, Flr, Fld,
                   a, b, c, d, Tnew, hh, Qtop, Qbot, Qleft, Qright, tol):
    """Igual que :func:`_step_slab` pero con ensamble/TDMA/error en ``prange``."""
    nx, ny = k.shape
    sx = dx * E * _SIGMA; sy = dy * E * _SIGMA
    iters = 0
    while True:
        iters += 1
        for cidx in range(n_cav):
            ci1 = cav_i1[cidx]; ci2 = cav_i2[cidx]
            tup = 0.0; tdn = 0.0; tlf = 0.0; trt = 0.0
            for i in range(ci1, ci2):
                tup += T[i, cj1 - 1]
                tdn += T[i, cj2]
            nud = ci2 - ci1
            for j in range(cj1, cj2):
                tlf += T[ci1 - 1, j]
                trt += T[ci2, j]
            nlr = cj2 - cj1
            tup /= nud; tdn /= nud; tlf /= nlr; trt /= nlr
            hh[cidx] = _slab_hh(tup, tdn, e22, beta, kair, gr, beta_exp, nu, alphaair)
            Tu = tup + 273.15; Td = tdn + 273.15
            Tl = tlf + 273.15; Tr = trt + 273.15
            Tu4 = Tu * Tu * Tu * Tu; Td4 = Td * Td * Td * Td
            Tl4 = Tl * Tl * Tl * Tl; Tr4 = Tr * Tr * Tr * Tr
            Qud = sx * (Tu4 - Td4) * Fud; Qul = sx * (Tu4 - Tl4) * Ful
            Qur = sx * (Tu4 - Tr4) * Fur
            Qru = sy * (Tr4 - Tu4) * Fru; Qrd = sy * (Tr4 - Td4) * Frd
            Qrl = sy * (Tr4 - Tl4) * Frl
            Qdu = sx * (Td4 - Tu4) * Fdu; Qdr = sx * (Td4 - Tr4) * Fdr
            Qdl = sx * (Td4 - Tl4) * Fdl
            Qlu = sy * (Tl4 - Tu4) * Flu; Qlr = sy * (Tl4 - Tr4) * Flr
            Qld = sy * (Tl4 - Td4) * Fld
            Qtop[cidx] = Qur + Qud + Qul
            Qbot[cidx] = Qdl + Qdu + Qdr
            Qleft[cidx] = Qlu + Qlr + Qld
            Qright[cidx] = Qrd + Qrl + Qru
        for j in prange(ny):
            for i in range(nx):
                nt = NT[i, j]
                kij = k[i, j]
                apo = rhoc[i, j] * dx * dy / dt
                aN = aS = aE = aW = 0.0
                if j > 0:
                    kk = k[i, j - 1]; aN = 2.0 * kk * kij / (kk + kij) * dx / dy
                if j < ny - 1:
                    kk = k[i, j + 1]; aS = 2.0 * kk * kij / (kk + kij) * dx / dy
                if i < nx - 1:
                    kk = k[i + 1, j]; aE = 2.0 * kk * kij / (kk + kij) * dy / dx
                if i > 0:
                    kk = k[i - 1, j]; aW = 2.0 * kk * kij / (kk + kij) * dy / dx
                if nt == 0:
                    cc = cav_of[i, j]
                    a[i, j] = 1.0; b[i, j] = 0.0; c[i, j] = 0.0; d[i, j] = Th[cc]
                elif nt == 9:
                    cc = cav_of[i, j]; hc = hh[cc]
                    a[i, j] = apo + aN + hc * dx + aE + aW
                    d[i, j] = aN * T[i, j - 1] + hc * dx * Th[cc] + apo * To[i, j] - Qtop[cc]
                    b[i, j] = aE; c[i, j] = aW
                elif nt == 10:
                    cc = cav_of[i, j]; hc = hh[cc]
                    a[i, j] = apo + hc * dx + aS + aE + aW
                    d[i, j] = aS * T[i, j + 1] + hc * dx * Th[cc] + apo * To[i, j] - Qbot[cc]
                    b[i, j] = aE; c[i, j] = aW
                elif nt == 11:
                    cc = cav_of[i, j]; hc = hh[cc]
                    a[i, j] = apo + aN + aS + hc * dy + aW
                    d[i, j] = aN * T[i, j - 1] + aS * T[i, j + 1] + apo * To[i, j] + hc * dy * Th[cc] - Qleft[cc]
                    b[i, j] = 0.0; c[i, j] = aW
                elif nt == 12:
                    cc = cav_of[i, j]; hc = hh[cc]
                    a[i, j] = apo + aN + aS + aE + hc * dy
                    d[i, j] = aN * T[i, j - 1] + aS * T[i, j + 1] + apo * To[i, j] + hc * dy * Th[cc] - Qright[cc]
                    b[i, j] = aE; c[i, j] = 0.0
                else:
                    aP = apo + aN + aS + aE + aW
                    dd = apo * To[i, j]
                    if j > 0:
                        dd += aN * T[i, j - 1]
                    if j < ny - 1:
                        dd += aS * T[i, j + 1]
                    if j == 0:
                        aP += ho * dx; dd += ho * dx * Tsa
                    if j == ny - 1:
                        aP += hi * dx; dd += hi * dx * Tint
                    a[i, j] = aP; b[i, j] = aE; c[i, j] = aW; d[i, j] = dd
        for j in prange(ny):
            P = np.empty(nx); Q = np.empty(nx)
            P[0] = b[0, j] / a[0, j]; Q[0] = d[0, j] / a[0, j]
            for i in range(1, nx):
                denom = a[i, j] - c[i, j] * P[i - 1]
                P[i] = b[i, j] / denom
                Q[i] = (d[i, j] + c[i, j] * Q[i - 1]) / denom
            Tnew[nx - 1, j] = Q[nx - 1]
            for i in range(nx - 2, -1, -1):
                Tnew[i, j] = P[i] * Tnew[i + 1, j] + Q[i]
        error = 0.0
        for j in prange(ny):
            for i in range(nx):
                error += (T[i, j] - Tnew[i, j]) / T[i, j]
                T[i, j] = Tnew[i, j]
        error = error / nx / ny
        if abs(error) <= tol:
            break
    return iters


@njit(cache=True)
def solve_day_slab_prod_par(NT, k, rhoc, Tsa_arr, ho, hi, dt, dx, dy, La, X,
                            rhoair, cair, T0, cav_of, cav_i1, cav_i2, cj1, cj2,
                            cavity_width, e22, E, beta,
                            Fud, Ful, Fur, Fru, Frd, Frl, Fdl, Fdr, Fdu, Flu, Flr, Fld,
                            tol_inner, tol_day, max_days):
    """Parallel version of :func:`solve_day_slab_prod` (uses ``_step_slab_par``)."""
    nx, ny = k.shape
    nsteps = Tsa_arr.shape[0]
    n_cav = cav_i1.shape[0]
    T = np.empty((nx, ny)); To = np.empty((nx, ny)); Told = np.empty((nx, ny))
    for j in range(ny):
        for i in range(nx):
            T[i, j] = T0
    Tint = T0
    Th = np.empty(n_cav)
    for cidx in range(n_cav):
        Th[cidx] = T0
    a = np.empty((nx, ny)); b = np.empty((nx, ny)); c = np.empty((nx, ny))
    d = np.empty((nx, ny)); Tnew = np.empty((nx, ny))
    hh = np.empty(n_cav)
    Qtop = np.empty(n_cav); Qbot = np.empty(n_cav)
    Qleft = np.empty(n_cav); Qright = np.empty(n_cav)
    Ti_s = np.empty(nsteps); Tso_s = np.empty(nsteps)
    Tsi_s = np.empty(nsteps); Th_s = np.empty(nsteps)
    Cair = rhoair * cair * La * X
    Ch = rhoair * cair * cavity_width * e22
    alphaair = _K_AIR / rhoair / cair
    days = 0; C = 1.0e9; Qin = Qout = 0.0
    while C > tol_day and days < max_days:
        for j in range(ny):
            for i in range(nx):
                Told[i, j] = T[i, j]
        Qin = Qout = 0.0
        for s in range(nsteps):
            tso = 0.0
            for i in range(nx):
                tso += T[i, 0]
            Tso_s[s] = tso / (nx - 1)
            for j in range(ny):
                for i in range(nx):
                    To[i, j] = T[i, j]
            Ti_old = Tint
            _step_slab_par(k, rhoc, To, T, Tsa_arr[s], Ti_old, Th, ho, hi, dt, dx, dy,
                           NT, cav_of, cav_i1, cav_i2, cj1, cj2, n_cav, e22, E, beta,
                           _K_AIR, _GR, _BETA_EXP, _NU_AIR, alphaair,
                           Fud, Ful, Fur, Fru, Frd, Frl, Fdl, Fdr, Fdu, Flu, Flr, Fld,
                           a, b, c, d, Tnew, hh, Qtop, Qbot, Qleft, Qright, tol_inner)
            thsum = 0.0
            for cidx in range(n_cav):
                ci1 = cav_i1[cidx]; ci2 = cav_i2[cidx]; hc = hh[cidx]
                qh = 0.0
                for i in range(ci1, ci2):
                    qh += hc * dx * (T[i, cj1 - 1] - Th[cidx])
                    qh += hc * dx * (T[i, cj2] - Th[cidx])
                for j in range(cj1, cj2):
                    qh += hc * dy * (T[ci1 - 1, j] - Th[cidx])
                    qh += hc * dy * (T[ci2, j] - Th[cidx])
                Th[cidx] = Th[cidx] + dt * qh / Ch
                thsum += Th[cidx]
            Th_s[s] = thsum / n_cav
            flux = 0.0
            for i in range(nx):
                flux += hi * dx * (T[i, ny - 1] - Ti_old)
            Tint = Ti_old + dt * flux / Cair
            Ti_s[s] = Tint
            tsi = 0.0
            for i in range(nx):
                tsi += T[i, ny - 1]
            Tsi_s[s] = tsi / (nx - 1)
            e = flux * dt / X
            if e > 0.0:
                Qin += e
            else:
                Qout -= e
        C = 0.0
        for j in range(ny):
            for i in range(nx):
                C += abs(Told[i, j] - T[i, j])
        C = C / nx / ny
        days += 1
    return Ti_s, Tso_s, Tsi_s, Th_s, T, days, Qin, Qout


# =================================================================
#  Air conditioning (AC) — Phase 9
# =================================================================
#
# Holds the indoor air FIXED at a setpoint `Tset` (not integrated) and accumulates
# the load: per step, the net flux at the inner surface `e=(Σ hi·dx·(T−Tset))·dt/X`
# is added to `Qcool` if heat comes in (e>0) or to `Qheat` if it goes out (e<0). In
# the air variants (hueca/slab) the cavity air `Th`/`Th[c]` KEEPS floating (the AC
# only controls the indoor air). They reuse the same `_step_*` as the free-running.


@njit(cache=True)
def solve_day_2d_ac(NT, k, rhoc, Tsa_arr, ho, hi, dt, dx, dy, La, X,
                    rhoair, cair, T0, Tset, tol_inner=1e-10, tol_day=5e-4, max_days=60):
    """AC for pure conduction (SOLID). Returns
    (Ti_series(=Tset), Tso, Tsi, T_field, days, Qcool, Qheat)."""
    nx, ny = k.shape
    nsteps = Tsa_arr.shape[0]
    T = np.empty((nx, ny)); To = np.empty((nx, ny)); Told = np.empty((nx, ny))
    for j in range(ny):
        for i in range(nx):
            T[i, j] = T0
    a = np.empty((nx, ny)); b = np.empty((nx, ny))
    c = np.empty((nx, ny)); d = np.empty((nx, ny)); Tnew = np.empty((nx, ny))
    P = np.empty(nx); Q = np.empty(nx); Tn = np.empty(nx)
    Ti_s = np.empty(nsteps); Tso_s = np.empty(nsteps); Tsi_s = np.empty(nsteps)
    days = 0; Qcool = Qheat = 0.0; C = 1.0e9
    while C > tol_day and days < max_days:
        for j in range(ny):
            for i in range(nx):
                Told[i, j] = T[i, j]
        Qcool = Qheat = 0.0
        for s in range(nsteps):
            tso = 0.0
            for i in range(nx):
                tso += T[i, 0]
            Tso_s[s] = tso / (nx - 1)
            for j in range(ny):
                for i in range(nx):
                    To[i, j] = T[i, j]
            _step_inner(k, rhoc, To, T, Tsa_arr[s], Tset, ho, hi, dt, dx, dy,
                        a, b, c, d, P, Q, Tn, Tnew, tol_inner)
            flux = 0.0
            for i in range(nx):
                flux += hi * dx * (T[i, ny - 1] - Tset)
            e = flux * dt / X
            if e > 0.0:
                Qcool += e
            else:
                Qheat -= e
            Ti_s[s] = Tset
            tsi = 0.0
            for i in range(nx):
                tsi += T[i, ny - 1]
            Tsi_s[s] = tsi / (nx - 1)
        C = 0.0
        for j in range(ny):
            for i in range(nx):
                C += abs(Told[i, j] - T[i, j])
        C = C / nx / ny
        days += 1
    return Ti_s, Tso_s, Tsi_s, T, days, Qcool, Qheat


@njit(cache=True)
def solve_day_2d_ac_par(NT, k, rhoc, Tsa_arr, ho, hi, dt, dx, dy, La, X,
                        rhoair, cair, T0, Tset, tol_inner=1e-10, tol_day=5e-4, max_days=60):
    """Parallel version of :func:`solve_day_2d_ac` (uses ``_step_inner_par``)."""
    nx, ny = k.shape
    nsteps = Tsa_arr.shape[0]
    T = np.empty((nx, ny)); To = np.empty((nx, ny)); Told = np.empty((nx, ny))
    for j in range(ny):
        for i in range(nx):
            T[i, j] = T0
    a = np.empty((nx, ny)); b = np.empty((nx, ny))
    c = np.empty((nx, ny)); d = np.empty((nx, ny)); Tnew = np.empty((nx, ny))
    Ti_s = np.empty(nsteps); Tso_s = np.empty(nsteps); Tsi_s = np.empty(nsteps)
    days = 0; Qcool = Qheat = 0.0; C = 1.0e9
    while C > tol_day and days < max_days:
        for j in range(ny):
            for i in range(nx):
                Told[i, j] = T[i, j]
        Qcool = Qheat = 0.0
        for s in range(nsteps):
            tso = 0.0
            for i in range(nx):
                tso += T[i, 0]
            Tso_s[s] = tso / (nx - 1)
            for j in range(ny):
                for i in range(nx):
                    To[i, j] = T[i, j]
            _step_inner_par(k, rhoc, To, T, Tsa_arr[s], Tset, ho, hi, dt, dx, dy,
                            a, b, c, d, Tnew, tol_inner)
            flux = 0.0
            for i in range(nx):
                flux += hi * dx * (T[i, ny - 1] - Tset)
            e = flux * dt / X
            if e > 0.0:
                Qcool += e
            else:
                Qheat -= e
            Ti_s[s] = Tset
            tsi = 0.0
            for i in range(nx):
                tsi += T[i, ny - 1]
            Tsi_s[s] = tsi / (nx - 1)
        C = 0.0
        for j in range(ny):
            for i in range(nx):
                C += abs(Told[i, j] - T[i, j])
        C = C / nx / ny
        days += 1
    return Ti_s, Tso_s, Tsi_s, T, days, Qcool, Qheat


@njit(cache=True)
def solve_day_hueca_ac(NT, k, rhoc, Tsa_arr, ho, hi, dt, dx, dy, La, X,
                       rhoair, cair, T0, Tset, i1, j1, i2, j2, a21, e22, E,
                       Fud, Ful, Fur, Fru, Frd, Frl, Fdl, Fdr, Fdu, Flu, Flr, Fld,
                       tol_inner, tol_day, max_days):
    """AC for a wall with air cavity (Thueco floats, Tint=Tset fixed). Returns
    (Ti(=Tset), Tso, Tsi, Th, T_field, days, Qcool, Qheat)."""
    nx, ny = k.shape
    nsteps = Tsa_arr.shape[0]
    T = np.empty((nx, ny)); To = np.empty((nx, ny)); Told = np.empty((nx, ny))
    for j in range(ny):
        for i in range(nx):
            T[i, j] = T0
    Th = T0
    a = np.empty((nx, ny)); b = np.empty((nx, ny)); c = np.empty((nx, ny))
    d = np.empty((nx, ny)); Tnew = np.empty((nx, ny))
    P = np.empty(nx); Q = np.empty(nx)
    Ti_s = np.empty(nsteps); Tso_s = np.empty(nsteps)
    Tsi_s = np.empty(nsteps); Th_s = np.empty(nsteps)
    Ch = rhoair * cair * a21 * e22
    days = 0; C = 1.0e9; Qcool = Qheat = 0.0
    while C > tol_day and days < max_days:
        for j in range(ny):
            for i in range(nx):
                Told[i, j] = T[i, j]
        Qcool = Qheat = 0.0
        for s in range(nsteps):
            tso = 0.0
            for i in range(nx):
                tso += T[i, 0]
            Tso_s[s] = tso / (nx - 1)
            for j in range(ny):
                for i in range(nx):
                    To[i, j] = T[i, j]
            Th_old = Th
            _, hh = _step_hueca(k, rhoc, To, T, Tsa_arr[s], Tset, Th_old,
                                ho, hi, dt, dx, dy, i1, j1, i2, j2, e22, E,
                                Fud, Ful, Fur, Fru, Frd, Frl, Fdl, Fdr, Fdu, Flu, Flr, Fld,
                                a, b, c, d, P, Q, Tnew, tol_inner)
            qh = 0.0
            for i in range(i1, i2):
                qh += hh * dx * (T[i, j1 - 1] - Th_old)
                qh += hh * dx * (T[i, j2] - Th_old)
            for j in range(j1, j2):
                qh += hh * dy * (T[i1 - 1, j] - Th_old)
                qh += hh * dy * (T[i2, j] - Th_old)
            Th = Th_old + dt * qh / Ch
            flux = 0.0
            for i in range(nx):
                flux += hi * dx * (T[i, ny - 1] - Tset)
            e = flux * dt / X
            if e > 0.0:
                Qcool += e
            else:
                Qheat -= e
            Ti_s[s] = Tset; Th_s[s] = Th
            tsi = 0.0
            for i in range(nx):
                tsi += T[i, ny - 1]
            Tsi_s[s] = tsi / (nx - 1)
        C = 0.0
        for j in range(ny):
            for i in range(nx):
                C += abs(Told[i, j] - T[i, j])
        C = C / nx / ny
        days += 1
    return Ti_s, Tso_s, Tsi_s, Th_s, T, days, Qcool, Qheat


@njit(cache=True)
def solve_day_hueca_ac_par(NT, k, rhoc, Tsa_arr, ho, hi, dt, dx, dy, La, X,
                           rhoair, cair, T0, Tset, i1, j1, i2, j2, a21, e22, E,
                           Fud, Ful, Fur, Fru, Frd, Frl, Fdl, Fdr, Fdu, Flu, Flr, Fld,
                           tol_inner, tol_day, max_days):
    """Parallel version of :func:`solve_day_hueca_ac` (uses ``_step_hueca_par``)."""
    nx, ny = k.shape
    nsteps = Tsa_arr.shape[0]
    T = np.empty((nx, ny)); To = np.empty((nx, ny)); Told = np.empty((nx, ny))
    for j in range(ny):
        for i in range(nx):
            T[i, j] = T0
    Th = T0
    a = np.empty((nx, ny)); b = np.empty((nx, ny)); c = np.empty((nx, ny))
    d = np.empty((nx, ny)); Tnew = np.empty((nx, ny))
    Ti_s = np.empty(nsteps); Tso_s = np.empty(nsteps)
    Tsi_s = np.empty(nsteps); Th_s = np.empty(nsteps)
    Ch = rhoair * cair * a21 * e22
    days = 0; C = 1.0e9; Qcool = Qheat = 0.0
    while C > tol_day and days < max_days:
        for j in range(ny):
            for i in range(nx):
                Told[i, j] = T[i, j]
        Qcool = Qheat = 0.0
        for s in range(nsteps):
            tso = 0.0
            for i in range(nx):
                tso += T[i, 0]
            Tso_s[s] = tso / (nx - 1)
            for j in range(ny):
                for i in range(nx):
                    To[i, j] = T[i, j]
            Th_old = Th
            _, hh = _step_hueca_par(k, rhoc, To, T, Tsa_arr[s], Tset, Th_old,
                                    ho, hi, dt, dx, dy, i1, j1, i2, j2, e22, E,
                                    Fud, Ful, Fur, Fru, Frd, Frl, Fdl, Fdr, Fdu, Flu, Flr, Fld,
                                    a, b, c, d, Tnew, tol_inner)
            qh = 0.0
            for i in range(i1, i2):
                qh += hh * dx * (T[i, j1 - 1] - Th_old)
                qh += hh * dx * (T[i, j2] - Th_old)
            for j in range(j1, j2):
                qh += hh * dy * (T[i1 - 1, j] - Th_old)
                qh += hh * dy * (T[i2, j] - Th_old)
            Th = Th_old + dt * qh / Ch
            flux = 0.0
            for i in range(nx):
                flux += hi * dx * (T[i, ny - 1] - Tset)
            e = flux * dt / X
            if e > 0.0:
                Qcool += e
            else:
                Qheat -= e
            Ti_s[s] = Tset; Th_s[s] = Th
            tsi = 0.0
            for i in range(nx):
                tsi += T[i, ny - 1]
            Tsi_s[s] = tsi / (nx - 1)
        C = 0.0
        for j in range(ny):
            for i in range(nx):
                C += abs(Told[i, j] - T[i, j])
        C = C / nx / ny
        days += 1
    return Ti_s, Tso_s, Tsi_s, Th_s, T, days, Qcool, Qheat


@njit(cache=True)
def solve_day_slab_ac(NT, k, rhoc, Tsa_arr, ho, hi, dt, dx, dy, La, X,
                      rhoair, cair, T0, Tset, cav_of, cav_i1, cav_i2, cj1, cj2,
                      cavity_width, e22, E, beta,
                      Fud, Ful, Fur, Fru, Frd, Frl, Fdl, Fdr, Fdu, Flu, Flr, Fld,
                      tol_inner, tol_day, max_days):
    """AC for an N-cavity roof (Thueco[c] float, Tint=Tset fixed). Returns
    (Ti(=Tset), Tso, Tsi, Th_mean, T_field, days, Qcool, Qheat)."""
    nx, ny = k.shape
    nsteps = Tsa_arr.shape[0]
    n_cav = cav_i1.shape[0]
    T = np.empty((nx, ny)); To = np.empty((nx, ny)); Told = np.empty((nx, ny))
    for j in range(ny):
        for i in range(nx):
            T[i, j] = T0
    Th = np.empty(n_cav)
    for cidx in range(n_cav):
        Th[cidx] = T0
    a = np.empty((nx, ny)); b = np.empty((nx, ny)); c = np.empty((nx, ny))
    d = np.empty((nx, ny)); Tnew = np.empty((nx, ny))
    P = np.empty(nx); Q = np.empty(nx)
    hh = np.empty(n_cav); Qtop = np.empty(n_cav); Qbot = np.empty(n_cav)
    Qleft = np.empty(n_cav); Qright = np.empty(n_cav)
    Ti_s = np.empty(nsteps); Tso_s = np.empty(nsteps)
    Tsi_s = np.empty(nsteps); Th_s = np.empty(nsteps)
    Ch = rhoair * cair * cavity_width * e22
    alphaair = _K_AIR / rhoair / cair
    days = 0; C = 1.0e9; Qcool = Qheat = 0.0
    while C > tol_day and days < max_days:
        for j in range(ny):
            for i in range(nx):
                Told[i, j] = T[i, j]
        Qcool = Qheat = 0.0
        for s in range(nsteps):
            tso = 0.0
            for i in range(nx):
                tso += T[i, 0]
            Tso_s[s] = tso / (nx - 1)
            for j in range(ny):
                for i in range(nx):
                    To[i, j] = T[i, j]
            _step_slab(k, rhoc, To, T, Tsa_arr[s], Tset, Th, ho, hi, dt, dx, dy,
                       NT, cav_of, cav_i1, cav_i2, cj1, cj2, n_cav, e22, E, beta,
                       _K_AIR, _GR, _BETA_EXP, _NU_AIR, alphaair,
                       Fud, Ful, Fur, Fru, Frd, Frl, Fdl, Fdr, Fdu, Flu, Flr, Fld,
                       a, b, c, d, P, Q, Tnew, hh, Qtop, Qbot, Qleft, Qright, tol_inner)
            thsum = 0.0
            for cidx in range(n_cav):
                ci1 = cav_i1[cidx]; ci2 = cav_i2[cidx]; hc = hh[cidx]
                qh = 0.0
                for i in range(ci1, ci2):
                    qh += hc * dx * (T[i, cj1 - 1] - Th[cidx])
                    qh += hc * dx * (T[i, cj2] - Th[cidx])
                for j in range(cj1, cj2):
                    qh += hc * dy * (T[ci1 - 1, j] - Th[cidx])
                    qh += hc * dy * (T[ci2, j] - Th[cidx])
                Th[cidx] = Th[cidx] + dt * qh / Ch
                thsum += Th[cidx]
            Th_s[s] = thsum / n_cav
            flux = 0.0
            for i in range(nx):
                flux += hi * dx * (T[i, ny - 1] - Tset)
            e = flux * dt / X
            if e > 0.0:
                Qcool += e
            else:
                Qheat -= e
            Ti_s[s] = Tset
            tsi = 0.0
            for i in range(nx):
                tsi += T[i, ny - 1]
            Tsi_s[s] = tsi / (nx - 1)
        C = 0.0
        for j in range(ny):
            for i in range(nx):
                C += abs(Told[i, j] - T[i, j])
        C = C / nx / ny
        days += 1
    return Ti_s, Tso_s, Tsi_s, Th_s, T, days, Qcool, Qheat


@njit(cache=True)
def solve_day_slab_ac_par(NT, k, rhoc, Tsa_arr, ho, hi, dt, dx, dy, La, X,
                          rhoair, cair, T0, Tset, cav_of, cav_i1, cav_i2, cj1, cj2,
                          cavity_width, e22, E, beta,
                          Fud, Ful, Fur, Fru, Frd, Frl, Fdl, Fdr, Fdu, Flu, Flr, Fld,
                          tol_inner, tol_day, max_days):
    """Parallel version of :func:`solve_day_slab_ac` (uses ``_step_slab_par``)."""
    nx, ny = k.shape
    nsteps = Tsa_arr.shape[0]
    n_cav = cav_i1.shape[0]
    T = np.empty((nx, ny)); To = np.empty((nx, ny)); Told = np.empty((nx, ny))
    for j in range(ny):
        for i in range(nx):
            T[i, j] = T0
    Th = np.empty(n_cav)
    for cidx in range(n_cav):
        Th[cidx] = T0
    a = np.empty((nx, ny)); b = np.empty((nx, ny)); c = np.empty((nx, ny))
    d = np.empty((nx, ny)); Tnew = np.empty((nx, ny))
    hh = np.empty(n_cav); Qtop = np.empty(n_cav); Qbot = np.empty(n_cav)
    Qleft = np.empty(n_cav); Qright = np.empty(n_cav)
    Ti_s = np.empty(nsteps); Tso_s = np.empty(nsteps)
    Tsi_s = np.empty(nsteps); Th_s = np.empty(nsteps)
    Ch = rhoair * cair * cavity_width * e22
    alphaair = _K_AIR / rhoair / cair
    days = 0; C = 1.0e9; Qcool = Qheat = 0.0
    while C > tol_day and days < max_days:
        for j in range(ny):
            for i in range(nx):
                Told[i, j] = T[i, j]
        Qcool = Qheat = 0.0
        for s in range(nsteps):
            tso = 0.0
            for i in range(nx):
                tso += T[i, 0]
            Tso_s[s] = tso / (nx - 1)
            for j in range(ny):
                for i in range(nx):
                    To[i, j] = T[i, j]
            _step_slab_par(k, rhoc, To, T, Tsa_arr[s], Tset, Th, ho, hi, dt, dx, dy,
                           NT, cav_of, cav_i1, cav_i2, cj1, cj2, n_cav, e22, E, beta,
                           _K_AIR, _GR, _BETA_EXP, _NU_AIR, alphaair,
                           Fud, Ful, Fur, Fru, Frd, Frl, Fdl, Fdr, Fdu, Flu, Flr, Fld,
                           a, b, c, d, Tnew, hh, Qtop, Qbot, Qleft, Qright, tol_inner)
            thsum = 0.0
            for cidx in range(n_cav):
                ci1 = cav_i1[cidx]; ci2 = cav_i2[cidx]; hc = hh[cidx]
                qh = 0.0
                for i in range(ci1, ci2):
                    qh += hc * dx * (T[i, cj1 - 1] - Th[cidx])
                    qh += hc * dx * (T[i, cj2] - Th[cidx])
                for j in range(cj1, cj2):
                    qh += hc * dy * (T[ci1 - 1, j] - Th[cidx])
                    qh += hc * dy * (T[ci2, j] - Th[cidx])
                Th[cidx] = Th[cidx] + dt * qh / Ch
                thsum += Th[cidx]
            Th_s[s] = thsum / n_cav
            flux = 0.0
            for i in range(nx):
                flux += hi * dx * (T[i, ny - 1] - Tset)
            e = flux * dt / X
            if e > 0.0:
                Qcool += e
            else:
                Qheat -= e
            Ti_s[s] = Tset
            tsi = 0.0
            for i in range(nx):
                tsi += T[i, ny - 1]
            Tsi_s[s] = tsi / (nx - 1)
        C = 0.0
        for j in range(ny):
            for i in range(nx):
                C += abs(Told[i, j] - T[i, j])
        C = C / nx / ny
        days += 1
    return Ti_s, Tso_s, Tsi_s, Th_s, T, days, Qcool, Qheat
