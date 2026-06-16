"""
=================================================================
 ehtools2d — Kernels 2D (Fase 2: ensamble de coeficientes)
=================================================================

Port fiel de ``calculate_coefficients`` del C (`legacy_eh/2dTfree/tools.h`) para
**bovedilla rellena** (``tipo 2``: nodos ``NT`` 1-8 y 13). Los nodos de cámara de
aire (0, 9-12) llegan en la Fase 6.

Esquema de volumen finito implícito, 5 puntos, conductividad por **media armónica**
``kh(a,b)=2ab/(a+b)``. Para cada nodo ``P=(i,j)``:

    aP·T_P = aE·T_E + aW·T_W + aN·T_N + aS·T_S + apo·T_P° + S_b

con vecinos en y (N=j-1 exterior, S=j+1 interior) **diferidos** al término ``d``
(Gauss-Seidel por líneas; el TDMA implícito solo actúa en x). El solver arma
``a=aP``, ``b=aE``, ``c=aW`` y ``d = aN·T_N + aS·T_S + apo·To + S_b``.

Convención de malla idéntica a `eh2d`: matrices ``(nx,ny)`` indexadas ``[i,j]``,
``j=0`` exterior (``Tsa,ho``), ``j=ny-1`` interior (``Tint,hi``); laterales
``i=0``/``i=nx-1`` adiabáticos.
"""

import numpy as np
from numba import njit, prange

# Tipos de nodo válidos para bovedilla rellena.
_RELLENA_TYPES = {1, 2, 3, 4, 5, 6, 7, 8, 13}


def _harmonic_faces(k, dx, dy):
    """
    Conductancias de cara por media armónica, ya escaladas por la geometría.
    Devuelve ``(aN, aS, aE, aW)`` con forma ``(nx,ny)``; cada una es 0 en la
    frontera donde el vecino no existe (lo que coincide con los coeficientes que
    el C pone a cero en esas fronteras).
    """
    aN = np.zeros_like(k)
    aS = np.zeros_like(k)
    aE = np.zeros_like(k)
    aW = np.zeros_like(k)

    def kh(a, b):
        return 2.0 * a * b / (a + b)

    # Norte (j-1, hacia exterior) y Sur (j+1, hacia interior): factor dx/dy.
    aN[:, 1:] = kh(k[:, 1:], k[:, :-1]) * dx / dy
    aS[:, :-1] = kh(k[:, :-1], k[:, 1:]) * dx / dy
    # Este (i+1) y Oeste (i-1): factor dy/dx.
    aE[:-1, :] = kh(k[:-1, :], k[1:, :]) * dy / dx
    aW[1:, :] = kh(k[1:, :], k[:-1, :]) * dy / dx
    return aN, aS, aE, aW


def calculate_coefficients_2d(NT, k, rhoc, To, T, Tsa, Tint, ho, hi, dt, dx, dy):
    """
    Ensamble de ``a,b,c,d`` para bovedilla rellena (``tipo 2``).

    Args:
        NT (int array (nx,ny)): tipos de nodo (debe ser ⊆ {1-8,13}).
        k, rhoc (float array (nx,ny)): conductividad y capacidad térmica por nodo.
        To (float array (nx,ny)): campo del paso de tiempo anterior (``T°``).
        T  (float array (nx,ny)): campo actual (para los vecinos diferidos N/S).
        Tsa, Tint (float): temperatura sol-aire (exterior) y del aire interior.
        ho, hi (float): coef. convectivos exterior e interior.
        dt, dx, dy (float): paso temporal y de malla.

    Returns:
        (a, b, c, d) arrays (nx,ny). ``a=aP``, ``b=aE``, ``c=aW``,
        ``d = aN·T_N + aS·T_S + apo·To + S_b``.
    """
    NT = np.asarray(NT)
    bad = set(np.unique(NT).tolist()) - _RELLENA_TYPES
    if bad:
        raise NotImplementedError(
            f"calculate_coefficients_2d solo soporta bovedilla rellena "
            f"(NT ⊆ {sorted(_RELLENA_TYPES)}); aparecieron {sorted(bad)}.")

    k = np.asarray(k, dtype=np.float64)
    rhoc = np.asarray(rhoc, dtype=np.float64)
    To = np.asarray(To, dtype=np.float64)
    T = np.asarray(T, dtype=np.float64)
    nx, ny = k.shape

    apo = rhoc * dx * dy / dt
    aN, aS, aE, aW = _harmonic_faces(k, dx, dy)

    # Vecinos diferidos en y (en las fronteras el coef. ya es 0 → término nulo).
    TN = np.zeros_like(T)
    TS = np.zeros_like(T)
    TN[:, 1:] = T[:, :-1]   # T_N = T[i, j-1]
    TS[:, :-1] = T[:, 1:]   # T_S = T[i, j+1]

    # aP = apo + (conducciones activas) ; las inactivas ya valen 0.
    a = apo + aN + aS + aE + aW
    d = apo * To + aN * TN + aS * TS

    # Fronteras convectivas: exterior en j=0 (Tsa, ho), interior en j=ny-1 (Tint, hi).
    a[:, 0] += ho * dx
    d[:, 0] += ho * dx * Tsa
    a[:, ny - 1] += hi * dx
    d[:, ny - 1] += hi * dx * Tint

    b = aE.copy()   # = aE (0 en i=nx-1, fronteras der. adiabáticas/types 2,4,7)
    c = aW.copy()   # = aW (0 en i=0,    fronteras izq. adiabáticas/types 1,3,6)

    return a, b, c, d


def _tdma_rows(a, b, c, d):
    """
    TDMA (Thomas) implícito en x para **todas** las filas ``j`` a la vez.
    ``a·T_i = b·T_{i+1} + c·T_{i-1} + d`` (b=aE, c=aW). Vectorizado sobre j;
    el barrido en i replica exactamente el orden de operaciones del C.
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
    Port de ``solve_PQ`` para bovedilla rellena: **un** paso de tiempo.

    Lazo interno (Gauss-Seidel por líneas): recalcula coeficientes con la `T` más
    reciente, resuelve cada fila con TDMA en x (vecinos en y diferidos), y repite
    hasta ``|error| <= tol`` con ``error = Σ (T-Tn)/T /(nx·ny)`` (con signo, como el
    C). Luego actualiza el aire interior `Tint` (nodo lumped).

    Args:
        To: campo del paso anterior; también es la condición inicial del paso
            (en el C, al iniciar el paso ``To == T``).
        Tint: temperatura del aire interior al inicio del paso.
        La, X, rhoair, cair: parámetros del nodo de aire interior.

    Returns:
        dict con ``T`` (campo resuelto), ``Tint`` (actualizado), ``iters``,
        ``Qin``, ``error``.
    """
    To = np.asarray(To, dtype=np.float64)
    nx, ny = To.shape
    Ti = float(Tint)            # Tint "viejo", constante durante el lazo interno
    T = To.copy()               # condición inicial = campo previo

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

    # Flujo convectivo interior y actualización del aire interior (nodo lumped).
    Tsurf = T[:, ny - 1]
    Qh = float(np.sum(hi * dt * dx * (Tsurf - Ti)))
    Qin = float(np.sum(np.where(Tsurf > Ti, hi * dx * (Tsurf - Ti), 0.0)))
    Cair = rhoair * cair * La * X
    Tint_new = (Qh + (Cair / dt) * Ti) * dt / Cair

    return {"T": T, "Tint": Tint_new, "iters": iters, "Qin": Qin, "error": error}


# =================================================================
#  Motor de producción (JIT numba) — Fase 5
# =================================================================
#
# Solver del día completo con convergencia día-a-día, kernels compilados.
# Diferencias CONSCIENTES respecto al port fiel (Fase 4), documentadas y
# respaldadas por la regresión de la Fase 4:
#   1. Actualización del aire interior con UN solo dt (físicamente correcta):
#        Tint += dt·Σ_i hi·dx·(Tsurf_i − Ti)/(ρc·La·X)
#      (el C tenía un dt² latente que solo coincide a dt=1). Así el 2D reduce
#      EXACTAMENTE al 1D del paquete a cualquier dt.
#   2. Promedios de superficie a /(nx-1) (medios nodos), para Tso y Tsi.
#   3. ho, hi se toman tal cual (sin el override de muro del C); la elección de
#      hi la hace la capa de API según corresponda.


@njit(cache=True)
def _step_inner(k, rhoc, To, T, Tsa, Tint, ho, hi, dt, dx, dy,
                a, b, c, d, P, Q, Tn, Tnew, tol):
    """Lazo interno (Gauss-Seidel por líneas) de un paso; actualiza ``T`` in situ.

    Buffers de trabajo (preasignados): ``a,b,c,d`` (nx,ny), ``P,Q,Tn`` (nx),
    ``Tnew`` (nx,ny). Asume NT ⊆ {1-8,13} (laterales adiabáticos, fronteras
    convectivas en j=0/j=ny-1) — válido para bovedilla rellena.
    """
    nx, ny = k.shape
    iters = 0
    while True:
        iters += 1
        # ---- ensamble a,b,c,d ----
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
        # ---- TDMA en x por cada fila j -> Tnew ----
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
        # ---- error con signo y commit T <- Tnew ----
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
    Motor de producción: corre el día completo repitiéndolo hasta régimen
    periódico (``mean|T_día−T_día_previo| < tol_day``).

    Returns:
        (Ti_series, Tso_series, Tsi_series, T_field, days, Qin, Qout)
        Series con forma ``(nsteps,)`` (un valor por paso de ``Tsa_arr``).
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
            # Tso: superficie exterior antes de resolver, /(nx-1)
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
            # actualización del aire interior (UN dt, físicamente correcto)
            flux = 0.0
            for i in range(nx):
                flux += hi * dx * (T[i, ny - 1] - Ti_old)
            Tint = Ti_old + dt * flux / Cair
            Ti_series[s] = Tint
            # Tsi: superficie interior tras resolver, /(nx-1)
            tsi = 0.0
            for i in range(nx):
                tsi += T[i, ny - 1]
            Tsi_series[s] = tsi / (nx - 1)
            # energía (por unidad de área interior)
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
#  Bovedilla con cámara de aire (tipo 1) — Fase 6
# =================================================================
#
# Física de cavidad: paredes del hueco (NT 9-12) convectan al aire del hueco
# (coef. hh por Nusselt) y radian entre sí (Stefan-Boltzmann + factores de vista);
# el aire del hueco es un nodo lumped Thueco (NT 0 fija T=Thueco). Port fiel de
# la rama tipo==1 de solve_PQ + casos 0,9-12 de calculate_coefficients.
# Sólo muro (beta=90): hh = 0.4005·|ΔT|^0.3033 / e22^0.0901.

_SIGMA = 5.6704e-8


@njit(cache=True)
def _step_hueca(k, rhoc, To, T, Tsa, Tint, Th, ho, hi, dt, dx, dy,
                i1, j1, i2, j2, e22, E,
                Fud, Ful, Fur, Fru, Frd, Frl, Fdl, Fdr, Fdu, Flu, Flr, Fld,
                a, b, c, d, P, Q, Tnew, tol):
    """Lazo interno de un paso para bovedilla con cámara de aire (muro)."""
    nx, ny = k.shape
    iters = 0
    hh = 1.0
    while True:
        iters += 1
        # --- temperaturas medias de las 4 paredes (T actual) ---
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
        # --- Nusselt (muro) ---
        hh = 0.4005 * (abs(tup - tdn) ** 0.3033) / (e22 ** 0.0901)
        # --- radiación entre paredes (factores de vista) ---
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
        # --- ensamble a,b,c,d ---
        for j in range(ny):
            for i in range(nx):
                nt_ij = 0
                # tipo de nodo deducido de las coordenadas del hueco
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
                    # aire del hueco: T = Th
                    a[i, j] = 1.0; b[i, j] = 0.0; c[i, j] = 0.0; d[i, j] = Th
                elif in_hole_cols and j == j1 - 1:
                    # pared superior (NT 9): hueco al sur
                    aP = apo + aN + hh * dx + aE + aW
                    dd = aN * T[i, j - 1] + hh * dx * Th + apo * To[i, j] - Qur - Qud - Qul
                    a[i, j] = aP; b[i, j] = aE; c[i, j] = aW; d[i, j] = dd
                elif in_hole_cols and j == j2:
                    # pared inferior (NT 10): hueco al norte
                    aP = apo + hh * dx + aS + aE + aW
                    dd = aS * T[i, j + 1] + hh * dx * Th + apo * To[i, j] - Qdl - Qdu - Qdr
                    a[i, j] = aP; b[i, j] = aE; c[i, j] = aW; d[i, j] = dd
                elif in_hole_rows and i == i1 - 1:
                    # pared izquierda (NT 11): hueco al este
                    aP = apo + aN + aS + hh * dy + aW
                    dd = aN * T[i, j - 1] + aS * T[i, j + 1] + apo * To[i, j] + hh * dy * Th - Qlu - Qlr - Qld
                    a[i, j] = aP; b[i, j] = 0.0; c[i, j] = aW; d[i, j] = dd
                elif in_hole_rows and i == i2:
                    # pared derecha (NT 12): hueco al oeste
                    aP = apo + aN + aS + aE + hh * dy
                    dd = aN * T[i, j - 1] + aS * T[i, j + 1] + apo * To[i, j] + hh * dy * Th - Qrd - Qrl - Qru
                    a[i, j] = aP; b[i, j] = aE; c[i, j] = 0.0; d[i, j] = dd
                else:
                    # nodo estándar (1-8,13)
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
        # --- TDMA en x por fila ---
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
        # --- error y commit ---
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
    """Factores de vista de la cavidad (h=e22 alto, l=a21 ancho), como el C."""
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
    Un paso de tiempo para bovedilla con cámara de aire (muro, beta=90).
    Devuelve dict con ``T, Tint, Thueco, hh, iters``.
    """
    if beta != 90.0:
        raise NotImplementedError("solve_step_hueca: sólo muro (beta=90) por ahora.")
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
    # aire del hueco (lumped, un solo dt)
    Qh_hole = 0.0
    for i in range(i1, i2):
        Qh_hole += hh * dx * (T[i, j1 - 1] - Th)
        Qh_hole += hh * dx * (T[i, j2] - Th)
    for j in range(j1, j2):
        Qh_hole += hh * dy * (T[i1 - 1, j] - Th)
        Qh_hole += hh * dy * (T[i2, j] - Th)
    Ch = rhoair * cair * a21 * e22
    Th_new = (Qh_hole + (Ch / dt) * Th) * dt / Ch
    # aire interior (idéntico al port fiel tipo 2: dt²)
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
    Día completo con convergencia día-a-día para bovedilla con cámara de aire.
    Réplica fiel del C (Tso /nx, Tsi /(nx-1), Tint con dt², Thueco lumped).
    El aire del hueco ``Th`` y el interior ``Tint`` marchan paso a paso.

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
            Tso_s[s] = tso / nx          # Tsout del C: /nx (réplica fiel)
            for j in range(ny):
                for i in range(nx):
                    To[i, j] = T[i, j]
            Ti_old = Tint
            Th_old = Th
            _, hh = _step_hueca(k, rhoc, To, T, Tsa_arr[s], Ti_old, Th_old,
                                ho, hi, dt, dx, dy, i1, j1, i2, j2, e22, E,
                                Fud, Ful, Fur, Fru, Frd, Frl, Fdl, Fdr, Fdu, Flu, Flr, Fld,
                                a, b, c, d, P, Q, Tnew, tol_inner)
            # aire del hueco
            qh = 0.0
            for i in range(i1, i2):
                qh += hh * dx * (T[i, j1 - 1] - Th_old)
                qh += hh * dx * (T[i, j2] - Th_old)
            for j in range(j1, j2):
                qh += hh * dy * (T[i1 - 1, j] - Th_old)
                qh += hh * dy * (T[i2, j] - Th_old)
            Th = (qh + (Ch / dt) * Th_old) * dt / Ch
            # aire interior (dt², fiel)
            qi = 0.0
            for i in range(nx):
                qi += hi * dt * dx * (T[i, ny - 1] - Ti_old)
            Tint = (qi + (Cair / dt) * Ti_old) * dt / Cair
            Ti_s[s] = Tint
            Th_s[s] = Th
            tsi = 0.0
            for i in range(nx):
                tsi += T[i, ny - 1]
            Tsi_s[s] = tsi / (nx - 1)    # max_min del C: /(nx-1)
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
    Versión de **producción** del día con cámara de aire (bloque hueco / bovedilla
    con aire). Igual que :func:`solve_day_hueca` pero con las correcciones
    conscientes de Fase 5: aire interior con **un solo `dt`** y superficies a
    `/(nx-1)`. La física del hueco (radiación + Nusselt + `Thueco`) es idéntica.
    Devuelve además `Qin, Qout` (energía por unidad de área interior, último día).

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
            Tso_s[s] = tso / (nx - 1)        # producción: /(nx-1)
            for j in range(ny):
                for i in range(nx):
                    To[i, j] = T[i, j]
            Ti_old = Tint
            Th_old = Th
            _, hh = _step_hueca(k, rhoc, To, T, Tsa_arr[s], Ti_old, Th_old,
                                ho, hi, dt, dx, dy, i1, j1, i2, j2, e22, E,
                                Fud, Ful, Fur, Fru, Frd, Frl, Fdl, Fdr, Fdu, Flu, Flr, Fld,
                                a, b, c, d, P, Q, Tnew, tol_inner)
            # aire del hueco (un solo dt, igual que la fiel — ya era correcto)
            qh = 0.0
            for i in range(i1, i2):
                qh += hh * dx * (T[i, j1 - 1] - Th_old)
                qh += hh * dx * (T[i, j2] - Th_old)
            for j in range(j1, j2):
                qh += hh * dy * (T[i1 - 1, j] - Th_old)
                qh += hh * dy * (T[i2, j] - Th_old)
            Th = (qh + (Ch / dt) * Th_old) * dt / Ch
            # aire interior (un solo dt, físicamente correcto)
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
#  Variante paralela (numba prange) — Fase 7
# =================================================================
#
# El método ya es Jacobi por líneas (las filas usan un único snapshot de T por
# iteración interna y se resuelven a un buffer Tnew antes de T<-Tnew), así que las
# filas son INDEPENDIENTES: paralelizar el barrido sobre j con prange no cambia el
# algoritmo ni la convergencia, solo reparte filas entre hilos. Portable: numba
# prange con su threading layer (workqueue por defecto, sin libs externas).


@njit(parallel=True, cache=True)
def _step_inner_par(k, rhoc, To, T, Tsa, Tint, ho, hi, dt, dx, dy,
                    a, b, c, d, Tnew, tol):
    """Igual que ``_step_inner`` pero con el barrido de filas en ``prange``.

    Cada fila j del TDMA usa buffers P,Q locales (un arreglo por iteración del
    prange) para evitar carreras. El error es una reducción de suma sobre j.
    """
    nx, ny = k.shape
    iters = 0
    while True:
        iters += 1
        # ---- ensamble (filas independientes) ----
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
        # ---- TDMA por fila (P,Q locales por hilo) ----
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
        # ---- error (reducción) y commit ----
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
    """Versión paralela de :func:`solve_day_2d`: el paralelismo (``prange`` sobre
    filas) vive en ``_step_inner_par``; el lazo día/paso es secuencial."""
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
