"""
Fase 4 — Integración temporal + convergencia día-a-día (bovedilla rellena).

Corre el lazo temporal completo (``t∈[0,86400]``, ``dt=1 s``) alimentado por la
``Tsa`` sinusoidal del C (``tests/c_boundary.py``), para UN día desde condición
inicial uniforme (igual que el golden de la Fase 0, ``DAYS=1``), y compara:

  - serie ``Tsa, Tso, Tsi, Tint`` (cada 600 s) vs ``golden/2d/gbv_5_1.csv``
  - índices vs ``golden/2d/indice_gbv_5_1.csv``

Aceptación "del orden" (PLAN-2D.md): serie ``atol 0.1 °C`` / ``rtol 1%``;
índices ~1-2 %.

La corrida del día es cara (~30 min en numpy puro). El driver guarda su resultado
en ``golden/2d/py_day_5_1.csv`` + ``py_indice_5_1.csv``; las pruebas comparan ese
caché contra el golden (rápido). Regenerar el caché:
    .venv/bin/python tests/test_eh2d_fullday.py
"""

import os

import numpy as np

from enerhabitat.ehtools2d import solve_step_2d

from test_eh2d_geometry import read_inp, read_field, section_from_inp, INP, GOLDEN
from c_boundary import Boundary

PY_SERIES = os.path.join(GOLDEN, "py_day_5_1.csv")
PY_INDEX = os.path.join(GOLDEN, "py_indice_5_1.csv")

RHOAIR = 1.1797660470258469
CAIR = 1005.458757


# --- driver del día ------------------------------------------------------------

def solve_day(sec, b, La, dt=1.0, t_max=86400.0, progress=False):
    """
    Replica el lazo temporal del `main` del C para bovedilla rellena, un día.

    Returns:
        (series, indices): ``series`` lista de filas
        ``[t_h, Is, Tsa, Ta, Tso, Tsi, Tint, Tc, DtaT]`` cada 600 s; ``indices``
        dict con las 10 columnas del índice del C.
    """
    nx, ny = sec.nx, sec.ny
    NT, k, rhoc = sec.NT, sec.kfield, sec.rhocfield
    dx, dy, X = sec.mesh.dx, sec.mesh.dy, sec.mesh.X
    ho = b.ho
    Tc, DtaT = b.Tc, b.DtaT
    hi = 8.1 if b.beta > 45.0 else 6.6   # convective_coefficients (muro/techo)

    Tint = Tc + DtaT
    T = np.full((nx, ny), Tint, dtype=np.float64)

    Tsamax = Tintmax = -100.0
    Tsamin = Tintmin = 100.0
    t_Tintmax = t_Tsamax = 0.0
    NumHot = NumCold = DenHot = DenCold = 0.0
    DDHhot = DDHcold = 0.0
    Qin_total = 0.0
    Tintaverage = 0.0

    series = []
    nsteps = int(round(t_max / dt)) + 1
    for s in range(nsteps):
        t = s * dt
        Tsa, Ta_inst, Is = b.tsa(t)
        Tsa0, _, _ = b.tsa(t, a=0.0)
        Tsa1, _, _ = b.tsa(t, a=1.0)

        # Tso: superficie exterior ANTES de resolver, /nx (Tsout del C).
        Tso = float(np.sum(T[:, 0]) / nx)

        res = solve_step_2d(NT, k, rhoc, T, Tsa, Tint, ho, hi, dt, dx, dy,
                            La, X, RHOAIR, CAIR)
        T = res["T"]
        Tint = res["Tint"]
        Qin_total += res["Qin"]

        # Tsi: superficie interior DESPUÉS de resolver, /(nx-1) (max_min del C).
        Tsi = float(np.sum(T[:, ny - 1]) / (nx - 1))

        # max_min
        if Tsa > Tsamax:
            Tsamax, t_Tsamax = Tsa, t
        if Tsa < Tsamin:
            Tsamin = Tsa
        if Tint < Tintmin:
            Tintmin = Tint
        if Tint > Tintmax:
            Tintmax, t_Tintmax = Tint, t

        # discomfort + TPI
        if Tint < Tc:
            DDHcold += (Tc - Tint) * dt / 3600.0
        if Tint > Tc:
            DDHhot += (Tint - Tc) * dt / 3600.0
        if Tint < Tc:
            NumCold += Tc - Tint
        if Tsa0 < Tc:
            DenCold += Tc - Tsa0
        if Tint > Tc:
            NumHot += Tint - Tc
        if Tsa1 > Tc:
            DenHot += Tsa1 - Tc

        Tintaverage += Tint

        if (t % 600.0) < dt:
            series.append([t / 3600.0, Is, Tsa, Ta_inst, Tso, Tsi, Tint, Tc, DtaT])

        if progress and s % 7200 == 0:
            print(f"  t={t/3600:5.1f} h  Tint={Tint:6.3f}  iters={res['iters']}")

    contador = 86400.0 / dt
    indices = {
        "Qin": Qin_total / 3600.0 / X,
        "decremento": (Tintmax - Tintmin) / (Tsamax - Tsamin),
        "retardo_h": (t_Tintmax - t_Tsamax) / 3600.0,
        "Tint_media": Tintaverage / contador,
        "Tintmin": Tintmin,
        "Tintmax": Tintmax,
        "TPIhot": (1.0 - NumHot / DenHot) * 100.0,
        "TPIcold": (1.0 - NumCold / DenCold) * 100.0,
        "DDHhot": DDHhot,
        "DDHcold": DDHcold,
    }
    return series, indices


# --- caché / IO ----------------------------------------------------------------

_IDX_ORDER = ["Qin", "decremento", "retardo_h", "Tint_media", "Tintmin",
              "Tintmax", "TPIhot", "TPIcold", "DDHhot", "DDHcold"]


def run_and_cache():
    params = read_inp(INP)
    sec = section_from_inp(params)
    b = Boundary(params)
    series, indices = solve_day(sec, b, La=float(params["La"]), progress=True)
    np.savetxt(PY_SERIES, np.array(series), delimiter="\t",
               header="t\tIs\tTsa\tTa\tTso\tTsi\tTint\tTc\tDtaT")
    with open(PY_INDEX, "w") as f:
        f.write("\t".join(f"{indices[k]:.6f}" for k in _IDX_ORDER) + "\n")
    return series, indices


def _load_series(path, skip):
    rows = []
    with open(path) as f:
        for i, line in enumerate(f):
            if i < skip or line.startswith("#"):
                continue
            rows.append([float(x) for x in line.split()])
    return np.array(rows)


def _load_py():
    if not os.path.exists(PY_SERIES):
        run_and_cache()
    series = _load_series(PY_SERIES, skip=0)   # header empieza con '#', se salta
    idx = _load_series(PY_INDEX, skip=0)[0]
    return series, idx


# --- pruebas -------------------------------------------------------------------

def test_series_vs_c():
    py = _load_py()[0]
    c = _load_series(os.path.join(GOLDEN, "gbv_5_1.csv"), skip=2)
    n = min(len(py), len(c))
    # columnas: 0 t,1 Is,2 Tsa,3 Ta,4 Tso,5 Tsi,6 Tint,7 Tc,8 DtaT
    for col, name in [(2, "Tsa"), (4, "Tso"), (5, "Tsi"), (6, "Tint")]:
        d = np.abs(py[:n, col] - c[:n, col])
        # el golden está a 2 decimales -> sumamos ese margen al atol
        assert np.all(d <= 0.1 + 5e-3), \
            f"serie {name}: max|Δ|={d.max():.3f} °C"


def test_indices_vs_c():
    idx = _load_py()[1]
    c = _load_series(os.path.join(GOLDEN, "indice_gbv_5_1.csv"), skip=0)[0]
    for j, name in enumerate(_IDX_ORDER):
        a, bvar = idx[j], c[j]
        tol = max(0.05, abs(bvar) * 0.02)   # ~2 % o el redondeo del golden
        assert abs(a - bvar) <= tol, f"índice {name}: py={a:.3f} c={bvar:.3f}"


def _demo():
    py, idx = _load_py()
    c = _load_series(os.path.join(GOLDEN, "gbv_5_1.csv"), skip=2)
    cidx = _load_series(os.path.join(GOLDEN, "indice_gbv_5_1.csv"), skip=0)[0]
    n = min(len(py), len(c))
    bar = "═" * 70
    print(bar)
    print("  FASE 4 · Día completo (RELLENA, tipo 2)   Python vs C golden")
    print(bar)

    print(f"\n  Serie — max |Δ| Python vs C (de {n} muestras a 600 s):")
    for col, name in [(2, "Tsa"), (4, "Tso"), (5, "Tsi"), (6, "Tint")]:
        d = np.abs(py[:n, col] - c[:n, col])
        print(f"    {name:<5} max|Δ| = {d.max():.3f} °C   media = {d.mean():.3f} °C")

    print(f"\n  Índices Python vs C:")
    print(f"    {'índice':<14}{'Python':>10}{'C':>10}")
    for j, name in enumerate(_IDX_ORDER):
        print(f"    {name:<14}{idx[j]:>10.3f}{cidx[j]:>10.3f}")

    # Curva diaria de Tint y Tsa (ASCII).
    print(f"\n{bar}\n  Tint(t) y Tsa(t) a lo largo del día (Python)")
    print(bar)
    Tsa = py[:n, 2]
    Tint = py[:n, 6]
    lo = min(Tsa.min(), Tint.min())
    hi = max(Tsa.max(), Tint.max())
    W = 56
    for r in range(n):
        h = py[r, 0]
        def pos(v):
            return int((v - lo) / (hi - lo) * (W - 1))
        line = [" "] * W
        line[pos(Tsa[r])] = "s"
        line[pos(Tint[r])] = "i"
        print(f"  {h:4.1f}h |{''.join(line)}|")
    print(f"  rango {lo:.1f}..{hi:.1f} °C   's'=Tsa(sol-aire)  'i'=Tint(interior)")

    print(f"\n{bar}")
    print("  Fase 4: la serie e índices reproducen el C 'del orden' ✅")
    print(bar)


if __name__ == "__main__":
    run_and_cache()
    test_series_vs_c()
    test_indices_vs_c()
    _demo()
