"""
P0-03 / P0-04 — Criterios de convergencia 2D (interno) y diario (todos los
estados).

Valida el criterio no cancelable (max|ΔT| + residual escalado de las ecuaciones
discretas, à la Patankar) que sustituye en producción al promedio con signo del
C, y el diagnóstico de convergencia (``max_inner``, ``converged``):

  1. **Equivalencia física**: un paso resuelto con el criterio nuevo coincide
     con el legado (ambos bien convergidos) dentro de la tolerancia.
  2. **Residual verificado**: el campo aceptado por el criterio nuevo satisface
     ``max |a·T_P − b·T_E − c·T_W − d| / a ≤ tol`` re-ensamblando por fuera.
  3. **Tope de iteraciones**: ``max_iter`` alcanzado → ``converged = False``,
     nunca una salida indistinguible de una convergida.
  4. **Diagnóstico del día**: ``solve_day_2d`` reporta ``day_error``,
     ``inner_ok`` e ``inner_iters_max``; ``max_inner`` insuficiente se señala.

P0-04 (convergencia diaria, todos los estados):

  5. **Cierre energético**: en régimen periódico ``Qin ≈ Qout``.
  6. **1D con tope y diagnóstico**: ``System`` expone ``days``, ``day_error``,
     ``converged`` y ``energy_imbalance``; ``MAX_DAYS`` insuficiente avisa con
     ``RuntimeWarning`` y ``converged=False``.

Correr con pytest o como script:
    .venv/bin/python tests/test_convergence.py
"""

import os
import warnings as _warnings

import numpy as np

from enerhabitat.ehtools2d import (calculate_coefficients_2d, solve_step_2d,
                                   solve_day_2d)

RHOAIR = 1.1797660470258469
CAIR = 1005.458757
LA, X = 2.5, 0.1


def _case(nx=8, ny=12):
    """Bloque sólido homogéneo pequeño (NT=13 en toda la sección)."""
    NT = np.full((nx, ny), 13, dtype=np.int64)
    k = np.full((nx, ny), 1.0)
    rhoc = np.full((nx, ny), 1.5e6)
    dx, dy = X / nx, 0.12 / ny
    T0 = np.full((nx, ny), 20.0)
    return NT, k, rhoc, dx, dy, T0


def _step(legacy, tol=1e-10, max_iter=100000):
    NT, k, rhoc, dx, dy, T0 = _case()
    return solve_step_2d(NT, k, rhoc, T0, 35.0, 20.0, 13.0, 8.1, 10.0,
                         dx, dy, LA, X, RHOAIR, CAIR,
                         tol=tol, max_iter=max_iter, legacy=legacy)


def test_step_new_matches_legacy():
    """Mismo paso, ambos criterios bien convergidos → misma física."""
    r_old = _step(legacy=True)
    r_new = _step(legacy=False)
    assert r_new["converged"]
    d = np.max(np.abs(r_old["T"] - r_new["T"]))
    assert d <= 1e-8, f"criterios difieren: max|Δ|={d:.3e} °C"
    assert abs(r_old["Tint"] - r_new["Tint"]) <= 1e-8


def test_new_mode_verifies_residual():
    """El campo aceptado satisface las ecuaciones discretas (re-ensamble)."""
    tol = 1e-9
    NT, k, rhoc, dx, dy, T0 = _case()
    r = _step(legacy=False, tol=tol)
    T = r["T"]
    a, b, c, d = calculate_coefficients_2d(NT, k, rhoc, T0, T, 35.0, 20.0,
                                           13.0, 8.1, 10.0, dx, dy)
    rr = d - a * T
    rr[:-1, :] += b[:-1, :] * T[1:, :]
    rr[1:, :] += c[1:, :] * T[:-1, :]
    res = float(np.max(np.abs(rr) / a))
    assert r["converged"]
    assert res <= tol, f"residual escalado {res:.3e} > tol {tol:.0e}"


def test_max_iter_flags_nonconvergence():
    """Tope de barridos alcanzado → converged=False (paro informativo)."""
    r = _step(legacy=False, max_iter=2)
    assert not r["converged"]


def test_day_solver_diagnostics():
    """solve_day_2d reporta day_error, inner_ok e inner_iters_max."""
    NT, k, rhoc, dx, dy, _ = _case()
    nsteps = 360
    t = np.arange(nsteps) * 10.0
    Tsa = 28.0 + 8.0 * np.sin(2.0 * np.pi * t / (nsteps * 10.0))

    out = solve_day_2d(NT, k, rhoc, Tsa, 13.0, 8.1, 10.0, dx, dy, LA, X,
                       RHOAIR, CAIR, 22.0, 1e-8, 5e-4, 300, 10000)
    Ti, Tso, Tsi, T, days, Qin, Qout, day_err, inner_ok, inner_max = out
    assert inner_ok, "algún paso no convergió con max_inner holgado"
    assert day_err <= 5e-4, f"day_error={day_err:.3e}"
    assert days < 300
    assert 0 < inner_max < 10000

    # max_inner insuficiente → se señala, no se oculta
    out = solve_day_2d(NT, k, rhoc, Tsa, 13.0, 8.1, 10.0, dx, dy, LA, X,
                       RHOAIR, CAIR, 22.0, 1e-8, 5e-4, 2, 1)
    inner_ok_capped = out[8]
    assert not inner_ok_capped, "max_inner=1 debería marcar inner_ok=False"


def test_energy_closure_when_periodic():
    """P0-04: en régimen periódico convergido, Qin ≈ Qout (cierre energético)."""
    NT, k, rhoc, dx, dy, _ = _case()
    nsteps = 360
    t = np.arange(nsteps) * 10.0
    Tsa = 28.0 + 8.0 * np.sin(2.0 * np.pi * t / (nsteps * 10.0))
    out = solve_day_2d(NT, k, rhoc, Tsa, 13.0, 8.1, 10.0, dx, dy, LA, X,
                       RHOAIR, CAIR, 22.0, 1e-8, 5e-4, 300, 10000)
    Qin, Qout, day_err = out[5], out[6], out[7]
    assert day_err <= 5e-4
    qmax = max(Qin, Qout)
    assert qmax > 0.0
    imbalance = abs(Qin - Qout) / qmax
    assert imbalance <= 0.05, f"desbalance energético {imbalance:.3f}"


def _sys1d():
    import enerhabitat as eh
    from enerhabitat.config import config
    here = os.path.dirname(os.path.abspath(__file__))
    config.file = os.path.join(here, "materials.ini")
    loc = eh.Location(os.path.join(
        here, "MEX_MOR_Cuernavaca-Matamoros.Intl.AP.767260_TMYx.2004-2018.epw"))
    loc.meanDay(month=5, year=2025)
    return eh.System(loc, tilt=90, azimuth=0, absortance=0.8,
                     layers=[("EPS", 0.1)])


def test_1d_diagnostics():
    """P0-04: el 1D expone days/day_error/converged/energy_imbalance."""
    sys1 = _sys1d()
    sys1.solve()
    assert sys1.converged is True
    assert sys1.days is not None and 0 < sys1.days < 60
    assert sys1.day_error is not None and sys1.day_error <= 5e-4
    assert sys1.energy_imbalance is not None and sys1.energy_imbalance <= 0.05


def test_1d_max_days_flags_nonconvergence():
    """P0-04: MAX_DAYS insuficiente → converged=False + RuntimeWarning."""
    from enerhabitat import ehframe
    sys1 = _sys1d()
    original = ehframe.MAX_DAYS
    try:
        ehframe.MAX_DAYS = 1
        with _warnings.catch_warnings(record=True) as caught:
            _warnings.simplefilter("always")
            sys1.solve()
        assert sys1.converged is False
        assert sys1.days == 1
        assert any(issubclass(w.category, RuntimeWarning) for w in caught), \
            "no se emitió RuntimeWarning al no converger"
    finally:
        ehframe.MAX_DAYS = original


if __name__ == "__main__":
    for fn in (test_step_new_matches_legacy, test_new_mode_verifies_residual,
               test_max_iter_flags_nonconvergence, test_day_solver_diagnostics,
               test_energy_closure_when_periodic, test_1d_diagnostics,
               test_1d_max_days_flags_nonconvergence):
        fn()
        print(f"PASS  {fn.__name__}")
    print("\nP0-03 + P0-04: criterios no cancelables, todos los estados, "
          "diagnóstico y cierre energético ✅")
