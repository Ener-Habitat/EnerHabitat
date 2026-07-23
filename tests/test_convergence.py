"""
P0-03 — Criterio de convergencia interna 2D.

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

Correr con pytest o como script:
    .venv/bin/python tests/test_convergence.py
"""

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


if __name__ == "__main__":
    for fn in (test_step_new_matches_legacy, test_new_mode_verifies_residual,
               test_max_iter_flags_nonconvergence, test_day_solver_diagnostics):
        fn()
        print(f"PASS  {fn.__name__}")
    print("\nP0-03: criterio no cancelable + residual + diagnóstico ✅")
