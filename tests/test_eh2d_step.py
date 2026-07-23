"""
Fase 3 — Solver de un paso 2D (bovedilla rellena, `tipo 2`).

Compara el campo ``T`` y ``Tint`` tras UN paso de ``solve_step_2d`` contra el dump
del C (``tests/golden/2d/dump_step_*.dat``, generado con ``make DUMPSTEP=1``),
partiendo del mismo campo inicial y los mismos escalares.

Tolerancia: ``atol 1e-6 °C`` en T y Tint; el nº de iteraciones internas debe
coincidir (ver PLAN-2D.md).

Correr con pytest o como script:
    .venv/bin/python tests/test_eh2d_step.py
"""

import os

import numpy as np

from enerhabitat.ehtools2d import solve_step_2d

from test_eh2d_geometry import (read_inp, read_meta, read_field,
                                section_from_inp, INP, GOLDEN,
                                HAS_LEGACY, LEGACY_REASON)

try:
    import pytest
    pytestmark = pytest.mark.skipif(not HAS_LEGACY, reason=LEGACY_REASON)
except ImportError:
    pass

ATOL = 1e-6


def _build():
    params = read_inp(INP)
    sec = section_from_inp(params)
    meta = read_meta(os.path.join(GOLDEN, "dump_step_meta.dat"))
    nx, ny = sec.nx, sec.ny
    T0 = read_field(os.path.join(GOLDEN, "dump_step_T0.dat"), nx, ny, float)
    T_c = read_field(os.path.join(GOLDEN, "dump_step_T.dat"), nx, ny, float)

    res = solve_step_2d(
        sec.NT, sec.kfield, sec.rhocfield, T0,
        Tsa=meta["Tsa"], Tint=meta["Tint_in"], ho=meta["ho"], hi=meta["hi"],
        dt=meta["dt"], dx=sec.mesh.dx, dy=sec.mesh.dy,
        La=meta["La"], X=meta["X"], rhoair=meta["rhoair"], cair=meta["cair"])
    return sec, res, T_c, meta, T0


def test_field_T():
    _, res, T_c, _, _ = _build()
    assert np.allclose(res["T"], T_c, atol=ATOL, rtol=0), \
        f"T difiere; max |Δ| = {np.max(np.abs(res['T']-T_c)):.2e} °C"


def test_Tint():
    _, res, _, meta, _ = _build()
    assert abs(res["Tint"] - meta["Tint_out"]) <= ATOL, \
        f"Tint py={res['Tint']:.9f} c={meta['Tint_out']:.9f}"


def test_inner_iters():
    _, res, _, meta, _ = _build()
    assert res["iters"] == int(meta["inner_iters"]), \
        f"iteraciones internas py={res['iters']} c={int(meta['inner_iters'])}"


def _demo():
    sec, res, T_c, meta, T0 = _build()
    nx, ny = sec.nx, sec.ny
    dT = res["T"] - T_c
    bar = "═" * 70
    print(bar)
    print("  FASE 3 · Solver de un paso (RELLENA, tipo 2)")
    print(f"  Python (solve_step_2d)  vs  C legacy   malla {nx}×{ny}")
    print(f"  Tsa={meta['Tsa']:g}  Tint_in={meta['Tint_in']:g}  "
          f"ho={meta['ho']:g}  hi={meta['hi']:g}  dt={meta['dt']:g}")
    print(bar)

    iters_ok = res["iters"] == int(meta["inner_iters"])
    tint_d = abs(res["Tint"] - meta["Tint_out"])
    print(f"\n  iteraciones internas : py={res['iters']}  "
          f"c={int(meta['inner_iters'])}   {'✓' if iters_ok else '✗'}")
    print(f"  campo T  max |Δ|     : {np.max(np.abs(dT)):.3e} °C   "
          f"{'✓' if np.max(np.abs(dT)) <= ATOL else '✗'}  (atol {ATOL:g})")
    print(f"  Tint  py={res['Tint']:.9f}  c={meta['Tint_out']:.9f}   "
          f"|Δ|={tint_d:.2e}   {'✓' if tint_d <= ATOL else '✗'}")
    print(f"  Tint cambió {meta['Tint_in']:.4f} -> {res['Tint']:.6f} °C "
          f"(aire interior se enfría hacia el muro)")

    # Perfil de temperatura en una columna (exterior->interior): se ve la difusión.
    icol = nx // 4   # columna en la vigueta (concreto), fuera del relleno
    print(f"\n{bar}\n  Perfil T en la columna i={icol} (exterior j=0 → interior j={ny-1})")
    print(bar)
    js = [0, ny // 8, ny // 4, ny // 2, 3 * ny // 4, 7 * ny // 8, ny - 1]
    print(f"  {'j':>5}{'T0 inicial':>14}{'T paso (PY)':>14}{'T paso (C)':>14}{'Δ vs C':>12}")
    for j in js:
        print(f"  {j:>5}{T0[icol,j]:>14.5f}{res['T'][icol,j]:>14.5f}"
              f"{T_c[icol,j]:>14.5f}{res['T'][icol,j]-T_c[icol,j]:>12.1e}")

    # Mapa de la magnitud del cambio |T - T0| (dónde se movió la temperatura).
    print(f"\n{bar}\n  Dónde cambió T en el paso  |T_paso - T_inicial|  (exterior arriba)")
    print("   ' ' ~0   '·' pequeño   '+' medio   '#' grande")
    print(bar)
    change = np.abs(res["T"] - T0)
    cmax = change.max() or 1.0
    glyphs = " ·+#"
    cols, rows = 60, 24
    istep = max(1, nx // cols)
    jstep = max(1, ny // rows)
    for j in range(0, ny, jstep):
        line = ""
        for i in range(0, nx, istep):
            level = int(min(3, change[i, j] / cmax * 4))
            line += glyphs[level]
        print("  " + line)
    print(f"  (cambio máximo en el paso: {cmax:.3f} °C, junto a la frontera exterior)")

    print(f"\n{bar}")
    print("  Fase 3: el paso de Python reproduce el del C (atol 1e-6, mismas iters) ✅")
    print(bar)


if __name__ == "__main__":
    if not HAS_LEGACY:
        raise SystemExit(f"SKIP: {LEGACY_REASON}")
    for fn in (test_field_T, test_Tint, test_inner_iters):
        fn()
    _demo()
