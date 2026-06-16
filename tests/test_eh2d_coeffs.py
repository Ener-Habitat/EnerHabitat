"""
Fase 2 — Ensamble de coeficientes 2D (bovedilla rellena, `tipo 2`).

Compara ``a,b,c,d`` de ``enerhabitat.ehtools2d.calculate_coefficients_2d`` contra
el dump del C (``tests/golden/2d/dump_coef_*.dat``, generado con
``make DUMPCOEF=1``), usando el MISMO campo determinista ``T,To`` y los mismos
``Tsa,Tint,ho,hi,dt`` que volcó el C.

Tolerancia: rtol 1e-10 (ver PLAN-2D.md).

Correr con pytest o como script:
    .venv/bin/python tests/test_eh2d_coeffs.py
"""

import os

import numpy as np

from enerhabitat.eh2d import Section2D, Bovedilla
from enerhabitat.ehtools2d import calculate_coefficients_2d

from test_eh2d_geometry import read_inp, read_meta, read_field, section_from_inp, INP, GOLDEN

RTOL = 1e-10


def _build():
    params = read_inp(INP)
    sec = section_from_inp(params)
    meta = read_meta(os.path.join(GOLDEN, "dump_coef_meta.dat"))
    nx, ny = sec.nx, sec.ny
    T = read_field(os.path.join(GOLDEN, "dump_coef_T.dat"), nx, ny, float)
    To = read_field(os.path.join(GOLDEN, "dump_coef_To.dat"), nx, ny, float)

    a, b, c, d = calculate_coefficients_2d(
        sec.NT, sec.kfield, sec.rhocfield, To, T,
        Tsa=meta["Tsa"], Tint=meta["Tint"], ho=meta["ho"], hi=meta["hi"],
        dt=meta["dt"], dx=sec.mesh.dx, dy=sec.mesh.dy)

    golden = {name: read_field(os.path.join(GOLDEN, f"dump_coef_{name}.dat"), nx, ny, float)
              for name in ("a", "b", "c", "d")}
    return sec, {"a": a, "b": b, "c": c, "d": d}, golden, meta, (T, To)


def _check(name, py, golden):
    assert np.allclose(py, golden, rtol=RTOL, atol=0), \
        f"{name}: difiere; max rel err = {np.nanmax(np.abs(py-golden)/np.abs(golden)):.2e}"


def test_coef_a():
    _, P, G, _, _ = _build()
    _check("a", P["a"], G["a"])


def test_coef_b():
    _, P, G, _, _ = _build()
    _check("b", P["b"], G["b"])


def test_coef_c():
    _, P, G, _, _ = _build()
    _check("c", P["c"], G["c"])


def test_coef_d():
    _, P, G, _, _ = _build()
    _check("d", P["d"], G["d"])


def _demo():
    sec, P, G, meta, (T, To) = _build()
    nx, ny = sec.nx, sec.ny
    bar = "═" * 70
    print(bar)
    print("  FASE 2 · Ensamble de coeficientes a,b,c,d (RELLENA, tipo 2)")
    print(f"  Python (ehtools2d)  vs  C legacy   malla {nx}×{ny}")
    print(f"  Campo T/To determinista; Tsa={meta['Tsa']:g} Tint={meta['Tint']:g} "
          f"ho={meta['ho']:g} hi={meta['hi']:g} dt={meta['dt']:g}")
    print(bar)

    print(f"\n  {'coef':<6}{'max |Δ|':>14}{'max rel err':>16}   estado")
    print("  " + "-" * 50)
    for name in ("a", "b", "c", "d"):
        py, g = P[name], G[name]
        amax = float(np.max(np.abs(py - g)))
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = np.abs(py - g) / np.abs(g)
            rmax = float(np.nanmax(np.where(np.isfinite(rel), rel, 0.0)))
        ok = np.allclose(py, g, rtol=RTOL, atol=0)
        print(f"  {name:<6}{amax:>14.3e}{rmax:>16.3e}   {'✓' if ok else '✗ DIFIERE'}")

    # Vista de un nodo de cada tipo, para que se entienda qué se ensambló.
    m = sec.mesh
    samples = [
        ("1  esq sup-izq (ext)", 0, 0),
        ("5  borde exterior",    nx // 2, 0),
        ("8  borde interior",    nx // 2, ny - 1),
        ("6  lateral izq",       0, ny // 2),
        ("13 interior (vigueta)", max(1, m.i1 // 2), ny // 2),
        ("13 relleno (bovedilla)", (m.i1 + m.i2) // 2, (m.j1 + m.j2) // 2),
    ]
    print(f"\n{bar}\n  Coeficientes Python vs C en nodos representativos")
    print(bar)
    print(f"  {'nodo':<24}{'':>6}{'a':>12}{'b(aE)':>10}{'c(aW)':>10}{'d':>12}")
    for label, i, j in samples:
        print(f"  {label:<24}{'PY':>6}{P['a'][i,j]:>12.4f}{P['b'][i,j]:>10.4f}"
              f"{P['c'][i,j]:>10.4f}{P['d'][i,j]:>12.4f}")
        print(f"  {'':<24}{'C':>6}{G['a'][i,j]:>12.4f}{G['b'][i,j]:>10.4f}"
              f"{G['c'][i,j]:>10.4f}{G['d'][i,j]:>12.4f}")

    print(f"\n{bar}")
    print("  Fase 2: el ensamble de Python reproduce el del C (rtol 1e-10) ✅")
    print(bar)


if __name__ == "__main__":
    for fn in (test_coef_a, test_coef_b, test_coef_c, test_coef_d):
        fn()
    _demo()
