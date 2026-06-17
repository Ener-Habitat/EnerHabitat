"""
Fase 7 — Desempeño y paralelización (bovedilla rellena).

El solver ya es **Jacobi por líneas** (filas independientes por iteración interna),
así que la variante paralela (`solve_day_2d_par`, `prange` sobre filas) ejecuta el
MISMO algoritmo repartido en hilos: sin costo en convergencia.

Mide:
  - **correctitud**: paralelo reproduce serial (atol 1e-6, en la práctica ~máquina);
  - **speedup** vs nº de hilos (numba `set_num_threads`);
  - **escalado** con el tamaño de malla.

Portable: numba `prange` con su threading layer por defecto (sin libs externas).

    .venv/bin/python tests/test_eh2d_perf.py
"""

import os
import time

import numpy as np
import numba

from enerhabitat.eh2d import Section2D, Fill
from enerhabitat.ehtools2d import solve_day_2d, solve_day_2d_par

RHOAIR = 1.1797660470258469
CAIR = 1005.458757
A = {"a11": 0.02, "a12": 0.03, "a13": 0.02, "a14": 0.0,
     "a21": 0.16, "a22": 0.16, "a23": 0.0}
E = {"e21": 0.02, "e22": 0.08, "e23": 0.02}


def _geom(nx, ny):
    sec = Section2D(nx=nx, ny=ny, L=[0.12, 0, 0, 0, 0, 0, 0],
                    k=[1.35, 0, 0, 0, 0, 0, 0], rhoc=[1.8e6, 0, 0, 0, 0, 0, 0],
                    kr=0.026, rhocr=64000, a=A, e=E, layer=1,
                    fill_type=Fill.SOLID).build()
    return sec


def _tsa(nsteps):
    t = np.arange(nsteps) * 1.0
    return 28.0 + 8.0 * np.sin(2 * np.pi * t / 86400.0)


def _run(fn, sec, Tsa, T0=28.0):
    m = sec.mesh
    return fn(sec.NT, sec.kfield, sec.rhocfield, Tsa, 13.0, 8.1, 1.0,
              m.dx, m.dy, 2.5, m.X, RHOAIR, CAIR, T0,
              1e-10, 5e-4, 1)   # max_days=1: una pasada de nsteps


def _time(fn, sec, Tsa):
    _run(fn, sec, Tsa)                      # warmup (compila)
    t0 = time.perf_counter()
    out = _run(fn, sec, Tsa)
    return time.perf_counter() - t0, out


def test_parallel_matches_serial():
    sec = _geom(80, 80)
    Tsa = _tsa(120)
    s = _run(solve_day_2d, sec, Tsa)
    p = _run(solve_day_2d_par, sec, Tsa)
    assert np.max(np.abs(s[0] - p[0])) <= 1e-6           # Ti series
    assert np.max(np.abs(s[3] - p[3])) <= 1e-6           # campo final


def _demo():
    bar = "═" * 70
    ncores = numba.config.NUMBA_DEFAULT_NUM_THREADS
    print(bar)
    print("  FASE 7 · Desempeño y paralelización (Jacobi por líneas, prange)")
    print(f"  Núcleos disponibles: {ncores}   (numba threading layer portable)")
    print(bar)

    # --- correctitud ---
    sec = _geom(80, 80); Tsa = _tsa(120)
    s = _run(solve_day_2d, sec, Tsa)
    p = _run(solve_day_2d_par, sec, Tsa)
    dTi = np.max(np.abs(s[0] - p[0])); dF = np.max(np.abs(s[3] - p[3]))
    print(f"\n  Correctitud paralelo vs serial (malla 80×80, 120 pasos):")
    print(f"    max|ΔTi| = {dTi:.2e} °C   max|ΔT_campo| = {dF:.2e} °C   "
          f"{'✓ idénticos' if max(dTi, dF) <= 1e-6 else '✗'}")

    # --- speedup vs hilos ---
    sec = _geom(160, 160); Tsa = _tsa(120)
    ts, _ = _time(solve_day_2d, sec, Tsa)
    print(f"\n  Speedup vs nº de hilos (malla 160×160, 120 pasos):")
    print(f"    {'hilos':>6}{'t [s]':>12}{'speedup':>10}{'efic.':>9}")
    print(f"    {'serial':>6}{ts:>12.3f}{'1.00×':>10}{'—':>9}")
    best = (1, 1.0)
    for nth in [1, 2, 4, 8]:
        if nth > ncores:
            continue
        numba.set_num_threads(nth)
        tp, _ = _time(solve_day_2d_par, sec, Tsa)
        sp = ts / tp
        print(f"    {nth:>6}{tp:>12.3f}{sp:>9.2f}×{sp/nth*100:>8.0f}%")
        if sp > best[1]:
            best = (nth, sp)
    numba.set_num_threads(ncores)

    # --- escalado de malla (serial) ---
    print(f"\n  Escalado de malla (serial, 60 pasos, ms por paso):")
    print(f"    {'malla':>10}{'nodos':>10}{'ms/paso':>12}")
    Tsa = _tsa(60)
    for n in [40, 80, 160, 240]:
        sec = _geom(n, n)
        tt, out = _time(solve_day_2d, sec, Tsa)
        print(f"    {f'{n}×{n}':>10}{n*n:>10}{tt/len(Tsa)*1000:>12.1f}")

    # --- barrido de dt (esquema implícito: menos pasos para la misma ventana) ---
    print(f"\n  Barrido de dt — ventana física de 4 h, malla 80×80 "
          f"(implícito: incondicionalmente estable):")
    print(f"    {'dt [s]':>7}{'pasos':>8}{'t [s]':>10}{'Tint final':>12}{'Δ vs dt=1':>12}")
    sec = _geom(80, 80)
    window = 4 * 3600.0
    ref = None
    for dt in [1.0, 10.0, 60.0, 300.0]:
        nsteps = int(window / dt)
        t = np.arange(nsteps) * dt
        Tsa = 28.0 + 8.0 * np.sin(2 * np.pi * t / 86400.0)
        m = sec.mesh
        _run_dt = lambda: solve_day_2d(sec.NT, sec.kfield, sec.rhocfield, Tsa,
                                       13.0, 8.1, dt, m.dx, m.dy, 2.5, m.X,
                                       RHOAIR, CAIR, 28.0, 1e-10, 5e-4, 1)
        _run_dt()                       # warmup
        t0 = time.perf_counter(); out = _run_dt(); el = time.perf_counter() - t0
        Tfin = out[0][-1]
        if ref is None:
            ref = Tfin
        print(f"    {dt:>7.0f}{nsteps:>8}{el:>10.3f}{Tfin:>12.4f}"
              f"{abs(Tfin-ref):>12.4f}")

    print(f"\n{bar}")
    print(f"  Recomendación:")
    print(f"   · Paralelizar con numba prange es PORTABLE (threading layer interno,")
    print(f"     sin libs externas) y SIN costo de convergencia (ya es Jacobi por")
    print(f"     líneas) → reproduce el serial al bit. Speedup modesto:")
    print(f"     {best[1]:.2f}× con {best[0]} hilos (granularidad fina: filas cortas).")
    print(f"   · La palanca mayor es subir dt (esquema implícito): la ventana se")
    print(f"     resuelve con ~1/dt de los pasos, con error 'del orden' acotado.")
    print(bar)


if __name__ == "__main__":
    test_parallel_matches_serial()
    print("PASS  test_parallel_matches_serial\n")
    _demo()
