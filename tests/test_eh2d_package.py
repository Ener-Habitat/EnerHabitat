"""
Fase 5 — Integración al paquete (API 2D de producción).

Valida la API ``System2D`` (reúsa EPW+pvlib vía el ``System`` 1D, motor JIT
``solve_day_2d``, ``config2d``):

  1. **Reducción al 1D**: una sección homogénea — un ``HollowBlock`` RELLENO de
     su mismo material, uniforme en x — reproduce la ``Ti`` del ``System`` 1D
     del paquete (``atol 0.1 °C``).
  2. **Periodicidad**: la convergencia día-a-día cierra antes de ``max_days``.
  3. **Balance de energía**: en régimen ``Qin ≈ Qout``.

Correr con pytest o como script:
    .venv/bin/python tests/test_eh2d_package.py
"""

import os

import numpy as np

from enerhabitat import Fill, HollowBlock, Location, System, config
from enerhabitat.config import config2d
from enerhabitat.eh2d import System2D

HERE = os.path.dirname(os.path.abspath(__file__))
EPW = os.path.join(HERE, "MEX_MOR_Cuernavaca-Matamoros.Intl.AP.767260_TMYx.2004-2018.epw")
MATERIALS = os.path.join(HERE, "materials.ini")

LAYERS = [("EPS", 0.1)]
# Same 0.1 m of EPS, built as the one supported homogeneous 2D section: a
# HollowBlock whose cavity is FILLED with the shell material (0.03+0.04+0.03).
GEOM = {"web": 0.02, "block_width": 0.16,
        "cover_top": 0.03, "cavity": 0.04, "cover_bottom": 0.03}


def _setup():
    config.file = MATERIALS
    config2d.ny = config.Nx       # mismo nº de nodos en el espesor que el 1D
    config2d.nx = 5               # ancho arbitrario (homogéneo en x -> cancela)
    loc = Location(EPW)
    return loc


_CACHE = {}


def _run():
    if "r" in _CACHE:
        return _CACHE["r"]
    loc = _setup()
    sys1 = System(loc, tilt=90, azimuth=0, absortance=0.8, layers=LAYERS)
    Ti1 = sys1.solve().to_numpy()

    block = HollowBlock("EPS", fill_type=Fill.SOLID, fill_material="EPS",
                        geometry=GEOM)
    sys2 = System2D(loc, layers=[block], tilt=90, azimuth=0, absortance=0.8)
    Ti2 = sys2.solve()
    _CACHE["r"] = (Ti1, Ti2, sys1, sys2)
    return _CACHE["r"]


def test_reduces_to_1d():
    Ti1, Ti2, _, _ = _run()
    d = np.abs(Ti1 - Ti2)
    assert d.max() <= 0.1, f"2D no reduce al 1D: max|Δ|={d.max():.3e} °C"


def test_periodicity():
    _, _, _, sys2 = _run()
    assert 0 < sys2.days < config2d.max_days, \
        f"no convergió día-a-día (days={sys2.days})"


def test_energy_balance():
    _, _, _, sys2 = _run()
    rel = sys2.energy_imbalance
    assert rel <= 0.02, f"Qin≠Qout en régimen: desbalance = {rel:.1%}"


def _demo():
    Ti1, Ti2, sys1, sys2 = _run()
    d = np.abs(Ti1 - Ti2)
    n = len(Ti1)
    bar = "═" * 70
    print(bar)
    print("  FASE 5 · API 2D de producción (System2D, EPW+pvlib, JIT)")
    print(bar)
    print(f"\n  Caso: capa homogénea {LAYERS} (sin bovedilla), "
          f"malla {config2d.nx}×{config2d.ny}")
    print(f"  Convergió en {sys2.days} días (tope {config2d.max_days})")

    print(f"\n  Reducción al 1D del paquete (Ti, {n} muestras):")
    print(f"    max|Ti_2D − Ti_1D| = {d.max():.3e} °C   media = {d.mean():.3e} °C"
          f"   {'✓' if d.max()<=0.1 else '✗'}  (atol 0.1)")

    qin, rel = sys2.energy_transfer, sys2.energy_imbalance
    print(f"\n  Balance de energía en régimen:")
    print(f"    Qin = {qin:.1f} (J/m²·día)   "
          f"desbalance = {rel:.2%}   {'✓' if rel<=0.02 else '✗'}")

    # Curva Ti(t): 1D vs 2D superpuestas (deben verse iguales).
    print(f"\n{bar}\n  Ti(t): '1'=1D  '2'=2D  (se superponen -> '#')")
    print(bar)
    lo, hi = min(Ti1.min(), Ti2.min()), max(Ti1.max(), Ti2.max())
    W = 56
    step = max(1, n // 36)
    for r in range(0, n, step):
        h = r * config.dt / 3600.0
        def pos(v):
            return int((v - lo) / (hi - lo) * (W - 1))
        line = [" "] * W
        p1, p2 = pos(Ti1[r]), pos(Ti2[r])
        line[p1] = "1"
        line[p2] = "#" if p2 == p1 else "2"
        print(f"  {h:4.1f}h |{''.join(line)}|")
    print(f"  rango {lo:.1f}..{hi:.1f} °C")

    print(f"\n{bar}")
    print("  Fase 5: System2D reduce al 1D, converge y balancea energía ✅")
    print(bar)


if __name__ == "__main__":
    for fn in (test_reduces_to_1d, test_periodicity, test_energy_balance):
        fn()
    _demo()
