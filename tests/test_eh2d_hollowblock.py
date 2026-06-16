"""
Fase 8a — Bloque hueco de concreto en muros (API `System2D` + `HollowBlock`).

Verifica el cableado de producción end-to-end, con la MISMA metodología que el
`System` 1D (config.file → Location → props → meanDay → Tsa → solve):

  - **metodología/flujo**: `solve()` devuelve una `pandas.Series` alineada a `Tsa()`;
  - **periodicidad**: converge día-a-día antes de `max_days`;
  - **balance de energía**: en régimen `Qin ≈ Qout`;
  - **orientación**: `HollowBlock` exige `tilt=90` (techo falla);
  - **capas antes/después**: el muro con capas a ambos lados del bloque resuelve y
    difiere de un muro con más aislante (sanidad).

Malla 2D reducida (`config2d`) para que la prueba corra rápido.

    .venv/bin/python tests/test_eh2d_hollowblock.py
"""

import os

import numpy as np
import pandas as pd

import enerhabitat as eh
from enerhabitat.config import config2d

HERE = os.path.dirname(os.path.abspath(__file__))
EPW = os.path.join(HERE, "MEX_MOR_Cuernavaca-Matamoros.Intl.AP.767260_TMYx.2004-2018.epw")
MATERIALS = os.path.join(HERE, "materials_2d.ini")

GEOM = {"web": 0.02, "block_width": 0.16,
        "cover_top": 0.02, "cavity": 0.08, "cover_bottom": 0.02}


def _setup():
    eh.config.file = MATERIALS
    config2d.nx, config2d.ny = 24, 60     # malla chica: prueba de API, no de precisión
    config2d.max_days = 30
    loc = eh.Location(EPW)
    loc.meanDay(month=5, year=2025)
    return loc


def _wall():
    loc = _setup()
    block = eh.HollowBlock("Concreto", emissivity=0.9, geometry=GEOM)
    wall = eh.System2D(location=loc)
    wall.tilt = 90
    wall.azimuth = 90
    wall.absortance = 0.6
    wall.layers = [("Mortero", 0.02), block, ("Yeso", 0.01)]
    return wall


_CACHE = {}


def _solved():
    if "w" not in _CACHE:
        w = _wall()
        ti = w.solve()
        _CACHE["w"] = (w, ti)
    return _CACHE["w"]


# --- pruebas -------------------------------------------------------------------

def test_methodology_returns_series():
    w, ti = _solved()
    assert isinstance(ti, pd.Series)
    assert len(ti) == len(w.Tsa())
    assert ti.index.equals(w.Tsa().index)
    assert np.isfinite(ti.to_numpy()).all()


def test_periodicity():
    w, _ = _solved()
    assert 0 < w.days < config2d.max_days


def test_energy_balance():
    w, _ = _solved()
    qin, qout = w.energy_transfer, w.Qout
    rel = abs(qin - qout) / max(abs(qin), 1e-9)
    assert rel <= 0.02, f"Qin={qin:.1f} Qout={qout:.1f} ({rel:.1%})"


def test_orientation_guard():
    loc = _setup()
    block = eh.HollowBlock("Concreto", emissivity=0.9, geometry=GEOM)
    roof = eh.System2D(location=loc)
    roof.tilt = 0                              # bloque hueco NO es para techo
    roof.layers = [("Mortero", 0.02), block, ("Yeso", 0.01)]
    try:
        roof.solve()
    except ValueError:
        pass
    else:
        raise AssertionError("se esperaba ValueError por tilt≠90")


def test_requires_one_element():
    loc = _setup()
    w = eh.System2D(location=loc, tilt=90)
    w.layers = [("Mortero", 0.02), ("Yeso", 0.01)]   # sin elemento 2D
    try:
        w.solve()
    except ValueError:
        pass
    else:
        raise AssertionError("se esperaba ValueError por falta de elemento 2D")


def _demo():
    w, ti = _solved()
    bar = "═" * 70
    print(bar)
    print("  FASE 8a · Muro con bloque hueco de concreto (System2D + HollowBlock)")
    print(bar)
    w.info()
    df = w.solve_dataframe
    Tsa = df["Tsa"].to_numpy()
    Ti = df["Ti"].to_numpy()
    Th = df["Thueco"].to_numpy()
    print(f"\n  Convergió en {w.days} días (tope {config2d.max_days})")
    qin, qout = w.energy_transfer, w.Qout
    rel = abs(qin - qout) / max(abs(qin), 1e-9)
    print(f"  Balance: Qin={qin:.1f}  Qout={qout:.1f}  desbalance={rel:.2%}")
    dec = (Ti.max() - Ti.min()) / (Tsa.max() - Tsa.min())
    print(f"  Tint rango {Ti.min():.2f}..{Ti.max():.2f} °C  ·  factor de decremento {dec:.3f}")
    print(f"  Thueco (aire del bloque) rango {Th.min():.2f}..{Th.max():.2f} °C")

    print(f"\n{bar}\n  Curva diaria: 's'=Tsa  'h'=Thueco(bloque)  'i'=Tint\n{bar}")
    n = len(Ti)
    lo, hi = min(Tsa.min(), Ti.min()), max(Tsa.max(), Ti.max())
    W = 54
    hrs = np.arange(n) * eh.config.dt / 3600.0
    for r in range(0, n, max(1, n // 24)):
        line = [" "] * W
        for v, ch in [(Tsa[r], "s"), (Th[r], "h"), (Ti[r], "i")]:
            line[int((v - lo) / (hi - lo) * (W - 1))] = ch
        print(f"  {hrs[r]:4.1f}h |{''.join(line)}|")
    print(f"  rango {lo:.1f}..{hi:.1f} °C")
    print(f"\n{bar}\n  Fase 8a: bloque hueco de concreto en muro, end-to-end ✅\n{bar}")


if __name__ == "__main__":
    for fn in (test_methodology_returns_series, test_periodicity, test_energy_balance,
               test_orientation_guard, test_requires_one_element):
        fn()
        print(f"PASS  {fn.__name__}")
    print()
    _demo()
