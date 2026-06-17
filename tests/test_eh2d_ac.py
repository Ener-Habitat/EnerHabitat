"""
Fase 9 — Aire acondicionado (AC) en los sistemas 2D (`System2D.solveAC`).

Espejo del `System.solveAC()` 1D: mantiene el aire interior en un setpoint y
calcula la energía de enfriamiento/calentamiento. Se verifica con la misma
metodología que el resto del 2D (config.file → Location → meanDay → Tsa → solveAC),
NO contra el C:

  - **reduce al 1D**: un techo `Slab` RELLENA de un solo material (homogéneo en x)
    reproduce `cooling_energy`/`heating_energy` del `System.solveAC` 1D equivalente;
  - **metodología/sanidad**: `Ti` constante = setpoint; `Qcool≥0`, `Qheat≥0`;
  - **periodicidad**: converge día-a-día antes de `max_days`;
  - **AIRE vs RELLENA**: la cámara de aire cambia la carga (más decremento → más
    `Qcool` que un relleno aislante en clima cálido).

Malla 2D reducida (`config2d`), motor serial.

    .venv/bin/python tests/test_eh2d_ac.py
"""

import os

import numpy as np
import pandas as pd

import enerhabitat as eh
from enerhabitat import System
from enerhabitat.config import config, config2d

HERE = os.path.dirname(os.path.abspath(__file__))
EPW = os.path.join(HERE, "MEX_MOR_Cuernavaca-Matamoros.Intl.AP.767260_TMYx.2004-2018.epw")
MATERIALS = os.path.join(HERE, "materials_2d.ini")

GEOM_WALL = {"web": 0.02, "block_width": 0.16,
             "cover_top": 0.02, "cavity": 0.08, "cover_bottom": 0.02}
GEOM_ROOF = {"web": 0.025, "foot": 0.025, "shoulder": 0.050, "n_cavities": 3,
             "cavity_width": 0.103, "topping": 0.100, "topping_cap": 0.050,
             "cover_top": 0.030, "cavity": 0.040, "cover_bottom": 0.030}


def _setup():
    config.file = MATERIALS
    config2d.nx, config2d.ny = 24, 60
    config2d.max_days = 30
    loc = eh.Location(EPW)
    loc.meanDay(month=5, year=2025)
    return loc


_CACHE = {}


def _wall_hueco():
    if "wall" not in _CACHE:
        loc = _setup()
        block = eh.HollowBlock("Concreto", emissivity=0.9, geometry=GEOM_WALL)
        w = eh.System2D(loc, tilt=90, azimuth=90, absortance=0.6)
        w.layers = [("Mortero", 0.02), block, ("Yeso", 0.01)]
        w.solveAC()
        _CACHE["wall"] = w
    return _CACHE["wall"]


def _roof(key, bovedilla, fill_material=None):
    if key not in _CACHE:
        loc = _setup()
        slab = eh.Slab("Concreto", fill_type=bovedilla, block_material="Bovedilla",
                       topping_material="Concreto", fill_material=fill_material,
                       emissivity=0.9, geometry=GEOM_ROOF)
        r = eh.System2D(loc, tilt=0, azimuth=0, absortance=0.3)
        r.layers = [("Aplanado", 0.003), slab, ("Yeso", 0.015)]
        r.solveAC()
        _CACHE[key] = r
    return _CACHE[key]


# --- pruebas -------------------------------------------------------------------

def test_methodology_and_sanity():
    w = _wall_hueco()
    ti = w.solveAC()
    assert isinstance(ti, pd.Series)
    assert len(ti) == len(w.Tsa()) and ti.index.equals(w.Tsa().index)
    # Ti constante = setpoint (AC mantiene el recinto)
    assert float(ti.max() - ti.min()) < 1e-9
    assert np.isfinite(ti.to_numpy()).all()
    assert w.cooling_energy >= 0.0 and w.heating_energy >= 0.0
    assert w.energy_transfer is None


def test_periodicity():
    w = _wall_hueco()
    assert 0 < w.days < config2d.max_days


def test_reduces_to_1d():
    # Techo Slab RELLENA con UN solo material → homogéneo en x → debe igualar al 1D.
    loc = _setup()
    config2d.ny = config.Nx                      # mismo nº de nodos en el espesor que el 1D
    mat, Lt = "Concreto", 0.20
    geom = {"web": 0.02, "foot": 0.02, "shoulder": 0.05, "n_cavities": 2,
            "cavity_width": 0.10, "topping": 0.10, "topping_cap": 0.0,
            "cover_top": 0.03, "cavity": 0.04, "cover_bottom": 0.03}   # suma 0.20
    slab = eh.Slab(mat, fill_type=eh.Fill.SOLID, block_material=mat,
                   topping_material=mat, fill_material=mat, geometry=geom)
    r2 = eh.System2D(loc, tilt=0, azimuth=0, absortance=0.3, layers=[slab])
    r2.solveAC()

    s1 = System(loc, tilt=0, azimuth=0, absortance=0.3, layers=[(mat, Lt)])
    s1.Tsa()
    s1.solveAC()

    def rel(x, y):
        return abs(x - y) / max(abs(x), abs(y), 1e-9)
    rc = rel(r2.cooling_energy, s1.cooling_energy)
    rh = rel(r2.heating_energy, s1.heating_energy)
    assert rc <= 0.05, f"cooling 2D={r2.cooling_energy:.1f} vs 1D={s1.cooling_energy:.1f} ({rc:.1%})"
    assert rh <= 0.05, f"heating 2D={r2.heating_energy:.1f} vs 1D={s1.heating_energy:.1f} ({rh:.1%})"
    config2d.ny = 60


def test_air_vs_insulating_fill():
    ra = _roof("roof_aire", eh.Fill.AIR)
    rf = _roof("roof_eps", eh.Fill.SOLID, fill_material="EPS")
    # La cámara de aire transfiere más que el relleno aislante → mayor carga de enfriamiento.
    assert ra.cooling_energy > rf.cooling_energy, \
        f"Qcool aire={ra.cooling_energy:.1f} debe superar relleno EPS={rf.cooling_energy:.1f}"


def _demo():
    w = _wall_hueco()
    ra = _roof("roof_aire", eh.Fill.AIR)
    rf = _roof("roof_eps", eh.Fill.SOLID, fill_material="EPS")
    bar = "═" * 70
    print(bar)
    print("  FASE 9 · Aire acondicionado en sistemas 2D (System2D.solveAC)")
    print(bar)
    w.info()
    sp = float(w.solveAC().iloc[0])
    print(f"\n  Setpoint (Tn.mean): {sp:.2f} °C   ·   convergió en {w.days} días")
    print(f"  MURO bloque hueco:  Qcool={w.cooling_energy:.1f}  Qheat={w.heating_energy:.1f} (J/m²·día)")
    print(f"  TECHO bovedilla AIRE:    Qcool={ra.cooling_energy:.1f}  Qheat={ra.heating_energy:.1f}")
    print(f"  TECHO bovedilla EPS:     Qcool={rf.cooling_energy:.1f}  Qheat={rf.heating_energy:.1f}")
    print(f"  → AIRE transfiere más que EPS: Qcool {ra.cooling_energy:.1f} > {rf.cooling_energy:.1f}  "
          f"{'✓' if ra.cooling_energy>rf.cooling_energy else '✗'}")
    print(f"\n{bar}\n  Fase 9: solveAC 2D (muro y techo), espejo del 1D ✅\n{bar}")


if __name__ == "__main__":
    for fn in (test_methodology_and_sanity, test_periodicity,
               test_reduces_to_1d, test_air_vs_insulating_fill):
        fn()
        print(f"PASS  {fn.__name__}")
    print()
    _demo()
