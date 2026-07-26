"""
hi dependiente de la direccion del flujo en techos (NOM-020).

Con ``config.hi_flow = True`` (default) un techo (tilt < 60) selecciona en cada
paso de tiempo ``hi_down = 6.6`` (superficie interior mas caliente que el aire,
estratificacion estable) o ``hi_up = 9.4`` (flujo ascendente); los muros siempre
usan ``hi = 8.1``. Se verifica:

  - **muro invariante**: tilt=90 da lo mismo con hi_flow True/False;
  - **techo sensible**: tilt=0 cambia las energias AC entre True/False, en el
    sentido fisico esperado (menos cooling: de dia domina hi_down < hi);
  - **fijar hi**: hi_up = hi_down = hi reproduce hi_flow=False en techos;
  - **2D consistente**: un Slab (techo) 2D tambien cambia con hi_flow;
  - **setters**: hi_up/hi_down validan > 0.

    .venv/bin/python tests/test_hi_flow.py
"""

import os

import pytest

import enerhabitat as eh
from enerhabitat.config import config, config2d

HERE = os.path.dirname(os.path.abspath(__file__))
EPW = os.path.join(HERE, "MEX_MOR_Cuernavaca-Matamoros.Intl.AP.767260_TMYx.2004-2018.epw")
MATERIALS = os.path.join(HERE, "materials_2d.ini")

LAYERS = [("Mortero", 0.025), ("Concreto", 0.10)]


@pytest.fixture(autouse=True)
def _restore_config():
    yield
    config.reset()


def _solve_1d_ac(tilt):
    config.file = MATERIALS
    s = eh.System(eh.Location(EPW))
    s.tilt = tilt
    if tilt != 0:
        s.azimuth = 90
    s.absortance = 0.3
    s.layers = list(LAYERS)
    s.location.meanDay(month=5, year=2025)
    s.solveAC()
    return s.cooling_energy, s.heating_energy


def test_wall_unaffected():
    config.reset()
    config.hi_flow = True
    on = _solve_1d_ac(90)
    config.hi_flow = False
    off = _solve_1d_ac(90)
    assert on == off


def test_roof_flow_dependent():
    config.reset()
    config.hi_flow = True
    qc_on, qh_on = _solve_1d_ac(0)
    config.hi_flow = False
    qc_off, qh_off = _solve_1d_ac(0)
    # daytime: hot ceiling -> stable stratification (hi_down < hi) -> less cooling
    assert qc_on < qc_off
    # nighttime: cold ceiling -> enhanced convection (hi_up > hi) -> more heating
    assert qh_on > qh_off


def test_pinning_equals_fixed():
    config.reset()
    config.hi_up = config.hi_down = config.hi
    pinned = _solve_1d_ac(0)
    config.reset()
    config.hi_flow = False
    fixed = _solve_1d_ac(0)
    assert pinned == pytest.approx(fixed, rel=1e-12)


def test_roof_2d_flow_dependent():
    config.reset()
    config.file = MATERIALS
    config2d.nx, config2d.ny = 24, 60
    config2d.max_days = 30

    def solve():
        slab = eh.Slab(
            rib_material="Concreto", block_material="Concreto",
            topping_material="Concreto", fill_type=eh.Fill.SOLID,
            fill_material="Concreto",
            geometry={"web": 0.025, "foot": 0.025, "shoulder": 0.050,
                      "n_cavities": 3, "cavity_width": 0.103,
                      "topping": 0.100, "topping_cap": 0.050,
                      "cover_top": 0.030, "cavity": 0.040,
                      "cover_bottom": 0.030},
        )
        roof = eh.System2D(eh.Location(EPW), tilt=0, absortance=0.3)
        roof.layers = [slab]
        roof.location.meanDay(month=5, year=2025)
        roof.solveAC()
        return roof.cooling_energy, roof.heating_energy

    config.hi_flow = True
    qc_on, qh_on = solve()
    config.hi_flow = False
    qc_off, qh_off = solve()
    assert qc_on < qc_off
    assert qh_on > qh_off


def test_setter_validation():
    with pytest.raises(ValueError):
        config.hi_up = 0
    with pytest.raises(ValueError):
        config.hi_down = -1


if __name__ == "__main__":
    test_wall_unaffected()
    test_roof_flow_dependent()
    test_pinning_equals_fixed()
    test_roof_2d_flow_dependent()
    print("hi_flow OK")
