"""
Cache por firmas en System (1D) y System2D.

El 1D usa el mismo patron que el 2D: una firma (tupla con los inputs actuales)
se compara en cada llamada, con cache separado por modo. Se verifica:

  - hit: repetir solve()/Tsa() sin cambios devuelve el MISMO objeto;
  - alternancia: solve() -> solveAC() -> solve() reusa cada cache y RESTAURA
    los escalares del modo (energy_transfer/cooling_energy correctos);
  - invalidacion: cambiar cualquier input (tilt, azimuth, absortance, layers
    por asignacion Y por mutacion in situ, config, meanDay) recomputa;
  - 2D: la alternancia tambien restaura escalares (fix del hit sin restore).

    .venv/bin/python tests/test_cache.py
"""

import os

import pytest

import enerhabitat as eh
from enerhabitat.config import config, config2d

HERE = os.path.dirname(os.path.abspath(__file__))
EPW = os.path.join(HERE, "MEX_MOR_Cuernavaca-Matamoros.Intl.AP.767260_TMYx.2004-2018.epw")
MATERIALS = os.path.join(HERE, "materials_2d.ini")


@pytest.fixture(autouse=True)
def _restore_config():
    yield
    config.reset()


def _wall():
    config.file = MATERIALS
    w = eh.System(eh.Location(EPW))
    w.tilt = 90
    w.azimuth = 90
    w.absortance = 0.3
    w.layers = [("Mortero", 0.025), ("Concreto", 0.10)]
    w.location.meanDay(month=5, year=2025)
    return w


def test_hit_returns_same_object():
    w = _wall()
    tsa1 = w.Tsa()
    assert w.Tsa() is tsa1
    ti1 = w.solve()
    assert w.solve() is ti1


def test_alternating_modes_cached_and_scalars_restored():
    w = _wall()
    ti_free = w.solve()
    e_free = w.energy_transfer
    w.solveAC()
    qc = w.cooling_energy
    assert w.energy_transfer is None        # los escalares son del ultimo modo

    ti_free2 = w.solve()
    assert ti_free2 is ti_free              # hit: no recomputo
    assert w.energy_transfer == e_free      # ...y los escalares regresaron
    assert w.cooling_energy is None

    w.solveAC()
    assert w.cooling_energy == qc


def test_every_input_invalidates():
    w = _wall()

    ti = w.solve()
    w.tilt = 0
    assert w.solve() is not ti

    w.tilt = 90
    ti = w.solve()
    w.azimuth = 180
    assert w.solve() is not ti

    w.azimuth = 90
    ti = w.solve()
    w.absortance = 0.6
    assert w.solve() is not ti

    w.absortance = 0.3
    ti = w.solve()
    w.layers = [("Mortero", 0.025), ("Concreto", 0.12)]
    assert w.solve() is not ti

    ti = w.solve()
    w.layers.append(("Yeso", 0.01))         # mutacion in situ: antes era el footgun
    assert w.solve() is not ti

    ti = w.solve()
    config.ho = 12
    assert w.solve() is not ti

    config.reset()
    ti = w.solve()
    w.location.meanDay(month=6, year=2025)
    assert w.solve() is not ti


def test_2d_alternating_scalars_restored():
    config.file = MATERIALS
    config2d.nx, config2d.ny = 12, 30
    config2d.max_days = 30
    block = eh.HollowBlock("Concreto", fill_type=eh.Fill.SOLID,
                           fill_material="EPS",
                           geometry={"web": 0.02, "block_width": 0.16,
                                     "cover_top": 0.03, "cavity": 0.04,
                                     "cover_bottom": 0.03})
    w = eh.System2D(eh.Location(EPW), tilt=90, azimuth=90, absortance=0.3)
    w.layers = [block]
    w.location.meanDay(month=5, year=2025)

    ti = w.solve()
    e = w.energy_transfer
    w.solveAC()
    qc = w.cooling_energy
    assert w.energy_transfer is None

    ti2 = w.solve()
    assert ti2 is ti                        # hit
    assert w.energy_transfer == e           # escalares restaurados (fix)
    assert w.cooling_energy is None

    w.solveAC()
    assert w.cooling_energy == qc


if __name__ == "__main__":
    test_hit_returns_same_object()
    test_alternating_modes_cached_and_scalars_restored()
    test_every_input_invalidates()
    test_2d_alternating_scalars_restored()
    print("cache OK")
