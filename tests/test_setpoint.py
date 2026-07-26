"""
Setpoint de AC en el System 1D (simetria con System2D).

Con ``setpoint = None`` (default) solveAC() mantiene Ti en Tn.mean(); con un
valor, lo mantiene en ese valor. Se verifica:

  - default: Ti constante = Tn.mean();
  - override: Ti constante = setpoint, y las energias cambian;
  - cache: cambiar setpoint invalida el resultado anterior;
  - solve() (free-running) ignora el setpoint.

    .venv/bin/python tests/test_setpoint.py
"""

import os

import pytest

import enerhabitat as eh
from enerhabitat.config import config

HERE = os.path.dirname(os.path.abspath(__file__))
EPW = os.path.join(HERE, "MEX_MOR_Cuernavaca-Matamoros.Intl.AP.767260_TMYx.2004-2018.epw")
MATERIALS = os.path.join(HERE, "materials_2d.ini")


def _wall():
    config.file = MATERIALS
    w = eh.System(eh.Location(EPW))
    w.tilt = 90
    w.azimuth = 90
    w.absortance = 0.3
    w.layers = [("Mortero", 0.025), ("Concreto", 0.10)]
    w.location.meanDay(month=5, year=2025)
    return w


def test_default_holds_tn():
    w = _wall()
    ti = w.solveAC()
    tn = float(w.Tsa()["Tn"].mean())
    assert ti.nunique() == 1
    assert ti.iloc[0] == pytest.approx(tn, abs=1e-9)


def test_setpoint_overrides_and_invalidates_cache():
    w = _wall()
    ti_tn = w.solveAC()
    qc_tn = w.cooling_energy

    w.setpoint = 24.0
    ti_24 = w.solveAC()                 # el cambio de setpoint debe recomputar
    assert ti_24.nunique() == 1
    assert ti_24.iloc[0] == pytest.approx(24.0, abs=1e-9)
    assert w.cooling_energy != qc_tn    # otra demanda a otro setpoint

    w.setpoint = None
    ti_back = w.solveAC()               # regresar a None recomputa a Tn
    assert ti_back.iloc[0] == pytest.approx(ti_tn.iloc[0], abs=1e-9)


def test_free_running_ignores_setpoint():
    w = _wall()
    ti_free = w.solve()
    w2 = _wall()
    w2.setpoint = 24.0
    ti_free2 = w2.solve()
    assert (ti_free.to_numpy() == ti_free2.to_numpy()).all()


if __name__ == "__main__":
    test_default_holds_tn()
    test_setpoint_overrides_and_invalidates_cache()
    test_free_running_ignores_setpoint()
    print("setpoint 1D OK")
