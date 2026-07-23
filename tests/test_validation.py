"""
P1 — Validación dirigida de rangos (lo que corrompe la física en silencio).

Valida que los valores no-físicos lanzan ``ValueError``/``KeyError`` claros en
vez de aceptarse:

  1. ``absortance`` fuera de [0,1] (1D y 2D).
  2. ``config``: ``La/ho/hi/aire`` ≤ 0 y ``Nx`` no entero o < 3.
  3. ``config2d.validate()``: malla, tolerancias y topes fuera de rango.
  4. Capas: espesor ≤ 0 y material inexistente (mensaje con disponibles).
  5. Materiales con ``k/rho/c ≤ 0`` rechazados al cargar el ``.ini``
     (conservando los materiales previos — rollback).

Correr con pytest o como script:
    .venv/bin/python tests/test_validation.py
"""

import os
import tempfile

from enerhabitat.config import Config, Config2D
from enerhabitat.ehframe import System
from enerhabitat.eh2d import System2D
from enerhabitat.ehtools import set_construction

HERE = os.path.dirname(os.path.abspath(__file__))
MATERIALS = os.path.join(HERE, "materials.ini")


def _raises(exc, fn, *args, **kw):
    try:
        fn(*args, **kw)
    except exc:
        return True
    raise AssertionError(f"{fn} debía lanzar {exc.__name__}")


def test_absortance_range():
    s1 = System(None, layers=[])
    _raises(ValueError, setattr, s1, "absortance", 1.5)
    _raises(ValueError, setattr, s1, "absortance", -0.1)
    s1.absortance = 0.8   # válido
    s2 = System2D(None)
    _raises(ValueError, setattr, s2, "absortance", 1.5)
    s2.absortance = 1.0   # frontera válida


def test_config_setters():
    c = Config()
    for attr, bad in (("La", 0), ("La", -1), ("ho", 0), ("hi", -8),
                      ("Nx", 2), ("Nx", 10.5),
                      ("AIR_DENSITY", 0), ("AIR_HEAT_CAPACITY", -1)):
        _raises(ValueError, setattr, c, attr, bad)
    c.Nx = 200.0          # entero disfrazado de float: aceptado como int
    assert c.Nx == 200 and isinstance(c.Nx, int)


def test_config2d_validate():
    c2 = Config2D()
    c2.validate()          # defaults válidos
    for attr, bad in (("nx", 2), ("ny", 0), ("nx", 10.5),
                      ("tol_inner", 0), ("tol_day", -1e-4),
                      ("max_days", 0), ("max_inner", 0.5)):
        good = getattr(c2, attr)
        setattr(c2, attr, bad)
        _raises(ValueError, c2.validate)
        setattr(c2, attr, good)


def test_layers_validation():
    c = Config()
    c.file = MATERIALS
    _raises(ValueError, set_construction, c.materials, [("EPS", 0.0)])
    _raises(ValueError, set_construction, c.materials, [("EPS", -0.05)])
    try:
        set_construction(c.materials, [("NoExiste", 0.1)])
    except KeyError as e:
        assert "NoExiste" in str(e) and "EPS" in str(e)
    else:
        raise AssertionError("material inexistente debía lanzar KeyError")


def test_material_positivity_on_load():
    c = Config()
    c.file = MATERIALS
    before = dict(c.materials)
    with tempfile.NamedTemporaryFile("w", suffix=".ini", delete=False) as f:
        f.write("[Malo]\nk = 0\nrho = 100\nc = 1000\n")
        bad_ini = f.name
    try:
        _raises(ValueError, setattr, c, "file", bad_ini)
        assert c.materials == before          # rollback intacto
        assert c.file == MATERIALS
    finally:
        os.unlink(bad_ini)


if __name__ == "__main__":
    for fn in (test_absortance_range, test_config_setters,
               test_config2d_validate, test_layers_validation,
               test_material_positivity_on_load):
        fn()
        print(f"PASS  {fn.__name__}")
    print("\nP1: validación dirigida — rangos físicos protegidos ✅")
