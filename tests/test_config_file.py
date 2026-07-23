"""
P1 — Manejo explícito del archivo de materiales (`config.file`).

Valida el contrato corregido (antes: errores por `print`, atributos sin crear
y getter que podía crashar con AttributeError):

  1. Instanciar `Config` sin `materials.ini` en el cwd no imprime ni falla, y
     deja `materials == {}` (lo que la docs siempre prometió).
  2. Asignar una ruta inexistente lanza `FileNotFoundError` de inmediato.
  3. El fallo conserva intactos la ruta y los materiales previos (rollback).
  4. El getter de `file` nunca hace I/O ni devuelve None: siempre la ruta
     configurada.

Correr con pytest o como script:
    .venv/bin/python tests/test_config_file.py
"""

import io
import os
import contextlib

from enerhabitat.config import Config

HERE = os.path.dirname(os.path.abspath(__file__))
MATERIALS = os.path.join(HERE, "materials.ini")


def _fresh_config_without_default():
    """Config() instanciada en un cwd sin materials.ini, capturando stdout."""
    cwd = os.getcwd()
    tmp = os.path.join(HERE, "__marimo__")   # dir existente sin materials.ini
    os.chdir(tmp)
    try:
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            cfg = Config()
        return cfg, out.getvalue()
    finally:
        os.chdir(cwd)


def test_init_without_default_is_silent_and_empty():
    cfg, printed = _fresh_config_without_default()
    assert printed == "", f"no debe imprimir nada, imprimió: {printed!r}"
    assert cfg.materials == {}
    assert cfg.file == "materials.ini"   # getter sin I/O, sin None


def test_missing_file_raises():
    cfg, _ = _fresh_config_without_default()
    try:
        cfg.file = "no_existe.ini"
    except FileNotFoundError:
        pass
    else:
        raise AssertionError("ruta inexistente debe lanzar FileNotFoundError")


def test_failed_assignment_keeps_previous_materials():
    cfg, _ = _fresh_config_without_default()
    cfg.file = MATERIALS
    assert "EPS" in cfg.materials
    before_file = cfg.file
    before_mats = dict(cfg.materials)
    try:
        cfg.file = "no_existe.ini"
    except FileNotFoundError:
        pass
    assert cfg.file == before_file
    assert cfg.materials == before_mats


if __name__ == "__main__":
    for fn in (test_init_without_default_is_silent_and_empty,
               test_missing_file_raises,
               test_failed_assignment_keeps_previous_materials):
        fn()
        print(f"PASS  {fn.__name__}")
    print("\nP1: config.file falla explícito, rollback y getter sin I/O ✅")
