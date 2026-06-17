"""
Inspector a escala de `System2D`: `section()`, `preview()` (mpl + ASCII), `section_report()`.

Verifica que se pueda revisar la asignación de materiales sin resolver.

    .venv/bin/python tests/test_eh2d_inspect.py
"""

import os

os.environ.setdefault("MPLBACKEND", "Agg")   # headless para matplotlib

import numpy as np

import enerhabitat as eh
from enerhabitat.config import config2d
from enerhabitat.eh2d import _categorize

HERE = os.path.dirname(os.path.abspath(__file__))
EPW = os.path.join(HERE, "MEX_MOR_Cuernavaca-Matamoros.Intl.AP.767260_TMYx.2004-2018.epw")
MATERIALS = os.path.join(HERE, "materials_2d.ini")


def _wall():
    eh.config.file = MATERIALS
    config2d.nx, config2d.ny = 40, 100
    loc = eh.Location(EPW)
    block = eh.HollowBlock("Concreto", emissivity=0.9, geometry={
        "web": 0.02, "block_width": 0.16,
        "cover_top": 0.02, "cavity": 0.08, "cover_bottom": 0.02})
    wall = eh.System2D(location=loc)
    wall.tilt = 90
    wall.layers = [("Aplanado", 0.02), block, ("Yeso", 0.01)]
    return wall


def test_section_arrays():
    w = _wall()
    sec = w.section()                      # no resuelve
    assert sec.NT.shape == (config2d.nx, config2d.ny)
    assert sec.kfield.shape == sec.NT.shape
    assert 0 in np.unique(sec.NT)          # hay aire del hueco


def test_categories_include_air_and_materials():
    w = _wall()
    sec = w.section()
    _, cats = _categorize(sec, "materials", eh.config.materials)
    texts = [t for t, _ in cats]
    assert any("Air" in t for t in texts)
    assert any("Concreto" in t for t in texts)


def test_preview_ascii_runs():
    w = _wall()
    assert w.preview(field="materials", backend="ascii") is None


def test_preview_mpl_saves(tmp_path=None):
    w = _wall()
    out = "/tmp/_eh2d_inspect_test.png"
    fig, axs = w.preview(panels=["nodetype", "k", "rhoc"], backend="mpl", save=out)
    assert os.path.exists(out)
    assert len(axs) == 3


def test_section_report_runs():
    _wall().section_report()               # no debe lanzar


if __name__ == "__main__":
    for fn in (test_section_arrays, test_categories_include_air_and_materials,
               test_preview_ascii_runs, test_preview_mpl_saves, test_section_report_runs):
        fn()
        print(f"PASS  {fn.__name__}")
    print("\nInspector a escala: OK ✅")
