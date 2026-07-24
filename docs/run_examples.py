"""Pre-computes the 2D examples of the documentation (usage.qmd).

A 2D solve with the default mesh takes tens of minutes, so the 2D examples are
NOT executed when the site is rendered. This script runs them once and stores
the results in docs/data/results/ (one CSV per case + summary.json with the
energies, convergence days and runtimes); the usage.qmd chunks only read those
files and plot.

Usage (from the repo root; the three cases run in parallel processes):
    uv run python docs/run_examples.py
"""
import json
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
OUT = DATA / "results"
EPW = str(DATA / "MEX_MOR_Cuernavaca-Matamoros.Intl.AP.767260_TMYx.2004-2018.epw")
MONTH, YEAR = 5, 2025


def _hollow_wall():
    import enerhabitat as eh
    eh.config.file = str(DATA / "materials.ini")
    block = eh.HollowBlock(
        material="Concreto", emissivity=0.9,
        geometry={"web": 0.02, "block_width": 0.16,
                  "cover_top": 0.02, "cavity": 0.08, "cover_bottom": 0.02},
    )
    wall = eh.System2D(eh.Location(EPW))
    wall.tilt = 90
    wall.azimuth = 90
    wall.absortance = 0.6
    wall.layers = [("Mortero", 0.02), block, ("Yeso", 0.01)]
    wall.location.meanDay(month=MONTH, year=YEAR)
    wall.Tsa()
    return wall


def _slab_roof():
    import enerhabitat as eh
    eh.config.file = str(DATA / "materials.ini")
    slab = eh.Slab(
        rib_material="ConcretoAltaDensidad", block_material="Bovedilla",
        topping_material="Concreto", fill_type=eh.Fill.AIR, emissivity=0.9,
        geometry={"web": 0.025, "foot": 0.025, "shoulder": 0.050,
                  "n_cavities": 3, "cavity_width": 0.103,
                  "topping": 0.100, "topping_cap": 0.050,
                  "cover_top": 0.030, "cavity": 0.040, "cover_bottom": 0.030},
    )
    roof = eh.System2D(eh.Location(EPW))
    roof.tilt = 0
    roof.absortance = 0.3
    roof.layers = [("Impermeabilizante", 0.003), slab, ("Yeso", 0.015)]
    roof.location.meanDay(month=MONTH, year=YEAR)
    roof.Tsa()
    return roof


def run_case(name):
    import numpy as np
    t0 = time.perf_counter()
    sys2d = _hollow_wall() if name.startswith("hollow") else _slab_roof()
    if name.endswith("_ac"):
        sys2d.solveAC()
        extra = {"cooling_energy": float(sys2d.cooling_energy),
                 "heating_energy": float(sys2d.heating_energy)}
    else:
        sys2d.solve()
        extra = {"energy_transfer": float(sys2d.energy_transfer)}
    runtime = time.perf_counter() - t0
    OUT.mkdir(parents=True, exist_ok=True)
    # full results: Ti, Tso, Tsi, Thueco + the whole Tsa grid, and the (nx, ny)
    # temperature field of the last (converged) day, with the mesh extent in mm
    sys2d.solve_dataframe.to_csv(OUT / f"{name}.csv")
    np.save(OUT / f"{name}_Tfield.npy", sys2d.Tfield)
    mesh = sys2d.section().mesh
    days = sys2d.days
    return name, {"runtime_s": round(runtime, 1),
                  "days": None if days is None else int(days),
                  "X_mm": round(mesh.X * 1000.0, 1),
                  "Y_mm": round(mesh.Y * 1000.0, 1), **extra}


if __name__ == "__main__":
    cases = ["hollow_free", "hollow_ac", "slab_free"]
    summary = {}
    t0 = time.perf_counter()
    with ProcessPoolExecutor(max_workers=len(cases)) as pool:
        for name, info in pool.map(run_case, cases):
            summary[name] = info
            print(name, info, flush=True)
    with open(OUT / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"done in {time.perf_counter() - t0:.0f} s → {OUT}", flush=True)
