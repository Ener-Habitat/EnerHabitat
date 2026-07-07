"""Solve the paper's joist-and-block roof (Listing 2 materials) in free-running
and air-conditioned modes; store the series and a summary for Figures/text.

Run from the repo root (takes ~15 min per case):

    uv run python softwareX/EnerHabitat/figures/run_slab_cases.py
"""
import json
import time

import pandas as pd

import enerhabitat as eh

EPW = "docs/data/MEX_MOR_Cuernavaca-Matamoros.Intl.AP.767260_TMYx.2004-2018.epw"
OUT = "softwareX/EnerHabitat/figures/data"


def build_roof():
    eh.config.file = "softwareX/EnerHabitat/figures/materials.ini"
    slab = eh.Slab(rib_material="High-density concrete", block_material="Filler block",
                   topping_material="Aerated concrete", fill_type=eh.Fill.AIR,
                   emissivity=0.9,
                   geometry={"web": 0.025, "foot": 0.025, "shoulder": 0.050,
                             "n_cavities": 3, "cavity_width": 0.103,
                             "topping": 0.100, "topping_cap": 0.050,
                             "cover_top": 0.030, "cavity": 0.040, "cover_bottom": 0.030})
    roof = eh.System2D(eh.Location(EPW))
    roof.tilt = 0
    roof.absortance = 0.3
    roof.layers = [("Waterproofing", 0.003), slab, ("Gypsum plaster", 0.015)]
    roof.location.meanDay(month="5")
    roof.Tsa()
    return roof


import os
os.makedirs(OUT, exist_ok=True)
summary = {}

roof = build_roof()
t0 = time.perf_counter()
ti = roof.solve()
summary["slab_free"] = {"runtime_s": round(time.perf_counter() - t0, 1),
                        "days": int(roof.days),
                        "energy_transfer": float(roof.energy_transfer)}
pd.concat([ti, roof.Tsa()], axis=1).to_csv(f"{OUT}/slab_free.csv")
print("free:", summary["slab_free"], flush=True)

roof = build_roof()
t0 = time.perf_counter()
ti = roof.solveAC()
summary["slab_ac"] = {"runtime_s": round(time.perf_counter() - t0, 1),
                      "days": int(roof.days),
                      "cooling_energy": float(roof.cooling_energy),
                      "heating_energy": float(roof.heating_energy)}
pd.concat([ti, roof.Tsa()], axis=1).to_csv(f"{OUT}/slab_ac.csv")
print("ac:", summary["slab_ac"], flush=True)

with open(f"{OUT}/summary.json", "w") as f:
    json.dump(summary, f, indent=2)
print("done", flush=True)
