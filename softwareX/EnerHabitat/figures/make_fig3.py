"""Figure 3 — Joist-and-block roof section drawn to scale by the package's
inspector, preview(field="materials") (paper: EnerHabitat, SoftwareX).

Materials as in Listing 2: high-density concrete joist, aerated-concrete
compression topping, bovedilla filler blocks, air cavities.

Run from the repo root:

    uv run --with matplotlib python softwareX/EnerHabitat/figures/make_fig3.py
"""
import matplotlib
matplotlib.use("Agg")

import enerhabitat as eh

EPW = "docs/data/MEX_MOR_Cuernavaca-Matamoros.Intl.AP.767260_TMYx.2004-2018.epw"
OUT = "softwareX/EnerHabitat/figures"

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

fig, axes = roof.preview(field="materials")
ax = axes[0] if hasattr(axes, "__len__") else axes
fig.suptitle("")
ax.set_title("")
ax.set_xlabel("width (mm)", fontsize=9)
ax.set_ylabel("thickness (mm), out → in", fontsize=9)
ax.tick_params(labelsize=8)
leg = ax.get_legend() or (fig.legends[0] if fig.legends else None)
if leg is not None:
    handles = getattr(leg, "legend_handles", None) or leg.legendHandles
    labels = [t.get_text() for t in leg.get_texts()]
    leg.remove()
    ax.legend(handles, labels, loc="center left", bbox_to_anchor=(1.02, 0.5),
              ncol=1, fontsize=8, frameon=False)
fig.set_size_inches(6.4, 2.7)
fig.savefig(f"{OUT}/Figure_3.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/Figure_3.png", dpi=300, bbox_inches="tight")
print("Figure_3 saved")
