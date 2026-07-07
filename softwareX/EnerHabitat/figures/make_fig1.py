"""Figure 1 — EnerHabitat ecosystem architecture (paper: EnerHabitat, SoftwareX).

Pure-matplotlib box diagram: inputs -> enerhabitat package (Location, System,
System2D, Config) -> outputs; the Shiny webapp calls the package; the validation
repository verifies it.

Run from the repo root:

    uv run --with matplotlib python softwareX/EnerHabitat/figures/make_fig1.py
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT = "softwareX/EnerHabitat/figures"
INK, MUTED = "#0b0b0b", "#52514e"
BLUE, AQUA, ORANGE = "#2a78d6", "#1baf7a", "#eb6834"
F_BLUE, F_AQUA, F_ORANGE, F_GRAY = "#eef4fc", "#e9f7f1", "#fdefe8", "#f4f3f0"

fig, ax = plt.subplots(figsize=(6.2, 3.4), dpi=300)
ax.set_xlim(0, 100)
ax.set_ylim(0, 64)
ax.axis("off")


def box(x, y, w, h, fc, ec, lw=1.1, style="round,pad=0.6"):
    b = FancyBboxPatch((x, y), w, h, boxstyle=style, facecolor=fc,
                       edgecolor=ec, lw=lw, mutation_scale=1)
    ax.add_patch(b)
    return b


def text(x, y, s, size=8, color=INK, weight="normal", ha="center", va="center"):
    ax.text(x, y, s, fontsize=size, color=color, ha=ha, va=va, weight=weight)


def arrow(p0, p1, label=None, lx=0, ly=0, style="-|>", color=MUTED):
    a = FancyArrowPatch(p0, p1, arrowstyle=style, color=color, lw=1.1,
                        mutation_scale=11, shrinkA=2, shrinkB=2)
    ax.add_patch(a)
    if label:
        text((p0[0] + p1[0]) / 2 + lx, (p0[1] + p1[1]) / 2 + ly, label,
             size=7.5, color=MUTED)


# --- web application (top) -----------------------------------------------
box(28, 52, 44, 9, F_AQUA, AQUA)
text(50, 58.4, "Web application — Shiny for Python", weight="bold")
text(50, 55.0, "enerhabitat.unam.mx · GUI for the 1D model · no programming required", size=7.5, color=MUTED)

# --- package (center) ------------------------------------------------------
box(22, 20, 56, 25, F_BLUE, BLUE)
text(50, 41.8, "enerhabitat Python package (PyPI)", weight="bold")

box(25.5, 23, 14.5, 14, "white", MUTED, lw=0.8)
text(32.7, 33.3, "Location", weight="bold", size=7.5)
text(32.7, 28.7, "EPW parsing\nmean day\nsolar (pvlib)", size=7)

box(42.5, 23, 18.5, 14, "white", MUTED, lw=0.8)
text(51.7, 33.3, "System · System2D", weight="bold", size=7.2)
text(51.7, 28.2, "Tsa()\nsolve() · solveAC()\nHollowBlock, Slab", size=7)

box(63.5, 23, 11.5, 14, "white", MUTED, lw=0.8)
text(69.2, 33.3, "Config", weight="bold", size=7.5)
text(69.2, 28.7, "materials.ini\nmesh, dt\nhi, ho", size=7)

arrow((40.2, 30), (42.3, 30))
arrow((63.3, 30), (61.2, 30))

# --- inputs (left) ---------------------------------------------------------
box(2, 24, 14, 12, F_GRAY, MUTED)
text(9, 32.6, "Inputs", weight="bold", size=7.5)
text(9, 28.4, "EPW file\nlayers +\n2D element", size=7)
arrow((16.4, 30), (21.6, 30))

# --- outputs (right) -------------------------------------------------------
box(84, 24, 14, 12, F_GRAY, MUTED)
text(91, 32.6, "Outputs", weight="bold", size=7.5)
text(91, 28.2, "Ta, Tsa, Ti\nenergies\nFD, time lag", size=7)
arrow((78.4, 30), (83.6, 30))

# --- validation (bottom) ---------------------------------------------------
box(28, 4, 44, 9, F_ORANGE, ORANGE)
text(50, 10.4, "Validation repository — eh_validation", weight="bold")
text(50, 7.0, "notebooks vs. EnergyPlus (1D) and guarded hot-plate data (2D)", size=7.5, color=MUTED)

arrow((50, 51.6), (50, 45.6), label="calls", lx=6, ly=0)
arrow((50, 13.4), (50, 19.4), label="verifies", lx=7.5, ly=0)

fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
fig.savefig(f"{OUT}/Figure_1.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/Figure_1.png", dpi=300, bbox_inches="tight")
print("Figure_1 saved")
