"""Figure 2 — Mean-day 1D comparison (paper: EnerHabitat, SoftwareX).

Reproduces the system of Listing 1 (Concreto 15 cm + EPS 5 cm, south-facing,
absortance 0.7, Cuernavaca, May) and compares it against the same wall without
insulation. Run from the repo root:

    uv run python softwareX/EnerHabitat/figures/make_fig2.py
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import enerhabitat as eh

EPW = "docs/data/MEX_MOR_Cuernavaca-Matamoros.Intl.AP.767260_TMYx.2004-2018.epw"
OUT = "softwareX/EnerHabitat/figures"

# Palette (dataviz-validated, light surface)
BLUE, AQUA, ORANGE = "#2a78d6", "#1baf7a", "#eb6834"
GRAY_LINE, INK, BAND = "#52514e", "#0b0b0b", "#e9e8e4"

eh.config.file = "softwareX/EnerHabitat/figures/materials.ini"
loc = eh.Location(EPW)
loc.meanDay(month="5")

systems = {
    "Concrete (12 cm)": [("Concrete", 0.12)],
    "EPS + concrete (1 in + 12 cm)": [("EPS", 0.0254), ("Concrete", 0.12)],
}

results = {}
tsa_df = None
for name, layers in systems.items():
    wall = eh.System(location=loc, tilt=90, azimuth=180)
    wall.absortance = 0.7
    wall.layers = layers
    df = wall.Tsa()
    if tsa_df is None:
        tsa_df = df
    results[name] = wall.solve()

hours = (tsa_df.index - tsa_df.index[0]).total_seconds() / 3600.0

fig, ax = plt.subplots(figsize=(5.5, 3.2), dpi=300)
plt.rcParams.update({"font.size": 9})

# Comfort zone band (neutral)
tn, dtn = tsa_df["Tn"], tsa_df["DeltaTn"]
ax.fill_between(hours, tn - dtn, tn + dtn, color=BAND, zorder=0)

series = [
    ("Ta",  tsa_df["Ta"],  GRAY_LINE, 1.3, (0, (4, 2))),
    ("Tsa", tsa_df["Tsa"], ORANGE,    1.5, "solid"),
    ("Ti concrete", results["Concrete (12 cm)"], BLUE, 1.7, "solid"),
    ("Ti EPS+concrete", results["EPS + concrete (1 in + 12 cm)"], AQUA, 1.7, "solid"),
]
for label, y, color, lw, ls in series:
    ax.plot(hours, y, color=color, lw=lw, ls=ls, zorder=2)

# Direct labels at the right edge, staggered to avoid collisions
MIN_GAP = 1.6  # °C between label anchors
ends = sorted(((s[1].iloc[-1], s[0]) for s in series), reverse=True)
placed = []
for y_end, label in ends:
    y_lab = y_end
    if placed and placed[-1][0] - y_lab < MIN_GAP:
        y_lab = placed[-1][0] - MIN_GAP
    placed.append((y_lab, y_end, label))
for y_lab, y_end, label in placed:
    ax.annotate(label, xy=(24, y_end), xytext=(24.35, y_lab),
                textcoords="data", va="center", ha="left",
                fontsize=8, color=INK, annotation_clip=False)

ax.annotate("comfort zone", xy=(0.6, float((tn + dtn).iloc[0]) - 0.4),
            fontsize=8, color=GRAY_LINE, va="top")

ax.set_xlim(0, 24)
ax.set_xticks(range(0, 25, 3))
ax.set_xlabel("Hour of the mean day (May)")
ax.set_ylabel("Temperature (°C)")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.grid(axis="y", color="#dddcd8", lw=0.6, zorder=1)
ax.tick_params(colors=GRAY_LINE)

fig.subplots_adjust(left=0.09, right=0.78, top=0.97, bottom=0.14)
fig.savefig(f"{OUT}/Figure_2.pdf")
fig.savefig(f"{OUT}/Figure_2.png")
print("saved Figure_2.pdf/.png")
for name, ti in results.items():
    print(f"{name}: Ti min {ti.min():.1f}, max {ti.max():.1f} °C")
print(f"Ta: {tsa_df['Ta'].min():.1f}–{tsa_df['Ta'].max():.1f}, Tsa: {tsa_df['Tsa'].min():.1f}–{tsa_df['Tsa'].max():.1f}")
