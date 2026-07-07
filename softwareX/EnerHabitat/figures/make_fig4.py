"""Figure 4 — Mean-day result of the 2D joist-and-block roof of Listing 2,
free-running (paper: EnerHabitat, SoftwareX).

Reads the series stored by run_slab_cases.py (figures/data/slab_free.csv).
Run from the repo root:

    uv run --with matplotlib python softwareX/EnerHabitat/figures/make_fig4.py
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

OUT = "softwareX/EnerHabitat/figures"

BLUE, ORANGE = "#2a78d6", "#eb6834"
GRAY_LINE, INK, BAND = "#52514e", "#0b0b0b", "#e9e8e4"

df = pd.read_csv(f"{OUT}/data/slab_free.csv", index_col=0, parse_dates=True)
hours = (df.index - df.index[0]).total_seconds() / 3600.0

fig, ax = plt.subplots(figsize=(5.5, 3.2), dpi=300)
plt.rcParams.update({"font.size": 9})

tn, dtn = df["Tn"], df["DeltaTn"]
ax.fill_between(hours, tn - dtn, tn + dtn, color=BAND, zorder=0)

series = [
    ("Ta",  df["Ta"],  GRAY_LINE, 1.3, (0, (4, 2))),
    ("Tsa", df["Tsa"], ORANGE,    1.5, "solid"),
    ("Ti",  df["Ti"],  BLUE,      1.8, "solid"),
]
for label, y, color, lw, ls in series:
    ax.plot(hours, y, color=color, lw=lw, ls=ls, zorder=2)

MIN_GAP = 1.6
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

fig.subplots_adjust(left=0.09, right=0.85, top=0.97, bottom=0.14)
fig.savefig(f"{OUT}/Figure_4.pdf")
fig.savefig(f"{OUT}/Figure_4.png")
print(f"saved; Ti {df['Ti'].min():.1f}-{df['Ti'].max():.1f}, Ta {df['Ta'].min():.1f}-{df['Ta'].max():.1f}, Tsa {df['Tsa'].min():.1f}-{df['Tsa'].max():.1f}")
