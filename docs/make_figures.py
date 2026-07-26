"""
Genera las figuras del README y del sitio de documentación:
  docs/img/hollow_block.png  — bloque hueco de concreto (muro, a escala)
  docs/img/slab.png          — vigueta y bovedilla (techo, 3 cavidades, vigueta en L, a escala)
  docs/img/domain_1d.png     — dominio 1D multicapa con condiciones de frontera

Uso:  .venv/bin/python docs/make_figures.py
Requiere matplotlib (extra `viz`).
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch, Patch

HERE = os.path.dirname(os.path.abspath(__file__))
IMG = os.path.join(HERE, "img")
os.makedirs(IMG, exist_ok=True)

# Paleta común
C_CONCRETE = "#7d7d7d"   # concreto (vigueta + colado/topping + cáscaras del bloque)
C_TOPPING = "#b9b9b9"    # capa de compresión (concreto, tono claro para distinguir)
C_BLOCK = "#9c9c9c"      # bloque de bovedilla
C_AIR = "#ffffff"        # cavidad de aire
EDGE = "#222222"


def _rect(ax, x0, x1, y0, y1, color, hatch=None, lw=0.8):
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor=color,
                           edgecolor=EDGE, lw=lw, hatch=hatch))


def _hdim(ax, x0, x1, y, label):
    """Cota horizontal (flecha doble + etiqueta) debajo del dibujo.

    Con el eje y invertido (exterior arriba), "debajo de la flecha" en pantalla
    es y creciente: la etiqueta se ancla en y+4 para no pisar la línea.
    """
    ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle="<->",
                                 mutation_scale=8, color="#444", lw=0.8))
    ax.text((x0 + x1) / 2, y + 4, label, ha="center", va="top", fontsize=7.5, color="#333")


def _vlabel(ax, x, y0, y1, name, mm):
    ax.text(x, (y0 + y1) / 2, f"{name}\n{mm:g} mm", ha="left", va="center", fontsize=8)


# =================================================================
#  Bloque hueco de concreto (muro)
# =================================================================
def hollow_block():
    web, bw, ct, cav, cb = 20, 160, 20, 80, 20
    X = web + bw + web
    Y = ct + cav + cb
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    # cáscaras (todo el ancho) + almas + cavidad ; y=0 exterior arriba (eje invertido)
    _rect(ax, 0, X, 0, ct, C_CONCRETE)                 # cover_top
    _rect(ax, 0, X, ct + cav, Y, C_CONCRETE)           # cover_bottom
    _rect(ax, 0, web, ct, ct + cav, C_CONCRETE)        # alma izq
    _rect(ax, X - web, X, ct, ct + cav, C_CONCRETE)    # alma der
    _rect(ax, web, X - web, ct, ct + cav, C_AIR, lw=1.0)  # cavidad
    ax.text(X / 2, ct + cav / 2, "AIR\ncavity", ha="center", va="center",
            fontsize=9, color="#555")

    # etiquetas de banda (derecha)
    _vlabel(ax, X + 8, 0, ct, "cover_top", ct)
    _vlabel(ax, X + 8, ct, ct + cav, "cavity", cav)
    _vlabel(ax, X + 8, ct + cav, Y, "cover_bottom", cb)
    # cotas horizontales (abajo)
    _hdim(ax, 0, web, Y + 14, f"web {web}")
    _hdim(ax, web, X - web, Y + 14, f"block_width {bw}")
    _hdim(ax, X - web, X, Y + 14, f"web {web}")
    # EXT arriba; INT con guía horizontal desde la izquierda
    ax.annotate("EXT (Tsa, ho)", xy=(X / 2, 0), xytext=(X / 2, -16),
                ha="center", fontsize=8, color="#b3471f",
                arrowprops=dict(arrowstyle="->", color="#b3471f"))
    ax.annotate("INT (Tint, hi)", xy=(-4, Y), xytext=(-98, Y),
                ha="left", va="center", fontsize=8, color="#1f5fb3",
                arrowprops=dict(arrowstyle="->", color="#1f5fb3"))

    ax.set_xlim(-104, X + 90)
    ax.set_ylim(Y + 46, -34)            # y invertido (exterior arriba)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("HollowBlock — wall cross-section  (x = thickness, y = width)\n"
                 f"thickness = cover_top + cavity + cover_bottom = {Y} mm", fontsize=10)
    ax.legend(handles=[Patch(fc=C_CONCRETE, ec=EDGE, label="block material (concrete)"),
                       Patch(fc=C_AIR, ec=EDGE, label="air cavity")],
              loc="lower center", bbox_to_anchor=(0.42, -0.02), ncol=2, fontsize=8, frameon=False)
    fig.tight_layout()
    out = os.path.join(IMG, "hollow_block.png")
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print("saved", out)


# =================================================================
#  Vigueta y bovedilla (techo, 3 cavidades, vigueta en L)
# =================================================================
def slab():
    web, foot, sh, n, cw = 25, 25, 50, 3, 103
    topping, cap, ct, cav, cb = 100, 50, 30, 40, 30
    X = 2 * (web + foot) + (n + 1) * sh + n * cw
    Y = topping + ct + cav + cb
    # bandas en y (exterior arriba)
    y_top0, y_top1 = 0, topping                      # topping (L2+L3)
    y_ct1 = topping + ct                             # fin cover_top (= inicio cavity)
    y_cav1 = y_ct1 + cav                             # fin cavity (= inicio cover_bottom)
    rib = web + foot

    fig, ax = plt.subplots(figsize=(11, 4.8))
    # topping a todo el ancho
    _rect(ax, 0, X, y_top0, y_top1, C_TOPPING)
    # banda de bovedilla (cover_top + cavity + cover_bottom) a todo el ancho
    _rect(ax, 0, X, y_top1, Y, C_BLOCK)
    # vigueta en L: alma (web) de cap..base en ambos bordes
    _rect(ax, 0, web, cap, Y, C_CONCRETE)
    _rect(ax, X - web, X, cap, Y, C_CONCRETE)
    # pie (foot) solo en cover_bottom
    _rect(ax, web, web + foot, y_cav1, Y, C_CONCRETE)
    _rect(ax, X - rib, X - web, y_cav1, Y, C_CONCRETE)
    # cavidades de aire (centradas en la banda cavity)
    x0 = rib + sh
    for c in range(n):
        xa = x0 + c * (cw + sh)
        _rect(ax, xa, xa + cw, y_ct1, y_cav1, C_AIR, lw=1.0)
        ax.text(xa + cw / 2, (y_ct1 + y_cav1) / 2, "AIR", ha="center", va="center",
                fontsize=8, color="#555")

    # línea de topping_cap (hasta aquí sube la tapa de topping; el alma sube por debajo)
    ax.plot([0, X], [cap, cap], ls="--", color="#c0392b", lw=1.0, zorder=5)
    ax.text(rib + sh, cap - 5, "topping_cap (web rises below this line)",
            ha="left", va="bottom", fontsize=6.8, color="#c0392b")
    # etiquetas de banda (derecha)
    _vlabel(ax, X + 8, y_top0, y_top1, "topping", topping)
    _vlabel(ax, X + 8, y_top1, y_ct1, "cover_top", ct)
    _vlabel(ax, X + 8, y_ct1, y_cav1, "cavity", cav)
    _vlabel(ax, X + 8, y_cav1, Y, "cover_bottom", cb)
    # cotas horizontales (abajo): web foot shoulder cavity ...
    segs = [("web", web), ("foot", foot), ("shoulder", sh)]
    for c in range(n):
        segs += [("cavity_width", cw), ("shoulder", sh)]
    segs += [("foot", foot), ("web", web)]
    xx = 0
    for nm, w in segs:
        _hdim(ax, xx, xx + w, Y + 18, f"{nm}\n{w}")
        xx += w
    # EXT arriba; INT con guía horizontal desde la izquierda
    ax.annotate("EXT (Tsa, ho)", xy=(X / 2, 0), xytext=(X / 2, -20),
                ha="center", fontsize=8, color="#b3471f",
                arrowprops=dict(arrowstyle="->", color="#b3471f"))
    ax.annotate("INT (Tint, hi)", xy=(-4, Y), xytext=(-108, Y),
                ha="left", va="center", fontsize=8, color="#1f5fb3",
                arrowprops=dict(arrowstyle="->", color="#1f5fb3"))

    ax.set_xlim(-116, X + 95)
    ax.set_ylim(Y + 66, -42)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("Slab (joist-and-block roof) — cross-section  ·  3 cavities, L-shaped rib\n"
                 f"width = 2·(web+foot)+(n+1)·shoulder+n·cavity_width = {X} mm   ·   "
                 f"thickness = topping+cover_top+cavity+cover_bottom = {Y} mm", fontsize=9.5)
    ax.legend(handles=[
        Patch(fc=C_TOPPING, ec=EDGE, label="topping (concrete)"),
        Patch(fc=C_CONCRETE, ec=EDGE, label="rib / joist (concrete, L-shaped)"),
        Patch(fc=C_BLOCK, ec=EDGE, label="filler block"),
        Patch(fc=C_AIR, ec=EDGE, label="air cavity (or solid fill if RELLENA)")],
        loc="lower center", bbox_to_anchor=(0.5, -0.16), ncol=2, fontsize=8, frameon=False)
    fig.tight_layout()
    out = os.path.join(IMG, "slab.png")
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print("saved", out)


# =================================================================
#  Dominio 1D multicapa con condiciones de frontera (model-1d.qmd)
# =================================================================
def domain_1d():
    widths = [22, 110, 18]                 # capas de ejemplo (solo esquema)
    labels = ["layer 1", "layer 2", "layer N"]
    colors = ["#b9b9b9", "#8f8f8f", "#d4d4d4"]
    H = 100
    X = sum(widths)
    AIRW = 85                              # aire interior (esquemático, no a escala)

    fig, ax = plt.subplots(figsize=(8.8, 4.2))
    x0 = 0
    centers = []
    for w, col in zip(widths, colors):
        _rect(ax, x0, x0 + w, 0, H, col)
        centers.append(x0 + w / 2)
        x0 += w

    # the wide middle layer holds its label inside; the thin side layers are
    # labelled from above with a leader line (their boxes are too narrow)
    ax.text(centers[1], H / 2, "layer 2\n$k_j,\\ \\rho_j,\\ c_j$",
            ha="center", va="center", fontsize=8.5)
    for c, lab in [(centers[0], "layer 1"), (centers[2], "layer N")]:
        ax.annotate(lab, xy=(c, 8), xytext=(c, -13),
                    ha="center", va="bottom", fontsize=8.5, color="#222",
                    arrowprops=dict(arrowstyle="->", color="#444", lw=0.7))

    # aire interior (nodo agrupado, punteado)
    ax.add_patch(Rectangle((X, 0), AIRW, H, facecolor="#eef4fb",
                           edgecolor="#1f5fb3", ls="--", lw=1.1))
    ax.text(X + AIRW / 2, H / 2, "indoor air\n$\\rho_a,\\ c_a,\\ L_a$\n$T_i(t)$",
            ha="center", va="center", fontsize=8.5, color="#1f5fb3")

    # condiciones de frontera
    ax.annotate("$T_{sa}(t)$,  $h_o$", xy=(0, H * 0.5), xytext=(-46, H * 0.5),
                ha="center", va="center", fontsize=9, color="#b3471f",
                arrowprops=dict(arrowstyle="->", color="#b3471f"))
    ax.annotate("$h_i$", xy=(X, H * 0.62), xytext=(X + 18, H * 0.94),
                ha="center", va="center", fontsize=9, color="#1f5fb3",
                arrowprops=dict(arrowstyle="->", color="#1f5fb3"))
    # continuidad de flujo en las juntas (arriba, lejos del eje x)
    ax.annotate("flux continuity at layer joints", xy=(widths[0], 14),
                xytext=(centers[1], -27), ha="center", va="bottom",
                fontsize=7.5, color="#444",
                arrowprops=dict(arrowstyle="->", color="#444", lw=0.8))

    # ejes x=0, x=L
    ax.annotate("", xy=(X + AIRW + 12, H + 18), xytext=(-12, H + 18),
                arrowprops=dict(arrowstyle="->", color="#333", lw=0.9))
    ax.text(X + AIRW + 14, H + 18, "$x$", ha="left", va="center", fontsize=9)
    for xpos, lab in [(0, "$x=0$\n(outside)"), (X, "$x=L$\n(inside)")]:
        ax.plot([xpos, xpos], [H + 14, H + 22], color="#333", lw=0.9)
        ax.text(xpos, H + 26, lab, ha="center", va="top", fontsize=8)

    ax.set_xlim(-70, X + AIRW + 30)
    ax.set_ylim(H + 46, -42)               # y invertido para dejar cotas abajo
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("1D domain — multilayer wall/roof (outside → inside)", fontsize=10)
    fig.tight_layout()
    out = os.path.join(IMG, "domain_1d.png")
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print("saved", out)


# =================================================================
#  Dominio 2D: una celda repetitiva con sus condiciones de frontera
#  (model-2d.qmd, junto a eq-bc2d)
# =================================================================
def domain_2d():
    web, bw = 20, 160                     # celda del bloque hueco de la doc
    ct, cav, cb = 20, 80, 20
    W = 2 * web + bw                      # 200 mm
    L = ct + cav + cb                     # 120 mm
    AIRH = 55                             # aire interior (esquemático)

    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    _rect(ax, 0, W, 0, L, "#c9c9c9")
    ax.add_patch(Rectangle((web, ct), bw, cav, facecolor="#dcefff",
                           edgecolor="k", ls="--", lw=0.9))
    ax.text(W / 2, ct + cav / 2,
            "cavity:  $h_c(\\Delta T)$ + radiation,  $T_h(t)$",
            ha="center", va="center", fontsize=8.5, color="#374151")

    # frontera exterior (arriba): sol-aire
    ax.annotate("$T_{sa}(t)$,  $h_o$", xy=(W * 0.5, 0), xytext=(W * 0.5, -32),
                ha="center", va="center", fontsize=9.5, color="#b3471f",
                arrowprops=dict(arrowstyle="->", color="#b3471f"))
    # aire interior (abajo, nodo agrupado)
    ax.add_patch(Rectangle((0, L), W, AIRH, facecolor="#eef4fb",
                           edgecolor="#1f5fb3", ls="--", lw=1.1))
    ax.text(W / 2, L + AIRH / 2, "indoor air   $T_i(t)$",
            ha="center", va="center", fontsize=9, color="#1f5fb3")
    ax.annotate("$h_i$", xy=(W * 0.22, L), xytext=(W * 0.22, L + AIRH * 0.62),
                ha="center", va="top", fontsize=9.5, color="#1f5fb3",
                arrowprops=dict(arrowstyle="->", color="#1f5fb3"))

    # costados adiabáticos (planos de simetría especular)
    for x, off in ((0, -7), (W, 7)):
        ax.text(x + off, L * 0.5, "$\\partial T/\\partial y = 0$  (adiabatic)",
                ha="center", va="center", fontsize=7.5, color="#374151",
                rotation=90)

    # caras y ejes (x a través del espesor, como en el 1D; y a lo ancho)
    ax.text(W + 16, 0, "$x = 0$ (outside)", ha="left", va="center",
            fontsize=8, color="#444")
    ax.text(W + 16, L, "$x = L$ (inside)", ha="left", va="center",
            fontsize=8, color="#444")
    ax.annotate("", xy=(-26, L), xytext=(-26, 0),
                arrowprops=dict(arrowstyle="<->", color="#444", lw=0.9))
    ax.text(-32, L / 2, "$L$", ha="right", va="center", fontsize=9)
    yW = L + AIRH + 16
    ax.annotate("", xy=(W, yW), xytext=(0, yW),
                arrowprops=dict(arrowstyle="<->", color="#444", lw=0.9))
    ax.text(W / 2, yW + 6, "$W$", ha="center", va="top", fontsize=9)

    ax.set_xlim(-50, W + 78)
    ax.set_ylim(yW + 22, -44)             # y invertido: exterior arriba
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("2D domain — one repeating cell with its boundary conditions",
                 fontsize=10)
    fig.tight_layout()
    out = os.path.join(IMG, "domain_2d.png")
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print("saved", out)


if __name__ == "__main__":
    hollow_block()
    slab()
    domain_1d()
    domain_2d()
