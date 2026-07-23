"""
=================================================================
 eh2d — 2D geometry of a joist-and-block slab (Phase 1)
=================================================================

Faithful port of the geometry/topology of the C solver `legacy_eh/2dTfree/`
(see ``PLAN-2D.md``). This phase builds, for a **solid filler block**
(``Fill.SOLID``, C ``tipo 2``):

    NT[i][j]        node-type mesh (1-8 boundaries/corners, 13 interior)
    k[i][j]         conductivity per node
    rhoc[i][j]      heat capacity per node
    X, Y, dx, dy    cell size and discretisation
    i1, j1, i2, j2  bounds of the fill block (filler block)

Mesh convention (identical to the C):
    i = 0 .. nx-1   width     X   (i=0 left, i=nx-1 right; adiabatic sides)
    j = 0 .. ny-1   thickness Y   (j=0 outside with Tsa/ho, j=ny-1 inside with Tint/hi)

Arrays are indexed ``A[i, j]`` with shape ``(nx, ny)``, same as the C dump
(``dump_NT.dat`` etc. iterate ``for i: for j: print i j A[i][j]``).

The physics (coefficient assembly, solver, time loop) arrives in later phases;
here only the geometry is reproduced, to validate it node by node.
"""

from dataclasses import dataclass, field
from enum import Enum

import numpy as np

from .config import config, config2d
from .ehtools2d import (solve_day_2d, solve_day_hueca_prod, solve_day_slab_prod,
                        solve_day_2d_ac, solve_day_hueca_ac, solve_day_slab_ac,
                        _view_factors)


class Fill(Enum):
    """State of the filler block (what the C numeric ``tipo`` controls)."""
    SOLID = "solid"            # tipo 2: solid block (fill kr, rhocr)
    AIR = "air"                  # tipo 1: air cavity (radiation + Nusselt)  [Phase 6]
    SOLID_SYMMETRIC = "solid_sym"  # tipo 4: symmetric half cell, solid    [later]


# Mapping to the C `tipo` integers, to read .inp files and legacy golden masters.
TIPO_C = {
    Fill.AIR: 1,
    Fill.SOLID: 2,
    Fill.SOLID_SYMMETRIC: 4,
}


@dataclass
class Mesh2D:
    """Dimensions and bounds of the fill block (result of the mesh computation)."""
    nx: int
    ny: int
    X: float
    Y: float
    dx: float
    dy: float
    i1: int
    j1: int
    i2: int
    j2: int


def compute_mesh(nx, ny, L, layer, a, e):
    """
    Literally reproduces the mesh computation of the C ``main`` (solid filler
    block / single cavity, ``a14 == 0``).

    Args:
        nx, ny (int): number of nodes in width and thickness.
        L (sequence[float]): layer thicknesses ``[L1..L7]`` (outside to inside).
        layer (int): layer number (from outside, 1-based) where the filler block sits.
        a (dict): horizontal geometry with keys ``a11,a12,a13,a14,a21,a22,a23``.
        e (dict): filler-block thicknesses with keys ``e21,e22,e23``.

    Returns:
        Mesh2D
    """
    # YY[0]=0, YY[1..7]=L1..L7; layer `layer` is replaced by the filler-block height.
    YY = [0.0] + [float(L[i]) if i < len(L) else 0.0 for i in range(7)]
    YY[layer] = e["e21"] + e["e22"] + e["e23"]

    # Cell width. a14==0 -> "single cavity" (the solid/hollow filler-block case).
    if a.get("a14", 0.0) == 0.0:
        X = a["a21"] + a["a11"] + a["a12"] / 2.0
    else:
        X = (a["a11"] + a["a21"] + a["a12"] + a["a22"]
             + a["a13"] + a["a23"] + a["a14"])

    dx = X / nx
    Y = sum(YY)
    dy = Y / ny

    # y offset of the filler-block layer (same integer truncation as the C).
    y1 = sum(YY[:layer])
    y1 = y1 / dy + 0.5

    i1 = int(a["a11"] / dx + 0.5)
    j1 = int(e["e21"] / dy + int(y1))
    i2 = int((a["a11"] + a["a21"]) / dx + 0.5)
    j2 = int((e["e21"] + e["e22"]) / dy + int(y1))

    return Mesh2D(nx=nx, ny=ny, X=X, Y=Y, dx=dx, dy=dy,
                  i1=i1, j1=j1, i2=i2, j2=j2)


def draw_rellena(nx, ny, i1, j1, i2, j2):
    """
    Port of ``draw_viguetabovedilla2rellena``: node-type mesh ``NT`` for a solid
    filler block. Marks the fill region with ``14`` (temporary); later
    :func:`set_krhoc_rellena` turns it into ``13`` when it assigns ``kr,rhocr``.

    Returns:
        np.ndarray (int) of shape (nx, ny).
    """
    NT = np.zeros((nx, ny), dtype=np.int64)

    # Corners.
    NT[0, 0] = 1            # top-left      (outside, adiabatic left)
    NT[nx - 1, 0] = 2       # top-right     (outside, adiabatic right)
    NT[0, ny - 1] = 3       # bottom-left   (inside, adiabatic left)
    NT[nx - 1, ny - 1] = 4  # bottom-right  (inside, adiabatic right)

    # Adiabatic sides.
    NT[0, 1:ny - 1] = 6
    NT[nx - 1, 1:ny - 1] = 7

    # Convective boundaries: outside (j=0) and inside (j=ny-1).
    NT[1:nx - 1, 0] = 5
    NT[1:nx - 1, ny - 1] = 8

    # Interior nodes.
    NT[1:nx - 1, 1:ny - 1] = 13

    # (In the C the cavity edges are rewritten to 13; here they are already 13, no-op.)

    # Fill region: temporary mark 14.
    NT[i1:i2, j1:j2] = 14

    return NT


def set_krhoc_rellena(nx, ny, dx, dy, L, k, rhoc, kr, rhocr, NT):
    """
    Port of ``set_krhocrelleno``: assigns ``k`` and ``rhoc`` per layer (j thresholds
    from cumulative thicknesses ``L1, L1+L2, ...``) and then overwrites the fill
    region (``NT==14``) with ``kr, rhocr``, marking it as interior node ``13``.

    Mutates ``NT`` in place (14 -> 13) and returns ``(k_field, rhoc_field)``.
    """
    k = [float(k[i]) if i < len(k) else 0.0 for i in range(7)]
    rhoc = [float(rhoc[i]) if i < len(rhoc) else 0.0 for i in range(7)]
    L = [float(L[i]) if i < len(L) else 0.0 for i in range(7)]

    kf = np.zeros((nx, ny), dtype=np.float64)
    rf = np.zeros((nx, ny), dtype=np.float64)

    # Cumulative thresholds (identical to the C `for (; j < (L1+...+Ln)/dy; ++j)`).
    thr = np.cumsum(L) / dy  # thr[n] = (L1+..+L_{n+1})/dy

    for i in range(nx):
        j = 0
        for n in range(7):
            while j < thr[n] and j < ny:
                kf[i, j] = k[n]
                rf[i, j] = rhoc[n]
                j += 1

    # Fill (filler block): kr, rhocr and NT 14 -> 13.
    fill = NT == 14
    kf[fill] = kr
    rf[fill] = rhocr
    NT[fill] = 13

    return kf, rf


def draw_hueca(nx, ny, i1, j1, i2, j2):
    """
    Port of ``draw_viguetabovedilla2hueca``: ``NT`` mesh for a filler block with an
    **air cavity** (``tipo 1``). Same frame as :func:`draw_rellena`, but the cavity
    region are air nodes ``0`` surrounded by walls:
    ``9`` (top, j=j1-1), ``10`` (bottom, j=j2), ``11`` (left, i=i1-1),
    ``12`` (right, i=i2).
    """
    NT = np.zeros((nx, ny), dtype=np.int64)
    NT[0, 0] = 1
    NT[nx - 1, 0] = 2
    NT[0, ny - 1] = 3
    NT[nx - 1, ny - 1] = 4
    NT[0, 1:ny - 1] = 6
    NT[nx - 1, 1:ny - 1] = 7
    NT[1:nx - 1, 0] = 5
    NT[1:nx - 1, ny - 1] = 8
    NT[1:nx - 1, 1:ny - 1] = 13
    NT[i1:i2, j1 - 1] = 9       # top wall of the cavity
    NT[i1:i2, j2] = 10          # bottom wall
    NT[i1 - 1, j1:j2] = 11      # left wall
    NT[i2, j1:j2] = 12          # right wall
    NT[i1:i2, j1:j2] = 0        # cavity air
    return NT


def set_krhoc_hueca(nx, ny, dx, dy, L, k, rhoc):
    """
    Port of ``set_krhoc`` (tipo 1): fills ``k,rhoc`` per layer in y (cumulative
    thresholds), **without** a fill override. Air/wall nodes keep the material of
    their layer (irrelevant in the air: the NT 0 case fixes T=Thueco).
    """
    k = [float(k[i]) if i < len(k) else 0.0 for i in range(7)]
    rhoc = [float(rhoc[i]) if i < len(rhoc) else 0.0 for i in range(7)]
    L = [float(L[i]) if i < len(L) else 0.0 for i in range(7)]
    kf = np.zeros((nx, ny))
    rf = np.zeros((nx, ny))
    thr = np.cumsum(L) / dy
    for i in range(nx):
        j = 0
        for n in range(7):
            while j < thr[n] and j < ny:
                kf[i, j] = k[n]
                rf[i, j] = rhoc[n]
                j += 1
    return kf, rf


# =================================================================
#  Roof geometry: joist and filler block, N cavities, 3 solids (Phase 8b)
# =================================================================
#  Three solid materials (topping, L-shaped joist, filler block) + N equal air
#  cavities. The joist is an **L**: web `web` (width d1, rises through the whole
#  element) + foot `foot` (adds d2, only in the bottom band `cover_bottom`),
#  forming the ledge the filler block rests on. The topping (compression layer)
#  spans the top band at full width. L1/finishes are NOT part of the element.


def compute_mesh_slab(nx, ny, L, layer, web, foot, shoulder, n_cav, cavity_width,
                      topping, cover_top, cavity, cover_bottom, topping_cap=0.0):
    """
    Mesh for the N-cavity roof slab. Returns ``(mesh, info)`` where ``info`` holds
    the element's internal integer bounds and the x-bounds of each cavity. Same
    integer truncation as :func:`compute_mesh`.
    """
    e_thick = topping + cover_top + cavity + cover_bottom
    YY = [0.0] + [float(L[i]) if i < len(L) else 0.0 for i in range(7)]
    YY[layer] = e_thick

    X = 2.0 * (web + foot) + (n_cav + 1) * shoulder + n_cav * cavity_width
    dx = X / nx
    Y = sum(YY)
    dy = Y / ny

    y1 = sum(YY[:layer]) / dy + 0.5
    base = int(y1)

    e21 = topping + cover_top          # "cap" above the cavity (topping + filler block)
    cj1 = int(e21 / dy + base)        # top-wall row = cj1-1; air from cj1 on
    cj2 = int((e21 + cavity) / dy + base)
    jet = base                        # top of the element
    jcap = int(topping_cap / dy + base)  # base of the topping cap (L2): the web stops here
    jcol = int(topping / dy + base)    # topping/filler-block boundary
    jeb = int((e21 + cavity + cover_bottom) / dy + base)  # base of the element

    niw = int(web / dx + 0.5)
    nif = int((web + foot) / dx + 0.5)

    cav_i1 = np.empty(n_cav, dtype=np.int64)
    cav_i2 = np.empty(n_cav, dtype=np.int64)
    x0 = (web + foot) + shoulder
    for c in range(n_cav):
        xs = x0 + c * (cavity_width + shoulder)
        cav_i1[c] = int(xs / dx + 0.5)
        cav_i2[c] = int((xs + cavity_width) / dx + 0.5)

    mesh = Mesh2D(nx=nx, ny=ny, X=X, Y=Y, dx=dx, dy=dy,
                  i1=int(cav_i1[0]), j1=cj1, i2=int(cav_i2[-1]), j2=cj2)
    info = {"jet": jet, "jcap": jcap, "jcol": jcol, "cj1": cj1, "cj2": cj2,
            "jeb": jeb, "niw": niw, "nif": nif, "cav_i1": cav_i1, "cav_i2": cav_i2}
    return mesh, info


def draw_slab_multi(nx, ny, cav_i1, cav_i2, cj1, cj2, hollow):
    """
    ``NT`` mesh of the N-cavity roof slab. ``hollow=True`` → each cavity is air
    (0) surrounded by walls 9/10/11/12; ``hollow=False`` (SOLID) → the cavity
    stays an interior node (13, fill material). Returns ``(NT, cav_of)`` with
    ``cav_of[i,j]`` = cavity index of the air/wall nodes (−1 elsewhere).
    """
    NT = np.zeros((nx, ny), dtype=np.int64)
    NT[0, 0] = 1; NT[nx - 1, 0] = 2; NT[0, ny - 1] = 3; NT[nx - 1, ny - 1] = 4
    NT[0, 1:ny - 1] = 6; NT[nx - 1, 1:ny - 1] = 7
    NT[1:nx - 1, 0] = 5; NT[1:nx - 1, ny - 1] = 8
    NT[1:nx - 1, 1:ny - 1] = 13
    cav_of = np.full((nx, ny), -1, dtype=np.int64)
    if hollow:
        for c in range(len(cav_i1)):
            i1, i2 = int(cav_i1[c]), int(cav_i2[c])
            NT[i1:i2, cj1 - 1] = 9
            NT[i1:i2, cj2] = 10
            NT[i1 - 1, cj1:cj2] = 11
            NT[i2, cj1:cj2] = 12
            NT[i1:i2, cj1:cj2] = 0
            cav_of[i1:i2, cj1 - 1] = c
            cav_of[i1:i2, cj2] = c
            cav_of[i1 - 1, cj1:cj2] = c
            cav_of[i2, cj1:cj2] = c
            cav_of[i1:i2, cj1:cj2] = c
    # SOLID: the cavity region stays 13 (the material is filled in set_krhoc).
    return NT, cav_of


def set_krhoc_slab(nx, ny, dx, dy, L, k, rhoc, layer, info,
                   k_topping, rc_topping, k_rib, rc_rib, k_block, rc_block,
                   k_fill, rc_fill, cav_i1, cav_i2, hollow):
    """
    Assigns ``k``/``rhoc`` per node for the roof slab with **three** solids:
    homogeneous layers (by y thresholds), and inside the element — topping (top
    band), filler block (the rest) and the **L**-shaped joist (web + foot). If
    ``hollow`` is false, fills the cavity with ``k_fill/rc_fill``.
    """
    kk = [float(k[i]) if i < len(k) else 0.0 for i in range(7)]
    rr = [float(rhoc[i]) if i < len(rhoc) else 0.0 for i in range(7)]
    LL = [float(L[i]) if i < len(L) else 0.0 for i in range(7)]
    kf = np.zeros((nx, ny)); rf = np.zeros((nx, ny))
    thr = np.cumsum(LL) / dy
    for i in range(nx):
        j = 0
        for n in range(7):
            while j < thr[n] and j < ny:
                kf[i, j] = kk[n]; rf[i, j] = rr[n]; j += 1

    jet, jcap, jcol = info["jet"], info["jcap"], info["jcol"]
    cj1, cj2, jeb = info["cj1"], info["cj2"], info["jeb"]
    niw, nif = info["niw"], info["nif"]
    for i in range(nx):
        is_alma = (i < niw) or (i >= nx - niw)
        is_foot = (i < nif) or (i >= nx - nif)
        for j in range(jet, jeb):
            if j < jcol:
                km, rm = k_topping, rc_topping      # topping band
            else:
                km, rm = k_block, rc_block         # filler block
            if is_alma and j >= jcap:
                km, rm = k_rib, rc_rib             # L web: from jcap (below the L2 cap) to the base
            elif is_foot and j >= cj2:
                km, rm = k_rib, rc_rib             # L foot (only cover_bottom)
            kf[i, j] = km; rf[i, j] = rm

    if not hollow:
        for c in range(len(cav_i1)):
            i1, i2 = int(cav_i1[c]), int(cav_i2[c])
            kf[i1:i2, cj1:cj2] = k_fill
            rf[i1:i2, cj1:cj2] = rc_fill
    return kf, rf


@dataclass
class SlabSection:
    """Roof joist-and-block section (N cavities, 3 solids). Exposes
    ``NT/kfield/rhocfield/mesh`` like :class:`Section2D` (for the inspector) plus
    the cavity arrays that the ``solve_day_slab_prod`` engine needs."""
    nx: int
    ny: int
    L: list
    k: list
    rhoc: list
    layer: int
    geom: dict            # web,foot,shoulder,n_cav,cavity_width,topping,cover_top,cavity,cover_bottom
    k_topping: float
    rc_topping: float
    k_rib: float
    rc_rib: float
    k_block: float
    rc_block: float
    k_fill: float = 0.0
    rc_fill: float = 0.0
    emissivity: float = 0.9
    beta: float = 0.0
    hollow: bool = True

    mesh: Mesh2D = field(init=False, default=None)
    NT: np.ndarray = field(init=False, default=None)
    kfield: np.ndarray = field(init=False, default=None)
    rhocfield: np.ndarray = field(init=False, default=None)
    cav_of: np.ndarray = field(init=False, default=None)
    cav_i1: np.ndarray = field(init=False, default=None)
    cav_i2: np.ndarray = field(init=False, default=None)
    info: dict = field(init=False, default=None)

    def build(self):
        g = self.geom
        mesh, info = compute_mesh_slab(
            self.nx, self.ny, self.L, self.layer,
            g["web"], g["foot"], g["shoulder"], g["n_cav"], g["cavity_width"],
            g["topping"], g["cover_top"], g["cavity"], g["cover_bottom"],
            g.get("topping_cap", 0.0))
        NT, cav_of = draw_slab_multi(self.nx, self.ny, info["cav_i1"],
                                     info["cav_i2"], info["cj1"], info["cj2"],
                                     self.hollow)
        kf, rf = set_krhoc_slab(
            self.nx, self.ny, mesh.dx, mesh.dy, self.L, self.k, self.rhoc,
            self.layer, info, self.k_topping, self.rc_topping, self.k_rib,
            self.rc_rib, self.k_block, self.rc_block, self.k_fill, self.rc_fill,
            info["cav_i1"], info["cav_i2"], self.hollow)
        self.mesh, self.NT, self.kfield, self.rhocfield = mesh, NT, kf, rf
        self.cav_of, self.cav_i1, self.cav_i2 = cav_of, info["cav_i1"], info["cav_i2"]
        self.info = info
        return self


@dataclass
class Section2D:
    """
    Description of a joist-and-block slab section.

    Layers L1..L7 (outside to inside) with ``k``/``rhoc`` per layer; filler-block
    fill with ``kr``/``rhocr``; horizontal geometry ``a*`` and filler-block
    thicknesses ``e2*``. For now only ``Fill.SOLID`` (``tipo 2``).
    """
    nx: int
    ny: int
    L: list           # [L1..L7]
    k: list           # [k1..k7]
    rhoc: list        # [rhoc1..rhoc7]
    kr: float
    rhocr: float
    a: dict           # a11,a12,a13,a14,a21,a22,a23
    e: dict           # e21,e22,e23
    layer: int = 1
    fill_type: Fill = Fill.SOLID

    mesh: Mesh2D = field(init=False, default=None)
    NT: np.ndarray = field(init=False, default=None)
    kfield: np.ndarray = field(init=False, default=None)
    rhocfield: np.ndarray = field(init=False, default=None)

    def build(self):
        """Builds ``mesh``, ``NT``, ``kfield``, ``rhocfield`` and returns self."""
        m = compute_mesh(self.nx, self.ny, self.L, self.layer, self.a, self.e)
        if self.fill_type is Fill.SOLID:
            NT = draw_rellena(m.nx, m.ny, m.i1, m.j1, m.i2, m.j2)
            kf, rf = set_krhoc_rellena(m.nx, m.ny, m.dx, m.dy,
                                       self.L, self.k, self.rhoc,
                                       self.kr, self.rhocr, NT)
        elif self.fill_type is Fill.AIR:
            NT = draw_hueca(m.nx, m.ny, m.i1, m.j1, m.i2, m.j2)
            kf, rf = set_krhoc_hueca(m.nx, m.ny, m.dx, m.dy,
                                     self.L, self.k, self.rhoc)
        else:
            raise NotImplementedError(
                f"{self.fill_type} (symmetric tipo 4) not ported yet.")
        self.mesh, self.NT, self.kfield, self.rhocfield = m, NT, kf, rf
        return self


# --- inspection -----------------------------------------------------------------

_NODE_GLYPH = {0: ".", 1: "1", 2: "2", 3: "3", 4: "4", 5: "5", 6: "6",
               7: "7", 8: "8", 9: "9", 10: "a", 11: "b", 12: "c",
               13: "·", 14: "#"}


def print_node_scheme(NT, max_cols=120, max_rows=60):
    """
    Prints the node-type map (j downwards = outside→inside, i across), like the
    diagram in ``PLAN-2D.md``. Subsamples if the mesh is large.
    """
    nx, ny = NT.shape
    istep = max(1, int(np.ceil(nx / max_cols)))
    jstep = max(1, int(np.ceil(ny / max_rows)))
    if istep > 1 or jstep > 1:
        print(f"# subsampling: every {istep} in i, {jstep} in j "
              f"(real mesh {nx}x{ny})")
    for j in range(0, ny, jstep):
        row = "".join(_NODE_GLYPH.get(int(NT[i, j]), "?")
                      for i in range(0, nx, istep))
        print(row)


def print_material_scheme(section, max_cols=120, max_rows=60):
    """
    Prints the material map (by ``k`` value) with a legend, so the **filler block**
    (fill) is distinguished from the **joist/layers** (hence the name "joist and
    filler block"). j downwards = outside→inside. Returns the legend
    ``{glyph: (k, rhoc)}``.
    """
    k = section.kfield
    rhoc = section.rhocfield
    nx, ny = k.shape

    # Glyph per material (rounded k value); most conductive = dense, fill = light.
    uniq = sorted({round(float(v), 12) for v in np.unique(k)})
    glyphs = "█▓▒░·:. "
    legend = {}
    val2glyph = {}
    for idx, kv in enumerate(sorted(uniq, reverse=True)):  # large k -> "dense" glyph
        gph = glyphs[idx] if idx < len(glyphs) else "?"
        val2glyph[kv] = gph
        # rhoc representative of that material.
        mask = np.isclose(k, kv)
        rc = float(rhoc[mask].flat[0]) if mask.any() else 0.0
        legend[gph] = (kv, rc)

    istep = max(1, int(np.ceil(nx / max_cols)))
    jstep = max(1, int(np.ceil(ny / max_rows)))
    if istep > 1 or jstep > 1:
        print(f"# subsampling: every {istep} in i, {jstep} in j "
              f"(real mesh {nx}x{ny})")
    for j in range(0, ny, jstep):
        row = "".join(val2glyph[round(float(k[i, j]), 12)]
                      for i in range(0, nx, istep))
        print(row)
    print("legend:")
    for gph, (kv, rc) in legend.items():
        print(f"  '{gph}'  k={kv:g} W/mK   rhoc={rc:g} J/m³K")
    return legend


def plot_node_scheme(section, ax=None, field="NT"):
    """
    Draws (matplotlib ``imshow``) ``NT``, ``k`` or ``rhoc``. ``field`` ∈
    {"NT","k","rhoc"}. Imports matplotlib lazily (optional dependency).
    """
    import matplotlib.pyplot as plt

    data = {"NT": section.NT, "k": section.kfield, "rhoc": section.rhocfield}[field]
    if ax is None:
        _, ax = plt.subplots()
    # Transpose to [j, i] so j (thickness) goes on the vertical axis (outside on top).
    im = ax.imshow(data.T, origin="upper", aspect="auto", interpolation="nearest")
    ax.set_xlabel("i (width)")
    ax.set_ylabel("j (thickness: outside→inside)")
    ax.set_title(field)
    plt.colorbar(im, ax=ax)
    return ax


# =================================================================
#  2D production API — reuses EPW+pvlib via System (1D)
# =================================================================
#  2D elements that can go as ONE layer inside `System2D.layers`:
#  - HollowBlock (Phase 8a): concrete hollow block, WALLS only (tilt=90).
#  - Slab        (Phase 8b): joist and filler block, ROOFS only (tilt=0).


def _geom_pick(g, friendly, raw, default=None):
    if friendly in g:
        return g[friendly]
    if raw in g:
        return g[raw]
    if default is not None:
        return default
    raise KeyError(f"geometry: missing '{friendly}' (or the raw '{raw}')")


class HollowBlock:
    """
    Concrete hollow block for **walls** (`tilt=90`). A shell of one material with
    one cell that is either an **air cavity** (``Fill.AIR``: wall Nusselt
    convection + radiation between the cavity walls) or **filled** with a solid
    material (``Fill.SOLID``: e.g. an insulating core), ``fill_material``.

    Args:
        material (str): shell/block material (e.g. "Concreto"), from ``config``.
        fill_type (Fill): ``AIR`` (air cavity) or ``SOLID`` (solid fill).
        fill_material (str|None): cavity fill material; required if ``SOLID``.
        emissivity (float): emissivity of the cavity walls (radiation, ``AIR``).
        geometry (dict): cell measures; friendly keys
            ``web``(=a11), ``block_width``(=a21), ``cover_top``(=e21),
            ``cavity``(=e22), ``cover_bottom``(=e23); the raw ``a11..e23`` are
            accepted too. By symmetry ``a12 = 2·web`` unless ``a12`` is given.
    """

    required_tilt = 90

    def __init__(self, material, fill_type=Fill.AIR, fill_material=None,
                 emissivity=0.9, geometry=None):
        self.material = material
        self.fill_type = fill_type
        self.fill_material = fill_material
        self.emissivity = emissivity
        self.geometry = dict(geometry or {})
        if fill_type is Fill.SOLID and not fill_material:
            raise ValueError("HollowBlock SOLID requires fill_material.")

    @property
    def material_main(self):
        return self.material

    def _ae(self):
        g = self.geometry
        web = _geom_pick(g, "web", "a11")
        a21 = _geom_pick(g, "block_width", "a21")
        e21 = _geom_pick(g, "cover_top", "e21")
        e22 = _geom_pick(g, "cavity", "e22")
        e23 = _geom_pick(g, "cover_bottom", "e23")
        a12 = g.get("a12", 2.0 * web)
        a = {"a11": web, "a12": a12, "a13": 0.0, "a14": 0.0,
             "a21": a21, "a22": 0.0, "a23": 0.0}
        e = {"e21": e21, "e22": e22, "e23": e23}
        return a, e

    @property
    def thickness(self):
        _, e = self._ae()
        return e["e21"] + e["e22"] + e["e23"]

    def signature(self):
        return ("HollowBlock", self.material, self.fill_type.value,
                self.fill_material, self.emissivity,
                tuple(sorted(self.geometry.items())))


class Slab:
    """
    Joist and filler block for **roofs** (`tilt=0`). Three solid materials —
    ``rib_material`` (joist, **L**-shaped), ``block_material`` (filler block, which
    surrounds the cavities) and ``topping_material`` (compression layer) — plus N
    equal cavities of air (``Fill.AIR``) or fill (``Fill.SOLID``). The
    cavity is horizontal → roof Nusselt (Rayleigh). L1/finishes are NOT part of the
    element: they go as homogeneous layers in ``System2D.layers``.

    Args:
        rib_material (str): joist material (web + foot).
        fill_type (Fill): ``AIR`` (cavity) or ``SOLID`` (solid).
        block_material (str): filler-block material (surrounds the cavity);
            defaults to the same as ``rib_material``.
        topping_material (str): topping material; defaults to ``rib_material``.
        fill_material (str|None): cavity fill material if ``SOLID``.
        emissivity (float): emissivity of the cavity walls (radiation) if AIR.
        geometry (dict): friendly keys ``web``(=d1), ``foot``(=d2),
            ``shoulder``(=d3), ``n_cavities``, ``cavity_width``(=d4), ``topping``
            (=L2+L3, total topping thickness), ``topping_cap`` (=L2, the full-width
            topping cap above the web; the L web rises only up to its base, leaving
            height L3+cover_top+cavity+cover_bottom; default 0 → web through the
            full height), ``cover_top`` (filler block above the cavity, =L4),
            ``cavity`` (=L5, cavity height), ``cover_bottom`` (=L6). Raw aliases
            ``d1..d4``.
    """

    required_tilt = 0

    def __init__(self, rib_material, fill_type=Fill.AIR, block_material=None,
                 topping_material=None, fill_material=None, emissivity=0.9,
                 geometry=None):
        self.rib_material = rib_material
        self.block_material = block_material or rib_material
        self.topping_material = topping_material or rib_material
        self.fill_type = fill_type
        self.fill_material = fill_material
        self.emissivity = emissivity
        self.geometry = dict(geometry or {})
        if fill_type is Fill.SOLID and not fill_material:
            raise ValueError("Slab SOLID requires fill_material.")

    @property
    def material_main(self):
        return self.topping_material

    def _geom(self):
        g = self.geometry
        web = _geom_pick(g, "web", "d1")
        foot = _geom_pick(g, "foot", "d2")
        shoulder = _geom_pick(g, "shoulder", "d3")
        cavity_width = _geom_pick(g, "cavity_width", "d4")
        n_cav = int(g.get("n_cavities", g.get("n_cav", 1)))
        cover_top = g.get("cover_top", 0.0)
        return {"web": web, "foot": foot, "shoulder": shoulder, "n_cav": n_cav,
                "cavity_width": cavity_width, "topping": g["topping"],
                "topping_cap": g.get("topping_cap", 0.0),
                "cover_top": cover_top, "cavity": _geom_pick(g, "cavity", "e22"),
                "cover_bottom": _geom_pick(g, "cover_bottom", "e23")}

    @property
    def thickness(self):
        g = self._geom()
        return g["topping"] + g["cover_top"] + g["cavity"] + g["cover_bottom"]

    def signature(self):
        return ("Slab", self.rib_material, self.block_material, self.topping_material,
                self.fill_type.value, self.fill_material, self.emissivity,
                tuple(sorted(self.geometry.items())))


# Types recognised as a "2D element" inside layers.
_ELEMENT_TYPES = (HollowBlock, Slab)


class System2D:
    """
    2D production constructive system (same methodology as the 1D ``System``).

    ``layers`` is the outside→inside stack of homogeneous layers ``(material, L)``
    plus **one** 2D element (``HollowBlock``/``Slab``); its position in the list is
    its order in the stack. ``Tsa(t)`` is reused from an internal 1D ``System``
    (EPW+pvlib, at the ``config.dt`` step); the geometry is built with ``Section2D``
    and solved with the JIT engine matching the element type.

    Usage (identical to the 1D one except for the 2D element)::

        block = eh.HollowBlock("Concreto", emissivity=0.9, geometry={...})
        wall = eh.System2D(location=loc); wall.tilt = 90
        wall.layers = [("Aplanado", 0.02), block, ("Yeso", 0.01)]
        loc.meanDay(month=5, year=2025); wall.Tsa()
        ti = wall.solve()
    """

    def __init__(self, location, tilt=90, azimuth=0, absortance=0.8, layers=None):
        self.location = location
        self.tilt = tilt
        self.azimuth = azimuth
        self.absortance = absortance
        self.layers = list(layers) if layers else []
        self.setpoint = None          # AC: if None → Tn.mean() (like the 1D)
        self._sys1 = None
        self._solve_df = None
        self._solve_sig = None
        self._ac_df = None
        self._ac_sig = None
        self.energy_transfer = None
        self.Qout = None
        self.cooling_energy = None
        self.heating_energy = None
        self.days = None
        self.day_error = None
        self.converged = None
        self.inner_iterations = None
        self.solve_dataframe = None

    # --- weather/solar: the 1D chain is reused ---
    def _system1d(self):
        from .ehframe import System
        if self._sys1 is None or self._sys1.location is not self.location:
            self._sys1 = System(self.location, tilt=self.tilt, azimuth=self.azimuth,
                                absortance=self.absortance, layers=[])
        else:
            self._sys1.tilt = self.tilt
            self._sys1.azimuth = self.azimuth
            self._sys1.absortance = self.absortance
        return self._sys1

    def Tsa(self):
        """``Tsa(t)`` DataFrame (reuses EPW+pvlib from the 1D ``System``)."""
        return self._system1d().Tsa()

    # --- validation and geometry ---
    def _element(self):
        elems = [l for l in self.layers if isinstance(l, _ELEMENT_TYPES)]
        if len(elems) != 1:
            raise ValueError(
                "layers must contain exactly one 2D element "
                "(HollowBlock or Slab).")
        return elems[0]

    def _validate(self):
        if len(self.layers) > 7:
            raise ValueError("at most 7 layers (including the 2D element).")
        elem = self._element()
        rt = getattr(elem, "required_tilt", None)
        if rt is not None and self.tilt != rt:
            raise ValueError(
                f"{type(elem).__name__} is only for tilt={rt}° "
                f"(this system has tilt={self.tilt}°).")

    def _build_section(self):
        layers = self.layers
        idx = next(i for i, l in enumerate(layers) if isinstance(l, _ELEMENT_TYPES))
        elem = layers[idx]
        mats = config.materials
        L = [0.0] * 7
        k = [0.0] * 7
        rhoc = [0.0] * 7
        for p, l in enumerate(layers):
            if isinstance(l, _ELEMENT_TYPES):
                m = mats[l.material_main]
                L[p] = elem.thickness
            else:
                name, Lv = l
                m = mats[name]
                L[p] = float(Lv)
            k[p] = m.k
            rhoc[p] = m.rho * m.c

        if isinstance(elem, Slab):
            def kr(name):
                mm = mats[name]
                return mm.k, mm.rho * mm.c
            k_rib, rc_rib = kr(elem.rib_material)
            k_block, rc_block = kr(elem.block_material)
            k_col, rc_col = kr(elem.topping_material)
            k_fill = rc_fill = 0.0
            hollow = elem.fill_type is Fill.AIR
            if not hollow:
                k_fill, rc_fill = kr(elem.fill_material)
            sec = SlabSection(
                nx=config2d.nx, ny=config2d.ny, L=L, k=k, rhoc=rhoc, layer=idx + 1,
                geom=elem._geom(), k_topping=k_col, rc_topping=rc_col, k_rib=k_rib,
                rc_rib=rc_rib, k_block=k_block, rc_block=rc_block, k_fill=k_fill,
                rc_fill=rc_fill, emissivity=elem.emissivity, beta=float(self.tilt),
                hollow=hollow).build()
            return sec, elem

        a, e = elem._ae()
        kr = rcr = 0.0
        if elem.fill_type is Fill.SOLID:
            fm = mats[elem.fill_material]
            kr, rcr = fm.k, fm.rho * fm.c
        sec = Section2D(nx=config2d.nx, ny=config2d.ny, L=L, k=k, rhoc=rhoc,
                        kr=kr, rhocr=rcr, a=a, e=e, layer=idx + 1,
                        fill_type=elem.fill_type).build()
        return sec, elem

    # --- solution ---
    def _signature(self):
        elem_sigs = tuple(l.signature() if isinstance(l, _ELEMENT_TYPES) else tuple(l)
                          for l in self.layers)
        return (id(self.location), self.location.flag().get("date"),
                self.tilt, self.azimuth, self.absortance, elem_sigs,
                config.version, tuple(sorted(config2d.to_dict().items())))

    def solve(self):
        """
        Runs the day to periodic steady state and returns ``Ti`` as a
        ``pandas.Series`` (aligned to the ``Tsa()`` grid). Stores
        ``energy_transfer`` (= Qin), ``Qout``, ``days`` and ``solve_dataframe``
        (with columns ``Ti, Tso, Tsi, Thueco``).
        """
        self._validate()
        df = self.Tsa()
        sig = self._signature()
        if self._solve_df is not None and self._solve_sig == sig:
            return self._solve_df["Ti"]

        import numpy as _np
        Tsa_arr = df["Tsa"].to_numpy(dtype=_np.float64)
        T0 = float(df["Tn"].mean())
        sec, elem = self._build_section()
        m = sec.mesh
        ho, hi, dt = config.ho, config.hi, float(config.dt)
        La = config.La
        rhoair, cair = config.AIR_DENSITY, config.AIR_HEAT_CAPACITY

        if isinstance(elem, Slab) and elem.fill_type is Fill.AIR:
            g = elem._geom()
            vf = _view_factors(g["cavity_width"], g["cavity"])
            out = solve_day_slab_prod(
                sec.NT, sec.kfield, sec.rhocfield, Tsa_arr, ho, hi, dt,
                m.dx, m.dy, La, m.X, rhoair, cair, T0,
                sec.cav_of, sec.cav_i1, sec.cav_i2, sec.info["cj1"], sec.info["cj2"],
                g["cavity_width"], g["cavity"], elem.emissivity, float(self.tilt), *vf,
                config2d.tol_inner, config2d.tol_day, config2d.max_days,
                config2d.max_inner)
            Ti, Tso, Tsi, Th, Tfield, days, Qin, Qout, day_err, inner_ok, inner_max = out
        elif elem.fill_type is Fill.AIR:   # HollowBlock (wall)
            a, e = elem._ae()
            a21, e22 = a["a21"], e["e22"]
            vf = _view_factors(a21, e22)
            out = solve_day_hueca_prod(
                sec.NT, sec.kfield, sec.rhocfield, Tsa_arr, ho, hi, dt,
                m.dx, m.dy, La, m.X, rhoair, cair, T0,
                m.i1, m.j1, m.i2, m.j2, a21, e22, elem.emissivity, *vf,
                config2d.tol_inner, config2d.tol_day, config2d.max_days,
                config2d.max_inner)
            Ti, Tso, Tsi, Th, Tfield, days, Qin, Qout, day_err, inner_ok, inner_max = out
        else:  # SOLID (HollowBlock or Slab): pure conduction
            out = solve_day_2d(
                sec.NT, sec.kfield, sec.rhocfield, Tsa_arr, ho, hi, dt,
                m.dx, m.dy, La, m.X, rhoair, cair, T0,
                config2d.tol_inner, config2d.tol_day, config2d.max_days,
                config2d.max_inner)
            Ti, Tso, Tsi, Tfield, days, Qin, Qout, day_err, inner_ok, inner_max = out
            Th = _np.full_like(Ti, _np.nan)

        res = df.copy()
        res["Ti"] = Ti
        res["Tso"] = Tso
        res["Tsi"] = Tsi
        res["Thueco"] = Th
        self.solve_dataframe = res
        self._solve_df, self._solve_sig = res, sig
        self.energy_transfer, self.Qout = Qin, Qout
        self.cooling_energy = self.heating_energy = None
        self.days, self.Tfield = days, Tfield
        self._store_convergence(day_err, inner_ok, inner_max, "solve")
        return res["Ti"]

    def _store_convergence(self, day_err, inner_ok, inner_max, who):
        """Store the convergence diagnostics and warn if the solve did not
        converge (day cycle or inner sweeps)."""
        self.day_error = float(day_err)
        self.inner_iterations = int(inner_max)
        self.converged = bool(inner_ok) and (day_err <= config2d.tol_day)
        if not self.converged:
            import warnings
            warnings.warn(
                f"System2D.{who}: not converged "
                f"(day_error={day_err:.3e} °C, tol_day={config2d.tol_day}, "
                f"inner_ok={bool(inner_ok)}, inner_iterations={inner_max}); "
                f"results may not be periodic. Increase config2d.max_days / "
                f"max_inner or relax the tolerances.", RuntimeWarning)

    def solveAC(self):
        """
        Runs the day with **air conditioning**: holds the indoor air fixed at the
        setpoint (``self.setpoint`` or ``Tn.mean()`` by default, like the 1D) and
        returns ``Ti`` as a **constant** ``pandas.Series`` (= setpoint). Stores
        ``cooling_energy`` (Qcool) and ``heating_energy`` (Qheat); ``energy_transfer``
        is left ``None``. Mirror of the 1D ``System.solveAC``; cache separate from
        ``solve()``. The cavity air (AIR) keeps floating: the AC only controls the
        indoor air.
        """
        self._validate()
        df = self.Tsa()
        sig = self._signature() + ("ac", self.setpoint)
        if self._ac_df is not None and self._ac_sig == sig:
            return self._ac_df["Ti"]

        import numpy as _np
        Tsa_arr = df["Tsa"].to_numpy(dtype=_np.float64)
        T0 = float(df["Tn"].mean())
        Tset = float(self.setpoint) if self.setpoint is not None else T0
        sec, elem = self._build_section()
        m = sec.mesh
        ho, hi, dt = config.ho, config.hi, float(config.dt)
        La = config.La
        rhoair, cair = config.AIR_DENSITY, config.AIR_HEAT_CAPACITY

        if isinstance(elem, Slab) and elem.fill_type is Fill.AIR:
            g = elem._geom()
            vf = _view_factors(g["cavity_width"], g["cavity"])
            out = solve_day_slab_ac(
                sec.NT, sec.kfield, sec.rhocfield, Tsa_arr, ho, hi, dt,
                m.dx, m.dy, La, m.X, rhoair, cair, T0, Tset,
                sec.cav_of, sec.cav_i1, sec.cav_i2, sec.info["cj1"], sec.info["cj2"],
                g["cavity_width"], g["cavity"], elem.emissivity, float(self.tilt), *vf,
                config2d.tol_inner, config2d.tol_day, config2d.max_days,
                config2d.max_inner)
            Ti, Tso, Tsi, Th, Tfield, days, Qcool, Qheat, day_err, inner_ok, inner_max = out
        elif elem.fill_type is Fill.AIR:   # HollowBlock (wall)
            a, e = elem._ae()
            a21, e22 = a["a21"], e["e22"]
            vf = _view_factors(a21, e22)
            out = solve_day_hueca_ac(
                sec.NT, sec.kfield, sec.rhocfield, Tsa_arr, ho, hi, dt,
                m.dx, m.dy, La, m.X, rhoair, cair, T0, Tset,
                m.i1, m.j1, m.i2, m.j2, a21, e22, elem.emissivity, *vf,
                config2d.tol_inner, config2d.tol_day, config2d.max_days,
                config2d.max_inner)
            Ti, Tso, Tsi, Th, Tfield, days, Qcool, Qheat, day_err, inner_ok, inner_max = out
        else:  # SOLID (HollowBlock or Slab): pure conduction
            out = solve_day_2d_ac(
                sec.NT, sec.kfield, sec.rhocfield, Tsa_arr, ho, hi, dt,
                m.dx, m.dy, La, m.X, rhoair, cair, T0, Tset,
                config2d.tol_inner, config2d.tol_day, config2d.max_days,
                config2d.max_inner)
            Ti, Tso, Tsi, Tfield, days, Qcool, Qheat, day_err, inner_ok, inner_max = out
            Th = _np.full_like(Ti, _np.nan)

        res = df.copy()
        res["Ti"] = Ti
        res["Tso"] = Tso
        res["Tsi"] = Tsi
        res["Thueco"] = Th
        self.solve_dataframe = res
        self._ac_df, self._ac_sig = res, sig
        self.cooling_energy, self.heating_energy = Qcool, Qheat
        self.energy_transfer = None
        self.days, self.Tfield = days, Tfield
        self._store_convergence(day_err, inner_ok, inner_max, "solveAC")
        return res["Ti"]

    # --- utilities mirroring System ---
    def add_layer(self, material, width):
        self.layers.append((material, width))
        self._solve_df = self._ac_df = None
        return self.layers

    def remove_layer(self, index):
        del self.layers[index]
        self._solve_df = self._ac_df = None
        return self.layers

    def copy(self):
        return System2D(self.location, tilt=self.tilt, azimuth=self.azimuth,
                        absortance=self.absortance, layers=list(self.layers))

    def info(self):
        print("<class 'enerhabitat.System2D'>")
        print(f"Location: {self.location.city}")
        print(f"Tilt: {self.tilt}°   Azimuth: {self.azimuth}°   "
              f"Absortance: {self.absortance}")
        print(f"2D mesh: {config2d.nx}×{config2d.ny}")
        print(f"Energy transfer (Qin): {self.energy_transfer}")
        print(f"Cooling energy: {self.cooling_energy}   "
              f"Heating energy: {self.heating_energy}")
        print("Layers (outside→inside):")
        for i, l in enumerate(self.layers):
            if isinstance(l, _ELEMENT_TYPES):
                print(f"\t{i+1}: {type(l).__name__}({l.material_main}, "
                      f"{l.fill_type.value}, {l.thickness:.3f} m)")
            else:
                print(f"\t{i+1}: {l[0]}, {l[1]} m")

    # --- geometry/material inspection (without solving) ---
    def section(self):
        """Returns the built ``Section2D`` (arrays ``NT,kfield,rhocfield`` + ``mesh``)."""
        self._validate()
        sec, _ = self._build_section()
        return sec

    def _layer_bounds_mm(self):
        """Internal layer boundaries (cumulative thicknesses) in mm, for the drawing."""
        acc, bounds = 0.0, []
        for l in self.layers:
            t = l.thickness if isinstance(l, _ELEMENT_TYPES) else l[1]
            acc += t
            bounds.append(acc * 1000.0)
        return bounds[:-1]

    def preview(self, field=None, panels=None, backend="auto", save=None):
        """
        Draws the section cut **to scale** (real ``X``×``Y`` ratio, outside on
        top) to review the material assignment, without solving.

        Args:
            field: a single field (``"nodetype"|"materials"|"k"|"rhoc"``).
            panels: list of fields; defaults to ``["nodetype","k","rhoc"]``.
            backend: ``"auto"`` (matplotlib if available, else ASCII), ``"mpl"`` or ``"ascii"``.
            save: path to save the figure (matplotlib only).

        Returns:
            matplotlib ``(fig, axes)``, or ``None`` if ASCII was used.
        """
        self._validate()
        sec, _ = self._build_section()
        mats = config.materials
        if panels is None:
            panels = [field] if field else ["nodetype", "k", "rhoc"]
        bounds = self._layer_bounds_mm()
        use_mpl = backend == "mpl" or (backend == "auto" and _has_mpl())
        if use_mpl:
            return plot_section(sec, panels, bounds, mats, save=save)
        for p in panels:
            print(ascii_section(sec, p, mats))
            print()
        return None

    def section_report(self):
        """Verification table (text): mesh, node types and assigned materials."""
        self._validate()
        sec, _ = self._build_section()
        NT, k, rc, m = sec.NT, sec.kfield, sec.rhocfield, sec.mesh
        nx, ny = NT.shape
        print(f"2D section  {nx}×{ny}   X={m.X*1000:.1f} mm  Y={m.Y*1000:.1f} mm   "
              f"dx={m.dx*1000:.3f}  dy={m.dy*1000:.3f} mm")
        print(f"Block/cavity: i∈[{m.i1},{m.i2})  j∈[{m.j1},{m.j2})")
        print("Node types (NT):")
        vals, counts = np.unique(NT, return_counts=True)
        for v, c in zip(vals, counts):
            print(f"    NT {int(v):>2}  {_NT_NAMES.get(int(v), '?'):<11} {int(c):>8} nodes")
        print("Assigned materials (by k, ρc):")
        lab, cats = _categorize(sec, "materials", config.materials)
        for idx, (txt, is_air) in enumerate(cats):
            mask = lab == idx
            cnt = int(mask.sum())
            jj = np.where(mask.any(axis=0))[0]
            y0, y1 = jj.min() * m.dy * 1000.0, (jj.max() + 1) * m.dy * 1000.0
            extra = "" if is_air else f"  (k={k[mask][0]:g}, ρc={rc[mask][0]:g})"
            print(f"    {txt:<26} {cnt:>8} nodes   y∈[{y0:.1f},{y1:.1f}] mm{extra}")


# =================================================================
#  To-scale section inspection (matplotlib + ASCII) — utility
# =================================================================

_NT_NAMES = {0: "air", 1: "corner", 2: "corner", 3: "corner", 4: "corner",
             5: "outer edge", 6: "side", 7: "side", 8: "inner edge",
             9: "top wall", 10: "bottom wall", 11: "left wall", 12: "right wall",
             13: "interior", 14: "fill"}


def _has_mpl():
    try:
        import matplotlib  # noqa: F401
        return True
    except Exception:
        return False


def _categorize(section, field, materials_dict):
    """
    Labels each node by category according to ``field`` and returns ``(lab, cats)``:
    ``lab`` (nx,ny int) = category index; ``cats`` = list of ``(text, is_air)``.
    """
    NT, k, rc = section.NT, section.kfield, section.rhocfield
    nx, ny = NT.shape
    lab = np.zeros((nx, ny), dtype=np.int64)
    if field == "nodetype":
        uniq = sorted(int(v) for v in np.unique(NT))
        idxmap = {v: i for i, v in enumerate(uniq)}
        for v, i in idxmap.items():
            lab[NT == v] = i
        cats = [(f"NT {v}: {_NT_NAMES.get(v, '?')}", v == 0) for v in uniq]
        return lab, cats
    # categories by material (materials/k/rhoc); the cavity air goes separately
    rev = {(round(mm.k, 9), round(mm.rho * mm.c, 9)): name
           for name, mm in materials_dict.items()}
    keys, kidx = [], {}
    for i in range(nx):
        for j in range(ny):
            key = ("air",) if NT[i, j] == 0 else (round(float(k[i, j]), 9),
                                                  round(float(rc[i, j]), 9))
            if key not in kidx:
                kidx[key] = len(keys)
                keys.append(key)
            lab[i, j] = kidx[key]
    cats = []
    for key in keys:
        if key == ("air",):
            cats.append(("Air (cavity)", True))
            continue
        kk, rr = key
        name = rev.get(key, f"material (k={kk:g})")
        if field == "k":
            txt = f"{name}: k={kk:g} W/mK"
        elif field == "rhoc":
            txt = f"{name}: ρc={rr:g} J/m³K"
        else:
            txt = name
        cats.append((txt, False))
    return lab, cats


def plot_section(section, fields, layer_bounds_mm, materials_dict, save=None):
    """Draws (matplotlib) one or several section fields to scale. Axes in mm."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch, Rectangle

    m = section.mesh
    Xmm, Ymm = m.X * 1000.0, m.Y * 1000.0
    palette = list(plt.cm.tab20.colors)
    air_color = (0.62, 0.80, 0.99)
    fields = list(fields)
    fig, axs = plt.subplots(1, len(fields), figsize=(4.4 * len(fields), 4.2))
    if len(fields) == 1:
        axs = [axs]
    for ax, field in zip(axs, fields):
        lab, cats = _categorize(section, field, materials_dict)
        colors, ci = [], 0
        for _, is_air in cats:
            if is_air:
                colors.append(air_color)
            else:
                colors.append(palette[ci % len(palette)])
                ci += 1
        ax.imshow(lab.T, extent=[0, Xmm, Ymm, 0], aspect="equal",
                  cmap=ListedColormap(colors), vmin=0, vmax=len(cats) - 1,
                  interpolation="nearest")
        for yb in layer_bounds_mm:
            ax.axhline(yb, color="k", lw=0.7, ls=":")
        x0, x1 = m.i1 * m.dx * 1000.0, m.i2 * m.dx * 1000.0
        y0, y1 = m.j1 * m.dy * 1000.0, m.j2 * m.dy * 1000.0
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False,
                               edgecolor="k", lw=1.2, ls="--"))
        ax.set_title(field)
        ax.set_xlabel("width [mm]")
        ax.set_ylabel("thickness [mm] (out→in)")
        ax.legend(handles=[Patch(color=c, label=t) for (t, _), c in zip(cats, colors)],
                  fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.13))
    fig.suptitle("2D section to scale — material assignment")
    fig.tight_layout()
    if save is not None:
        fig.savefig(save, dpi=130, bbox_inches="tight")
    return fig, axs


def ascii_section(section, field, materials_dict, target_cols=60):
    """**To-scale** ASCII drawing of a section field (fallback without matplotlib)."""
    lab, cats = _categorize(section, field, materials_dict)
    nx, ny = lab.shape
    m = section.mesh
    cols = min(target_cols, nx)
    # characters ~2:1 (height:width); adjust rows to respect the X:Y ratio.
    rows = max(1, min(ny, int(round(cols * (m.Y / m.X) / 2.0))))
    glyphs = "█▓▒░+=*o#x%@O="
    gmap, gi = [], 0
    for _, is_air in cats:
        if is_air:
            gmap.append("·")
        else:
            gmap.append(glyphs[gi % len(glyphs)])
            gi += 1
    out = [f"  [{field}]  to-scale cut  X={m.X*1000:.0f}mm × Y={m.Y*1000:.0f}mm "
           f"(outside on top)"]
    for r in range(rows):
        j = min(ny - 1, int((r + 0.5) * ny / rows))
        row = "".join(gmap[lab[min(nx - 1, int((c + 0.5) * nx / cols)), j]]
                      for c in range(cols))
        out.append("  " + row)
    out.append("  legend:")
    for (txt, _), g in zip(cats, gmap):
        out.append(f"    '{g}'  {txt}")
    return "\n".join(out)
