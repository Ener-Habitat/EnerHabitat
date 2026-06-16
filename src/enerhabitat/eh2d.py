"""
=================================================================
 eh2d — Geometría 2D de losa de vigueta y bovedilla (Fase 1)
=================================================================

Port fiel de la geometría/topología del solver C `legacy_eh/2dTfree/`
(ver ``PLAN-2D.md``). Esta fase construye, para una **bovedilla rellena**
(``Bovedilla.RELLENA``, ``tipo 2`` del C):

    NT[i][j]        malla de tipos de nodo (1-8 fronteras/esquinas, 13 interior)
    k[i][j]         conductividad por nodo
    rhoc[i][j]      capacidad térmica por nodo
    X, Y, dx, dy    tamaño de la celda y discretización
    i1, j1, i2, j2  límites del bloque de relleno (bovedilla)

Convención de malla (idéntica al C):
    i = 0 .. nx-1   ancho   X   (i=0 izquierda, i=nx-1 derecha; laterales adiabáticos)
    j = 0 .. ny-1   espesor Y   (j=0 exterior con Tsa/ho, j=ny-1 interior con Tint/hi)

Las matrices se indexan ``A[i, j]`` con forma ``(nx, ny)``, igual que el dump del C
(``dump_NT.dat`` etc. recorren ``for i: for j: print i j A[i][j]``).

La física (ensamble de coeficientes, solver, lazo temporal) llega en fases
posteriores; aquí solo se reproduce la geometría para validarla nodo a nodo.
"""

from dataclasses import dataclass, field
from enum import Enum

import numpy as np

from .config import config, config2d
from .ehtools2d import (solve_day_2d, solve_day_2d_par, solve_day_hueca_prod,
                        solve_day_hueca_prod_par, solve_day_slab_prod,
                        solve_day_slab_prod_par, _view_factors)


class Bovedilla(Enum):
    """Estado de la bovedilla (lo que el ``tipo`` numérico del C controla)."""
    RELLENA = "rellena"            # tipo 2: bloque sólido (relleno kr, rhocr)
    AIRE = "aire"                  # tipo 1: cámara de aire (radiación + Nusselt)  [Fase 6]
    RELLENA_SIMETRICA = "rellena_sim"  # tipo 4: media celda simétrica, rellena    [posterior]


# Mapeo a los enteros `tipo` del C, para leer .inp y golden masters legacy.
TIPO_C = {
    Bovedilla.AIRE: 1,
    Bovedilla.RELLENA: 2,
    Bovedilla.RELLENA_SIMETRICA: 4,
}


@dataclass
class Mesh2D:
    """Dimensiones y límites del bloque de relleno (resultado del cálculo de malla)."""
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
    Reproduce literalmente el cálculo de malla del ``main`` del C (bovedilla
    rellena / una cámara, ``a14 == 0``).

    Args:
        nx, ny (int): número de nodos en ancho y espesor.
        L (sequence[float]): espesores de capa ``[L1..L7]`` (de afuera hacia adentro).
        layer (int): número de capa (desde afuera, 1-based) donde va la bovedilla.
        a (dict): geometría horizontal con claves ``a11,a12,a13,a14,a21,a22,a23``.
        e (dict): espesores de la bovedilla con claves ``e21,e22,e23``.

    Returns:
        Mesh2D
    """
    # YY[0]=0, YY[1..7]=L1..L7; la capa `layer` se sustituye por el alto de la bovedilla.
    YY = [0.0] + [float(L[i]) if i < len(L) else 0.0 for i in range(7)]
    YY[layer] = e["e21"] + e["e22"] + e["e23"]

    # Ancho de celda. a14==0 -> "una cámara" (caso de la bovedilla rellena/hueca).
    if a.get("a14", 0.0) == 0.0:
        X = a["a21"] + a["a11"] + a["a12"] / 2.0
    else:
        X = (a["a11"] + a["a21"] + a["a12"] + a["a22"]
             + a["a13"] + a["a23"] + a["a14"])

    dx = X / nx
    Y = sum(YY)
    dy = Y / ny

    # Desplazamiento en y de la capa de la bovedilla (mismo truncamiento entero que el C).
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
    Port de ``draw_viguetabovedilla2rellena``: malla de tipos de nodo ``NT`` para
    bovedilla rellena. Marca la zona del relleno con ``14`` (temporal); luego
    :func:`set_krhoc_rellena` la convierte en ``13`` al asignarle ``kr,rhocr``.

    Returns:
        np.ndarray (int) de forma (nx, ny).
    """
    NT = np.zeros((nx, ny), dtype=np.int64)

    # Esquinas.
    NT[0, 0] = 1            # sup-izq  (exterior, adiabático izq)
    NT[nx - 1, 0] = 2       # sup-der  (exterior, adiabático der)
    NT[0, ny - 1] = 3       # inf-izq  (interior, adiabático izq)
    NT[nx - 1, ny - 1] = 4  # inf-der  (interior, adiabático der)

    # Laterales adiabáticos.
    NT[0, 1:ny - 1] = 6
    NT[nx - 1, 1:ny - 1] = 7

    # Fronteras convectivas exterior (j=0) e interior (j=ny-1).
    NT[1:nx - 1, 0] = 5
    NT[1:nx - 1, ny - 1] = 8

    # Nodos interiores.
    NT[1:nx - 1, 1:ny - 1] = 13

    # (En el C se reescriben a 13 los bordes del hueco; aquí ya son 13, es no-op.)

    # Zona de relleno: marca temporal 14.
    NT[i1:i2, j1:j2] = 14

    return NT


def set_krhoc_rellena(nx, ny, dx, dy, L, k, rhoc, kr, rhocr, NT):
    """
    Port de ``set_krhocrelleno``: asigna ``k`` y ``rhoc`` por capa (umbrales en j
    por espesores acumulados ``L1, L1+L2, ...``) y luego sobrescribe la zona del
    relleno (``NT==14``) con ``kr, rhocr``, marcándola como nodo interior ``13``.

    Modifica ``NT`` in situ (14 -> 13) y devuelve ``(k_field, rhoc_field)``.
    """
    k = [float(k[i]) if i < len(k) else 0.0 for i in range(7)]
    rhoc = [float(rhoc[i]) if i < len(rhoc) else 0.0 for i in range(7)]
    L = [float(L[i]) if i < len(L) else 0.0 for i in range(7)]

    kf = np.zeros((nx, ny), dtype=np.float64)
    rf = np.zeros((nx, ny), dtype=np.float64)

    # Umbrales acumulados (idénticos a los `for (; j < (L1+...+Ln)/dy; ++j)` del C).
    thr = np.cumsum(L) / dy  # thr[n] = (L1+..+L_{n+1})/dy

    for i in range(nx):
        j = 0
        for n in range(7):
            while j < thr[n] and j < ny:
                kf[i, j] = k[n]
                rf[i, j] = rhoc[n]
                j += 1

    # Relleno (bovedilla): kr, rhocr y NT 14 -> 13.
    fill = NT == 14
    kf[fill] = kr
    rf[fill] = rhocr
    NT[fill] = 13

    return kf, rf


def draw_hueca(nx, ny, i1, j1, i2, j2):
    """
    Port de ``draw_viguetabovedilla2hueca``: malla ``NT`` para bovedilla con
    **cámara de aire** (``tipo 1``). Igual al marco de :func:`draw_rellena`, pero
    la zona del hueco son nodos de aire ``0`` rodeados por paredes:
    ``9`` (superior, j=j1-1), ``10`` (inferior, j=j2), ``11`` (izquierda, i=i1-1),
    ``12`` (derecha, i=i2).
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
    NT[i1:i2, j1 - 1] = 9       # pared superior del hueco
    NT[i1:i2, j2] = 10          # pared inferior
    NT[i1 - 1, j1:j2] = 11      # pared izquierda
    NT[i2, j1:j2] = 12          # pared derecha
    NT[i1:i2, j1:j2] = 0        # aire del hueco
    return NT


def set_krhoc_hueca(nx, ny, dx, dy, L, k, rhoc):
    """
    Port de ``set_krhoc`` (tipo 1): llena ``k,rhoc`` por capas en y (umbrales
    acumulados), **sin** override de relleno. Los nodos de aire/paredes conservan
    el material de su capa (irrelevante en el aire: el caso NT 0 fija T=Thueco).
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
#  Geometría de techo: vigueta y bovedilla, N cavidades, 3 sólidos (Fase 8b)
# =================================================================
#  Tres materiales sólidos (colado, vigueta en L, bovedilla) + N cavidades de
#  aire iguales. La vigueta es una **L**: alma `web` (ancho d1, sube por todo el
#  elemento) + pie `foot` (suma d2, solo en la banda inferior `cover_bottom`),
#  formando la repisa donde se apoya la bovedilla. El colado (capa de compresión)
#  ocupa la banda superior a todo el ancho. L1/acabados NO son parte del elemento.


def compute_mesh_slab(nx, ny, L, layer, web, foot, shoulder, n_cav, cavity_width,
                      colado, cover_top, cavity, cover_bottom, colado_cap=0.0):
    """
    Malla para la losa de techo de N cavidades. Devuelve ``(mesh, info)`` donde
    ``info`` trae los límites enteros internos del elemento y los bounds en x de
    cada cavidad. Mismo truncamiento entero que :func:`compute_mesh`.
    """
    e_thick = colado + cover_top + cavity + cover_bottom
    YY = [0.0] + [float(L[i]) if i < len(L) else 0.0 for i in range(7)]
    YY[layer] = e_thick

    X = 2.0 * (web + foot) + (n_cav + 1) * shoulder + n_cav * cavity_width
    dx = X / nx
    Y = sum(YY)
    dy = Y / ny

    y1 = sum(YY[:layer]) / dy + 0.5
    base = int(y1)

    e21 = colado + cover_top          # "tapa" sobre la cavidad (colado + bovedilla)
    cj1 = int(e21 / dy + base)        # fila de la pared superior = cj1-1; aire desde cj1
    cj2 = int((e21 + cavity) / dy + base)
    jet = base                        # tope del elemento
    jcap = int(colado_cap / dy + base)  # base de la tapa de colado (L2): el alma no sube de aquí
    jcol = int(colado / dy + base)    # frontera colado/bovedilla
    jeb = int((e21 + cavity + cover_bottom) / dy + base)  # base del elemento

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
    Malla ``NT`` de la losa de techo con N cavidades. ``hollow=True`` → cada
    cavidad es aire (0) rodeada de paredes 9/10/11/12; ``hollow=False`` (RELLENA)
    → la cavidad queda como nodo interior (13, material de relleno). Devuelve
    ``(NT, cav_of)`` con ``cav_of[i,j]`` = índice de cavidad de los nodos aire/pared
    (−1 en el resto).
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
    # RELLENA: la zona de cavidad queda 13 (se rellena el material en set_krhoc).
    return NT, cav_of


def set_krhoc_slab(nx, ny, dx, dy, L, k, rhoc, layer, info,
                   k_colado, rc_colado, k_rib, rc_rib, k_block, rc_block,
                   k_fill, rc_fill, cav_i1, cav_i2, hollow):
    """
    Asigna ``k``/``rhoc`` por nodo para la losa de techo con **tres** sólidos:
    capas homogéneas (por umbrales en y), y dentro del elemento — colado (banda
    superior), bovedilla (resto) y vigueta en **L** (alma + pie). Si ``hollow`` es
    falso, rellena la cavidad con ``k_fill/rc_fill``.
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
                km, rm = k_colado, rc_colado      # banda de colado
            else:
                km, rm = k_block, rc_block         # bovedilla
            if is_alma and j >= jcap:
                km, rm = k_rib, rc_rib             # alma de la L: de jcap (bajo la tapa L2) a la base
            elif is_foot and j >= cj2:
                km, rm = k_rib, rc_rib             # pie de la L (solo cover_bottom)
            kf[i, j] = km; rf[i, j] = rm

    if not hollow:
        for c in range(len(cav_i1)):
            i1, i2 = int(cav_i1[c]), int(cav_i2[c])
            kf[i1:i2, cj1:cj2] = k_fill
            rf[i1:i2, cj1:cj2] = rc_fill
    return kf, rf


@dataclass
class SlabSection:
    """Sección de techo de vigueta y bovedilla (N cavidades, 3 sólidos). Expone
    ``NT/kfield/rhocfield/mesh`` como :class:`Section2D` (para el inspector) más los
    arreglos de cavidad que necesita el motor ``solve_day_slab_prod``."""
    nx: int
    ny: int
    L: list
    k: list
    rhoc: list
    layer: int
    geom: dict            # web,foot,shoulder,n_cav,cavity_width,colado,cover_top,cavity,cover_bottom
    k_colado: float
    rc_colado: float
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
            g["colado"], g["cover_top"], g["cavity"], g["cover_bottom"],
            g.get("colado_cap", 0.0))
        NT, cav_of = draw_slab_multi(self.nx, self.ny, info["cav_i1"],
                                     info["cav_i2"], info["cj1"], info["cj2"],
                                     self.hollow)
        kf, rf = set_krhoc_slab(
            self.nx, self.ny, mesh.dx, mesh.dy, self.L, self.k, self.rhoc,
            self.layer, info, self.k_colado, self.rc_colado, self.k_rib,
            self.rc_rib, self.k_block, self.rc_block, self.k_fill, self.rc_fill,
            info["cav_i1"], info["cav_i2"], self.hollow)
        self.mesh, self.NT, self.kfield, self.rhocfield = mesh, NT, kf, rf
        self.cav_of, self.cav_i1, self.cav_i2 = cav_of, info["cav_i1"], info["cav_i2"]
        self.info = info
        return self


@dataclass
class Section2D:
    """
    Descripción de una sección de losa de vigueta y bovedilla.

    Capas L1..L7 (de afuera hacia adentro) con ``k``/``rhoc`` por capa; relleno de
    la bovedilla con ``kr``/``rhocr``; geometría horizontal ``a*`` y espesores de la
    bovedilla ``e2*``. Por ahora solo ``Bovedilla.RELLENA`` (``tipo 2``).
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
    bovedilla: Bovedilla = Bovedilla.RELLENA

    mesh: Mesh2D = field(init=False, default=None)
    NT: np.ndarray = field(init=False, default=None)
    kfield: np.ndarray = field(init=False, default=None)
    rhocfield: np.ndarray = field(init=False, default=None)

    def build(self):
        """Construye ``mesh``, ``NT``, ``kfield``, ``rhocfield`` y devuelve self."""
        m = compute_mesh(self.nx, self.ny, self.L, self.layer, self.a, self.e)
        if self.bovedilla is Bovedilla.RELLENA:
            NT = draw_rellena(m.nx, m.ny, m.i1, m.j1, m.i2, m.j2)
            kf, rf = set_krhoc_rellena(m.nx, m.ny, m.dx, m.dy,
                                       self.L, self.k, self.rhoc,
                                       self.kr, self.rhocr, NT)
        elif self.bovedilla is Bovedilla.AIRE:
            NT = draw_hueca(m.nx, m.ny, m.i1, m.j1, m.i2, m.j2)
            kf, rf = set_krhoc_hueca(m.nx, m.ny, m.dx, m.dy,
                                     self.L, self.k, self.rhoc)
        else:
            raise NotImplementedError(
                f"{self.bovedilla} (tipo 4 simétrico) aún no portado.")
        self.mesh, self.NT, self.kfield, self.rhocfield = m, NT, kf, rf
        return self


# --- inspección ----------------------------------------------------------------

_NODE_GLYPH = {0: ".", 1: "1", 2: "2", 3: "3", 4: "4", 5: "5", 6: "6",
               7: "7", 8: "8", 9: "9", 10: "a", 11: "b", 12: "c",
               13: "·", 14: "#"}


def print_node_scheme(NT, max_cols=120, max_rows=60):
    """
    Imprime el mapa de tipos de nodo (j hacia abajo = exterior→interior, i a lo
    ancho), como el diagrama de ``PLAN-2D.md``. Submuestrea si la malla es grande.
    """
    nx, ny = NT.shape
    istep = max(1, int(np.ceil(nx / max_cols)))
    jstep = max(1, int(np.ceil(ny / max_rows)))
    if istep > 1 or jstep > 1:
        print(f"# submuestreo: cada {istep} en i, {jstep} en j "
              f"(malla real {nx}x{ny})")
    for j in range(0, ny, jstep):
        row = "".join(_NODE_GLYPH.get(int(NT[i, j]), "?")
                      for i in range(0, nx, istep))
        print(row)


def print_material_scheme(section, max_cols=120, max_rows=60):
    """
    Imprime el mapa de materiales (por valor de ``k``) con una leyenda, de modo que
    la **bovedilla** (relleno) se distinga de la **vigueta/capas** (de ahí el nombre
    "vigueta y bovedilla"). j hacia abajo = exterior→interior. Devuelve la leyenda
    ``{glifo: (k, rhoc)}``.
    """
    k = section.kfield
    rhoc = section.rhocfield
    nx, ny = k.shape

    # Glifo por material (valor de k redondeado); el más conductor = denso, relleno = ligero.
    uniq = sorted({round(float(v), 12) for v in np.unique(k)})
    glyphs = "█▓▒░·:. "
    legend = {}
    val2glyph = {}
    for idx, kv in enumerate(sorted(uniq, reverse=True)):  # k grande -> glifo "denso"
        gph = glyphs[idx] if idx < len(glyphs) else "?"
        val2glyph[kv] = gph
        # rhoc representativo de ese material.
        mask = np.isclose(k, kv)
        rc = float(rhoc[mask].flat[0]) if mask.any() else 0.0
        legend[gph] = (kv, rc)

    istep = max(1, int(np.ceil(nx / max_cols)))
    jstep = max(1, int(np.ceil(ny / max_rows)))
    if istep > 1 or jstep > 1:
        print(f"# submuestreo: cada {istep} en i, {jstep} en j "
              f"(malla real {nx}x{ny})")
    for j in range(0, ny, jstep):
        row = "".join(val2glyph[round(float(k[i, j]), 12)]
                      for i in range(0, nx, istep))
        print(row)
    print("leyenda:")
    for gph, (kv, rc) in legend.items():
        print(f"  '{gph}'  k={kv:g} W/mK   rhoc={rc:g} J/m³K")
    return legend


def plot_node_scheme(section, ax=None, field="NT"):
    """
    Dibuja (matplotlib ``imshow``) ``NT``, ``k`` o ``rhoc``. ``field`` ∈
    {"NT","k","rhoc"}. Importa matplotlib de forma perezosa (dependencia opcional).
    """
    import matplotlib.pyplot as plt

    data = {"NT": section.NT, "k": section.kfield, "rhoc": section.rhocfield}[field]
    if ax is None:
        _, ax = plt.subplots()
    # Transponer a [j, i] para que j (espesor) vaya en el eje vertical (exterior arriba).
    im = ax.imshow(data.T, origin="upper", aspect="auto", interpolation="nearest")
    ax.set_xlabel("i (ancho)")
    ax.set_ylabel("j (espesor: exterior→interior)")
    ax.set_title(field)
    plt.colorbar(im, ax=ax)
    return ax


# =================================================================
#  API de producción 2D — reúsa EPW+pvlib vía System (1D)
# =================================================================
#  Elementos 2D que pueden ir como UNA capa dentro de `System2D.layers`:
#  - HollowBlock (Fase 8a): bloque hueco de concreto, solo MUROS (tilt=90).
#  - Slab        (Fase 8b): vigueta y bovedilla, solo TECHOS (tilt=0).  [pendiente]


def _geom_pick(g, friendly, raw, default=None):
    if friendly in g:
        return g[friendly]
    if raw in g:
        return g[raw]
    if default is not None:
        return default
    raise KeyError(f"geometry: falta '{friendly}' (o el crudo '{raw}')")


class HollowBlock:
    """
    Bloque hueco de concreto para **muros** (`tilt=90`). Un solo material con una
    cámara de aire (convección Nusselt de muro + radiación entre paredes). Es una
    bovedilla ``AIRE`` cuyo material de marco es el propio ``material``.

    Args:
        material (str): material del bloque (p.ej. "Concreto"), de ``config``.
        emissivity (float): emisividad de las paredes del hueco (radiación).
        geometry (dict): medidas de la celda; claves amistosas
            ``web``(=a11), ``block_width``(=a21), ``cover_top``(=e21),
            ``cavity``(=e22), ``cover_bottom``(=e23); se aceptan los crudos
            ``a11..e23``. Por simetría ``a12 = 2·web`` salvo que se dé ``a12``.
    """

    bovedilla = Bovedilla.AIRE
    required_tilt = 90

    def __init__(self, material, emissivity=0.9, geometry=None):
        self.material = material
        self.emissivity = emissivity
        self.geometry = dict(geometry or {})

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
        return ("HollowBlock", self.material, self.emissivity,
                tuple(sorted(self.geometry.items())))


class Slab:
    """
    Vigueta y bovedilla para **techos** (`tilt=0`). Tres materiales sólidos —
    ``rib_material`` (vigueta, en **L**), ``block_material`` (bovedilla, que rodea
    las cavidades) y ``colado_material`` (capa de compresión)— más N cavidades de
    aire iguales (``Bovedilla.AIRE``) o de relleno (``Bovedilla.RELLENA``). La
    cavidad es horizontal → Nusselt de techo (Rayleigh). L1/acabados NO son parte
    del elemento: van como capas homogéneas en ``System2D.layers``.

    Args:
        rib_material (str): material de la vigueta (alma + pie).
        bovedilla (Bovedilla): ``AIRE`` (cámara) o ``RELLENA`` (sólida).
        block_material (str): material del bloque de bovedilla (rodea la cavidad);
            por defecto el mismo que ``rib_material``.
        colado_material (str): material del colado; por defecto ``rib_material``.
        fill_material (str|None): material de relleno de la cavidad si ``RELLENA``.
        emissivity (float): emisividad de las paredes del hueco (radiación) si AIRE.
        geometry (dict): claves amistosas ``web``(=d1), ``foot``(=d2),
            ``shoulder``(=d3), ``n_cavities``, ``cavity_width``(=d4), ``colado``
            (=L2+L3, espesor total del colado), ``colado_cap`` (=L2, la tapa de
            colado a todo el ancho por encima del alma; el alma de la L sube solo
            hasta su base, dejando altura L3+cover_top+cavity+cover_bottom; default
            0 → alma a toda la altura), ``cover_top`` (bovedilla sobre la cavidad,
            =L4), ``cavity`` (=L5, alto del hueco), ``cover_bottom`` (=L6). Alias
            crudos ``d1..d4``.
    """

    required_tilt = 0

    def __init__(self, rib_material, bovedilla=Bovedilla.AIRE, block_material=None,
                 colado_material=None, fill_material=None, emissivity=0.9,
                 geometry=None):
        self.rib_material = rib_material
        self.block_material = block_material or rib_material
        self.colado_material = colado_material or rib_material
        self.bovedilla = bovedilla
        self.fill_material = fill_material
        self.emissivity = emissivity
        self.geometry = dict(geometry or {})
        if bovedilla is Bovedilla.RELLENA and not fill_material:
            raise ValueError("Slab RELLENA requiere fill_material.")

    @property
    def material_main(self):
        return self.colado_material

    def _geom(self):
        g = self.geometry
        web = _geom_pick(g, "web", "d1")
        foot = _geom_pick(g, "foot", "d2")
        shoulder = _geom_pick(g, "shoulder", "d3")
        cavity_width = _geom_pick(g, "cavity_width", "d4")
        n_cav = int(g.get("n_cavities", g.get("n_cav", 1)))
        cover_top = g.get("cover_top", 0.0)
        return {"web": web, "foot": foot, "shoulder": shoulder, "n_cav": n_cav,
                "cavity_width": cavity_width, "colado": g["colado"],
                "colado_cap": g.get("colado_cap", 0.0),
                "cover_top": cover_top, "cavity": _geom_pick(g, "cavity", "e22"),
                "cover_bottom": _geom_pick(g, "cover_bottom", "e23")}

    @property
    def thickness(self):
        g = self._geom()
        return g["colado"] + g["cover_top"] + g["cavity"] + g["cover_bottom"]

    def signature(self):
        return ("Slab", self.rib_material, self.block_material, self.colado_material,
                self.bovedilla.value, self.fill_material, self.emissivity,
                tuple(sorted(self.geometry.items())))


# Tipos reconocidos como "elemento 2D" dentro de layers.
_ELEMENT_TYPES = (HollowBlock, Slab)


class System2D:
    """
    Sistema constructivo 2D de producción (misma metodología que ``System`` 1D).

    ``layers`` es la pila afuera→adentro de capas homogéneas ``(material, L)`` y
    **un** elemento 2D (``HollowBlock``/``Slab``); su posición en la lista es su
    orden en la pila. La ``Tsa(t)`` se reúsa de un ``System`` 1D interno
    (EPW+pvlib, al paso ``config.dt``); la geometría se arma con ``Section2D`` y se
    resuelve con el motor JIT correspondiente al tipo de elemento.

    Uso (idéntico al 1D salvo el elemento 2D)::

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
        self._sys1 = None
        self._solve_df = None
        self._solve_sig = None
        self.energy_transfer = None
        self.Qout = None
        self.days = None
        self.solve_dataframe = None

    # --- clima/solar: se reúsa la cadena 1D ---
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
        """DataFrame de ``Tsa(t)`` (reúsa EPW+pvlib del ``System`` 1D)."""
        return self._system1d().Tsa()

    # --- validación y geometría ---
    def _element(self):
        elems = [l for l in self.layers if isinstance(l, _ELEMENT_TYPES)]
        if len(elems) != 1:
            raise ValueError(
                "layers debe contener exactamente un elemento 2D "
                "(HollowBlock o Slab).")
        return elems[0]

    def _validate(self):
        if len(self.layers) > 7:
            raise ValueError("máximo 7 capas (incluido el elemento 2D).")
        elem = self._element()
        rt = getattr(elem, "required_tilt", None)
        if rt is not None and self.tilt != rt:
            raise ValueError(
                f"{type(elem).__name__} es solo para tilt={rt}° "
                f"(este sistema tiene tilt={self.tilt}°).")

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
            k_col, rc_col = kr(elem.colado_material)
            k_fill = rc_fill = 0.0
            hollow = elem.bovedilla is Bovedilla.AIRE
            if not hollow:
                k_fill, rc_fill = kr(elem.fill_material)
            sec = SlabSection(
                nx=config2d.nx, ny=config2d.ny, L=L, k=k, rhoc=rhoc, layer=idx + 1,
                geom=elem._geom(), k_colado=k_col, rc_colado=rc_col, k_rib=k_rib,
                rc_rib=rc_rib, k_block=k_block, rc_block=rc_block, k_fill=k_fill,
                rc_fill=rc_fill, emissivity=elem.emissivity, beta=float(self.tilt),
                hollow=hollow).build()
            return sec, elem

        a, e = elem._ae()
        kr = rcr = 0.0
        if elem.bovedilla is Bovedilla.RELLENA:
            fm = mats[elem.fill_material]
            kr, rcr = fm.k, fm.rho * fm.c
        sec = Section2D(nx=config2d.nx, ny=config2d.ny, L=L, k=k, rhoc=rhoc,
                        kr=kr, rhocr=rcr, a=a, e=e, layer=idx + 1,
                        bovedilla=elem.bovedilla).build()
        return sec, elem

    # --- solución ---
    def _signature(self):
        elem_sigs = tuple(l.signature() if isinstance(l, _ELEMENT_TYPES) else tuple(l)
                          for l in self.layers)
        return (id(self.location), self.location.flag().get("date"),
                self.tilt, self.azimuth, self.absortance, elem_sigs,
                config.version, tuple(sorted(config2d.to_dict().items())))

    def solve(self):
        """
        Corre el día hasta régimen periódico y devuelve ``Ti`` como
        ``pandas.Series`` (alineada a la rejilla de ``Tsa()``). Guarda
        ``energy_transfer`` (= Qin), ``Qout``, ``days`` y ``solve_dataframe``
        (con columnas ``Ti, Tso, Tsi, Thueco``).
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

        par = config2d.parallel   # motor paralelo por default (numba prange sobre filas)

        if isinstance(elem, Slab) and elem.bovedilla is Bovedilla.AIRE:
            g = elem._geom()
            vf = _view_factors(g["cavity_width"], g["cavity"])
            slab_engine = solve_day_slab_prod_par if par else solve_day_slab_prod
            out = slab_engine(
                sec.NT, sec.kfield, sec.rhocfield, Tsa_arr, ho, hi, dt,
                m.dx, m.dy, La, m.X, rhoair, cair, T0,
                sec.cav_of, sec.cav_i1, sec.cav_i2, sec.info["cj1"], sec.info["cj2"],
                g["cavity_width"], g["cavity"], elem.emissivity, float(self.tilt), *vf,
                config2d.tol_inner, config2d.tol_day, config2d.max_days)
            Ti, Tso, Tsi, Th, Tfield, days, Qin, Qout = out
        elif elem.bovedilla is Bovedilla.AIRE:   # HollowBlock (muro)
            a, e = elem._ae()
            a21, e22 = a["a21"], e["e22"]
            vf = _view_factors(a21, e22)
            hueca_engine = solve_day_hueca_prod_par if par else solve_day_hueca_prod
            out = hueca_engine(
                sec.NT, sec.kfield, sec.rhocfield, Tsa_arr, ho, hi, dt,
                m.dx, m.dy, La, m.X, rhoair, cair, T0,
                m.i1, m.j1, m.i2, m.j2, a21, e22, elem.emissivity, *vf,
                config2d.tol_inner, config2d.tol_day, config2d.max_days)
            Ti, Tso, Tsi, Th, Tfield, days, Qin, Qout = out
        else:  # RELLENA (HollowBlock o Slab): conducción pura
            rellena_engine = solve_day_2d_par if par else solve_day_2d
            out = rellena_engine(
                sec.NT, sec.kfield, sec.rhocfield, Tsa_arr, ho, hi, dt,
                m.dx, m.dy, La, m.X, rhoair, cair, T0,
                config2d.tol_inner, config2d.tol_day, config2d.max_days)
            Ti, Tso, Tsi, Tfield, days, Qin, Qout = out
            Th = _np.full_like(Ti, _np.nan)

        res = df.copy()
        res["Ti"] = Ti
        res["Tso"] = Tso
        res["Tsi"] = Tsi
        res["Thueco"] = Th
        self.solve_dataframe = res
        self._solve_df, self._solve_sig = res, sig
        self.energy_transfer, self.Qout = Qin, Qout
        self.days, self.Tfield = days, Tfield
        return res["Ti"]

    # --- utilidades espejo de System ---
    def add_layer(self, material, width):
        self.layers.append((material, width))
        self._solve_df = None
        return self.layers

    def remove_layer(self, index):
        del self.layers[index]
        self._solve_df = None
        return self.layers

    def copy(self):
        return System2D(self.location, tilt=self.tilt, azimuth=self.azimuth,
                        absortance=self.absortance, layers=list(self.layers))

    def info(self):
        print("<class 'enerhabitat.System2D'>")
        print(f"Location: {self.location.city}")
        print(f"Tilt: {self.tilt}°   Azimuth: {self.azimuth}°   "
              f"Absortance: {self.absortance}")
        print(f"Malla 2D: {config2d.nx}×{config2d.ny}")
        print(f"Energy transfer (Qin): {self.energy_transfer}")
        print("Layers (afuera→adentro):")
        for i, l in enumerate(self.layers):
            if isinstance(l, _ELEMENT_TYPES):
                print(f"\t{i+1}: {type(l).__name__}({l.material_main}, "
                      f"{l.bovedilla.value}, {l.thickness:.3f} m)")
            else:
                print(f"\t{i+1}: {l[0]}, {l[1]} m")

    # --- inspección de la geometría/materiales (sin resolver) ---
    def section(self):
        """Devuelve el ``Section2D`` construido (arrays ``NT,kfield,rhocfield`` + ``mesh``)."""
        self._validate()
        sec, _ = self._build_section()
        return sec

    def _layer_bounds_mm(self):
        """Fronteras internas de capa (espesores acumulados) en mm, para el dibujo."""
        acc, bounds = 0.0, []
        for l in self.layers:
            t = l.thickness if isinstance(l, _ELEMENT_TYPES) else l[1]
            acc += t
            bounds.append(acc * 1000.0)
        return bounds[:-1]

    def preview(self, field=None, panels=None, backend="auto", save=None):
        """
        Dibuja el corte de la sección **a escala** (proporción real ``X``×``Y``,
        exterior arriba) para revisar la asignación de materiales, sin resolver.

        Args:
            field: un campo suelto (``"nodetype"|"materials"|"k"|"rhoc"``).
            panels: lista de campos; por defecto ``["nodetype","k","rhoc"]``.
            backend: ``"auto"`` (matplotlib si está, si no ASCII), ``"mpl"`` o ``"ascii"``.
            save: ruta para guardar la figura (solo matplotlib).

        Returns:
            ``(fig, axes)`` de matplotlib, o ``None`` si se usó ASCII.
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
        """Tabla de verificación (texto): malla, tipos de nodo y materiales asignados."""
        self._validate()
        sec, _ = self._build_section()
        NT, k, rc, m = sec.NT, sec.kfield, sec.rhocfield, sec.mesh
        nx, ny = NT.shape
        print(f"Sección 2D  {nx}×{ny}   X={m.X*1000:.1f} mm  Y={m.Y*1000:.1f} mm   "
              f"dx={m.dx*1000:.3f}  dy={m.dy*1000:.3f} mm")
        print(f"Bloque/cavidad: i∈[{m.i1},{m.i2})  j∈[{m.j1},{m.j2})")
        print("Tipos de nodo (NT):")
        vals, counts = np.unique(NT, return_counts=True)
        for v, c in zip(vals, counts):
            print(f"    NT {int(v):>2}  {_NT_NAMES.get(int(v), '?'):<11} {int(c):>8} nodos")
        print("Materiales asignados (por k, ρc):")
        lab, cats = _categorize(sec, "materials", config.materials)
        for idx, (txt, is_air) in enumerate(cats):
            mask = lab == idx
            cnt = int(mask.sum())
            jj = np.where(mask.any(axis=0))[0]
            y0, y1 = jj.min() * m.dy * 1000.0, (jj.max() + 1) * m.dy * 1000.0
            extra = "" if is_air else f"  (k={k[mask][0]:g}, ρc={rc[mask][0]:g})"
            print(f"    {txt:<26} {cnt:>8} nodos   y∈[{y0:.1f},{y1:.1f}] mm{extra}")


# =================================================================
#  Inspección a escala de la sección (matplotlib + ASCII) — utilidad
# =================================================================

_NT_NAMES = {0: "aire", 1: "esquina", 2: "esquina", 3: "esquina", 4: "esquina",
             5: "borde ext", 6: "lateral", 7: "lateral", 8: "borde int",
             9: "pared sup", 10: "pared inf", 11: "pared izq", 12: "pared der",
             13: "interior", 14: "relleno"}


def _has_mpl():
    try:
        import matplotlib  # noqa: F401
        return True
    except Exception:
        return False


def _categorize(section, field, materials_dict):
    """
    Etiqueta cada nodo por categoría según ``field`` y devuelve ``(lab, cats)``:
    ``lab`` (nx,ny int) = índice de categoría; ``cats`` = lista ``(texto, es_aire)``.
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
    # categorías por material (materials/k/rhoc); el aire del hueco va aparte
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
            cats.append(("Aire (hueco)", True))
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
    """Dibuja a escala (matplotlib) uno o varios campos de la sección. Ejes en mm."""
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
        ax.set_xlabel("ancho [mm]")
        ax.set_ylabel("espesor [mm] (ext→int)")
        ax.legend(handles=[Patch(color=c, label=t) for (t, _), c in zip(cats, colors)],
                  fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.13))
    fig.suptitle("Sección 2D a escala — asignación de materiales")
    fig.tight_layout()
    if save is not None:
        fig.savefig(save, dpi=130, bbox_inches="tight")
    return fig, axs


def ascii_section(section, field, materials_dict, target_cols=60):
    """Dibujo ASCII **a escala** de un campo de la sección (respaldo sin matplotlib)."""
    lab, cats = _categorize(section, field, materials_dict)
    nx, ny = lab.shape
    m = section.mesh
    cols = min(target_cols, nx)
    # caracteres ~2:1 (alto:ancho); ajustar filas para respetar la proporción X:Y.
    rows = max(1, min(ny, int(round(cols * (m.Y / m.X) / 2.0))))
    glyphs = "█▓▒░+=*o#x%@O="
    gmap, gi = [], 0
    for _, is_air in cats:
        if is_air:
            gmap.append("·")
        else:
            gmap.append(glyphs[gi % len(glyphs)])
            gi += 1
    out = [f"  [{field}]  corte a escala  X={m.X*1000:.0f}mm × Y={m.Y*1000:.0f}mm "
           f"(exterior arriba)"]
    for r in range(rows):
        j = min(ny - 1, int((r + 0.5) * ny / rows))
        row = "".join(gmap[lab[min(nx - 1, int((c + 0.5) * nx / cols)), j]]
                      for c in range(cols))
        out.append("  " + row)
    out.append("  leyenda:")
    for (txt, _), g in zip(cats, gmap):
        out.append(f"    '{g}'  {txt}")
    return "\n".join(out)
