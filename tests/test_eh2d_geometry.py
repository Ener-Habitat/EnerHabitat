"""
Fase 1 — Geometría/topología 2D (bovedilla rellena, `tipo 2`).

Compara la geometría construida por ``enerhabitat.eh2d`` contra el *golden master*
volcado por el C legacy (``tests/golden/2d/dump_*.dat``, generado en la Fase 0):

    NT          igualdad EXACTA nodo a nodo
    k, rhoc     rtol 1e-12
    X,Y,dx,dy   rtol 1e-12 ; i1,j1,i2,j2 exactos

Se puede correr con pytest o directamente como script:
    .venv/bin/python tests/test_eh2d_geometry.py
"""

import os

import numpy as np

from enerhabitat.eh2d import Section2D, Fill, TIPO_C

HERE = os.path.dirname(os.path.abspath(__file__))
GOLDEN = os.path.join(HERE, "golden", "2d")
INP = os.path.join(HERE, "..", "legacy_eh", "2dTfree", "standalone", "conduction.e.inp")

# Las fuentes C de referencia (legacy_eh/) no viven en main: están archivadas
# en el tag `archive/0.2.0-dev` y en el repositorio de validación. Para
# restaurarlas localmente:  git archive archive/0.2.0-dev legacy_eh | tar -x
HAS_LEGACY = os.path.isfile(INP)
LEGACY_REASON = ("requiere las fuentes C de referencia (legacy_eh/): "
                 "restaurar con `git archive archive/0.2.0-dev legacy_eh | tar -x`")
try:
    import pytest
    pytestmark = pytest.mark.skipif(not HAS_LEGACY, reason=LEGACY_REASON)
except ImportError:
    pass


# --- lectura de archivos legacy ------------------------------------------------

def read_inp(path):
    """Parsea el `.inp` del C (líneas `clave : valor # comentario`)."""
    params = {}
    with open(path, encoding="latin-1") as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if not line or ":" not in line:
                continue
            key, val = line.split(":", 1)
            params[key.strip()] = val.strip()
    return params


def read_meta(path):
    meta = {}
    with open(path) as f:
        for line in f:
            key, val = line.split()
            meta[key] = float(val) if "." in val or "e" in val.lower() else int(val)
    return meta


def read_field(path, nx, ny, dtype):
    """Lee un dump `i j valor` en una matriz (nx, ny)."""
    arr = np.zeros((nx, ny), dtype=dtype)
    with open(path) as f:
        for line in f:
            i, j, v = line.split()
            arr[int(i), int(j)] = dtype(v)
    return arr


def section_from_inp(params):
    """Construye un Section2D (bovedilla rellena) a partir del `.inp`."""
    g = lambda key: float(params[key])
    L = [g(f"L{n}") for n in range(1, 8)]
    k = [g(f"k{n}") for n in range(1, 8)]
    rhoc = [g(f"rhoc{n}") for n in range(1, 8)]
    a = {key: g(key) for key in ("a11", "a12", "a13", "a14", "a21", "a22", "a23")}
    e = {key: g(key) for key in ("e21", "e22", "e23")}
    sec = Section2D(
        nx=int(params["nx"]), ny=int(params["ny"]),
        L=L, k=k, rhoc=rhoc, kr=g("kr"), rhocr=g("rhocr"),
        a=a, e=e, layer=int(params["layer"]),
        fill_type=Fill.SOLID,
    )
    return sec.build()


# --- fixtures / build ----------------------------------------------------------

def _build():
    params = read_inp(INP)
    assert int(params["tipo"]) == TIPO_C[Fill.SOLID], \
        "el .inp golden debe ser tipo 2 (bovedilla rellena)"
    sec = section_from_inp(params)
    meta = read_meta(os.path.join(GOLDEN, "dump_meta.dat"))
    return sec, meta


# --- pruebas -------------------------------------------------------------------

def test_mesh_matches_c():
    sec, meta = _build()
    m = sec.mesh
    assert (m.nx, m.ny) == (int(meta["nx"]), int(meta["ny"]))
    assert (m.i1, m.j1, m.i2, m.j2) == (
        int(meta["i1"]), int(meta["j1"]), int(meta["i2"]), int(meta["j2"]))
    for name, got in [("X", m.X), ("Y", m.Y), ("dx", m.dx), ("dy", m.dy)]:
        assert np.isclose(got, meta[name], rtol=1e-12, atol=0), \
            f"{name}: py={got!r} c={meta[name]!r}"


def test_NT_exact():
    sec, meta = _build()
    NT_c = read_field(os.path.join(GOLDEN, "dump_NT.dat"), sec.nx, sec.ny, int)
    assert np.array_equal(sec.NT, NT_c), \
        f"NT difiere en {int((sec.NT != NT_c).sum())} nodos"


def test_k_rhoc():
    sec, meta = _build()
    k_c = read_field(os.path.join(GOLDEN, "dump_k.dat"), sec.nx, sec.ny, float)
    rhoc_c = read_field(os.path.join(GOLDEN, "dump_rhoc.dat"), sec.nx, sec.ny, float)
    assert np.allclose(sec.kfield, k_c, rtol=1e-12, atol=0)
    assert np.allclose(sec.rhocfield, rhoc_c, rtol=1e-12, atol=0)


def _demo():
    """Salida visual entendible al correr el módulo como script."""
    from enerhabitat.eh2d import print_node_scheme, print_material_scheme

    sec, meta = _build()
    NT_c = read_field(os.path.join(GOLDEN, "dump_NT.dat"), sec.nx, sec.ny, int)
    k_c = read_field(os.path.join(GOLDEN, "dump_k.dat"), sec.nx, sec.ny, float)
    rhoc_c = read_field(os.path.join(GOLDEN, "dump_rhoc.dat"), sec.nx, sec.ny, float)
    m = sec.mesh

    bar = "═" * 70
    print(bar)
    print(f"  FASE 1 · Geometría 2D vigueta y bovedilla (RELLENA, tipo 2)")
    print(f"  Python (eh2d)  vs  C legacy (golden master)   malla {sec.nx}×{sec.ny}")
    print(bar)

    # Tabla de parámetros de malla: Python vs C.
    print(f"\n  {'parámetro':<10}{'Python':>16}{'C (golden)':>16}   estado")
    print("  " + "-" * 56)
    rows = [
        ("X",  m.X,  meta["X"]),  ("Y",  m.Y,  meta["Y"]),
        ("dx", m.dx, meta["dx"]), ("dy", m.dy, meta["dy"]),
        ("i1", m.i1, int(meta["i1"])), ("j1", m.j1, int(meta["j1"])),
        ("i2", m.i2, int(meta["i2"])), ("j2", m.j2, int(meta["j2"])),
    ]
    for name, py, c in rows:
        ok = np.isclose(py, c, rtol=1e-12, atol=0)
        pys = f"{py:.8g}" if isinstance(py, float) else str(py)
        cs = f"{c:.8g}" if isinstance(c, float) else str(c)
        print(f"  {name:<10}{pys:>16}{cs:>16}   {'✓' if ok else '✗ DIFIERE'}")

    # Diferencias nodo a nodo.
    nd_NT = int((sec.NT != NT_c).sum())
    nd_k = int((~np.isclose(sec.kfield, k_c, rtol=1e-12, atol=0)).sum())
    nd_rc = int((~np.isclose(sec.rhocfield, rhoc_c, rtol=1e-12, atol=0)).sum())
    total = sec.nx * sec.ny
    print(f"\n  Nodos que difieren del C  (de {total}):")
    print(f"    NT   : {nd_NT:>6}   {'✓ idénticos' if nd_NT == 0 else '✗'}")
    print(f"    k    : {nd_k:>6}   {'✓ rtol 1e-12' if nd_k == 0 else '✗'}")
    print(f"    rhoc : {nd_rc:>6}   {'✓ rtol 1e-12' if nd_rc == 0 else '✗'}")

    # Esquema de tipos de nodo (malla chica para verlo entero).
    a = {"a11": 0.02, "a12": 0.03, "a13": 0.02, "a14": 0.0,
         "a21": 0.16, "a22": 0.16, "a23": 0.0}
    e = {"e21": 0.02, "e22": 0.08, "e23": 0.02}
    small = Section2D(nx=40, ny=40, L=[0.12, 0, 0, 0, 0, 0, 0],
                      k=[1.35, 0, 0, 0, 0, 0, 0], rhoc=[1.8e6, 0, 0, 0, 0, 0, 0],
                      kr=0.026, rhocr=64000, a=a, e=e, layer=1,
                      fill_type=Fill.SOLID).build()

    print(f"\n{bar}\n  Esquema de tipos de nodo NT (malla {small.nx}×{small.ny}, "
          f"exterior arriba)")
    print("  1-4 esquinas · 5 ext(Tsa) · 8 int(Tint) · 6/7 laterales · '·' interior")
    print(bar)
    print_node_scheme(small.NT)

    print(f"\n{bar}\n  Mapa de materiales — se ve la BOVEDILLA dentro de la vigueta")
    print(bar)
    print_material_scheme(small)

    print(f"\n{bar}")
    print("  Fase 1: la geometría de Python reproduce el golden del C ✅")
    print(bar)


if __name__ == "__main__":
    if not HAS_LEGACY:
        raise SystemExit(f"SKIP: {LEGACY_REASON}")
    for fn in (test_mesh_matches_c, test_NT_exact, test_k_rhoc):
        fn()
    _demo()
