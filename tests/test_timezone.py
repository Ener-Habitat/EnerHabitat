"""
P2 — Zona horaria del EPW: offsets fraccionarios preservados.

El parser truncaba el offset UTC decimal del encabezado EPW
(``int(datos[8].split('.')[0])``: +5.5 → +5, 30 min de error en hora solar) y
usaba ``Etc/GMT±N``, que ni siquiera existe para fracciones. Ahora se usa
``datetime.timezone(timedelta(minutes=round(offset·60)))`` (stdlib, sin pytz):

  1. Offsets fraccionarios (+5.5, +5.75, −3.5, −9.5) → utcoffset exacto.
  2. Offset entero (−6, control) → idéntico al comportamiento previo.
  3. meanDay funciona con un offset fraccionario (índice tz-aware sano).

Correr con pytest o como script:
    .venv/bin/python tests/test_timezone.py
"""

import os
import tempfile
from datetime import timedelta

import enerhabitat as eh

HERE = os.path.dirname(os.path.abspath(__file__))
EPW = os.path.join(HERE, "MEX_MOR_Cuernavaca-Matamoros.Intl.AP.767260_TMYx.2004-2018.epw")


def _epw_with_offset(offset_str):
    """Copia del EPW de Cuernavaca con el campo 8 (offset UTC) editado."""
    with open(EPW, "r") as f:
        lines = f.readlines()
    header = lines[0].split(",")
    header[8] = offset_str
    lines[0] = ",".join(header)
    tmp = tempfile.NamedTemporaryFile("w", suffix=".epw", delete=False)
    tmp.writelines(lines)
    tmp.close()
    return tmp.name


def _offset_of(epw_path):
    loc = eh.Location(epw_path)
    return loc.timezone.utcoffset(None)


def test_fractional_offsets_preserved():
    cases = {"5.5": timedelta(hours=5, minutes=30),
             "5.75": timedelta(hours=5, minutes=45),
             "-3.5": timedelta(hours=-3, minutes=-30),
             "-9.5": timedelta(hours=-9, minutes=-30)}
    for field, expected in cases.items():
        path = _epw_with_offset(field)
        try:
            got = _offset_of(path)
            assert got == expected, f"offset {field}: {got} != {expected}"
        finally:
            os.unlink(path)


def test_integer_offset_unchanged():
    got = _offset_of(EPW)          # Cuernavaca: -6.0
    assert got == timedelta(hours=-6)


def test_meanday_with_fractional_offset():
    path = _epw_with_offset("5.5")
    try:
        loc = eh.Location(path)
        df = loc.meanDay(month=5, year=2025)
        assert df.index.tz is not None
        assert df.index.tz.utcoffset(None) == timedelta(hours=5, minutes=30)
        assert len(df) > 0
    finally:
        os.unlink(path)


if __name__ == "__main__":
    for fn in (test_fractional_offsets_preserved, test_integer_offset_unchanged,
               test_meanday_with_fractional_offset):
        fn()
        print(f"PASS  {fn.__name__}")
    print("\nP2: offsets fraccionarios del EPW preservados ✅")
