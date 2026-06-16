# Fase 0 — Golden master del C legacy 2D (vigueta y bovedilla, flotación libre)

Build *standalone* del solver C de `legacy_eh/2dTfree/`, sin PostgreSQL ni envío de
correo, para usarlo como **golden master** de validación y **baseline de tiempo** del
port a Python (ver `PLAN-2D.md`).

## Qué se parchó respecto al original

Copias de `conduction.c`, `arrays.h`, `sol.h`, `new_input.h`, `tools.h`, `conduction.e.inp`
desde `legacy_eh/2dTfree/`, con estos cambios mínimos:

- **`conduction.c`**: se quitó `#include <libpq-fe.h>` y la llamada final
  `system("java … Mailparametros …")`. Se añadió un bloque de *debug dump* guardado por
  `#ifdef DUMP` (ver abajo). El físico/numérico queda intacto: `dt=1 s`, `nx=ny=160` del `.inp`.
- **`tools.h`**: `database_begin` / `database_end` reducidas a *stubs* vacíos (eran las únicas
  consumidoras de PostgreSQL/`PQ*`). Conservan su firma para no tocar los *call sites*.
- **`Makefile`**: nuevo, sin `-lpq` ni include de PostgreSQL. `readline`/`ncurses` se toman
  de Homebrew (`brew --prefix readline`). Se compila con `-std=gnu89` porque el código usa
  retornos *implicit-int* de K&R que clang moderno rechaza por defecto.

## Compilar y correr

```sh
make            # solver
make DUMP=1     # solver + dumps de geometría (dat/dump_*.dat)
./conduction.e conduction.e.inp
```

Salidas (todas en `./dat/`):

- `gbv_5_1.csv`  — serie temporal (`Is Tsa Ta Tparedext Tparedint Tint Tc DeltaT`), cada 600 s.
  El nombre es `<user>_<mes>_<cont>.csv`; con el `.inp` actual: `user=gbv`, `mes=5`, `cont=1`.
- `indice_gbv_5_1.csv` — índices (Qin, factor de decremento, retardo, Tint media/min/max,
  TPI calor/frío, DDH calor/frío).
- `dump_meta.dat`, `dump_NT.dat`, `dump_k.dat`, `dump_rhoc.dat` — solo con `DUMP=1`.

## Dumps de depuración (`-DDUMP`)

Volcados una sola vez, tras armar la geometría (antes del lazo temporal):

- `dump_meta.dat`: `nx ny X Y dx dy i1 j1 i2 j2 tipo` (clave‑valor).
- `dump_NT.dat`:  `i  j  NT[i][j]`  (tipo de nodo entero).
- `dump_k.dat`:   `i  j  k[i][j]`   (conductividad, `%.17g`).
- `dump_rhoc.dat`:`i  j  rhoc[i][j]`(capacidad térmica, `%.17g`).

Convención de malla (ver `PLAN-2D.md`): `i=0..nx-1` ancho, `j=0..ny-1` espesor
(`j=0` exterior, `j=ny-1` interior). Sirven de referencia nodo‑a‑nodo para las Fases 1–2.

### Verificación de la geometría (caso del `.inp`, `tipo 2` bovedilla rellena)

- `nx=ny=160`, `X=0.195 m`, `Y=0.12 m`, `dx=0.00121875`, `dy=0.00075`.
- `i1=16 j1=26 i2=148 j2=133`.
- Histograma `NT`: esquinas {1,2,3,4}×1, bordes {5,6,7,8}×158, interior 13×24964 → 25600 nodos.
  No aparecen nodos de cámara de aire (0, 9–12): correcto para bovedilla rellena.
- Relleno: exactamente `(i2-i1)·(j2-j1)=14124` celdas con `k=0.026`, `rhoc=64000`
  (resto `k=1.35`, `rhoc=1.8e6` de la capa de concreto).

## Entregables congelados en `tests/golden/2d/`

`dump_meta.dat`, `dump_NT.dat`, `dump_k.dat`, `dump_rhoc.dat` (geometría/material) y,
al terminar la corrida, `gbv_5_1.csv` + `indice_gbv_5_1.csv` (serie e índices) y el
tiempo de pared (baseline).
