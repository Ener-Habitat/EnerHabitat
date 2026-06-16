# Golden master 2D — Fase 0 (`PLAN-2D.md`)

Referencias congeladas del solver C legacy (`legacy_eh/2dTfree/standalone/`), caso del
`.inp`: bovedilla **rellena** (`tipo 2`), `nx=ny=160`, `dt=1 s`, mes=5.

## Archivos

| Archivo | Qué es |
|---------|--------|
| `dump_meta.dat` | `nx ny X Y dx dy i1 j1 i2 j2 tipo` (clave‑valor) |
| `dump_NT.dat`   | `i j NT[i][j]` — tipo de nodo (Fase 1) |
| `dump_k.dat`    | `i j k[i][j]` — conductividad (Fase 1) |
| `dump_rhoc.dat` | `i j rhoc[i][j]` — capacidad térmica (Fase 1) |
| `gbv_5_1.csv`   | serie temporal: `t[h] Is Tsa Ta Tparedext Tparedint Tint Tc DeltaT`, cada 600 s |
| `indice_gbv_5_1.csv` | índices: Qin, factor decremento, retardo, Tint media/min/max, TPIcalor/frío, DDHcalor/frío |

## Cómo se generó

```sh
cd legacy_eh/2dTfree/standalone
make DUMP=1 DAYS=1
./conduction.e conduction.e.inp     # salidas en ./dat/
```

**Alcance:** `DAYS=1` → un solo día desde condición inicial uniforme (`Tint=Tc+DtaT`),
**no** el régimen periódico convergido. La serie por tanto no coincide con el
`legacy_eh/2dTfree/dat/gbv_5_1.csv` committeado (ése es convergido, sin columna de tiempo).
Para el golden convergido completo: `make` sin `DAYS` (converge a `error<1e-5`).

## Baseline de tiempo (Fase 0)

- Máquina: macOS arm64 (Apple clang), `-O2 -std=gnu89`.
- **Un día, `nx=ny=160`, `dt=1 s`: 722.55 s de pared** (719.03 user, 2.36 sys), RSS ~4.2 MB.
- `error` día‑a‑día tras el día 1: `1.198e+00` (aún lejos del régimen; converger requiere
  varios días → varios×~12 min).

## Verificación de geometría (ya confirmada)

`tipo 2`; `X=0.195 Y=0.12 dx=0.00121875 dy=0.00075`; `i1,j1,i2,j2 = 16,26,148,133`.
`NT` ∈ {1–8,13}; 25600 nodos; relleno = `(i2-i1)·(j2-j1)=14124` celdas con `k=0.026`,
`rhoc=64000` (resto `k=1.35`, `rhoc=1.8e6`).
