# PLAN-2D — Bitácora de lo hecho (Fases 0–8a)

Detalle de las fases completadas del port 2D (objetivo, entregables, archivos y
resultados de prueba contra el C). El plan vivo (diseño + pendiente) está en
[`PLAN-2D.md`](PLAN-2D.md).

## Fases

Cada fase es un PR pequeño con su prueba. Una fase no se da por buena hasta que su prueba pasa.

### Fase 0 — Golden master del C legacy ✅ HECHA
**Objetivo:** compilar y correr el C en esta máquina, congelar referencias.
- ✅ `legacy_eh/2dTfree/standalone/` con el C parchado: stubs vacíos de
  `database_begin/database_end` (eran las únicas consumidoras de `libpq/PQ*`), eliminado el
  `system("java … Mailparametros")` y el `#include <libpq-fe.h>`. `Makefile` nuevo sin `-lpq`,
  con `readline/ncurses` de Homebrew y `-std=gnu89` (el C usa retornos *implicit-int* de K&R
  que clang moderno rechaza). `dt=1`, `nx=ny=160` del `.inp` intactos.
- ✅ Flag de *debug dump* (`make DUMP=1`, `-DDUMP`) que escribe `dat/dump_meta.dat` (nx,ny,X,Y,
  dx,dy,i1,j1,i2,j2,tipo), `dump_NT.dat`, `dump_k.dat`, `dump_rhoc.dat` (formato `i j valor`),
  una vez tras armar la geometría. (`a,b,c,d`/`T` en un paso fijo quedan para cuando los pida
  la Fase 2/3.)
- ✅ Corrido el caso del `.inp` (bovedilla rellena, `tipo 2`). Series (`gbv_5_1.csv` **con
  columna de tiempo `t[h]`**), índices (`indice_gbv_5_1.csv`) y dumps congelados en
  `tests/golden/2d/` (ver su `README.md`).
- **Resultado:** binario compila y corre. **Baseline de tiempo: 722.55 s/día**
  (`nx=ny=160`, `dt=1 s`, macOS arm64, `-O2`).

**Cambios de alcance acordados durante la Fase 0:**
1. **Columna de tiempo en el CSV.** Se activó (estaba en el `printf` comentado) la escritura
   de `t/3600.` como **primera columna**, en **horas decimales** (1.0 = 1 h, 1.5 = 1 h 30 min).
   Encabezado nuevo: `t Is Tsa Ta T_paredext Tparedint Tint T_n DeltaT_n` / `[h] [W/m2] [oC]…`.
   Antes esta columna se pegaba desde otro archivo.
2. **Corrida de un día (configurable).** Se añadió un tope `MAXDAYS` al lazo de convergencia
   día‑a‑día: `make DAYS=N` corre N días (`DAYS=1` = un día y para); sin `DAYS` (o `0`)
   converge a `error<1e-5` como el original. El golden de Fase 0 es de **un día** desde la
   condición inicial uniforme, **no** el régimen convergido → **no** coincide con el
   `dat/gbv_5_1.csv` committeado (ése es convergido, sin columna de tiempo). C y Python
   correrán el mismo único día para validar. El golden convergido completo queda disponible
   con `make` sin `DAYS` si se necesita como referencia de régimen.

### Fase 1 — Geometría/topología (sin física) ✅ HECHA
**Objetivo:** `eh2d.py`: construir `NT`, `k`, `rhoc`, `dx,dy,X,Y,i1,j1,i2,j2` para bovedilla
rellena (`tipo 2`).
- ✅ `src/enerhabitat/eh2d.py`: enum `Bovedilla` (+`TIPO_C`), `compute_mesh` (cálculo de
  malla del `main`: `X,Y,dx,dy,i1,j1,i2,j2` con el mismo truncamiento entero del C),
  `draw_rellena` (port de `draw_viguetabovedilla2rellena`), `set_krhoc_rellena` (port de
  `set_krhocrelleno`), dataclass `Section2D.build()`, e inspectores `print_node_scheme` /
  `plot_node_scheme`. Matrices con forma `(nx,ny)` indexadas `[i,j]` (igual que el dump).
- ✅ **Prueba** `tests/test_eh2d_geometry.py`: lee el `.inp` + dumps de la Fase 0 y compara
  `NT` (igualdad exacta), `k`/`rhoc` (`rtol 1e-12`) y `X,Y,dx,dy,i1,j1,i2,j2`. **Pasa** en la
  malla `nx=ny=160` del golden. Verificado además el esquema visual en malla chica.
  *(Nota: `pytest` no está instalado en `.venv`; la prueba es auto-ejecutable como script:
  `.venv/bin/python tests/test_eh2d_geometry.py`.)*
- *Pendiente menor:* el cálculo de **índices** del `main` (Qin, decremento, retardo, …) se
  porta junto con el lazo temporal en la Fase 4, donde se necesitan.

### Fase 2 — Ensamble de coeficientes ✅ HECHA
**Objetivo:** `ehtools2d.calculate_coefficients_2d` → `a,b,c,d` para un `T` dado (bovedilla
rellena: solo casos `NT` 1–8 y 13).
- ✅ `src/enerhabitat/ehtools2d.py`: `calculate_coefficients_2d(NT,k,rhoc,To,T,Tsa,Tint,ho,
  hi,dt,dx,dy)`. Port vectorizado (numpy) de `calculate_coefficients` para `NT ⊆ {1-8,13}`;
  media armónica en caras, vecinos N/S diferidos a `d`, `a=aP`, `b=aE`, `c=aW`. Lanza
  `NotImplementedError` si aparecen nodos de hueco (0,9-12). (JIT numba → Fase 5.)
- ✅ **Dump del C** (añadido al standalone): `make DUMPCOEF=1` llena `T,To` con un campo
  determinista (rampa + ruido entero, reproducible exacto en Python), `Tsa=30,Tint=24`, llama
  `calculate_coefficients` una vez y vuelca `dump_coef_{a,b,c,d,T,To}.dat` + `dump_coef_meta`
  (Tsa,Tint,ho,hi,dt), luego sale (sin lazo temporal). Congelados en `tests/golden/2d/`.
- ✅ **Prueba** `tests/test_eh2d_coeffs.py`: compara `a,b,c,d` con el dump usando el mismo
  campo. **Pasa** con error relativo ~1e-16 (≪ `rtol 1e-10`). Muestra nodos representativos
  (vigueta vs bovedilla) lado a lado Python/C.

### Fase 3 — Solver de un paso (línea-TDMA + aire interior) ✅ HECHA
**Objetivo:** `ehtools2d.solve_PQ_2d`: barrido Gauss-Seidel por líneas con lazo interno
`|error|>1e-10` + actualización lumped de `Tint`.
- ✅ `src/enerhabitat/ehtools2d.py`: `solve_step_2d(...)` + helper `_tdma_rows` (TDMA en x
  vectorizado sobre todas las filas, mismo orden de operaciones del C). Lazo interno
  `do…while(|error|>1e-10)` con `error = Σ(T-Tn)/T/(nx·ny)` (con signo), coeficientes
  recalculados cada iteración con la `T` más reciente (vecinos en y diferidos), y
  actualización del aire interior `Tint = (Qh + (ρc·La·X/dt)·Ti)·dt/(ρc·La·X)`.
- ✅ **Dump del C** (`make DUMPSTEP=1`): un `solve_PQ` desde un campo determinista, vuelca
  `dump_step_{T0,T}.dat` + meta (Tsa,Tint_in/out,ho,**hi efectivo**,dt,La,X,ρair,cair,Qin,
  inner_iters). Congelado en `tests/golden/2d/`.
- ✅ **Hallazgo:** `convective_coefficients` **sobrescribe `hi`** según `beta` (muro `beta=90`
  → `hi=8.1`, ignora el `hi=6.6` del `.inp`); el solver usa 8.1. Se replica. `interchange(A,B)`
  copia `A←B`.
- ✅ **Prueba** `tests/test_eh2d_step.py`: **iteraciones internas py=32 = c=32**; campo `T`
  max|Δ|=**1.4e-14 °C**; `Tint` coincide exacto — todo ≪ `atol 1e-6`. Muestra el perfil de
  temperatura exterior→interior.

### Fase 4 — Integración temporal + convergencia día-a-día (bovedilla rellena) ✅ HECHA
**Objetivo:** lazo `t∈[0,86400]` y `while(error>1e-5)`, alimentado por la Tsa sinusoidal
del C (portada al harness de prueba). Calcular series e índices.
- ✅ **Frontera del C** portada a `tests/c_boundary.py`: `sol.h` (día juliano, declinación,
  orto/ocaso, duración del día) + `Ta_Tc_DtaT` + `time_evolution_Tsa`. Validada contra las
  columnas `Tsa,Ta,Is` del CSV golden (max|Δ| ≤ redondeo a 2 decimales; `Tc=29.21`,
  `DtaT=1.25` exactos).
- ✅ **Driver del día** `solve_day` (en `tests/test_eh2d_fullday.py`): replica el `main` para
  bovedilla rellena, **un día** desde IC uniforme (igual que el golden `DAYS=1`). Detalles
  finos replicados: `Tso` pre-solve `/nx`, `Tsi` post-solve `/(nx-1)`, `hi=8.1` (muro),
  acumuladores TPI/DDH/Qin, índices con `contador=86400/dt`. Guarda caché
  `golden/2d/py_day_5_1.csv` + `py_indice_5_1.csv`.
- ✅ **Prueba** `tests/test_eh2d_fullday.py` (~28 min de corrida, cacheada): serie
  `Tsa,Tso,Tsi,Tint` **max|Δ|=0.005 °C** (= redondeo del golden) ≪ `atol 0.1`; **10 índices**
  dentro del redondeo (decremento 0.43, retardo 6.1 h, TPI 76/68, DDH 34.5/14.7). Visual:
  tabla de índices y curva `Tsa(t)`/`Tint(t)` (amortiguamiento + retardo).
- *Nota:* solo se corre 1 día (golden `DAYS=1`); la convergencia día-a-día completa y el
  costo se atacan con JIT en Fase 5/7.

### Fase 5 — Integración al paquete (API 2D, producción) ✅ HECHA
**Objetivo:** API de alto nivel reutilizando EPW+pvlib.
- ✅ `eh2d.py`: clase **`System2D`** que reúsa la cadena clima/solar del paquete (la `Tsa(t)`
  sale de un `System` 1D interno → EPW+pvlib, al paso `config.dt`), arma la geometría por
  capas (+ bovedilla rellena opcional) y corre el motor del día. Exportada en `__init__`
  junto a `Section2D`, `Bovedilla`, `config2d`.
- ✅ **`config2d`** (en `config.py`): `nx, ny, tol_inner, tol_day, max_days` sin tocar el 1D.
- ✅ **JIT numba** (`ehtools2d.solve_day_2d` + `_step_inner`, `@njit(cache=True)`): motor del
  día completo con convergencia día-a-día, kernels compilados (≈ del orden de segundos/día tras
  compilar, vs ~28 min en numpy puro).
- ✅ **Correcciones conscientes** (respaldadas por la regresión de Fase 4, documentadas en
  `ehtools2d`): (1) actualización del aire interior con **un solo `dt`** (el C tenía un `dt²`
  latente que solo coincide a `dt=1`) → el 2D reduce EXACTO al 1D a cualquier `dt`;
  (2) promedios de superficie `Tso`/`Tsi` a `/(nx-1)`.
- ✅ **Prueba** `tests/test_eh2d_package.py` sobre el EPW de `tests/`:
  **reduce al 1D** (capa homogénea, `max|Δ|=3.1e-3 °C` ≪ 0.1) · **periodicidad** (converge en
  3 días) · **balance de energía** (`Qin=Qout`, desbalance 0.00 %). Visual: curvas `Ti(t)`
  1D y 2D superpuestas.

### Fase 6 — Bovedilla con cámara de aire (`tipo 1`) ✅ HECHA
**Objetivo:** radiación entre paredes (factores de vista), `hh` (Nusselt muro/techo),
nodo de aire del hueco `Thueco`, y casos `NT` 0, 9–12.
- ✅ **Geometría** (`eh2d`): `draw_hueca`, `set_krhoc_hueca`, `Bovedilla.AIRE` en
  `Section2D.build`. Golden `tipo 1` regenerado (`.inp` `tipo 1`, `cont=2`) en
  `tests/golden/2d/hueca/` (geometría + un paso + día).
- ✅ **Física del hueco** (`ehtools2d`): kernel njit `_step_hueca` con factores de vista
  (`_view_factors`), radiación Stefan-Boltzmann entre las 4 paredes, Nusselt
  `hh=0.4005·|ΔTud|^0.3033/e22^0.0901` (muro), nodo de aire `Thueco` (lumped, un solo `dt`),
  y casos `NT` 0/9-12. Wrappers `solve_step_hueca` (un paso) y `solve_day_hueca` (día, njit).
- ✅ **Prueba** `tests/test_eh2d_hueca.py`:
  - geometría exacta (0 diffs; 14124 aire, paredes 9/10=132, 11/12=107);
  - **un paso a precisión de máquina**: `hh` Δ=3e-16, `Thueco`/`Tint` exactos, campo `T`
    1.8e-14 °C, iters 32=32;
  - **día completo** (1 día, casa con golden `DAYS=1`): serie `Tsa/Tso/Tsi/Tint` max|Δ|≤0.013 °C
    e índices dentro del redondeo (decremento 0.65, retardo 4.0 h, TPI 63.9/49.3, DDH 51.8/23.5).
- *Notas:* sólo muro (`beta=90`); el techo (`beta=0`, Rayleigh) queda pendiente. `System2D`
  aún enruta al motor de relleno; conectar `Bovedilla.AIRE` a `solve_day_hueca` en producción
  es un seguimiento menor.

### Fase 7 — Desempeño y paralelización ✅ HECHA
**Objetivo:** medir y analizar viabilidad.
- ✅ **Hallazgo clave:** el método del C **ya es Jacobi por líneas** (en cada iteración interna
  se ensambla con un único snapshot de `T` y se resuelven todas las filas a un buffer antes de
  `T←Tnew`). Las filas son **independientes** → paralelizar el barrido sobre `j` con `prange`
  ejecuta el MISMO algoritmo, **sin costo de convergencia**. *Red-black es innecesario* (no hay
  acoplamiento Gauss-Seidel que desenredar).
- ✅ **Variante paralela** (`ehtools2d`): `_step_inner_par` + `solve_day_2d_par`
  (`@njit(parallel=True)`, `prange` en ensamble/TDMA/reducción; buffers `P,Q` locales por hilo).
  Portable: threading layer interno de numba (workqueue), **sin libs externas**.
- ✅ **Benchmark** `tests/test_eh2d_perf.py` (esta máquina, 18 núcleos):
  - **correctitud**: paralelo reproduce serial **al bit** (Δ=0.0);
  - **speedup por hilos modesto**: 1.27× @ 4 hilos, no escala (filas cortas, granularidad fina);
  - **escalado de malla** (serial): 0.1/0.5/6.0/19.6 ms/paso para 40²/80²/160²/240²;
  - **barrido de `dt`** (esquema implícito): subir `dt` da mejora **sublineal** porque el lazo
    interno necesita más iteraciones por paso — `dt=10` → ~1.6× con Δ<0.001 °C (punto dulce);
    `dt=60` Δ≈0.005 °C; `dt=300` Δ≈0.02 °C.
- **Recomendación final:** la paralelización con `numba prange` es portable y correcta pero
  rinde poco (~1.3×) por el TDMA secuencial en filas cortas; la palanca práctica mayor es
  **`dt=10`** (~1.6×, error despreciable). Combinadas ~2×. Ganancias mayores exigirían cambio
  algorítmico (line-Gauss-Seidel real para reducir iteraciones internas) o GPU (no portable).

---

## Fase 8a — Bloque hueco de concreto (muros): entregado

- `ehtools2d.solve_day_hueca_prod`: versión de **producción** del día con cámara de aire
  (aire interior con un solo `dt`, superficies `/(nx-1)`, + `Qin/Qout`); no toca la fiel.
- `eh2d.HollowBlock(material, emissivity, geometry)`: elemento 2D de muro (bovedilla `AIRE`,
  `required_tilt=90`, geometría amistosa `web/block_width/cover_*` + alias `a*/e*`). Exportado.
- `eh2d.System2D` **reescrita** (patrón de `System`): `Tsa()` por composición de un `System`
  1D, `_build_section()` mapea `layers`(+elemento) → `Section2D`, `solve()` enruta por tipo
  (AIRE→`solve_day_hueca_prod`), devuelve `pandas.Series` `Ti` alineada a `Tsa()`, guarda
  `energy_transfer/Qout/days/solve_dataframe`. Validación elemento↔`tilt`. Se borró el
  `_build_fields` simplificado.
- Fixture `tests/materials_2d.ini` (Concreto, Mortero, Yeso, EPS).
- **Prueba** `tests/test_eh2d_hollowblock.py` (5/5 pasan, malla chica): metodología (flujo
  idéntico al 1D), periodicidad (4 días), balance de energía (`Qin=Qout`, 0.00 %), guardas de
  orientación (`tilt≠90` falla) y de "un elemento". Demo: factor de decremento 0.47.

## Inspector a escala de la asignación de materiales: entregado

Implementado en `eh2d.py` y `System2D`: `section()`, `preview(field/panels/backend/save)`
(matplotlib a escala — ejes en mm, `aspect="equal"`, fronteras de capa, contorno del bloque,
leyenda; ASCII a escala de respaldo), `section_report()` (tabla de NT + materiales con `k,ρc`
y rango en `y`). matplotlib como extra opcional `enerhabitat[viz]`. Prueba
`tests/test_eh2d_inspect.py` (5/5).

## Fase 8b — Vigueta y bovedilla (techos): entregado

Losa de techo (`tilt=0`) con **tres sólidos** (colado / vigueta en **L** / bovedilla) + **N
cavidades** iguales (aire o relleno), modelo derivado de la Fig. 2b del paper (revisado con el
usuario). NO validado bit-a-bit contra el C: criterio de aceptación = metodología + periodicidad
+ energía (como 8a).
- **Kernels** (`ehtools2d`): `_slab_hh` (Nusselt de muro `beta=90` y de techo Rayleigh `beta=0`,
  portado del C: `Ra=gr·β·(Tdn−Tup)·e22³/ν/α`, `hh=kair/e22·(1+1.44·dot11+dot22)`, estable →
  conducción), `_step_slab` (njit; N cavidades, ensamble guiado por `NT`+`cav_of`, radiación
  Stefan-Boltzmann + Nusselt por hueco), `solve_day_slab_prod` (día con **N nodos `Th`**, un solo
  `dt`, superficies `/(nx-1)`, `Qin/Qout`), wrapper `solve_step_slab`.
- **Geometría** (`eh2d`): `compute_mesh_slab` (`X=2·(d1+d2)+(n+1)·d3+n·d4`; bandas verticales
  colado/cover_top/cavity/cover_bottom; bounds x de cada cavidad), `draw_slab_multi` (N huecos
  + paredes 9-12 + `cav_of`, o relleno), `set_krhoc_slab` (vigueta en L: alma `web` de `jcap`
  (base de la tapa L2) a la base + pie `foot` solo en `cover_bottom`; colado banda superior;
  bovedilla resto; fill si RELLENA),
  dataclass `SlabSection` (expone `NT/kfield/rhocfield/mesh` para el inspector).
- **API** (`eh2d`): clase **`Slab`** (`required_tilt=0`) en `_ELEMENT_TYPES`; `_build_section`/
  `solve()` rutean (RELLENA→`solve_day_2d`, AIRE→`solve_day_slab_prod`, `beta=tilt`); export en
  `__init__`. Material `Bovedilla` añadido al fixture `tests/materials_2d.ini`.
- **Motor serial/paralelo (`config2d.parallel`, default serial):** se añadieron variantes
  `prange` de todos los motores (`solve_day_2d_par`, `solve_day_hueca_prod_par`,
  `solve_day_slab_prod_par`+`_step_slab_par`). `System2D.solve()` elige por la perilla. Default
  **serial** porque el barrido por líneas es de grano fino: el paralelo rinde ~1.3× solo en
  mallas finas y es más lento en chicas (overhead de hilos / busy-wait; test serial 95 s →
  paralelo-default 26 min). Verificado paralelo == serial al bit (prueba opt-in
  `EH_TEST_PARALLEL=1`). Para volumen, paralelizar a nivel de procesos (cada solve serial).
- **Prueba** `tests/test_eh2d_slab.py` (7/7 + 1 opt-in): metodología (Series alineada),
  periodicidad (~5 días), balance `Qin≈Qout` (0.00 %), **decremento AIRE 0.26 > RELLENA-EPS
  0.225**, guarda de orientación (`tilt≠0` falla), inspector (3 cavidades + vigueta en L + 3
  materiales), altura de la L (`test_l_shape_cap_height`), y paralelo==serial (opt-in).
- **Vigueta en L — altura (corregido con el usuario):** la L tiene altura **L3+L4+L5+L6** (sube
  por L3 pero NO por la tapa L2). Implementado con `colado_cap` (=L2, tapa de colado a todo el
  ancho sobre el alma); el alma ocupa los bordes de `jcap` (base de la tapa) a la base del
  elemento. Verificado en `test_l_shape_cap_height`. `colado_cap=0` → alma a toda la altura.


---

# Apéndice — Diseño y referencia (movido desde PLAN-2D)

Diseño, nomenclatura y referencias de verificación que guiaron las fases hechas (antes vivían en `PLAN-2D.md`).

## Decisiones tomadas

1. **Referencia de validación: recompilar el C legacy y resolver.** Se compila el C en
   esta máquina (quitando PostgreSQL y el envío de correo), se corre como *golden master*
   y se vuelcan campos internos para comparar. El C es además el **baseline de tiempo** del
   método original.
2. **Criterio de aceptación: "del orden".** El acuerdo Python↔C se exige a tolerancia
   razonable (mismo orden de magnitud, ver tabla de tolerancias por fase), **no** bit-exacto.
   Pasos internos alimentados con la *misma* Tsa pueden ser más estrictos; la comparación
   final con EPW+pvlib es solo de orden.
3. **Primer corte: bovedilla RELLENA (`tipo 2` del C).** Sin física de cavidad (`hh=1`,
   sin radiación interna). La bovedilla con CÁMARA DE AIRE (`tipo 1`: radiación + Nusselt +
   aire de hueco) y la variante simétrica (`tipo 4`) van en fases posteriores. Ver la
   sección [Nomenclatura de configuraciones](#nomenclatura-de-configuraciones-tipo).
   *(Estado: RELLENA y AIRE ya implementadas y validadas —Fases 2–6 y 8a—; falta `tipo 4`.)*
4. **Tsa inyectable.** El core 2D recibe un arreglo `Tsa(t)` arbitrario. Para validar contra
   el C, el harness usa el modelo sinusoidal del C (`time_evolution_Tsa`, `Ta_Tc_DtaT`);
   para producción se conecta al `Location.Tsa()` (EPW+pvlib) ya existente.
5. **Port fiel primero.** Se replican literalmente los detalles del C (incluido `dt=1 s`)
   para poder comparar; el barrido de `dt` va en una fase posterior, ya con regresión que lo
   respalde. **El `Tsi/(nx-1)` de `max_min` NO es un error**: es el promedio de superficie
   consistente con la discretización (medios nodos en las superficies, ver
   [Hallazgo de revisión](#hallazgo-de-revisión-promedios-de-superficie-nx-1)); se replica
   tal cual.
6. **Código nuevo en módulos propios:** kernels 2D en `src/enerhabitat/ehtools2d.py`,
   geometría/orquestación 2D en `src/enerhabitat/eh2d.py`. `ehframe.py`/`ehtools.py` (1D)
   quedan intactos.

## Mapa del método (qué se porta)

El C resuelve conducción 2D en una sección transversal de losa (`nx`×`ny`; `x`=ancho,
`y`=espesor, capas de afuera→adentro). Piezas:

- **Geometría / topología** (`tools.h`): malla de tipos de nodo `NT[i][j]` (esquinas 1–4,
  fronteras 5–8, paredes de hueco 9–12, interior 13, marca de relleno 14, aire de hueco 0)
  y campos `k[i][j]`, `rhoc[i][j]` por capa + relleno. Índices `i1,j1,i2,j2`, `X,Y,dx,dy`.
- **Ensamble de coeficientes** (`calculate_coefficients`): estencil de 5 puntos con
  conductividad por media armónica → `a,b,c,d` por nodo según `NT`.
- **Solver** (`solve_PQ`): **Gauss-Seidel por líneas** — TDMA implícito en `x` para cada
  fila `j`, con vecinos en `y` tomados del campo `T` que se actualiza in situ; itera
  `do…while(|error|>1e-10)` hasta convergencia interna del paso.
- **Aire interior** `Tint`: nodo de capacitancia concentrada integrado en el tiempo
  (igual patrón que el 1D del paquete).
- **Cavidad (solo tipo 1)**: radiación entre 4 paredes (factores de vista `Fur,Ful,…`),
  `hh` por Nusselt, y aire del hueco `Thueco`. *(fase 6)*
- **Lazo temporal**: `for t in [0,86400] paso dt=1`. **Lazo día-a-día**:
  `while(error>1e-5)` repite el día hasta régimen periódico.
- **Índices** (`abrefile_indice`): Qin, factor de decremento, retardo, Tint media/min/max,
  TPI frío/calor, DDH frío/calor.

La ruta de **clima/solar/Tsa** es física idéntica a la del paquete (pvlib) → se reutiliza
`Location.Tsa()` en producción; no se reimplementa.

---

## Nomenclatura de configuraciones (`tipo`)

En el C la geometría se selecciona con el entero `tipo`, cuyo significado está escondido en
qué pareja de funciones `draw_*` / `set_krhoc*` se llama. Una losa de **vigueta y bovedilla**
se compone de la *vigueta* (la nervadura/viga estructural, normalmente concreto) y la
*bovedilla* (el bloque de relleno entre viguetas). Lo que `tipo` realmente controla es **el
estado de la bovedilla**: rellena de material sólido, o con cámara de aire.

| `tipo` (C) | Funciones C | Qué es físicamente | Física de cavidad |
|------------|-------------|--------------------|-------------------|
| 2 | `draw_viguetabovedilla2rellena` + `set_krhocrelleno` | Bovedilla **rellena** de material (`kr`, `rhocr`). Bloque sólido. | Ninguna (`hh=1`, sin radiación). Nodos `NT` 1–8, 13. |
| 1 | `draw_viguetabovedilla2hueca` + `set_krhoc` | Bovedilla con **cámara de aire** (hueco). | Radiación entre 4 paredes (factores de vista) + convección Nusselt (`hh`) + aire del hueco `Thueco`. Nodos `NT` 0, 9–12. |
| 4 | `draw_viguetabovedillarellena` + `set_krhoc_viguetabovedillarellena` | Bovedilla rellena, **media celda simétrica** (usa `a11/2`). | Ninguna (como tipo 2, geometría simétrica). |

**Nomenclatura nueva propuesta** (sustituir el `tipo` numérico por un enum descriptivo y
alusivo). Sugerencia:

```python
from enum import Enum

class Bovedilla(Enum):
    RELLENA           = "rellena"      # tipo 2: bloque sólido (relleno kr, rhocr)
    AIRE              = "aire"         # tipo 1: cámara de aire (radiación + Nusselt)
    RELLENA_SIMETRICA = "rellena_sim" # tipo 4: media celda simétrica, rellena
```

El nombre alude al **estado de la bovedilla** (que es lo que cambia la física), no a la
vigueta. Los nombres exactos (`RELLENA`/`AIRE` vs. p.ej. `SOLIDA`/`CAMARA_AIRE`) quedan a tu
elección antes de la Fase 1; el resto del plan usa "bovedilla rellena" y "bovedilla con
cámara de aire". El mapeo `tipo→enum` se conserva documentado para poder leer los `.inp` y
los golden masters del C.

---

## Tipos de nodo, geometría y ecuación discretizada

Esta sección documenta **qué ecuación resuelve cada tipo de nodo `NT[i][j]`**, a qué posición
y geometría corresponde, para poder **verificar** la portación contra el C. Todo se deriva de
`calculate_coefficients` (`tools.h`).

### Convención de la malla

- `i = 0 … nx-1` recorre el **ancho** `X` (eje horizontal). `i=0` izquierda, `i=nx-1` derecha.
- `j = 0 … ny-1` recorre el **espesor** `Y` (eje vertical). **`j=0` = superficie EXTERIOR
  (arriba, recibe `Tsa` con `ho`); `j=ny-1` = superficie INTERIOR (abajo, intercambia con el
  aire del recinto `Tint` con `hi`).** Las capas L1…L7 se apilan de afuera (`j` chico) hacia
  adentro (`j` grande).
- Laterales `i=0` e `i=nx-1`: **adiabáticos** (simetría de la celda repetida).
- `Δx = X/nx`, `Δy = Y/ny`.

```
                 i=0                                   i=nx-1
                (lateral                              (lateral
                 adiabático)                           adiabático)
              ┌───────────────────────────────────────────────┐
   j=0        │  1   5   5   5   5   5   5   5   5   5   5  2 │  ← EXTERIOR  (Tsa, ho)
   (exterior) │                                               │
              │  6   13  13  13  13  13  13  13  13  13  13 7 │
   j ↓        │  6   13  13  ┌────────────────────┐  13  13 7 │
   (hacia     │  6   13  13  │  BOVEDILLA         │  13  13 7 │   rellena → 13 (con kr,ρcr)
    adentro)  │  6   13  13  │ i∈[i1,i2) j∈[j1,j2)│  13  13 7 │   aire    → 0 / 9-12
              │  6   13  13  └────────────────────┘  13  13 7 │
              │  6   13  13  13  13  13  13  13  13  13  13 7 │
   j=ny-1     │  3   8   8   8   8   8   8   8   8   8   8  4 │  ← INTERIOR (Tint, hi)
   (interior) └───────────────────────────────────────────────┘
```

### Forma general de la ecuación

Volumen finito implícito, estencil de 5 puntos. Para cada nodo `P=(i,j)`:

```
aP·T_P = aE·T_E + aW·T_W + aN·T_N + aS·T_S + apo·T_P° + S_b
```

- `T_E=(i+1,j)`, `T_W=(i-1,j)` (vecinos en x), `T_N=(i,j-1)` (hacia exterior),
  `T_S=(i,j+1)` (hacia interior). `T_P°` = T del paso de tiempo anterior (`To`).
- Coeficientes de conducción (media armónica de `k`):
  `aN = kh(k_{i,j-1},k_{i,j})·Δx/Δy`, `aS = kh(k_{i,j},k_{i,j+1})·Δx/Δy`,
  `aE = kh(k_{i,j},k_{i+1,j})·Δy/Δx`, `aW = kh(k_{i-1,j},k_{i,j})·Δy/Δx`,
  con `kh(a,b)=2ab/(a+b)`.
- Término temporal: `apo = ρc·Δx·Δy/Δt`.
- `S_b` = fuente de frontera (convección exterior `ho·Δx·Tsa`, interior `hi·Δx·Tint`,
  hueco `hh·Δx·Th` o `hh·Δy·Th`, y radiación del hueco `±Q`).
- `aP` = `apo` + suma de los coeficientes de conducción **activos** + términos convectivos
  activos.

El solver arma `a=aP`, `b=aE`, `c=aW`, y mete `aN·T_N + aS·T_S + apo·T_P° + S_b` en `d`
(los vecinos en `y` van "diferidos" → Gauss-Seidel por líneas; TDMA implícito solo en `x`).

### Catálogo de nodos (`tipo 2` / bovedilla rellena: usa 1–8 y 13)

| `NT` | Posición / geometría | Coef. activos | Frontera (en `d`) | Ecuación física |
|------|----------------------|---------------|-------------------|-----------------|
| **1** | Esquina sup-izq `(0,0)` | aS, aE | `+ho·Δx·Tsa` | Conv. exterior arriba + adiabático izq |
| **2** | Esquina sup-der `(nx-1,0)` | aS, aW | `+ho·Δx·Tsa` | Conv. exterior arriba + adiabático der |
| **3** | Esquina inf-izq `(0,ny-1)` | aN, aE | `+hi·Δx·Tint` | Conv. interior abajo + adiabático izq |
| **4** | Esquina inf-der `(nx-1,ny-1)` | aN, aW | `+hi·Δx·Tint` | Conv. interior abajo + adiabático der |
| **5** | Borde superior `(i,0)` | aS, aE, aW | `+ho·Δx·Tsa` | Frontera convectiva exterior |
| **6** | Borde izquierdo `(0,j)` | aN, aS, aE | — | Lateral adiabático izq (`aW=0`) |
| **7** | Borde derecho `(nx-1,j)` | aN, aS, aW | — | Lateral adiabático der (`aE=0`) |
| **8** | Borde inferior `(i,ny-1)` | aN, aE, aW | `+hi·Δx·Tint` | Frontera convectiva interior |
| **13** | Nodo interior | aN, aS, aE, aW | — | Conducción 2D pura |

Ejemplos explícitos (para comparar con el dump del C):

- **NT 5** (borde exterior): `aP = apo + ho·Δx + aS + aE + aW`,
  `d = ho·Δx·Tsa + aS·T_S + apo·To`. Flujo convectivo `ho·Δx·(Tsa − T_P)` en la cara superior.
- **NT 8** (borde interior): `aP = apo + aN + hi·Δx + aE + aW`,
  `d = aN·T_N + hi·Δx·Tint + apo·To`. Flujo `hi·Δx·(Tint − T_P)` en la cara inferior.
- **NT 13** (interior): `aP = apo + aN + aS + aE + aW`, `d = aN·T_N + aS·T_S + apo·To`.

### Catálogo de nodos de cámara de aire (`tipo 1`, Fase 6: añade 0 y 9–12)

La bovedilla con aire reemplaza la zona `i∈[i1,i2)`, `j∈[j1,j2)` por nodos de aire y sus
paredes. `Th` = temperatura del aire del hueco, `hh` = coef. convectivo del hueco (Nusselt),
`Q*` = intercambio radiativo entre paredes (factores de vista).

```
       i1-1  i1            i2-1  i2
        │    │              │    │
 j1-1   ·   9    9    9    9   ·       NT 9  : pared SUPERIOR del hueco (j=j1-1)
 j1     11   0    0    0    0   12      NT 11 : pared IZQUIERDA (i=i1-1)
        11   0    0    0    0   12      NT 0  : aire del hueco  → T = Th
 j2-1   11   0    0    0    0   12      NT 12 : pared DERECHA   (i=i2)
 j2     ·   10   10   10   10   ·       NT 10 : pared INFERIOR  (j=j2)
```

| `NT` | Posición | Coef. activos | Frontera (en `d`) | Ecuación física |
|------|----------|---------------|-------------------|-----------------|
| **9**  | Pared superior del hueco, `j=j1-1` | aN, aE, aW | `+hh·Δx·Th − Qur − Qud − Qul` | Conv. al aire del hueco (cara inferior) + radiación |
| **10** | Pared inferior del hueco, `j=j2` | aS, aE, aW | `+hh·Δx·Th − Qdl − Qdu − Qdr` | Conv. al aire del hueco (cara superior) + radiación |
| **11** | Pared izquierda, `i=i1-1` | aN, aS, aW | `+hh·Δy·Th − Qlu − Qlr − Qld` | Conv. al aire del hueco (cara derecha) + radiación |
| **12** | Pared derecha, `i=i2` | aN, aS, aE | `+hh·Δy·Th − Qrd − Qrl − Qru` | Conv. al aire del hueco (cara izquierda) + radiación |
| **0**  | Aire del hueco | — | — | `aP=1, aE=aW=0, d=Th` ⇒ `T=Th` (se fija al nodo lumped del hueco) |

> Nota: `NT 14` es **marca temporal** de "relleno" en el `draw_*`; `set_krhoc*` la convierte a
> `13` tras asignar `kr,ρcr`. No llega al solver.

### Cómo verificar (entregable de verificación)

- El C ya trae los `printf` del esquema `NT` comentados en cada `draw_*`. En la **Fase 0** se
  activan (vía `-DDUMP`) para volcar `NT`, `k`, `rhoc` a archivo.
- En la **Fase 1** `eh2d.py` expondrá un `print_node_scheme(NT)` que imprime el mismo mapa
  (como el diagrama de arriba) y un `plot` opcional (matplotlib `imshow`) coloreando por tipo
  de nodo y por material, para inspección visual.
- La prueba de Fase 1 compara el `NT` de Python contra el dump del C **nodo a nodo** (igualdad
  exacta), de modo que el esquema visual y el numérico coinciden.

---

## Hallazgo de revisión: promedios de superficie `(nx-1)`

Revisión pedida sobre el `Tsi/(nx-1)` de `max_min`. **No es un error.** Con `nx` nodos a lo
ancho hay `nx-1` celdas y los nodos laterales (`i=0`, `i=nx-1`, fronteras adiabáticas) son
medios volúmenes de control; el promedio de temperatura de superficie consistente con esa
discretización divide la suma de nodos entre `nx-1`, no entre `nx`. Es resultado de la
discretización en superficie, no descuido.

Lo que sí conviene anotar es una **inconsistencia interna del C** detectada al revisar:

- La columna impresa `Tparedint` (`Tsi`) sale de `max_min` → divide entre **`nx-1`** (sobre
  el campo `T` ya resuelto del paso).
- La columna impresa `Tparedext` (`Tso`) sale de `Tsout` → divide entre **`nx`** (sobre el
  campo `T` previo a resolver).
- El `Tsint(...)` de la línea 266 (que divide entre `nx`) queda **muerto**: `max_min` lo
  sobrescribe antes de imprimir.

Por el mismo argumento de medios nodos, **ambas** superficies (exterior `j=0` e interior
`j=ny-1`) deberían dividir entre `nx-1`. Decisión para el plan:

- **Fase 0–4 (port fiel):** replicar exactamente lo que hace cada *call site* del C
  (`Tsi`=`/(nx-1)`, `Tso`=`/nx`, timing pre/post-solve incluido) para que el golden master
  coincida.
- **Fase 5 (producción):** estandarizar **todos** los promedios de superficie a `/(nx-1)`
  (y, si se quiere ser riguroso, a la forma trapezoidal `[½T₀+T₁+…+T_{n-2}+½T_{n-1}]/(nx-1)`),
  documentando el cambio como mejora consciente respaldada por la regresión de la Fase 4.

---

## Fase 8 — Cableado de la API de producción (8a ✅ · 8b ✅)
> **8a (muros, bloque hueco)** y **8b (techos, vigueta y bovedilla)**: ambas hechas (resumen
> arriba en este archivo). La **Referencia de diseño** (modelo conceptual, esquemas, API
> `System2D`) de esta sección se realizó en 8a y es la base de 8b.

**Objetivo:** exponer el motor 2D (ya validado) detrás de una API de alto nivel ergonómica,
con la **misma metodología** que el `System` 1D, para **dos** sistemas constructivos con
elemento 2D, cada uno restringido a su orientación física:

- **Bloque hueco de concreto** → **muros** (`tilt=90`). Un solo material (concreto) con
  cámara(s) de aire. La cavidad es **vertical** → Nusselt de **muro** (`beta=90`), que **ya
  está validado** (Fase 6) → es lo más fácil de cablear.
- **Vigueta y bovedilla** → **techos/entrepisos** (`tilt=0`). Vigueta (concreto) + bovedilla
  (rellena de material, o cámara de aire). La cavidad es **horizontal** → Nusselt de **techo**
  (correlación de Rayleigh, `beta=0`), que está **pendiente** de portar.

**Orden de la fase (de fácil a difícil):**
- **Fase 8a — Bloque hueco (muros). ✅ HECHA.** `System2D` cablea el bloque hueco de concreto en
  muro reusando el kernel de cámara de aire validado (`beta=90`).
- **Fase 8b — Vigueta y bovedilla (techos). ⏳ PENDIENTE.** Cablear la losa de techo: bovedilla
  **rellena** (kernel `solve_day_2d`, ya validado) y bovedilla con **aire** (requiere portar el
  Nusselt de techo `beta=0` / Rayleigh en `_step_hueca`).

**Fase 8a — entregado y verificado:** `solve_day_hueca_prod` (producción), `HollowBlock`,
`System2D` reescrita (patrón de `System`, `_build_section`→`Section2D`, enrutado por tipo),
fixture `materials_2d.ini` y prueba `test_eh2d_hollowblock.py` (5/5). Detalle en
[`PLAN-2D-hecho.md`](PLAN-2D-hecho.md).

**Estado del paquete tras 8a:** la 8a unificó la API por el camino correcto (`System2D` →
`Section2D` → motor por tipo, reusando pvlib). Esto **resolvió** los huecos que tenía el
paquete antes de 8a: enrutado por tipo de bovedilla, geometría unificada en `Section2D`
(se borró el `_build_fields` simplificado), conexión `Section2D`+pvlib en `solve()`, y
convenciones de producción del hueco (`solve_day_hueca_prod`). **Queda** lo listado en
[Trabajo pendiente](#trabajo-pendiente-diseño-y-propuestas): Fase 8b (techo), extras
(serial/paralelo, `dt`, multi-hueco, tipo 4) e inspector a escala.

#### Referencia de diseño (implementado en 8a, base de 8b)

*(Esta sección documenta el modelo, los esquemas y la API ya implementados en la Fase 8a;
la Fase 8b reúsa todo esto y solo añade el techo. Los esquemas de bloque hueco y de
vigueta-bovedilla viven aquí.)*

##### Modelo conceptual: pila de capas con UN elemento 2D

Ambos sistemas son una **pila de capas en el espesor `y`** (afuera→adentro). La mayoría de
las capas son **homogéneas** (uniformes en `x`: aplanados, morteros, yeso). **Exactamente
una** capa es el **elemento 2D** (bloque hueco *o* vigueta-bovedilla): internamente alterna
en `x` material sólido y la cavidad/bovedilla. **Su posición en la lista `layers` es su orden**
en la pila (desplaza el bloque en `y`; lo maneja `compute_mesh` con `y1=Σ YY[:layer]`).

Convención de la celda (común a ambos): `x` = ancho de la celda repetida (laterales `i=0` e
`i=nx-1` **adiabáticos** por simetría), `y` = espesor (`j=0` exterior con `Tsa,ho`; `j=ny-1`
interior con `Tint,hi`). La celda vertical del elemento es `recubrimiento | banda media |
recubrimiento` = `e21 | e22 | e23`; la horizontal es `½sólido | hueco | ½sólido` = `a11 | a21
| a12/2` (un hueco).

##### (8a) Bloque hueco de concreto — MURO (`tilt=90`)
Un **solo material** (concreto) con una cámara de aire. Cavidad **vertical** → Nusselt de muro.

```
   BLOQUE HUECO DE CONCRETO  ·  MURO (tilt=90)
   x = a lo largo del muro   ·   y = espesor (exterior → interior)

         i=0      i1                          i2      i=nx-1
          │  a11   │            a21            │ a12/2  │      (a11 = a12/2 = ½ alma)
          ├────────┼───────────────────────────┼────────┤
   EXT →  │████████████████████████████████████████████│  e21   cáscara exterior  (concreto)
  (Tsa,ho)├────────┼───────────────────────────┼────────┤
          │████████│        A I R E            │████████│  e22   alma │ cavidad │ alma
          │  alma  │   (Thueco · radiación)    │  alma  │
          ├────────┼───────────────────────────┼────────┤
   INT →  │████████████████████████████████████████████│  e23   cáscara interior  (concreto)
  (Tint,hi)────────┴───────────────────────────┴────────┘
   █ = concreto (material único)   cavidad = AIRE
   alma a11 (= a12/2) ·  cavidad: ancho a21, alto e22 ·  cáscaras e21 (ext), e23 (int)
   espesor del bloque = e21 + e22 + e23   ·   ancho de celda  X = a11 + a21 + a12/2
```

##### (8b) Vigueta y bovedilla — TECHO (`tilt=0`) · modelo definitivo (Fig. 2b del paper)
**Tres sólidos + aire** (revisado con el usuario sobre la Fig. 2b): **colado** (capa de
compresión, todo el ancho), **vigueta en L** (alma vertical + pie horizontal = repisa de
apoyo), **bovedilla** (bloque que rodea las cavidades) y **N cavidades de aire**. Cavidad
**horizontal** → Nusselt de techo (Rayleigh). El elemento 2D (`Slab`) es la pila **L2–L6**;
**L1 y demás acabados/recubrimientos NO son parte del `Slab`** → entran como capas homogéneas
normales en `layers[]`.

```
   VIGUETA Y BOVEDILLA  ·  TECHO (tilt=0)  ·  Slab = L2..L6   (3 cavidades de ejemplo)
   x = ancho (d)   ·   y = espesor (exterior = arriba → interior = abajo)

        │ d1 │ d2 │ d3 │   d4   │ d3 │   d4   │ d3 │   d4   │ d3 │ d2 │ d1 │
   L2   │▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒│  COLADO (todo el ancho)
   L3   │██│▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒│██│  colado + almas que asoman
   ─────┼──┼──────────────────────────────────────────┼──┼─  línea L3/L4 (colado/bovedilla)
   L4   │██│░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░│██│  bovedilla (techo del hueco)
   L5   │██│░░┌────┐░░┌────┐░░┌────┐░░│██│  N cavidades de AIRE
   L6   │████│░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░│████│  pie de la L + bovedilla
   █ = vigueta en L (rib_material)   ▒ = colado (colado_material)
   ░ = bovedilla (block_material)    ┌┐ = AIRE (cavidad)
   vigueta en L = alma d1 (sube de L6 hasta asomar en L3) + pie d2 (solo L6) ← repisa de apoyo
   X = 2·(d1+d2) + (n+1)·d3 + n·d4    (n cavidades; n+1 segmentos d3 de bovedilla)
   ej. 2·50 + 4·50 + 3·103 = 609 mm   ·   L1/acabados → capas en layers[] (fuera del Slab)
```

##### Física de la cavidad (cuando hay AIRE, ambos sistemas)
El centro de la banda `e22` es un hueco rodeado por 4 paredes; el aire es un nodo lumped:

```
            ┌──── pared SUP (NT 9) ────┐
   pared    │  ·  ·  ·  ·  ·  ·  ·  ·   │  pared      · AIRE del hueco = NT 0 → T = Thueco
   IZQ   →  │  ·   aire del hueco       │  ← DER      · convección hh (Nusselt) pared ⇄ aire
   (NT 11)  │  ·   (nodo Thueco)        │  (NT 12)        muro→beta90 · techo→Rayleigh(beta0)
            └──── pared INF (NT 10) ────┘             · radiación entre las 4 paredes
                                                         (factores de vista + emisividad)
```

##### Nomenclatura (nombres amistosos ↔ medidas crudas del C)
| amistoso | crudo | significado |
|----------|-------|-------------|
| `rib_half` / `web` | `a11` | media nervadura/alma (lado izq.); por simetría `=a12/2` |
| `block_width` | `a21` | ancho del hueco/bovedilla |
| `cover_top` | `e21` | recubrimiento/cáscara exterior |
| `cavity` | `e22` | alto del hueco/bovedilla |
| `cover_bottom` | `e23` | recubrimiento/cáscara interior |
| (multi-hueco) | `a12,a13,a14,a22,a23` | varios huecos por celda (`a14≠0`); default = un hueco |

> **Nota (8b, `Slab`):** la vigueta y bovedilla del techo NO es un rectángulo de un material.
> Usa **tres sólidos** (`rib_material` vigueta, `block_material` bovedilla, `colado_material`
> colado) y la **vigueta en L** (`web=d1` alma + `foot=d2` pie en L6). Su geometría se nombra
> `web/foot/shoulder/n_cavities/cavity_width/colado/cavity/cover_bottom` (ver
> [(8b) modelo definitivo](#8b-vigueta-y-bovedilla--techo-tilt0--modelo-definitivo-fig-2b-del-paper)).
> La tabla de arriba aplica tal cual a `HollowBlock` (8a, un material, vigueta rectangular).

**La bovedilla puede ser:**
- **`RELLENA`** (`tipo 2`): bloque sólido de un material de relleno (`fill_material` → `kr,ρcr`).
  Solo **conducción**.
- **`AIRE`** (`tipo 1`): **cámara de aire**. Las 4 paredes del hueco intercambian con el aire
  del hueco (nodo `Thueco`) por **convección** (coef. `hh` por Nusselt: correlación de **muro**
  si `tilt=90`, de **techo/Rayleigh** si `tilt=0`) y entre sí por **radiación** (factores de
  vista + emisividad `emissivity`, Stefan-Boltzmann). Requiere `emissivity`.

> Los kernels **ya soportan** capas antes/después: la malla usa `L1..L7` + un índice `layer`,
> y `compute_mesh` desplaza el bloque con `y1=Σ YY[:layer]`. Lo que falta es **exponerlo**.

##### Principio rector: misma metodología que el `System` 1D

`System2D` se usa **igual** que `System`. La diferencia: `layers` admite, además de tuplas
homogéneas `(material, L)`, **un** objeto de elemento 2D. El patrón es **definir el elemento
primero y luego insertarlo** en `layers` (su posición = su orden en la pila). Hay dos clases
de elemento, una por sistema constructivo:

- **`eh.HollowBlock(...)`** — bloque hueco de concreto (solo **muros**, `tilt=90`).
- **`eh.Slab(...)`** — vigueta y bovedilla (solo **techos**, `tilt=0`).

**(8a) Muro con bloque hueco de concreto:**
```python
import enerhabitat as eh
eh.config.file = "./materials.ini"
loc = eh.Location("./epw/example.epw")

# 1) define el elemento 2D
block = eh.HollowBlock(
    material   = "Concreto",          # material único del bloque
    emissivity = 0.9,                 # radiación en la cámara de aire
    geometry   = {                    # medidas de la celda del bloque
        "web":          0.02,   # a11  alma (= a12/2 por simetría)
        "block_width":  0.16,   # a21  ancho de la cavidad
        "cover_top":    0.02,   # e21  cáscara exterior
        "cavity":       0.08,   # e22  alto de la cavidad
        "cover_bottom": 0.02,   # e23  cáscara interior
    },
)
# 2) insértalo en la pila (afuera → adentro)
wall = eh.System2D(location=loc)
wall.tilt = 90                        # muro (obligatorio para HollowBlock)
wall.azimuth = 90
wall.absortance = 0.6
wall.layers = [("Aplanado", 0.02), block, ("Yeso", 0.01)]

loc.meanDay(month=5, year=2025)
wall.Tsa()
ti = wall.solve()                     # → pandas Series Ti
```

**(8b) Techo de vigueta y bovedilla** (3 sólidos + aire, vigueta en L):
```python
# 1) define el elemento 2D (techo, tilt=0)
slab = eh.Slab(
    rib_material    = "Concreto",          # vigueta (alma + pie, en L)
    block_material  = "Bovedilla",         # bloque que rodea las cavidades
    colado_material = "Concreto",          # capa de compresión (L2+L3)
    bovedilla       = eh.Bovedilla.AIRE,   # o RELLENA (entonces fill_material)
    fill_material   = None,                # material del hueco si RELLENA
    emissivity      = 0.9,                 # requerido si AIRE
    geometry = {                           # cotas crudas del paper (mm→m)
        "web":          0.025,   # d1   alma de la vigueta (sube L6→asoma en L3)
        "foot":         0.025,   # d2   pie de la L (solo L6) → repisa de apoyo
        "shoulder":     0.050,   # d3   bovedilla entre vigueta/cavidades (n+1 segmentos)
        "n_cavities":   3,
        "cavity_width": 0.103,   # d4
        "colado":       0.100,   # L2+L3  capa de compresión (espesor total)
        "colado_cap":   0.050,   # L2     tapa de colado SOBRE el alma (el alma sube L3+L4+L5+L6)
        "cover_top":    0.030,   # L4     bovedilla sobre la cavidad
        "cavity":       0.040,   # L5     alto de la cavidad de aire
        "cover_bottom": 0.030,   # L6     parte baja (aloja el pie de la L)
    },
)
# 2) insértalo en la pila; L1/acabados son capas normales (FUERA del Slab)
roof = eh.System2D(location=loc)
roof.tilt = 0                              # techo (obligatorio para Slab)
roof.absortance = 0.3
roof.layers = [("Impermeabilizante", 0.003), slab, ("Yeso", 0.015)]

loc.meanDay(month=5, year=2025)
roof.Tsa()
ti = roof.solve()
```

El espesor del elemento se **deriva** de su geometría (`cover_top+cavity+cover_bottom`), así
que no se repite un `L`. Se aceptan las claves crudas `a11..e23` como alias (reproducir
`.inp`/golden del C). `System2D` **valida la orientación**: `HollowBlock` exige `tilt=90`,
`Slab` exige `tilt=0`.

##### Diseño de `System2D` (reescritura espejo de `System`)

La `System2D` actual (clase plana) se **reescribe** replicando el patrón de `System`
(`ehframe.py`): propiedades con setters que invalidan caché, `__flag` con
`tsa_date/solve_date/config`, recomputo por `config.version` **y** `config2d.version`.

0. **Clases de elemento 2D** (en `eh2d.py`, exportadas), con base común:
   - **`HollowBlock`** (muro): `material` (str, único), `emissivity` (float),
     `geometry` (`web/block_width/cover_top/cavity/cover_bottom` + alias `a*/e*`). Internamente
     es una bovedilla `AIRE` con `rib_material == fill-context == material`. Exige `tilt=90`.
   - **`Slab`** (techo, **3 sólidos + aire**): `rib_material` (vigueta), `block_material`
     (bloque que rodea la cavidad), `colado_material` (capa de compresión), `bovedilla`
     (`Bovedilla`), `fill_material` (str|None, requerido si RELLENA), `emissivity` (float,
     requerido si AIRE), `geometry` (vigueta en **L**: `web`=d1, `foot`=d2, `shoulder`=d3,
     `n_cavities`, `cavity_width`=d4, `colado`=L2+L3, `cavity`=L5, `cover_bottom`=L6 + alias
     `a*/d*`). Exige `tilt=0`. La vigueta es una **L** (alma `web` que sube de L6 hasta asomar
     en L3 + pie `foot` solo en L6 = repisa de apoyo de la bovedilla), **no** un rectángulo.
   - `HollowBlock`: `thickness = cover_top+cavity+cover_bottom`. `Slab`:
     `thickness = colado + cavity + cover_bottom` (espesor de L2–L6); **L1 y acabados quedan
     fuera** del elemento → van en `layers[]`. Validación de campos en ambas.
1. **Constructor:** `System2D(location, tilt=90, azimuth=0, absortance=0.8, layers=[])`.
   `layers` = lista afuera→adentro de tuplas `(material, L)` y **un** elemento 2D
   (`HollowBlock` o `Slab`). (Sin `bovedilla/geometry` sueltos: viven en el elemento.)
2. **Propiedades asignables** (con `__invalidate_cache`): `location, tilt, azimuth,
   absortance, layers` (+ `add_layer/remove_layer`). Valida coherencia elemento↔`tilt`.
3. **`Tsa()`** — **reúsa la cadena 1D**: internamente compone un `System` 1D
   (misma `location/tilt/azimuth/absortance`, y como `layers` las capas homogéneas + una
   capa equivalente para el `Slab` con su espesor) y devuelve su `Tsa()` (EPW+pvlib, al
   paso `config.dt`). Cero reimplementación de clima/solar. Con caché igual que el 1D.
4. **`solve()`** — flotación libre:
   - **mapea `layers` → `Section2D`**: las tuplas y el `Slab` (en su posición) dan
     `L[1..7]`, `k[1..7]`, `rhoc[1..7]` (vía `config.materials`); el índice del `Slab` en la
     lista (1-based) → `layer`; su `rib_material` → material de esa capa; su `geometry` →
     `a*/e*`; `fill_material` → `kr,ρcr`; `emissivity` → `E`. Construye `Section2D(...).build()`
     (camino validado: `compute_mesh` + `draw_rellena/hueca` + `set_krhoc_*`), `nx,ny` de
     `config2d`. **Soporta capas antes y después** del `Slab` por construcción (`y1=Σ YY[:layer]`);
   - obtiene `Tsa_arr = self.Tsa().Tsa.values` y `T0 = Tn.mean()`;
   - enruta por tipo: `RELLENA` → `solve_day_2d`; `AIRE` → `solve_day_hueca` (**versión de
     producción**, ver punto 6);
   - usa `config.ho`, `config.hi`, `config.dt`, `config.La`, aire de `config` (sin el
     override de muro del C, igual que el 1D);
   - devuelve `Ti` como **`pandas.Series`** alineada a la rejilla de `Tsa()` (mismo índice
     temporal que el 1D), y guarda `energy_transfer` (= `Qin`). Expone además
     `Tso/Tsi/Thueco` (p.ej. atributo `solve_dataframe` o método `series()`).
5. **Métodos espejo:** `info()`, `copy()`, `flag()`, read-only `energy_transfer`
   (`cooling_energy/heating_energy` quedan `None`; `solveAC` 2D se difiere y se anota).
6. **Convenciones de producción del hueco:** añadir variante corregida del solver de aire
   (un solo `dt` en `Thueco` y en `Tint`, superficies `/(nx-1)`), consistente con las
   correcciones de Fase 5. Implementar como `solve_day_hueca` con flag `faithful=False` o un
   `solve_day_hueca_prod` aparte; **no** tocar la versión fiel que valida el golden de Fase 6.
7. **Extras opcionales:** parámetro/propiedad para elegir motor **serial/paralelo**
   (`solve_day_2d` vs `solve_day_2d_par`) y para `dt` efectivo si se decide exponerlo en
   `config2d`. Techo `beta=0` (Rayleigh) de la cámara de aire: portar la rama del C
   (`hh` por Rayleigh) en `_step_hueca` para soportar `tilt=0` con aire.

##### Restricciones y notas del mapeo `layers`(+elemento 2D) → `Section2D`
- **Exactamente un elemento 2D** (`HollowBlock` o `Slab`) en `layers`; su posición (1-based)
  es el `layer` del C. Las capas homogéneas antes/después se apilan normal (el desplazamiento
  `y1=Σ YY[:layer]` ya lo maneja `compute_mesh`).
- El **espesor del elemento** (`cover_top+cavity+cover_bottom`) define `YY[layer]`; coherente
  con `e21+e22+e23`.
- `nx,ny` de `config2d`: `ny` debe ser suficiente para resolver las bandas `e21/e22/e23`
  (avisar si `cavity/dy` es muy chico). Multi-hueco (`a12,a13,a14,a22,a23`) queda soportado por
  `compute_mesh` (rama `a14≠0`) pero la prueba base usa un hueco.
- `System2D._build_section()` arma `Section2D(...).build()` y toma
  `NT, kfield, rhocfield, mesh(dx,dy,X,i1,j1,i2,j2)`. **Se elimina** el `_build_fields`
  simplificado actual (se borra), unificando en un solo camino de geometría.

---

### ✅ Fase 8b — Vigueta y bovedilla (techos) — HECHA

**Entregado y verificado** (prueba `tests/test_eh2d_slab.py`, 7/7 + 1 opt-in, malla chica):
- **Kernels** (`ehtools2d`): `_slab_hh` (Nusselt de **muro** `beta=90` y de **techo** Rayleigh
  `beta=0`, portado del C), `_step_slab` (njit; N cavidades, ensamble guiado por `NT`+`cav_of`,
  3 materiales por `k/rhoc`, radiación+Nusselt por hueco), `solve_day_slab_prod` (día con
  **N nodos `Th`**, un solo `dt`, superficies `/(nx-1)`, `Qin/Qout`), wrapper `solve_step_slab`.
  **Variantes paralelas** (`prange`): `_step_slab_par`+`solve_day_slab_prod_par` (techo) y
  `_step_hueca_par`+`solve_day_hueca_prod_par` (muro); ver extras (default **serial**).
- **Geometría** (`eh2d`): `compute_mesh_slab` (`X=2·(d1+d2)+(n+1)·d3+n·d4`, bandas verticales,
  bounds x de cada cavidad), `draw_slab_multi` (N huecos `0`+paredes `9-12` para AIRE / relleno
  `13` para RELLENA, + `cav_of`), `set_krhoc_slab` (**3 sólidos**: colado / vigueta en **L**
  —alma `web` de `jcap` (base de la tapa L2) a la base + pie `foot` solo en `cover_bottom`— /
  bovedilla; + fill si RELLENA),
  dataclass `SlabSection` (expone `NT/kfield/rhocfield/mesh` para el inspector).
- **API** (`eh2d`): clase **`Slab`** (`rib_material/block_material/colado_material/bovedilla/
  fill_material/emissivity/geometry`, `required_tilt=0`), sumada a `_ELEMENT_TYPES`;
  `System2D._build_section`/`solve()` rutean `Slab` (RELLENA→`solve_day_2d`,
  AIRE→`solve_day_slab_prod` con `beta=tilt`); exportada en `__init__`. Motor serial/paralelo
  según `config2d.parallel` (default serial).
- **Resultados** (techo del paper, 3 cavidades): converge en ~5 días, balance `Qin≈Qout` 0.00 %,
  **decremento AIRE 0.26 > RELLENA(EPS) 0.225** (el hueco transfiere más que el relleno
  aislante), inspector muestra las 3 cavidades + vigueta en L + 3 materiales a escala.
- **Material nuevo** en el fixture `tests/materials_2d.ini`: `Bovedilla` (bloque ligero).

**Vigueta en L — altura (corregido con el usuario):** el alma de la L tiene altura
**L3+L4+L5+L6** (sube por L3 pero **no** por la tapa L2 del colado). Se implementó con la clave
`colado_cap` (=L2, tapa de colado a todo el ancho por encima del alma): el alma ocupa los bordes
desde `jcap` (base de la tapa) hasta la base del elemento. Verificado en
`test_l_shape_cap_height` (con vigueta de `k` distinto al colado: la tapa L2 queda colado, el
alma de `jcap` a la base queda vigueta, altura ≈150 mm). `colado_cap=0` (default) → alma a toda
la altura. Multi-cavidad usa celda completa (sin simetría). Reúsa la Referencia de diseño
(`System2D`/`_build_section`,
solver de aire de producción). **Decisiones tomadas:** geometría de **N cavidades iguales** con
**física de cavidad por hueco** (enfoque A), **celda completa**, **tres materiales sólidos +
aire** y **vigueta en L** (revisado con el usuario sobre la Fig. 2b).

El diseño detallado abajo se conserva como referencia de implementación.

**Modelo definitivo (Fig. 2b).** Tres sólidos: **colado** (capa de compresión, todo el
ancho), **vigueta en L** (alma vertical en el borde + pie horizontal en L6 = repisa de apoyo)
y **bovedilla** (bloque que rodea las N cavidades de aire). El elemento `Slab` es la pila
**L2–L6**; **L1 y los acabados/recubrimientos NO son parte del `Slab`** → van como capas
homogéneas en `layers[]`. (Esquema en
[Referencia de diseño › (8b)](#8b-vigueta-y-bovedilla--techo-tilt0--modelo-definitivo-fig-2b-del-paper).)

```
        │ d1 │ d2 │ d3 │   d4   │ d3 │   d4   │ d3 │   d4   │ d3 │ d2 │ d1 │
   L2   │▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒│  COLADO (todo el ancho)
   L3   │██│▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒│██│  colado + almas que asoman
   L4   │██│░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░│██│  bovedilla (techo del hueco)
   L5   │██│░░┌──┐░░┌──┐░░┌──┐░░│██│  N cavidades de AIRE
   L6   │████│░░░░░░░░░░░░░░░░░░░░░░░░░░░░░│████│  pie de la L + bovedilla
   █=vigueta en L (rib_material)  ▒=colado (colado_material)  ░=bovedilla (block_material)
```
Mapeo a las medidas crudas del C / paper: vigueta `web=d1` (alma) + `foot=d2` (pie, solo L6);
`shoulder=d3` (bovedilla entre vigueta/cavidades, **n+1** segmentos); `cavity_width=d4`
(**n** cavidades). Ancho `X = 2·(d1+d2) + (n+1)·d3 + n·d4` (ej. `2·50+4·50+3·103 = 609 mm`).
Vertical: `colado` = L2+L3 (todo el ancho), `cavity` = L5 (alto del hueco), `cover_bottom`
= L6 (aloja el pie). El alma de la vigueta sube de L6 y **asoma sobre la línea L3/L4** dentro
del colado (detalle a honrar de la figura: el alma interrumpe el colado en L3 en los bordes).

**API objetivo (`Slab`):** ver el ejemplo completo en
[Referencia de diseño › (8b) Techo](#principio-rector-misma-metodología-que-el-system-1d)
(`rib_material/block_material/colado_material/bovedilla/fill_material/emissivity` +
`geometry` con `web/foot/shoulder/n_cavities/cavity_width/colado/cavity/cover_bottom`).

**Pasos:**
1. Clase **`Slab`** (`rib_material/block_material/colado_material/bovedilla/fill_material/
   emissivity/geometry`) + validación `Slab`↔`tilt=0`; reconocerla como elemento 2D
   (`_ELEMENT_TYPES`). `thickness = colado+cavity+cover_bottom` (L2–L6); L1/acabados quedan
   fuera, en `layers[]`.
2. **Geometría N cavidades + vigueta en L + 3 materiales**: `draw_slab_multi` (genera el
   colado a todo el ancho, la **vigueta en L** —alma `d1` que sube y asoma en L3 + pie `d2`
   en L6—, los N rectángulos de aire con sus paredes 9–12, las almas/hombros `d3` de bovedilla
   y la bovedilla que rodea los huecos) y un `set_krhoc_slab` que asigna **tres** `k/ρc`
   (vigueta, bovedilla, colado). `compute_mesh` para N huecos (`X=2·(d1+d2)+(n+1)·d3+n·d4`).
   Celda completa (sin media celda). **Nuevo tipo de nodo** para la pared izq/der de la L y la
   interfaz colado/bovedilla si hace falta distinguir materiales en `calculate_coefficients`.
3. **Física por hueco**: generalizar `_step_hueca`/`solve_day_hueca_prod` a **N nodos
   `Thueco`** (uno por cavidad), con radiación (factores de vista con `l=cavity_width`,
   `h=cavity`) y Nusselt **por cavidad**. Las almas/hombros (`d3`) y la vigueta conducen entre
   cavidades (puente térmico) de forma natural.
4. **Nusselt de techo** (`tilt=0`, correlación de Rayleigh) en la física de cavidad
   (hoy solo muro `beta=90`).
5. `System2D.solve()` enruta `Slab`: RELLENA→`solve_day_2d`, AIRE→solver de aire N-cavidades.
6. **Pruebas 8b**: techo RELLENA y AIRE end-to-end (metodología, periodicidad, energía);
   decremento(aire) > decremento(rellena); capas antes/después (L1/acabados); `inspector`
   muestra las N cavidades, la vigueta en L y los 3 materiales a escala.

## Tolerancias por fase (resumen)

| Fase | Qué se compara | Tolerancia |
|------|----------------|------------|
| 1 | `NT` / `k`,`rhoc` | exacto / `rtol 1e-12` |
| 2 | `a,b,c,d` | `rtol 1e-10` |
| 3 | `T`, `Tint` (1 paso) | `atol 1e-6 °C` |
| 4 | serie `Ti` / índices vs C | `atol 0.1 °C`,`rtol 1%` / ~1–2% |
| 5 | reducción al 1D del paquete | `atol 0.1 °C` |
| 6 | serie/índices bovedilla con aire vs C | "del orden" |
| 7 | paralelo vs serial (ya es Jacobi por líneas) | bit-exacto |
| 8a | muro bloque hueco end-to-end | periodicidad + `Qin≈Qout` |

## Estructura de archivos (actual)

```
legacy_eh/2dTfree/standalone/   # C parchado + dumps + Makefile (Fase 0; flags DUMP/DUMPCOEF/
                                #   DUMPSTEP/DAYS) + conduction_hueca.inp (tipo 1)
src/enerhabitat/eh2d.py         # geometría 2D + Section2D + HollowBlock + System2D + inspectores
src/enerhabitat/ehtools2d.py    # kernels 2D numba (rellena, aire fiel/producción, paralelo)
src/enerhabitat/config.py       # + config2d (nx, ny, tolerancias)
tests/golden/2d/                # referencias del C: rellena (raíz) + hueca/ (tipo 1)
tests/c_boundary.py             # frontera sol-aire del C (sol.h + Ta_Tc_DtaT + time_evolution_Tsa)
tests/materials_2d.ini          # fixture de materiales 2D (Concreto, Mortero, Yeso, EPS)
tests/test_eh2d_geometry.py     # Fase 1
tests/test_eh2d_coeffs.py       # Fase 2
tests/test_eh2d_step.py         # Fase 3
tests/test_eh2d_fullday.py      # Fase 4
tests/test_eh2d_package.py      # Fase 5 (reducción al 1D)
tests/test_eh2d_hueca.py        # Fase 6 (cámara de aire vs C)
tests/test_eh2d_perf.py         # Fase 7 (benchmark/paralelización)
tests/test_eh2d_hollowblock.py  # Fase 8a (muro bloque hueco end-to-end)
tests/test_eh2d_inspect.py      # Inspector a escala (section/preview/section_report)
tests/test_eh2d_slab.py         # Fase 8b (techo vigueta y bovedilla, N cavidades)
```

## Riesgos / notas

- **Convergencia del lazo interno:** el C usa `error` con signo y `T[i][j]` en el
  denominador; replicar el criterio exacto en Fase 3 para igualar nº de iteraciones.
- **`dt=1 s`** hace la corrida cara (86 400 pasos × varios días × iteración interna). Por eso
  Fase 7 explora `dt` mayor (el esquema es implícito) además de paralelizar.
- **Fuente de Tsa distinta** (C sinusoidal vs pvlib): la comparación final con el C es solo
  "del orden"; la validación estricta vive en las Fases 1–4 con Tsa idéntica.
- **Promedios de superficie:** `Tsi/(nx-1)` es correcto (medios nodos, ver hallazgo de
  revisión). ✅ Estandarizado a `/(nx-1)` en producción (Fase 5 rellena, 8a hueco), conservando
  la convención fiel (`Tso /nx`) solo en los drivers de validación. Pendiente menor: revisar
  unidades de los índices.
