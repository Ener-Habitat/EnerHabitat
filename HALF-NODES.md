# HALF-NODES: nodos en superficie con volúmenes de control fraccionarios

Especificación para la rama `half-nodes`. Migra el 1D y el 2D del esquema
actual (celdas completas centradas, película aplicada a la celda de frontera)
al esquema de nodos sobre las superficies con volúmenes fraccionarios. Regla
central (por tipo de nodo, con `w` el factor de volumen respecto a `ΔxΔy`):

| Nodo | Volumen `w` | Caras especiales |
|---|---|---|
| interior | 1 | — |
| borde del dominio (superficie exterior/interior, costado adiabático) | 1/2 | película sobre `w·Δ` |
| esquina del dominio | 1/4 | dos caras de `Δ/2` |
| pared de cavidad (NT 9–12) | 1/2 | convección+radiación con `T_h` sobre su segmento de cara |
| **esquina re-entrante de cavidad (NT 1–4)** | **3/4** (la cavidad ocupa un cuadrante) | **caras hacia la cavidad de `Δx/2` y `Δy/2`** |

Los tipos `NT` existentes ya clasifican todos estos casos (herencia del C, que
era de nodos-en-vértice); la implementación consiste en asignar a cada tipo su
`w` y sus áreas, no en inventar estructura nueva.

## Estado de partida (verificado en el código)

- 1D: `dx = L/Nx`, celdas completas; masa `ρc·dx/dt` uniforme; película sumada
  a la celda de frontera. Error de frontera de primer orden ~0.3 % con los
  defaults (validación Dirichlet-analítica: 0.58 % a `Nx = 200`, 0.13 % a 800).
- 2D: `dx = X/nx`, `dy = Y/ny`; `apo = ρc·dx·dy/Δt` uniforme para TODOS los
  nodos, incluidos corners y paredes de cavidad; el reporte `Tso`/`Tsi` usa la
  convención `Σ/(nx−1)` del C, sesgada `+T/(nx−1)` en la malla actual.
- Los kernels legacy (golden C) comparten parte del ensamblaje: hay que
  auditar y aislar (gate `wx=None → uniforme`, o duplicar) para que los golden
  no cambien ni un bit.

## Fase 0 — línea base de aceptación

- Congelar números de referencia: validación Dirichlet-analítica 1D (script de
  matrices de transferencia), los tres casos pre-computados de la doc, y la
  suite completa en verde.
- Criterios que gobiernan TODAS las fases: (i) el error Dirichlet 1D debe
  BAJAR a igual `Nx`; (ii) la reducción 2D→1D homogénea se mantiene (ambos
  esquemas cambian consistentemente); (iii) cierre de energía
  `energy_imbalance ≈ 0`; (iv) golden legacy intactos; (v) conservación exacta
  de masa térmica a cualquier malla.

## Fase 1 — 1D (pequeña, autocontenida)

- `dx = L/(Nx−1)`; nodos 0 y `Nx−1` sobre las superficies con `w = 1/2`
  (`mass_coeff[extremos] *= 0.5`); `Gf` integra tramos nodo-a-nodo.
- El mapeo material→celda recorta los volúmenes extremos a `[0, dx/2]` y
  `[L−dx/2, L]`; la conservación de `Σ ρc_j·L_j` sigue exacta por promedio
  ponderado.
- `T[0]`/`T[Nx−1]` pasan a ser temperaturas de superficie exactas (mejora
  `Tso`/`Tsi`, el criterio de `hi_flow` y el acoplamiento del aire).
- Aceptación: error Dirichlet < 0.58 % a `Nx = 200`; suite 1D.

## Fase 2 — 2D, fronteras externas

- HALLAZGO (fase 1, verificado en código): los constructores de malla
  (`compute_mesh`, `compute_mesh_slab`, `draw_*`, `set_krhoc_*`) son
  COMPARTIDOS con los caminos legacy de los golden tests. El gate de
  compatibilidad debe estar en la CONSTRUCCIÓN DE MALLA (parámetro
  `vertex=True/False` o funciones paralelas), no solo en los kernels:
  `dx = X/(nx−1)` cambia índices de cavidad (`i1..j2`), NT y campos k/ρc.
  Auditar primero qué constructores usan los golden (tests/golden,
  test_eh2d_hueca, test_eh2d_step, test_eh2d_coeffs) y aislarlos.
- Confirmado también: los NT 1–4 son esquinas DE DOMINIO y 6–7 los costados
  (la clasificación que esta fase necesita ya existe desde el C).
- `Δx = W/(nx−1)`, `Δy = L/(ny−1)`; arreglos de pesos `wx[nx]`, `wy[ny]`
  (1/2 en los extremos) pasados a los kernels.
- Masa `ρc·(wxΔx)(wyΔy)/Δt`; áreas de cara N/S = `wxΔx`, E/W = `wyΔy`;
  películas `h·(wxΔx)` en `j = 0, ny−1`; flujo al aire interior
  `Σ h·(wxΔx)(T−Ti)` con `Σ wxΔx = W` exacto.
- El *snapping* de la geometría (cavidades, interfaces del rib/topping) pasa a
  líneas de nodo: `i = round(pos/Δx)`.
- Reporte `Tso`/`Tsi` como promedio ponderado `Σ wT/(nx−1)` (desaparece el
  sesgo del C); el criterio de `hi_flow` usa el mismo promedio.
- Aceptación: reducción 2D→1D, energía, golden intactos.

## Fase 3 — cavidad con nodos en sus paredes (la parte delicada)

- Paredes de cavidad (NT 9–12): `w = 1/2`; intercambio convectivo+radiativo
  con `T_h` a través de su segmento de cara (segmentos de los extremos a
  `Δ/2`, de modo que `Σ segmentos = perímetro` exacto).
- **Esquinas re-entrantes (NT 1–4): `w = 3/4`; caras hacia la cavidad de
  `Δx/2` y `Δy/2`** — la consideración que motiva esta fase.
- Radiosidad: las temperaturas medias de cada pared (entrada de los factores
  de transferencia y del Nusselt) pasan a promedios ponderados con extremos a
  mitad de peso.
- Nodo `T_h`: la integral de perímetro usa los mismos segmentos ponderados.
- Aplica a `_step_hueca` y `_step_slab` (por cavidad) y a sus 4 drivers; los
  kernels legacy quedan aislados con el gate.
- Aceptación: golden intactos; hueca/slab full-day; comparación de energías
  contra fase 2 (corrimientos pequeños y documentados); Borbón (hot-box) como
  referencia física de que no nos alejamos del experimento.

## Fase 4 — documentación y regeneración

- `numerics.qmd`: la ecuación discreta gana los pesos
  (`a_P = ρc·w_xΔx·w_yΔy/Δt + …`), la sección de mapeo capa→malla y las
  definiciones de `Δ`; nota de que las `Tso`/`Tsi` reportadas son ahora
  temperaturas de superficie exactas (se retira la nota de convención C en
  `api.qmd`).
- Regenerar los tres pre-cómputos (~25 min) y re-renderizar; actualizar los
  números citados en las páginas.

## Estimación y riesgo

- Fase 1: una sesión. Fase 2: 1–2 sesiones (mecánica, muchos puntos de
  contacto). Fase 3: 2–3 sesiones (la única con física de ensamblaje nueva) +
  regeneración. Total: ~una semana de sesiones cuidadosas.
- Riesgo principal: romper la fidelidad de los golden legacy o introducir un
  error de área/volumen silencioso. Mitigación: gate de pesos, aceptación por
  fase, y el criterio (i) — el error Dirichlet analítico — que detecta
  cualquier inconsistencia de frontera de inmediato.
- Las fases 1–2 son valiosas por sí solas (superficies exactas, sesgo de
  reporte eliminado) y pueden mergearse sin la 3; la 3 puede diferirse si sus
  corrimientos exigen re-validar contra experimento con más calma.
