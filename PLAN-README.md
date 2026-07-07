# PLAN-README — Documentación del paquete (física + numérica + uso)

Plan para llevar la documentación de EnerHabitat al nivel de un paquete científico:
**problema físico, ecuaciones con condiciones de frontera (1D y 2D), método numérico y uso del
paquete**, integrando el paper de referencia
([`pdfs/eh_paper.pdf`](pdfs/eh_paper.pdf) — Barrios et al., *Solar Energy* 131 (2016) 296–304,
[doi:10.1016/j.solener.2015.12.017](http://dx.doi.org/10.1016/j.solener.2015.12.017)).

## Diagnóstico (dónde estamos)

El [`README.md`](README.md) actual ya cubre bien el **uso** (instalación, quickstart, workflow,
API 1D y 2D, `config`, materiales). Lo que falta es la **base científica**:

- La sección *Theoretical background* solo tiene la ecuación de calor 1D y `Tsa` en ASCII-art,
  **sin condiciones de frontera**, sin el balance del aire interior, sin el modelo 2D
  (cavidades: radiación + convección), sin el método numérico (volúmenes de control, TDMA,
  media armónica, criterios de convergencia) y sin validación.
- El paper describe la herramienta en línea (EH v2.2, 2016); el paquete **difiere en varios
  puntos** (clima por EPW, `dt`, coeficientes, tolerancias, geometrías 2D más generales).
  Documentar el paquete fielmente exige una tabla de divergencias paper ↔ paquete.
- No hay forma canónica de citar el paquete (sin `CITATION.cff`, sin sección *How to cite*).

## Decisiones propuestas (para tu revisión) 🔍

1. **Estructura: README + sitio Quarto en GitHub Pages, mismo repo** — el README se queda
   como puerta de entrada centrada en **uso** (con un *Theoretical background* compacto:
   ecuaciones clave + enlaces al sitio), y la teoría completa vive en tres documentos
   Quarto en `docs/`, que además son las páginas del sitio:
   - `docs/model-1d.qmd` — problema físico y modelo 1D
   - `docs/model-2d.qmd` — modelo 2D (bloque hueco y vigueta-bovedilla)
   - `docs/numerics.qmd` — método numérico, convergencia y validación

   El sitio se publica en **`https://ener-habitat.github.io/EnerHabitat/`** (GitHub Pages,
   gratis — mismo esquema que `python-windrose.github.io/windrose`). Se elige **Quarto**
   (el autor ya lo domina) por sus ventajas para docs científicas: **ecuaciones numeradas
   con referencias cruzadas** (`$$...$$ {#eq-heat}` → `@eq-heat`), **bibliografía BibTeX**
   (citar el paper, Patankar, Xamán, Hollands, Humphreys & Nicol…) y **chunks de Python
   ejecutables** — los ejemplos corren al renderizar, así que el build verifica los
   snippets. Ver Fase 5.

   **Mismo repo, sin repo aparte**: los fuentes (`docs/*.qmd`, `_quarto.yml`,
   `references.bib`) viven en `main` junto al paquete; el HTML compilado vive solo en la
   rama **`gh-pages`**, que `quarto publish` crea y mantiene. GitHub Pages sirve esa rama.
   El wheel de PyPI solo empaqueta `src/enerhabitat`, así que la doc no engorda el paquete.

   *Alternativa descartada:* meter todo al README (quedaría >1200 líneas y PyPI lo muestra
   completo como landing page). *Alternativas descartadas:* MkDocs Material / Sphinx —
   herramientas nuevas para el autor, sin numeración de ecuaciones ni BibTeX nativos.

   **Compatibilidad con PyPI** (el paquete ya está publicado): PyPI renderiza **solo** el
   README (`long_description`) — no hospeda subpáginas. Los enlaces del README hacia la
   documentación deben ser **URLs absolutas** (al sitio de GitHub Pages) — los relativos se
   rompen en PyPI (las imágenes ya usan `raw.githubusercontent.com`, mismo patrón). Además,
   `Documentation = "https://ener-habitat.github.io/EnerHabitat/"` en `[project.urls]` pone
   el enlace en el sidebar de PyPI (eso son las "subpáginas" que se ven en otros proyectos).
2. **Idioma: inglés** — consistente con el README, el código y PyPI (el paper y la audiencia
   internacional también). Los planes internos siguen en español.
3. **Matemáticas: LaTeX** — en el sitio Quarto con **numeración y referencias cruzadas**
   (`$$...$$ {#eq-heat}` → `@eq-heat`); en el README, LaTeX nativo de GitHub
   (`$...$`/`$$...$$`) solo para las 2–3 ecuaciones esenciales, reemplazando el ASCII-art.
   Nota: PyPI **no** renderiza math en el README, razón de más para que las ecuaciones
   densas vivan en el sitio.
4. **Citabilidad: `CITATION.cff`** en la raíz + sección *How to cite* en el README (paper de
   *Solar Energy* + el paquete con su versión).

## Estructura del sitio Quarto (cómo se construye)

Todo vive en **este mismo repo**. La rama `main` guarda los fuentes; la rama `gh-pages`
guarda solo el HTML compilado y es la que sirve GitHub Pages.

### Árbol de archivos (rama `main`)

```
EnerHabitat/
├── src/enerhabitat/            ← el paquete (lo único que empaqueta el wheel de PyPI)
├── README.md                   ← landing de GitHub y PyPI: uso + resumen teórico + enlaces al sitio
├── CITATION.cff
└── docs/                       ← proyecto Quarto (fuentes del sitio)
    ├── _quarto.yml             ← configuración: website, navbar, tema, crossref, bibliografía
    ├── index.qmd               ← Home: overview + instalación + quickstart
    ├── usage.qmd               ← Usage: workflow, 4 ejemplos (1D/2D × libre/AC), config, materiales
    ├── model-1d.qmd            ← Theory: problema físico y modelo 1D           (Fase 1)
    ├── model-2d.qmd            ← Theory: modelo 2D, cavidades                  (Fase 2)
    ├── numerics.qmd            ← Theory: método numérico, convergencia, validación (Fase 3)
    ├── api.qmd                 ← API: Location, System, System2D, HollowBlock, Slab, config
    ├── about.qmd               ← About: changelog, how to cite, autores, licencia
    ├── references.bib          ← bibliografía BibTeX (se cita con @barrios2016, @patankar1980…)
    ├── img/                    ← figuras (hollow_block.png, slab.png, dominio 1D…)
    ├── make_figures.py         ← genera las figuras (script, no es página del sitio)
    ├── data/                   ← materials.ini + EPW pequeño para los chunks ejecutables
    ├── _freeze/                ← resultados de los chunks (SÍ se versiona: CI no re-ejecuta)
    ├── _site/                  ← HTML generado (gitignored — solo vive en gh-pages)
    └── .quarto/                ← caché de Quarto (gitignored)
```

### `_quarto.yml` (esqueleto)

```yaml
project:
  type: website
  output-dir: _site

website:
  title: "EnerHabitat"
  navbar:
    left:
      - {href: index.qmd,  text: Home}
      - {href: usage.qmd,  text: Usage}
      - text: Theory
        menu:
          - {href: model-1d.qmd, text: "1D model"}
          - {href: model-2d.qmd, text: "2D model"}
          - {href: numerics.qmd, text: "Numerical method"}
      - {href: api.qmd,    text: API}
      - {href: about.qmd,  text: About}

format:
  html:
    theme: cosmo        # o flatly; decidir al montar
    toc: true

bibliography: references.bib

execute:
  freeze: auto          # los chunks solo se re-ejecutan si su .qmd cambia
```

### Flujo de construcción y publicación

```
 rama main                    build                         publicación
┌──────────────┐   quarto render docs   ┌───────────┐   push automático   ┌───────────────┐
│ docs/*.qmd   │ ─────────────────────► │ docs/_site│ ──────────────────► │ rama gh-pages │
│ _quarto.yml  │  (ejecuta los chunks,  │  (HTML)   │  (quarto publish /  │  (solo HTML)  │
│ references...│   numera ecuaciones,   └───────────┘   GitHub Action)    └───────┬───────┘
└──────────────┘   resuelve citas)                                                │
                                                                                  ▼
      trabajo local: `quarto preview docs`             GitHub Pages sirve
      primer deploy:  `quarto publish gh-pages`        https://ener-habitat.github.io/EnerHabitat/
```

- **Local**: se escribe en `docs/*.qmd` y se revisa con `quarto preview docs`.
- **Publicación (manual, decisión de Guillermo)**: `quarto publish gh-pages` desde
  `docs/` tras cada cambio de documentación. **No hay CI**: se evaluó un workflow de
  GitHub Actions y se descartó — publica solo el mantenedor, a mano. Regla práctica:
  commit de fuentes a `main` + `quarto publish`, para que el sitio no se desfase del
  repo.
- El README enlaza a las páginas del sitio con URL absoluta
  (`https://ener-habitat.github.io/EnerHabitat/model-1d.html`, etc.), que funcionan igual
  desde GitHub y desde PyPI.

## Fuente de cada contenido

| Contenido | Fuente |
|---|---|
| Ecuación de calor 1D, continuidad entre capas, BCs exterior/interior | Paper §2.1, Eqs. (1)–(5) |
| Temperatura sol-aire `Tsa`, absortancia, `RF` | Paper Eq. (4); código `ehframe.py` (`RF=-3.9` horizontal / `0` vertical) |
| Aire interior: `Tn` (AC) y balance libre | Paper Eqs. (6)–(7); código `ehtools.py:67` (`Tn = 13.5 + 0.54·T̄a`), `La = 2.5 m` |
| Amplitud de la zona de confort (setpoint AC) | Morillón (referencia a añadir; código `ehframe.py`) |
| Modelo 2D: extensión de la ec. de calor, BCs, simetría (∂T/∂y=0) | Paper §2.2; `PLAN-2D-hecho.md` §Tipos de nodo |
| Física de cavidades: radiación (gris-difusa, factores de vista) + Nusselt | Paper §2.2; código `ehtools2d.py:352-382` (muro), `:705-725` (techo, Rayleigh–Bénard) |
| Geometrías 2D: bloque hueco (muros), vigueta-bovedilla con N cavidades y vigueta en L (techos) | `PLAN-2D-hecho.md` Fases 8a/8b; figuras `docs/img/*.png` |
| Método numérico: volúmenes de control implícitos, TDMA, media armónica | Paper §2.3; Patankar (1980); `ehtools.py` |
| 2D: línea-por-línea (line-TDMA), tolerancias, convergencia día-a-día | Paper Eqs. (8)–(9); `ehtools2d.py`, `config2d` |
| Construcción del día promedio (EPW, pvlib, Chow & Levermore) | Paper §2.4 (adaptado: EPW en vez de Meteonorm); `ehframe.py` |
| Validación | Paper §3 (vs EnergyPlus) y §4 (bloque hueco vs experimento, dif. máx. 6.5%); golden masters del C legacy (`PLAN-2D-hecho.md` Fase 0); reducción 2D→1D |

## Tabla de divergencias paper ↔ paquete (a incluir en `docs/numerics.md`)

Puntos donde el paquete **no** es el EH del paper — documentar el paquete, señalando el cambio:

| Aspecto | Paper (EH online 2016) | Paquete |
|---|---|---|
| Clima | BD Meteonorm de 80 ciudades MX | **Archivo EPW** del usuario + pvlib (cualquier sitio) |
| Irradiancia del día promedio | Seno positivo (amplitud = máx. del mes); difusa lineal con la inclinación; geometría solar de Duffie & Beckman | **Medias horarias del mes** interpoladas; `Is` en superficie inclinada por **transposición pvlib** |
| `RF` (radiación de onda larga) | 3.9 °C horizontal → 0 vertical, **lineal con la inclinación** | **Binario**: 3.9 °C solo si `tilt==0`, si no 0 |
| Setpoint del AC | `Tn` (neutralidad, Humphreys & Nicol) | Igual (`Tn` constante); `DeltaTn` (Morillón) solo como columna de datos |
| Paso de tiempo | 1 s | **10 s, fijo** (estabilidad del nodo de aire explícito) |
| Malla | 1D `Nx=100`; 2D 160×160 | 1D `Nx=200`; 2D `nx=80 × ny=160` (`config2d`) |
| `hi` | 8.1 (muro), 9.4/6.6 (techo s/ dirección) | **8.6** W/m²K (NOM-020/008-ENER), configurable |
| Convergencia periódica | `Cs = 1e-5` °C | `tol_day = 5e-4` (y `tol_inner = 1e-10`, `max_days=60`) |
| Geometría 2D | bloque relleno / bloque hueco (Fig. 1) | bloque hueco **muros** + vigueta-bovedilla **techos** (N cavidades, vigueta en L), relleno `AIR`/`SOLID` |
| Materiales | BD integrada + BD de usuario | **solo** `materials.ini` del usuario (sin BD integrada) |

> Verificar cada valor contra el código al redactar (no confiar en esta tabla ni en el paper
> de memoria): `config.py`, `ehframe.py`, `ehtools.py`, `ehtools2d.py`.

## Fases

### Fase 0 — Esqueleto del sitio publicado, para revisión ✅ (render/publish: Guillermo)

Deja la estructura completa del sitio como páginas *stub* (título + propósito + secciones
vacías marcadas 🚧), de modo que se pueda **publicar y revisar la navegación** antes de
escribir contenido. Sin chunks ejecutables todavía (el render no requiere `enerhabitat`).

- [x] `docs/_quarto.yml` — website, navbar (Home/Usage/Theory/API/About), tema, `toc`,
      `bibliography`, `freeze`.
- [x] Stubs: `index.qmd`, `usage.qmd`, `model-1d.qmd`, `model-2d.qmd`, `numerics.qmd`,
      `api.qmd`, `about.qmd` — cada uno con el outline de secciones de su fase.
- [x] `docs/references.bib` — semilla con las referencias clave (se completa en Fase 5).
- [x] `.gitignore` — añadir `docs/_site/` y `docs/.quarto/`.
- [x] **Guillermo**: `quarto preview docs` para revisar local, y primer
      `quarto publish gh-pages` (desde `docs/`) + *Settings → Pages* si hace falta.
      **Publicado**: <https://ener-habitat.github.io/EnerHabitat/> ✅ (2026-07-06)
- [ ] Revisión conjunta de la estructura publicada → ajustes de navegación antes de Fase 1.

### Fase 1 — `docs/model-1d.qmd`: problema físico y modelo 1D ✅ (revisión: Guillermo)

El documento ancla: qué problema se resuelve y con qué modelo.

- [x] **Planteamiento físico**: componente único de la envolvente (muro/techo opaco, sin
      ventanas/ventilación/infiltración ni cargas térmicas al interior), día promedio del mes, régimen periódico.
      Esquema del dominio: capas 1..N, `x=0` exterior → `x=L` interior, aire interior a `La`.
- [x] **Ecuación de calor** por capa (Eq. 1) + continuidad de flujo en las juntas (Eq. 2).
- [x] **Condición de frontera exterior** (Eq. 3) con **temperatura sol-aire** (Eq. 4):
      absortancia, `ho`, `RF` — **verificado en código**: `Tsa = Ta + a·Is/ho − LWR` con
      `LWR = 3.9` solo si `tilt==0`, si no `0` (binario, sin interpolación con la
      inclinación); `Is` con pvlib (`get_total_irradiance`).
- [x] **Condición de frontera interior** (Eq. 5) y los **dos modos**:
      *AC* — `Ti = Tn` constante, `Tn = 0.54·T̄a + 13.5` (Humphreys & Nicol); **verificado**:
      `DeltaTn` (amplitud de Morillón, tabla por rangos de `ΔTa`) se entrega como columna
      pero `solveAC` NO la usa; *libre* — balance del aire interior con `La = 2.5 m`,
      `ρa=1.18`, `ca=1005`.
- [x] **Salidas y unidades**: `Ti`, `energy_transfer`, `cooling/heating_energy` en
      J/(m²·día), con las integrales explícitas de cada energía.
- [x] **Figura**: dominio 1D con BCs etiquetadas (`domain_1d()` en `docs/make_figures.py`
      → `docs/img/domain_1d.png`).
- [x] **Sección "The average day"**: modelo coseno de Chow & Levermore para `Ta` (mín. al
      amanecer, máx. a la hora media del máximo del mes), irradiancias = medias horarias
      del mes interpoladas, `Tn`/`DeltaTn`.
- [ ] **Guillermo**: revisar y publicar (`quarto publish gh-pages`). Verificar la
      referencia de Morillón usada (`morillon2004`: *Human bioclimatic atlas for Mexico*,
      Solar Energy 76, 781–792) — confirma que es la fuente correcta de la tabla de
      `DeltaTn`.

### Fase 2 — `docs/model-2d.qmd`: sistemas no homogéneos ✅ (revisión: Guillermo)

- [x] **Por qué 2D**: capa heterogénea a lo ancho; celda repetitiva con **simetría lateral**
      adiabática. Convención de ejes declarada explícitamente: `x` = ancho, `y` = espesor
      (la del paquete), con nota de que el paper usa la transpuesta.
- [x] **Ecuación de calor 2D** con propiedades por posición + BCs exterior/interior iguales
      al 1D y lados adiabáticos (`∂T/∂x = 0`), en un solo bloque de ecuaciones.
- [x] **Física de la cavidad de aire** (`Fill.AIR`), todo verificado en `ehtools2d.py`:
      (i) conducción; (ii) convección — **muro**: `hc = 0.4005·|ΔT|^0.3033/d^0.0901`
      (Xamán); **techo**: Rayleigh–Bénard con Hollands
      `Nu = 1 + 1.44[1−1708/Ra]⁺ + [(Ra/5830)^⅓−1]⁺`, estable → `kair/d`
      (`_slab_hh`, `ehtools2d.py:719-735`); (iii) **radiación** gris-difusa con emisividad
      `E` y factores de vista por cuerdas cruzadas (`_view_factors`, `ehtools2d.py:479`);
      (iv) nodo de cavidad agrupado `Th` (Euler explícito, `Ch = ρa·ca·w·d`).
- [x] **Relleno sólido** (`Fill.SOLID`): solo conducción; caso límite = 1D (enlazado a la
      validación).
- [x] **Geometrías**: `HollowBlock` y `Slab` con las figuras existentes y las fórmulas de
      ancho/espesor.
- [x] **Restricciones y salidas**: exactamente un elemento 2D, `tilt` obligado por tipo,
      espesor derivado de la geometría, resultados por unidad de área interior.
- [ ] **Guillermo**: revisar y publicar (`quarto publish gh-pages`).

### Fase 3 — `docs/numerics.qmd`: método numérico, convergencia y validación ✅ (revisión: Guillermo)

- [x] **Discretización**: volúmenes de control implícitos (Patankar 1980), ecuación
      discreta general, media armónica en las caras (`prepare_static_coefficients`,
      `ehtools.py:220`), malla (`Nx=200` 1D; `80×160` 2D).
- [x] **Solución**: TDMA directo y exacto por paso (1D); línea-por-línea con vecinos `y`
      retardados (Gauss–Seidel) + TDMA en `x` por fila, iterado a `tol_inner=1e-10`, con
      radiación/Nusselt re-evaluados **en cada barrido** (2D, `_step_hueca`/`_step_slab`).
- [x] **Nodo de aire interior (acople, modo libre)**: conducción implícita sin criterio;
      acople segregado con paso explícito del aire, `Fo < 2` / `Fo ≲ 1`, `Fo ≈ 0.03` con
      `dt=10 s` — razón del `dt` fijo. Nodos de cavidad `Th` con el mismo tratamiento.
      En AC no hay criterio (implícito puro).
- [x] **Convergencia al régimen periódico**: criterio día-a-día, `Cs = 5e-4` °C
      (**hardcoded en 1D**, `config2d.tol_day` en 2D), `max_days=60`, arranque en `T̄n`;
      no convergencia detectable vía `System2D.days == max_days` (documentado; el solver
      no emite warning).
- [x] **Tabla de divergencias paper ↔ paquete** (10 filas, cada una verificada en código).
- [x] **Validación**: 1D vs EnergyPlus (0.5 °C, DF 0.03, LT 0.4 h, 3–4% energías);
      2D relleno = 1D; 2D hueco vs hot-box (6.5%, Borbón et al. 2010 — añadido a la
      bibliografía); golden masters del C en `tests/golden/` + suite `test_eh2d_*.py`.
- [x] **Desempeño**: serial JIT numba (costo de compilación en la 1.ª llamada); `prange`
      medido ~1.06× y removido; procesos ~6× a 8 / ~10× a 16.
- [ ] **Guillermo**: revisar y publicar (`quarto publish gh-pages`).

### Fase 4 — README: reestructurar y enlazar ✅ (revisión: Guillermo)

**Propósito del README** (define qué se queda y qué se va al sitio): es la **landing page
doble** — la portada del repo en GitHub y la *única* página que PyPI renderiza. Su lector es
alguien decidiendo si el paquete le sirve y queriendo su primer resultado en minutos. Debe
responder: *qué es, cómo se instala, cómo corro el primer ejemplo, dónde está el detalle* —
y nada más. La profundidad (teoría completa, ejemplos extensos, API exhaustiva) vive en el
sitio Quarto; el README enlaza. Restricción técnica que refuerza el reparto: PyPI no
renderiza LaTeX ni resuelve enlaces relativos.

- [x] Reescribir *Theoretical background* como **resumen** (~30 líneas): ecuación de calor,
      `Tsa` y `Tn` en LaTeX, los dos modos de solución, un párrafo sobre el 2D, y enlaces a
      las páginas *Theory* del sitio. Sin ASCII-art.
- [x] **Adelgazar hacia el sitio**: los ejemplos extensos, la referencia API detallada y las
      notas técnicas (`dt` fijo, paralelización) pasan a `usage.qmd`/`api.qmd`/`numerics.qmd`;
      en el README queda quickstart (1D + 2D compacto) + tabla mínima de API + notas de una
      línea con enlace. **Regla anti-deriva**: cada contenido vive en un solo lugar; la única
      duplicación permitida es el quickstart (README ↔ `index.qmd`).
- [x] **Escribir `docs/usage.qmd`** (destino de la migración): workflow, estructura de
      carpetas, `materials.ini` completo, los 5 ejemplos (matriz 1D/2D × libre/AC + `Slab`),
      `Fill.SOLID`, inspector, `config` y `config2d` — por ahora como **listados estáticos**
      (los chunks ejecutables son de la Fase 5).
- [x] **Escribir `docs/api.qmd` — referencia de API completa** (manual): `Location`,
      `System`, `System2D` (incl. `setpoint` y `days`, verificados en `eh2d.py`),
      `HollowBlock` y `Slab` (tablas de `geometry`), `Fill`, inspector, `config`
      (incl. `AIR_DENSITY`/`AIR_HEAT_CAPACITY`) y `config2d`, con métodos verificados
      contra `config.py`. (`quartodoc` autodoc queda opcional/futuro.)
- [x] Todos los enlaces del README a la documentación con **URL absoluta al sitio**
      (`https://ener-habitat.github.io/EnerHabitat/...`) + badge *docs* junto a los de PyPI.
- [x] **Corregir imprecisiones detectadas en Fase 1**: setpoint AC = `Tn` a secas
      (`DeltaTn` solo columna de datos, mencionada como tal) y convención
      `Tsa = Ta + a·Is/ho − RF` con `RF=3.9` en techos.
- [x] Sección **How to cite** (cita del paper + BibTeX) enlazando `CITATION.cff`
      (creado — ver Fase 6).
- [ ] **Guillermo**: revisar el nuevo README (GitHub lo renderiza con LaTeX; en PyPI las
      ecuaciones se verán como texto — decidir si aceptable o quitar los `$$`) y publicar
      el sitio actualizado (`quarto publish gh-pages`).

### Fase 5 — Sitio de documentación: Quarto + GitHub Pages (mismo repo) ✅ (revisión: Guillermo)

Publica `docs/` como sitio Quarto en `https://ener-habitat.github.io/EnerHabitat/` (estilo
windrose). Los `.qmd` de las Fases 1–3 son las páginas.

- [x] **Proyecto Quarto en `docs/`** (hecho en Fase 0): `_quarto.yml` con navbar, tema
      `cosmo`, `toc`, `bibliography`, `freeze`. Sitio publicado en
      <https://ener-habitat.github.io/EnerHabitat/>.
- [x] **`docs/references.bib`**: Barrios 2016, Patankar 1980, Xamán 2005, Hollands 1975,
      Humphreys & Nicol 2000, Chow & Levermore 2007, Morillón 2004, ASHRAE 1997,
      Incropera 2007, Holmgren 2018 (pvlib), Borbón 2010. (Añadir Barrios 2011/2012 y
      NOM-008/020 si alguna página los cita al final.)
- [x] **Ejemplos 1D como chunks ejecutables** — 1D libre y 1D AC corren al renderizar,
      con gráfica (`Ti`/`Ta`/`Tsa` + banda de confort `Tn ± ΔTn`), paleta validada
      (CVD/contraste) y etiquetas al final de línea con anti-colisión. Datos reales:
      mayo, Cuernavaca.
- [x] **Ejemplos 2D pre-computados** — medido (malla por defecto, mayo/Cuernavaca):
      bloque hueco libre **21 min** (4 días), bloque hueco AC **10 min** (2 días),
      `Slab` libre **15 min** (5 días) → chunks vivos inviables (`freeze` es POR
      DOCUMENTO: cualquier edición de `usage.qmd` re-ejecutaría todo).
      **Implementado**: `docs/run_examples.py` (versionado, 3 casos en paralelo,
      ~21 min total) → CSVs + `summary.json` en `docs/data/results/`; los chunks de
      `usage.qmd` solo leen y grafican, con nota de cómo regenerar. La página muestra
      el listado "como lo escribe el usuario" + el resultado pre-computado (energías,
      días de convergencia, tiempo de cómputo y gráfica).
- [x] **Datos para los chunks**: `docs/data/` con `materials.ini` (todos los materiales
      de la doc) y el EPW de Cuernavaca (1.5 MB, copiado de `tests/`). `docs/_freeze/`
      generado — **debe versionarse**.
- [x] **Entorno de render**: grupo de dependencias `docs` (PEP 735: `jupyter`,
      `matplotlib`) vía `uv add --group docs` — no afecta al paquete publicado.
      Renderizar con `uv run quarto render docs`.
- [x] **Primer deploy** con `quarto publish gh-pages` + *Settings → Pages* (Fase 0).
- [x] ~~Deploy automático (CI)~~ — **descartado por decisión** (2026-07-06): la
      publicación es **manual** con `quarto publish gh-pages`; el workflow que se había
      escrito se eliminó antes de llegar al remoto.
- [x] **`.gitignore`**: `docs/_site/` y `docs/.quarto/` (Fase 0).
- [x] **`index.qmd` y `about.qmd`**: Home (overview + install + quickstart) y About
      (how to cite + changelog + autores + licencia) llenados.
- [ ] **Guillermo**: revisar, commit + push de fuentes (incluir `docs/_freeze/` y
      `docs/data/`) y `quarto publish gh-pages` cuando los ejemplos 2D estén integrados.

### Fase 6 — Figuras y citabilidad ⏳

- [x] `docs/make_figures.py`: figura del **dominio 1D con BCs** añadida (Fase 1).
      (Opcional pendiente: celda 2D con BCs/simetría anotadas sobre las figuras existentes.)
- [x] **`CITATION.cff`** (paper como `preferred-citation`, autores, DOI, versión 0.2.0,
      URL del sitio) — creado en Fase 4 para respaldar la sección *How to cite* del README.
- [ ] **`pyproject.toml`**: añadir `Documentation = "https://ener-habitat.github.io/EnerHabitat/"`
      a `[project.urls]` → enlace *Documentation* en el sidebar de PyPI (visible en el
      siguiente release).
- [ ] `CHANGELOG.md`: entrada de documentación; evaluar aquí el **bump de versión** pendiente
      de [`PLAN-2D.md`](PLAN-2D.md) (0.2.x con docs, o dejarlo para el siguiente release).

### Fase 7 — Verificación cruzada (control de calidad) ⏳

- [ ] **Ecuaciones vs código**: cada constante y ecuación documentada se contrasta con el
      código (`RF`, `Tn`, `ho/hi`, `dt`, `La`, tolerancias, correlaciones de Nusselt,
      emisividad). Ninguna cifra entra a docs solo por venir del paper.
- [ ] **Snippets ejecutables**: los del sitio los verifica el propio `quarto render`
      (chunks ejecutables — si un ejemplo se rompe, el build falla); los del README (pocos,
      quickstart) se corren a mano con el `materials.ini` de ejemplo.
- [ ] **Render**: revisar el sitio con `quarto preview` (ecuaciones numeradas, citas,
      figuras) y que el README siga legible en GitHub **y** en PyPI (donde el math no
      renderiza — por eso lleva solo las ecuaciones esenciales).
- [ ] Enlaces internos (README ↔ docs ↔ figuras) sin rotos.

## Criterio de terminado

Un usuario nuevo puede, sin abrir el paper ni el código: (1) entender qué problema físico se
resuelve y con qué ecuaciones/BCs en 1D y 2D, (2) saber cómo se resuelve numéricamente y con
qué criterios de convergencia y validación, (3) correr una simulación 1D y 2D copiando los
snippets, y (4) citar el paquete y el paper correctamente.
