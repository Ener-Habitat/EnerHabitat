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
├── .github/workflows/docs.yml  ← CI: re-publica el sitio en cada push a main que toque docs/
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
- **Primer deploy (manual, una vez)**: `quarto publish gh-pages` desde `docs/` — crea la
  rama `gh-pages` y configura Pages.
- **Después (automático)**: el workflow `docs.yml` re-renderiza y publica en cada push a
  `main` que toque `docs/`. Gracias a `freeze`, el CI no necesita ejecutar Python.
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
- [ ] **Guillermo**: `quarto preview docs` para revisar local, y primer
      `quarto publish gh-pages` (desde `docs/`) + *Settings → Pages* si hace falta.
- [ ] Revisión conjunta de la estructura publicada → ajustes de navegación antes de Fase 1.

### Fase 1 — `docs/model-1d.qmd`: problema físico y modelo 1D ⏳

El documento ancla: qué problema se resuelve y con qué modelo.

- [ ] **Planteamiento físico**: componente único de la envolvente (muro/techo opaco, sin
      ventanas/ventilación/infiltración ni cargas térmicas al interior), día promedio del mes, régimen periódico.
      Esquema del dominio: capas 1..N, `x=0` exterior → `x=L` interior, aire interior a `La`.
- [ ] **Ecuación de calor** por capa (Eq. 1) + continuidad de flujo en las juntas (Eq. 2).
- [ ] **Condición de frontera exterior** (Eq. 3) con **temperatura sol-aire** (Eq. 4):
      absortancia, `ho`, `RF` (−3.9 °C horizontal, 0 vertical, interpolación con la
      inclinación), y cómo se construye `Is` sobre superficie inclinada (pvlib).
- [ ] **Condición de frontera interior** (Eq. 5) y los **dos modos**:
      *aire acondicionado* — `Ti = Tn` constante, `Tn = 0.54·T̄o + 13.5` (Humphreys & Nicol) y
      amplitud de confort de Morillón; *libre* — balance del aire interior (Eq. 7) con
      `H = La = 2.5 m`.
- [ ] **Salidas y unidades**: `Ti`, `energy_transfer`, `cooling/heating_energy` en J/(m²·día)
      (mover/duplicar la nota de unidades que hoy está en el README).
- [ ] **Figura**: dominio 1D con BCs etiquetadas (añadir a `docs/make_figures.py`).

### Fase 2 — `docs/model-2d.qmd`: sistemas no homogéneos ⏳

- [ ] **Por qué 2D**: capa heterogénea a lo ancho; celda repetitiva con **simetría lateral**
      (adiabática, ∂T/∂y = 0 en el paper; en el paquete `x` = ancho, `y` = espesor — unificar
      la convención de ejes y decirlo explícitamente).
- [ ] **Ecuación de calor 2D** (extensión de Eq. 1) + BCs exterior/interior iguales al 1D.
- [ ] **Física de la cavidad de aire** (`Fill.AIR`): (i) conducción en los sólidos,
      (ii) convección natural — coeficiente dependiente de temperatura vía correlaciones de
      Nusselt: **muro** (cavidad vertical) y **techo** (Rayleigh–Bénard; estable → solo
      conducción), (iii) **radiación** entre las 4 superficies (grises-difusas, factores de
      vista, aire no participante), (iv) nodo de aire de cavidad agrupado (balance tipo Eq. 7
      con las 4 superficies). Citar Xamán et al. (2005) y Hollands et al. (1975) como el paper.
- [ ] **Relleno sólido** (`Fill.SOLID`): solo conducción; caso límite = capa homogénea 1D.
- [ ] **Geometrías**: `HollowBlock` (muro) y `Slab` (techo: 3 sólidos, vigueta en L,
      N cavidades) — reutilizar las figuras existentes `docs/img/hollow_block.png` y
      `docs/img/slab.png` y las fórmulas de ancho/espesor que ya están en el README.
- [ ] **Restricciones**: exactamente un elemento 2D en la pila, `tilt` obligado por tipo.

### Fase 3 — `docs/numerics.qmd`: método numérico, convergencia y validación ⏳

- [ ] **Discretización**: volúmenes de control implícitos (Patankar 1980), media armónica de
      conductividad en las caras, malla (`Nx` 1D; `nx×ny` 2D).
- [ ] **Solución**: TDMA directo (1D); línea-por-línea TDMA + iteración interna (2D) con
      `tol_inner` (Eq. 8) y el acople no lineal (radiación/Nusselt re-evaluados por iteración).
- [ ] **Nodo de aire interior (acople, modo libre)**: dejar clara la distinción — la
      conducción es **implícita** (incondicionalmente estable, sin criterio de `dt`), pero
      el acople muro↔aire está **segregado**: el TDMA usa la `Tint` vieja y luego el aire
      avanza **explícito** (forward-Euler, `ehtools.py:320`,
      `Tint += Fo·(T_surf − Tint)` con `Fo = hi·dt/(ρa·ca·La)`, `ehframe.py:601`). Ese
      paso escalar exige `Fo < 2` (estabilidad) y `Fo ≲ 1` (sin oscilación); con `dt=10 s`,
      `Fo ≈ 0.03` — es la razón por la que `dt` está fijo (mover la nota del README aquí y
      dejar allá una línea). En modo **AC** (`Tint` = setpoint constante) no hay criterio:
      el esquema es implícito puro.
- [ ] **Convergencia al régimen periódico**: criterio día-a-día (Eq. 9) con `tol_day` y
      `max_days`; qué pasa si no converge.
- [ ] **Tabla de divergencias paper ↔ paquete** (la de arriba, verificada contra código).
- [ ] **Validación** (resumen con números): 1D vs EnergyPlus (ΔT ≤ 0.5 °C, DF ≤ 0.03,
      LT ≤ 0.4 h; energías ≤ 3–4%); 2D relleno = 1D; 2D hueco vs experimento hot-box
      (dif. máx. 6.5%); paquete vs golden masters del C legacy y reducción 2D→1D
      (enlazar `tests/`).
- [ ] **Desempeño**: solver serial JIT (numba); paralelizar por **procesos** para barridos
      (~6× a 8 procesos) — condensar la nota que ya está en el README.

### Fase 4 — README: reestructurar y enlazar ⏳

**Propósito del README** (define qué se queda y qué se va al sitio): es la **landing page
doble** — la portada del repo en GitHub y la *única* página que PyPI renderiza. Su lector es
alguien decidiendo si el paquete le sirve y queriendo su primer resultado en minutos. Debe
responder: *qué es, cómo se instala, cómo corro el primer ejemplo, dónde está el detalle* —
y nada más. La profundidad (teoría completa, ejemplos extensos, API exhaustiva) vive en el
sitio Quarto; el README enlaza. Restricción técnica que refuerza el reparto: PyPI no
renderiza LaTeX ni resuelve enlaces relativos.

- [ ] Reescribir *Theoretical background* como **resumen** (~30 líneas): ecuación de calor y
      `Tsa` en LaTeX, los dos modos de solución, una línea sobre el 2D, y enlaces a las
      páginas *Theory* del sitio. Sin ASCII-art.
- [ ] **Adelgazar hacia el sitio**: los ejemplos extensos, la referencia API detallada y las
      notas técnicas (`dt` fijo, paralelización) pasan a `usage.qmd`/`api.qmd`/`numerics.qmd`;
      en el README queda quickstart + tabla mínima de API + notas de una línea con enlace.
      **Regla anti-deriva**: cada contenido vive en un solo lugar; la única duplicación
      permitida es el quickstart (README ↔ `index.qmd`).
- [ ] Todos los enlaces del README a la documentación con **URL absoluta al sitio**
      (`https://ener-habitat.github.io/EnerHabitat/...`) — PyPI no resuelve enlaces
      relativos; seguir el patrón de las imágenes (`raw.githubusercontent.com`).
- [ ] Añadir la ecuación de `Tn` (hoy solo se nombra Humphreys & Nicol) o enlazar a
      `docs/model-1d.md`.
- [ ] Sección **How to cite** (paper + paquete) enlazando `CITATION.cff`.
- [ ] Revisar TOC y que ninguna sección de uso pierda contenido (el uso ya está bien; solo se
      mueve teoría). Mantener las notas breves de `dt` fijo y de paralelización por procesos
      con enlace al detalle en `docs/numerics.md`.

### Fase 5 — Sitio de documentación: Quarto + GitHub Pages (mismo repo) ⏳

Publica `docs/` como sitio Quarto en `https://ener-habitat.github.io/EnerHabitat/` (estilo
windrose). Los `.qmd` de las Fases 1–3 son las páginas.

- [ ] **Proyecto Quarto en `docs/`**: `docs/_quarto.yml` (`project: type: website`), con
      `crossref` para ecuaciones, `bibliography: references.bib` y tema HTML (p.ej.
      `cosmo`/`flatly` + `toc`). El árbol de archivos, la navegación y el `_quarto.yml` de
      arranque están en [Estructura del sitio Quarto](#estructura-del-sitio-quarto-cómo-se-construye).
      La API se documenta manual al inicio (`api.qmd`); `quartodoc` para autodoc queda
      opcional/futuro.
- [ ] **`docs/references.bib`**: Barrios et al. 2016, Barrios et al. 2011/2012, Patankar
      1980, Xamán 2005, Hollands 1975, Humphreys & Nicol 2000, Chow & Levermore 2007,
      Morillón, NOM-008/020-ENER, ASHRAE 1997 (sacar del PDF y de `PLAN-2D-hecho.md`).
- [ ] **Ejemplos de `usage.qmd`** — matriz **sistema × modo**, los 4 como chunks
      ejecutables con su gráfica (`Ti` vs `Ta`/`Tsa`/`Tn`; energías en el caso AC):
      1. **1D libre** — muro de dos capas, `solve()` + `energy_transfer`
      2. **1D con AC** — mismo muro, `solveAC()` + `cooling_energy`/`heating_energy`
      3. **2D bloque hueco libre** — `HollowBlock` + `System2D.solve()` (+ `preview()` de
         la sección)
      4. **2D bloque hueco con AC** — mismo muro, `solveAC()`
      (Opcional: un 5.º ejemplo con `Slab` para techos, hoy presente en el README.)
- [ ] **Chunks ejecutables**: los ejemplos usan `enerhabitat` instalado + un `materials.ini`
      y un EPW de ejemplo dentro de `docs/data/` (decidir un EPW pequeño versionable).
      Activar `freeze: auto` para no re-ejecutar en cada render/CI y versionar
      `docs/_freeze/`.
- [ ] **Deploy**: al inicio manual con `quarto publish gh-pages` (crea la rama `gh-pages`);
      después, workflow `.github/workflows/docs.yml` con `quarto-dev/quarto-actions` en cada
      push a `main` que toque `docs/`.
- [ ] **Único paso manual en GitHub (una vez)**: *Settings → Pages* → source: rama
      `gh-pages` (`quarto publish` lo suele configurar solo; verificar).
- [ ] **`.gitignore`**: `docs/_site/` y `docs/.quarto/` (el HTML solo vive en `gh-pages`).
- [ ] Verificar el render local con `quarto preview docs` antes del primer publish
      (ecuaciones numeradas, citas, figuras, chunks).
- [ ] Cuidado: `docs/make_figures.py` y `docs/img/` conviven con el sitio — Quarto solo
      renderiza lo listado/`.qmd`; si molesta, mover el script a `tools/` o excluirlo con
      `render: [...]` en `_quarto.yml`.

### Fase 6 — Figuras y citabilidad ⏳

- [ ] `docs/make_figures.py`: añadir la figura del **dominio 1D con BCs** (y opcional: celda
      2D con BCs/simetría anotadas sobre las figuras existentes).
- [ ] **`CITATION.cff`** (paper como `preferred-citation`, autores, DOI, versión del paquete).
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
