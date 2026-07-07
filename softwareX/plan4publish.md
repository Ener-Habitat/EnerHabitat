# Plan de publicación — EnerHabitat en SoftwareX

**Objetivo:** llevar `softwareX/EnerHabitat/enerhabitat.tex` a un manuscrito listo para
enviarse a SoftwareX como *Original Software Publication* (OSP), cumpliendo la guía de
autores al 100 %, con hechos verificados contra el estado real del software.

**Archivo a editar:** `softwareX/EnerHabitat/enerhabitat.tex` (LaTeX puro, sin paquetes
ni trucos sofisticados; `main.tex` queda solo como cantera de texto — tiene mejores
párrafos de motivación/landscape que se pueden aprovechar).

**Template oficial:** `softwareX/softwarex-osp-template.tex` (versión 2026, descargado
2026-07-07). Manda sobre cualquier otra referencia. Diferencias clave vs el template
2023 con el que se armó el borrador: (i) tabla de código **C1–C8, sin fila
"Reproducible Capsule"**; (ii) encabezados `Required Metadata` + `Current code
version` antes de la tabla; (iii) la información del ejecutable va en la sección
`Current executable software version` **al final, después de las referencias**;
(iv) límite 4,000 palabras / 6 páginas de texto principal.

**Referencias de estilo:** los dos ejemplos publicados en `softwareX/examples/`
(EnTiSe, SoftwareX 2026, 102784; AprioriAllPkg, SoftwareX 2026, 102806) — útiles para
tono y contenido, pero se publicaron con el template viejo: **en formato gana el
template 2026**.

---

## 1. Estado actual verificado (hechos, no lo que dice el manuscrito)

| Elemento | Realidad verificada |
|---|---|
| Paquete | `enerhabitat` **v0.2.1** en PyPI, Python ≥ 3.10, deps: `numba>=0.62.1`, `pvlib>=0.13.1`; extra `[viz]` con matplotlib |
| Repo del paquete | `https://github.com/Ener-Habitat/EnerHabitat` — README completo, LICENSE (MIT), código en `src/` ✓, CITATION.cff ✓, CHANGELOG ✓, tests ✓ |
| Tags de release | Existen v0.1.8, v0.1.10, v0.2.0 — **falta el tag v0.2.1** |
| API real 1D | `Location`, `System`, `config`: `loc.meanDay(month=...)`, `sys.Tsa()`, `sys.solve()` (free-running), `sys.solveAC()`, atributos `energy_transfer`, `cooling_energy`, `heating_energy`. El Listing 1 del manuscrito es consistente con la API real (verificar ejecutándolo) |
| **API 2D (¡ya existe!)** | `System2D`, `HollowBlock` (bloque hueco de concreto, cavidad de aire o relleno sólido), `Slab` (vigueta y bovedilla), `Fill`, `Section2D`, `config2d`. Física de cavidad: radiación con factores de vista + convección Nusselt dependiente de T. Inspector de sección a escala (`preview()`, `section_report()`) |
| Documentación | `https://ener-habitat.github.io/EnerHabitat/` (Quarto): usage con ejemplos ejecutables 1D/2D, teoría 1D, teoría 2D, método numérico, API reference |
| Validación | `https://github.com/Ener-Habitat/eh_validation` — notebooks: `000` (1D vs EnergyPlus sin AC), `001` (1D vs EnergyPlus con AC), `002` (2D sin AC), `003` (R estacionario 2D vs datos experimentales de Borbón et al. 2010, placa caliente ASTM C177). **README.md vacío (0 bytes), sin LICENSE** |
| Webapp (código) | `https://github.com/Ener-Habitat/EnerHabitat-webapp` — versión **2.9.0** en pyproject (`eh-shiny`), deps: `enerhabitat>=0.1.7`, `shiny>=1.3.0`, `plotly>=6.0.1`, `shinywidgets>=0.5.2`. **README.md vacío (0 bytes), sin LICENSE** |
| Webapp (servicio) | `https://enerhabitat.unam.mx/` — modelo 1D, para usuarios sin experiencia en Python |
| Palabras del borrador | ~1,840 (límite 4,000) → hay espacio de sobra para incorporar el 2D |
| Figuras del borrador | **0** (la Fig. de arquitectura está comentada, pero `Figure~\ref{fig:architecture}` sigue citada en el texto → produce "??") |

## 2. Discrepancias del manuscrito a corregir (fact-check)

1. **C1 versión del código:** dice v0.1.7 → debe ser **v0.2.1** (y crear el tag).
2. **C2 repo:** dice `github.com/AltamarMx/EnerHabitat` → debe ser
   `https://github.com/Ener-Habitat/EnerHabitat` (idealmente el permalink al tag:
   `.../EnerHabitat/tree/v0.2.1`, como hace EnTiSe).
3. **Fila de documentación (C7 en el template 2026):** dice el repo → debe ser
   `https://ener-habitat.github.io/EnerHabitat/`. Además la tabla actual usa la
   numeración vieja C1–C9 con "Reproducible Capsule": renumerar a C1–C8 y mover la
   tabla S al final del documento (sección `Current executable software version`).
4. **Repo de validación:** dice `eh_validation_eplus` → es `eh_validation`.
5. **URL webapp:** dice `www.enerhabitat.unam.mx` → confirmar canónica (`https://enerhabitat.unam.mx/`).
6. **S1 versión webapp:** dice "v0.9.0 (Beta)" → pyproject dice **2.9.0**. Decidir cuál reportar.
7. **S6:** actualizar deps (`enerhabitat>=0.2.1` si la webapp se actualiza, `shinywidgets` falta).
8. **El 2D aparece como "trabajo futuro"** en Conclusions → **ya está implementado y
   validado**. Es EL diferenciador científico (el paper de 2016 se distinguía justamente
   por sistemas no homogéneos). Hay que promoverlo a funcionalidad central del artículo.
9. **`\ref{fig:architecture}` roto** (figura comentada). O se crea la figura o se quita la referencia.
10. **Autores — RESUELTO (sección 7):** nuevo orden Cruz Salas, Rodríguez Calderón,
    Huelsz, Rojas, Barrios (corresponding). Rehacer frontmatter y CRediT para los 5.
11. **Afiliación IPN incompleta:** la guía exige dirección postal completa — pendiente
    de Fernando; dejar `%% TODO` en el .tex.
12. **Validación descrita:** el texto habla de 5 materiales en Temixco/mayo — cotejar con lo
    que realmente hacen los notebooks 000/001 y añadir la validación 2D (002/003, Borbón).
13. **Claim "one million evaluations"**: mantener solo si es defendible (aparecía en el sitio
    original); en Impact citarlo como historial de la herramienta Java 2005–202X.
14. **Ejemplo 2 (validación) referencia repo equivocado** y no menciona notebooks concretos.

## 3. Checklist de cumplimiento — guía de autores SoftwareX

### Estructura obligatoria (el paper se devuelve si falta)
- [ ] Las **5 secciones del template**: 1. Motivation and significance · 2. Software
      description (2.1 Architecture, 2.2 Functionalities, 2.3 snippets opcional) ·
      3. Illustrative examples · 4. Impact · 5. Conclusions. *(Ya están; conservar numeración.)*
- [ ] **Tabla de metadatos de código C1–C8 (template 2026)** bajo los encabezados
      `Required Metadata` / `Current code version`: C1 versión · C2 link GitHub
      permanente (obligatorio) · C3 licencia · C4 versionado · C5 lenguajes/herramientas ·
      C6 requisitos y dependencias · C7 documentación · C8 email de soporte.
      **Sin fila "Reproducible Capsule"** (eliminada en el template 2026).
- [ ] Info del ejecutable (webapp) en la sección **`Current executable software
      version` al final del documento** (después de referencias), como tabla S1–S8
      — opcional, la usamos para la webapp.
- [ ] Máx **4,000 palabras** (abstract + texto + captions + footnotes) — presupuesto abajo.
- [ ] Máx **6 figuras**.
- [ ] Abstract ≈ **100 palabras** (el actual tiene ~120 y hay que reescribirlo con el 2D).
- [ ] **1–7 keywords** (hay 6; revisar; evitar multi-palabra donde se pueda).
- [ ] Título formato "SoftwareName: short title" ✓ (revisar longitud/claridad).
- [ ] Secciones numeradas y citadas por número (no "see the text").

### Repositorio (la guía es explícita: el paper se devuelve si falta)
- [ ] Repo GitHub público ✓
- [ ] `README.md` bien documentado ✓ (paquete) / **✗ webapp y validación (vacíos)**
- [ ] `LICENSE.txt` — el paquete tiene `LICENSE`; **añadir copia `LICENSE.txt`** (costo
      cero, evita una devolución administrativa) y añadir licencia a webapp y validación.
- [ ] Código fuente en `repo/src` ✓ (paquete).
- [ ] Licencia aprobada OSI: MIT ✓

### Declaraciones y secciones finales
- [ ] **CRediT** completo y consistente con la lista de autores (incluye a todos).
- [ ] **Declaration of competing interest** ✓ (ya está).
- [ ] **Declaración de IA generativa** ✓ (ya está; requerida; va antes de referencias).
- [ ] **Funding — RESUELTO:** sin financiamiento del desarrollo nuevo → frase estándar
      "This research did not receive any specific grant from funding agencies in the
      public, commercial, or not-for-profit sectors." (CONACYT-SENER 118665 del original
      se menciona solo en Acknowledgements.)
- [ ] Acknowledgements en sección separada antes de referencias ✓.
- [ ] **Data statement**: los EPW y datos de validación están en `eh_validation` → citarlo.

### Referencias (estilo numérico Elsevier)
- [ ] `elsarticle-num` ✓; números en corchetes; orden de aparición.
- [ ] **Cita de software con formato FORCE11**: citar el propio EnerHabitat
      (versión + PID). Recomendado: **archivar release v0.2.1 en Zenodo → DOI** y añadir
      `@misc` tipo software (como la ref. [7] modelo de la guía: "Version 0.2.1
      [software]. Zenodo; 2026. https://doi.org/10.5281/zenodo.XXXXX").
- [ ] Citar también: pvlib (Anderson 2023 — actualizar de Holmgren 2018), numba, Shiny
      for Python, Plotly, EnergyPlus, Patankar (TDMA), Humphreys & Nicol (confort
      adaptativo), Borbón et al. 2010 (validación experimental 2D), Barrios 2011/2012/2016,
      NOM-020-ENER. Añadir URL/DOI y fecha de acceso a las referencias web.
- [ ] Sin referencias huérfanas (todas citadas ↔ todas listadas).

### Envío (Editorial Manager)
- [ ] Article Type: **Original Software Publication**, https://www.editorialmanager.com/softx
- [ ] Subir **PDF compilado** + **zip con fuentes** (.tex, .bib/.bbl, figuras, .cls/.bst).
- [ ] **Highlights**: archivo aparte, 3–5 bullets ≤ 85 caracteres c/u.
- [ ] Graphical abstract (opcional): 531×1328 px min.
- [ ] Declaración de conflictos vía el "declarations tool" (.docx).
- [ ] Figuras además como archivos separados (`Figure_1.pdf`, ...): vector PDF/EPS
      o PNG ≥ 300 dpi (1063 px ancho de columna mínimo).

## 4. Estructura y contenido propuestos (presupuesto ≈ 3,300–3,800 palabras)

**Título:** "EnerHabitat: A Python package and web application for the dynamic thermal
evaluation of homogeneous and non-homogeneous building envelope components" *(u opción
más corta; el guiño "homogeneous and non-homogeneous" conecta con el paper 2016)*.

**Abstract (~100 palabras):** qué es (paquete + webapp), qué resuelve (conducción
1D multicapa y **2D para bloque hueco / vigueta y bovedilla con cavidades de aire**),
cómo (volúmenes de control implícitos + TDMA + numba; EPW vía pvlib; día promedio
mensual), modos free-running/AC, validación (EnergyPlus + experimental), sucesor del
Ener-Habitat Java (>1M evaluaciones).

**Metadata C1–C8 (template 2026, al frente):** v0.2.1 · permalink al tag GitHub · MIT ·
git · Python/numba/pvlib/pandas/numpy · Python ≥ 3.10 multiplataforma · docs en
GitHub Pages · gbv@ier.unam.mx. (Ya no existe fila de Reproducible Capsule; el
Binder/Zenodo, si se hace, se cita en References o README.)
**Sección final `Current executable software version` (webapp, S1–S8):** **v3.0**
(BLOQUEO: liberar tras integrar la retro de Guadalupe y coautores — ver sección 7.3) ·
https://enerhabitat.unam.mx/ · MIT · web (cualquier navegador) · deps
shiny/plotly/enerhabitat · manual = misma URL. Va **después de References**, según
el template 2026.

1. **Motivation and significance (~750 pal.)** — base: párrafos de `main.tex` (mejor
   argumentados que los de enerhabitat.tex): U-value estacionario engaña en clima cálido
   y vivienda free-running; historia Ener-Habitat Java (validado, >1M evaluaciones,
   descontinuado por infraestructura); landscape (EnergyPlus/TRNSYS = curva de
   aprendizaje; calculadoras R = sin masa térmica; THERM/Heat2 = sin métricas dinámicas
   de envolvente); el nicho: herramienta ligera centrada en el componente, open source,
   1D **y 2D**, con webapp para no programadores. Workflow del usuario en un párrafo.
2. **Software description (~1,000 pal.)**
   - 2.1 Arquitectura: ecosistema de 3 repos (paquete / webapp / validación) +
     clases `Location`, `System`, `System2D` (+ `HollowBlock`, `Slab`), `Config`;
     numba JIT, pvlib; webapp Shiny como capa fina sin lógica propia. → **Figura 1**.
     **Obligatorio aclarar (decisión 2026-07-07):** (i) la tabla S refiere a la
     webapp = forma ejecutable online del software para no programadores, y el texto
     debe decir explícitamente que la webapp expone el **modelo 1D**; (ii) enlazar el
     repo del código de la webapp (`EnerHabitat-webapp`) en esta sección.
   - 2.2 Funcionalidades: día promedio desde EPW; Tsa (tilt/azimuth/absortancia, RF);
     `solve()` free-running vs `solveAC()` (Tn adaptativa Humphreys–Nicol);
     multicapa hasta N capas; 2D con cavidades (radiación + Nusselt); métricas
     (energía, FD, TR); webapp: hasta 5 sistemas, 8 orientaciones, export CSV.
   - 2.3 Snippet: Listing 1 (1D, el actual **verificado ejecutándolo**) y Listing 2
     corto (2D `Slab`, vigueta y bovedilla) — los listings no cuentan como figuras.
3. **Illustrative examples (~800 pal.)**
   - Ej. 1: comparación de 2–3 sistemas constructivos 1D free-running para Temixco
     (código + **Figura 2**: Ta/Tsa/Ti del día promedio con banda de confort).
   - Ej. 2: sistema 2D — **techo de vigueta y bovedilla** (`Slab`), el caso más
     representativo de la construcción mexicana (**Figura 3**: dibujo del sistema
     constructivo + sección del solver; **Figura 4**: Ti o energías comparadas).
   - Ej. 3: validación — notebooks 000/001 (1D vs EnergyPlus, con/sin AC) y 003
     (R estacionario vs Borbón exp.); números de acuerdo concretos (0.5 °C, etc.,
     verificados contra los notebooks). Referencia al repo `eh_validation`.
   - Ej. 4 (breve): workflow webapp (**Figura 5**: captura de pantalla) — el camino
     para no programadores.
4. **Impact (~700 pal.)** — la sección que más pesan los revisores:
   preguntas nuevas (barridos paramétricos programáticos, indicadores regulatorios
   más allá del U-value, acoplamiento con LCA/optimización); mejora de flujos
   existentes (reemplaza el servicio Java cerrado; reproducible; pip install);
   práctica diaria (docencia IER-UNAM, cribado temprano de envolventes, NOM-020);
   adopción — historial >1M evaluaciones del original con la redacción transparente
   acordada en 7.6 ("internal usage records... no longer archived"), PyPI, webapp
   pública UNAM; comercial (MIT permite integración en consultorías). Citar papers
   que usaron Ener-Habitat [2011, 2012, 2016 y los que confirme Guillermo].
5. **Conclusions (~200 pal.)** — qué es, nicho, y trabajo futuro REAL (ya no el 2D):
   más geometrías no homogéneas, métricas de confort/sobrecalentamiento, contenedores,
   acoplamiento con simuladores de edificio completo.

Secciones finales: CRediT (todos los autores) · Competing interests · Declaración IA ·
Funding · Acknowledgements · References.

## 5. Figuras (máx 6; usar ≤ 5 para margen)

| # | Contenido | Fuente |
|---|---|---|
| 1 | Arquitectura del ecosistema (paquete–webapp–validación, flujo EPW→meanDay→Tsa→solve) | crear (diagrama simple, vector) |
| 2 | Día promedio 1D: Ta, Tsa, Ti de 2–3 sistemas + zona de confort | script con el paquete (matplotlib) |
| 3 | **Vigueta y bovedilla**: dibujo del sistema constructivo (esquema 3D/isométrico o corte anotado: vigueta L, bovedilla, capa de compresión, cavidades) + sección a escala del solver (`preview()`) — puede ser figura compuesta (a)/(b) | crear dibujo + `preview()`/docs `make_figures.py` |
| 4 | Resultado 2D del techo vigueta y bovedilla (Ti free-running y/o energías AC; o panel de validación vs EnergyPlus/Borbón) | notebooks de `eh_validation` |
| 5 | Captura webapp | screenshot https://enerhabitat.unam.mx/ |

Todas citadas en el texto, numeradas por orden, caption con título breve + descripción,
archivos separados con nombres `Figure_1`, `Figure_2`, ...

## 6. Acciones en los repos (fuera del manuscrito, pero bloquean el envío)

1. ~~Tag v0.2.1~~ **HECHO 2026-07-07** (`git tag -a v0.2.1` en `f01f96b`, pusheado;
   el permalink de C2 ya resuelve). Opcional: crear también el "Release" en GitHub
   con notas del changelog.
2. **README.md de `eh_validation`** (hoy 0 bytes): qué valida, estructura (epw/, idf/,
   notebooks/, pdfs/), cómo reproducir, tabla resumen de resultados.
3. **README.md de `EnerHabitat-webapp`** (hoy 0 bytes): qué es, URL del servicio,
   cómo correrla localmente, relación con el paquete.
4. **Licencia en `eh_validation` y `EnerHabitat-webapp`** (hoy no tienen): MIT.
5. **`LICENSE.txt`** junto a `LICENSE` en el repo principal (la guía lo pide literal).
6. *(Pendiente — 7.7)* **Zenodo DOI** del release v0.2.1 para la cita de software FORCE11.
7. **Liberar webapp v3.0** — bloqueado por integrar la retro de Guadalupe y coautores
   (ver 7.3); alinear pyproject (hoy 2.9.0) con la versión reportada en S1.
8. **Crear repo `eh_usage`** (Ener-Habitat/eh_usage) con variaciones de uso del paquete:
   recetas cortas (barrido de absortancia, comparación de N sistemas, barrido de meses,
   materiales custom, uso de config, 2D básico), cada una como notebook o script
   autocontenido + README índice. Si llega antes del envío, se cita en el paper.

## 7. Decisiones tomadas (2026-07-07)

1. **Autores (orden definitivo):**
   1. Miriam Cruz Salas (IER-UNAM)
   2. Fernando Rodríguez Calderón (IPN)
   3. Guadalupe Huelsz (IER-UNAM)
   4. Jorge Rojas (IER-UNAM)
   5. **Guillermo Barrios (IER-UNAM) — corresponding author** (gbv@ier.unam.mx)

   CRediT debe cubrir a los 5. Roles propuestos (confirmar con cada quien antes del envío):
   - M. Cruz Salas: Validation, Investigation, Writing – review & editing
   - F. Rodríguez Calderón: Software, Validation, Writing – review & editing
   - G. Huelsz: Conceptualization, Methodology, Writing – review & editing
   - J. Rojas: Conceptualization, Methodology, Writing – review & editing
   - G. Barrios: Conceptualization, Methodology, Software, Validation, Supervision,
     Writing – original draft
2. **Afiliación de Fernando (IPN): PENDIENTE** — dejar placeholder
   `\address[ipn]{Instituto Politécnico Nacional, [UNIDAD Y DIRECCIÓN PENDIENTES], Mexico}`
   y marcar con `%% TODO` en el .tex. Bloquea el envío, no la redacción.
3. **Versión webapp (S1): 3.0** — reportar "v3.0" cuando se libere.
   > **NOTA-BLOQUEO:** la 3.0 no existe aún; requiere que Guillermo **reciba la
   > retroalimentación de Guadalupe y demás coautores y la integre** a la webapp.
   > Hasta entonces el .tex lleva `%% TODO v3.0` en S1/S6 y no se envía el paper.
4. **Funding:** ninguno para el desarrollo nuevo → usar la frase estándar de la guía:
   "This research did not receive any specific grant from funding agencies in the
   public, commercial, or not-for-profit sectors." (El CONACYT-SENER 118665 del
   Ener-Habitat original se queda solo en Acknowledgements, como contexto histórico.)
5. **Ejemplo 2D: vigueta y bovedilla (`Slab`, techo)** — es el caso más complejo y
   luce mejor. Incluir **dibujo del sistema constructivo** (ver figuras, sección 5).
6. **Claim ">1M evaluaciones":** SÍ se sostiene, pero las estadísticas eran auditorías
   internas de la plataforma Java que se perdieron. Redacción honesta y defendible:
   "Internal usage records of the original platform registered over one million
   thermal evaluations during its ~10 years of operation; these logs are no longer
   archived." — un claim histórico transparente sobre su procedencia; si un revisor
   lo objeta, se degrada a "widely used" citando las tesis/papers que la usaron.
7. **Zenodo DOI: PENDIENTE** — mantener en la lista de la Fase 0; la cita de software
   en References queda con el permalink del tag de GitHub mientras tanto.
8. **Repo de "usage" (nuevo):** crear un repositorio con variaciones de uso del
   paquete (recetas/casos). **Nombre recomendado: `eh_usage`** — consistente con
   `eh_validation` (misma convención `eh_*`, mismo rol de repo compañero) y con la
   página "Usage" de la documentación. Alternativas si se prefiere otro matiz:
   `eh_examples` (más obvio para quien busca ejemplos), `eh_cookbook` (connotación
   de recetas paso a paso), `EnerHabitat-examples` (convención de la webapp).
   Si se crea antes del envío, citarlo en Illustrative examples junto a `eh_validation`.

## 8. Fases de ejecución

- **Fase 0 — Decisiones tomadas (sección 7)** ✓ 2026-07-07. Quedan pendientes que
  bloquean el ENVÍO pero no la redacción: afiliación de Fernando (7.2), liberación
  webapp 3.0 tras integrar retro de coautores (7.3), Zenodo (7.7), confirmación de
  roles CRediT con cada coautor. Acciones de repos (sección 6) pueden ir en paralelo.
- **Fase 1 — Correcciones factuales** ✓ **HECHA 2026-07-07** (cambios 1–7 autorizados
  uno a uno): frontmatter con orden nuevo de autores (Cruz Salas → ... → Barrios corr.),
  tabla C1–C8 template 2026 (v0.2.1, permalink al tag creado, docs GitHub Pages),
  tabla S1–S7 movida al final (`Current executable software version`, v3.0 con
  TODO-bloqueo), `\ref` roto eliminado (TODO Fase 3), URLs corregidas (eh_validation,
  enerhabitat.unam.mx), CRediT de 5 autores, sección Funding. Compila limpio,
  0 refs rotas, 1,875 palabras. TODOs vivos en el .tex: afiliación IPN, roles CRediT
  por confirmar, webapp v3.0, reescritura 2D en Conclusions (Fase 2).
- **Fase 2 — Contenido** ✓ **HECHA 2026-07-07** (cambios 8–15 autorizados uno a uno):
  título "time-dependent ... homogeneous and non-homogeneous"; abstract ~105 pal. con
  1D+2D (ejemplo: vigueta y bovedilla); keywords nuevas; Motivation fusionada con
  main.tex (U-value, landscape, 1M honesto, workflow); 2.1 ecosistema 3 repos +
  System2D + volúmenes de control + link webapp repo + convención exterior→interior
  (como EnergyPlus); 2.2 con funcionalidades 2D y unidades corregidas (J/m²·día,
  verificado ejecutando); sección 3 con 4 ejemplos (Listing 1D corregido y ejecutado;
  Listing 2D idéntico al script de docs, energía 27,041 J/m²·día, ~15 min; validación
  según notebooks reales 000–003 incl. Borbón; webapp); Impact reestructurado por
  criterios de la guía; Conclusions con 2D como presente y future work real; refs:
  pvlib→Anderson 2023, autocita de software FORCE11 (tag v0.2.1), Borbón añadida,
  Incropera eliminada (huérfana). Compila limpio, 0 refs rotas, **2,791 palabras**.
  Pendiente Fase 4: verificar/añadir DOIs faltantes (barrios2011/2012, crawley2008).
- **Fase 3 — Figuras** ✓ **HECHA 2026-07-07** (cambios 16–19 + screenshot): 5 figuras
  insertadas, todas con script reproducible en `softwareX/EnerHabitat/figures/`
  (`make_fig1..4.py`, `run_slab_cases.py`, `materials.ini` del paper en inglés):
  Fig 1 arquitectura del ecosistema · Fig 2 día promedio 1D del Listing 1 (EPS 1 in +
  concreto 12 cm vs concreto solo; oscilación 14→2 °C) · Fig 3 sección a escala del
  `preview()` con leyenda = strings del Listing 2 (High-density concrete / Aerated
  concrete / Filler block) · Fig 4 resultado 2D free-running (Ti amortiguada ~4 °C,
  Tsa nocturna < Ta por radiación al cielo; 11,539 J/m²·día, 11 días, ~21 min) ·
  Fig 5 screenshot webapp 2560×1493 px (caso AC, 3 sistemas, muro este, Cuernavaca/mayo;
  interfaz en español — aclarado en caption y texto). Datos AC de respaldo en
  `figures/data/summary.json` (cooling 259,959 / heating 151,361 J/m²·día).
  Compila limpio, 0 refs rotas, **3,222 palabras**. TODOs vivos en el .tex: solo los
  3 bloqueos de envío (afiliación IPN, roles CRediT, webapp v3.0).
- **Fase 4 — Verificación** ✓ **HECHA 2026-07-07** (cambio 20 incluido):
  - compila limpio (0 errores, 0 warnings bibtex, 0 refs indefinidas, 0 "??");
  - Listing 1 **ejecutado verbatim** (archivos nombrados como en el código) contra
    v0.2.1 — corre sin errores, unidades correctas; Listing 2 verificado idéntico
    parámetro a parámetro al script que produjo los números publicados;
  - 3,222 palabras ≤ 4,000 · abstract 103 ≈ 100 · 6 keywords · 5 figuras ≤ 6;
  - validación cotejada contra notebooks reales (Temixco_2018CST.epw, este, mayo;
    5 materiales; 2D consistencia; Borbón ASTM C177);
  - DOIs verificados vía Crossref y añadidos (Barrios 2011/2012, Crawley 2008,
    Borbón 2010 + number=6); refs web con fecha de acceso;
  - checklist sección 3: todo ✓ salvo los bloqueos de envío listados abajo.
- **Fase 5 — Paquete de envío** (parcial, 2026-07-07): `highlights.txt` creado
  (5 bullets, todos ≤ 85 caracteres, validado por script) · figura renombrada
  `Figure_5.png` (consistente con Figure_1..4) · PDF compilado limpio.
  **Pospuesto por decisión de Guillermo:** zip de fuentes (se arma al momento del
  envío) y commit (tras su revisión). **Al enviar:** declaración de conflictos se
  genera en el "declarations tool" de Editorial Manager (.docx) · graphical abstract
  opcional (decidir) · Article Type = "Original Software Publication" en
  https://www.editorialmanager.com/softx · subir PDF + zip de fuentes (.tex, .bib,
  .bst, .cls, figures/) + highlights.txt como archivo aparte.

---
*Notas de estilo: mantener LaTeX plano (elsarticle + amssymb/hyperref/listings/graphicx
ya cargados, nada más). Inglés, lenguaje inclusivo, sin "we believe" — hechos citables.*
