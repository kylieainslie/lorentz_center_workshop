# Intro to Infectious Disease Modelling — for survival analysts

Quarto reveal.js deck. ~40 slides, paced for 45–60 minutes.

## Render

```bash
# from this folder:
quarto render intro_idm_for_survival.qmd            # produces .html
quarto preview intro_idm_for_survival.qmd           # live reload while editing
```

## Required R packages

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "deSolve", "scales", "svglite"))
```

`svglite` is used as the figure device so the deck builds on stock R / macOS
without needing XQuartz (X11) or Cairo. To switch to PNG instead, change
`dev: "svglite"` to `dev: "ragg_png"` in the YAML header (install the
`ragg` package first).

Mermaid diagrams render natively in Quarto — no extra packages needed.

## Structure

1. Hazard ↔ force of infection (incl. Bernoulli 1760, frailty)
2. SIR / SEIR as multistate models (incl. Reed–Frost 1928)
3. R₀, Rₜ, renewal equations (incl. final-size relation, 1918 Philadelphia)
4. Censoring & truncation in outbreak data (incl. John Snow 1854)
5. Synthesis (dictionary, software bridges, take-aways)

Each section opens with a coloured divider slide. Most content slides
are two-column with prose on the left and a generated figure (R or
Mermaid) on the right.

## Editing tips

- All figures are reproducible from the inline code chunks — tweak
  parameters (β, γ, σ, frailty multipliers, …) directly in the .qmd.
- Speaker notes live in `::: {.notes}` blocks (only the title slide
  has one so far — easy to add more).
- To time the talk, use `quarto preview` and the `s` key in reveal.js
  to open speaker view with a built-in timer.
- For a darker theme, swap `theme: simple` for `theme: dark` or
  `theme: blood` in the YAML.

## Potential additions (if you want to grow it)

- A live data demo — pull recent COVID-19 incidence and run `EpiNow2`
  in an appendix slide.
- A worked example connecting Cox PH to an exposure-stratified FoI
  model on a single dataset.
- A phylogenetic / coalescent slide showing how trees give independent
  windows on the same hazard process.
