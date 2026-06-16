# How to serial interval?

Slides for the Lorentz Center workshop *Connecting Survival Analysis and
Infectious Disease Modeling* (Leiden, 16 June 2026).

A 20-minute talk on estimating the serial interval: what it is, why it is
really a censored and truncated delay-estimation problem, a within-host
informed approach, generic alternatives, lessons from applied work, and
in-flight tooling for composing the distributions involved.

## Build

```bash
task          # render to _site/index.html
task preview  # live preview
```

Rendered slides are deployed to GitHub Pages by the
[`render-and-deploy`](.github/workflows/render-and-deploy.yml) workflow on
every push to `main`.

## Sources and attribution

Figures are reused from the projects and papers they originate in; each slide
links to its source. See `references.bib` for full citations.
