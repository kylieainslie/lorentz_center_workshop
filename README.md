# Connecting Survival Analysis and Infectious Disease Modeling

**Lorentz Center Workshop | 15–19 June 2026 | Leiden, the Netherlands**

Website: <https://kylieainslie.github.io/lorentz_center_workshop>

---

## About

This repository contains the materials, case study data, analysis scripts, and vignettes for the SA–IDM workshop. The workshop brings together researchers from survival analysis (SA) and infectious disease modelling (IDM) to make methodological connections explicit and develop shared tools for cross-disciplinary analysis.

Over five days, participants work in mixed groups on a shared simulated epidemic dataset and collectively produce vignettes, code examples, and a "Ten Simple Rules" guide for SA–IDM integration.

## Repository Structure

```
.
├── _quarto.yml              # Quarto website configuration
├── index.qmd                # Home page
├── case-studies.qmd         # Case study overview and group questions
├── preparation.qmd          # Pre-workshop reading and setup
├── day1.qmd – day5.qmd      # Daily materials
├── organizers.qmd           # Scientific organizers
├── assets/
│   ├── presentations/       # Workshop presentations (.qmd)
│   └── case_study/          # Case study documents (PDF, DOCX)
├── data/                    # Simulated epidemic datasets (.rds)
├── analysis/
│   ├── index.qmd            # Contributor guide for group scripts
│   └── group-*/             # Group analysis scripts (one folder per group)
├── vignettes/
│   ├── index.qmd            # Vignette listing
│   ├── resources.qmd        # Papers, software, and links
│   └── vignette-group*.qmd  # Cross-disciplinary worked examples
└── docs/                    # Rendered website output (do not edit directly)
```

## Groups and Questions

| Group | Topic |
|-------|-------|
| A | Key Epidemiological Parameters |
| B | Serial Interval Distribution |
| C | Reporting Delays and Underreporting |
| D | Vaccine Effectiveness |
| E | Excess Mortality |

Each group has a dedicated branch (`group-A` through `group-E`) for analysis scripts and vignettes.

## Working with Branches

The website lives on `main`. Group work lives on branch `group-X`:

```bash
git checkout group-D          # switch to your group's branch
git pull origin group-D       # get latest changes
```

To pick up updates from `main`:

```bash
git merge main
```

## Scientific Organizers

| Name | Affiliation |
|------|-------------|
| Liesbeth de Wreede | Leiden University Medical Center (LUMC) |
| Hein Putter | Leiden University Medical Center (LUMC) |
| Kylie Ainslie | University of Melbourne & University of Hong Kong |
| Don Klinkenberg | RIVM & Wageningen University and Research |
| Jacco Wallinga | RIVM & Leiden University Medical Center (LUMC) |
| Steven Abrams | UHasselt & University of Antwerp |

## License

Workshop materials are shared for educational use. Please contact the organizers before reusing or adapting content.
