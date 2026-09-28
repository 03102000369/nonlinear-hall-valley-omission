# Hall-silent valleys in anomalous-velocity nonlinear Hall transport

A conducting valley can have zero direct Berry-curvature Hall weight while still changing the Hall response through collisions. This repository accompanies **When Can Hall-Silent Valleys Be Neglected in Anomalous-Velocity Nonlinear Hall Transport?**

## Main result

Exact sector elimination and an ordered physical mode hierarchy separate active current, active valley distortion, passive current and passive valley polarization. Contact disorder closes on this space. Finite-range leakage has an observable-dependent effect. On the specified 3,510-point positive-coupling map, four modes reproduce the Hall response within 0.36%; safe omission of the entire passive sector is a separate, tolerance-dependent decision.

## Repository contents

- `src/nonlinear_hall/`: model, independent numerical checks, physical projection and figures.
- `tests/`, `scripts/`, `config/`: automated tests, portable entry points and explicit final settings.
- `data/processed/`, `validation/`: one canonical copy of each retained table and concise validation summaries.
- `figures/main/`, `figures/supplement/`: final vector figures and previews.
- `manuscript/`, `supplement/`: LaTeX sources and compiled PDFs.
- `docs/`: claims, methods, provenance and release verification.

Read the [manuscript](manuscript/main.pdf) and [supplement](supplement/supplement.pdf).

## Installation

Use Python 3.12 for the pinned environment. Other Python versions compatible with the package may produce small numerical differences.

```sh
git clone https://github.com/03102000369/nonlinear-hall-valley-omission.git
cd nonlinear-hall-valley-omission
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements-lock.txt
python -m pip install -e '.[test]'
```

LaTeX builds require `pdflatex`, `bibtex`, and the packages listed in the document preambles (TeX Live with the recommended fonts is suitable).

## Quick validation

```sh
make test
make validate
```

Validation verifies dataset hashes, runs independent physics checks, and recomputes a finite-range representative at the production angular resolution.

## Reproduce main results

```sh
make reproduce-core
```

This recomputes the 3,600-point principal map, two-kernel ablation, invariant and spectral residuals, fixed-budget decisions, four representative mechanism coefficients and finite-mass controls. It compares all 13 central tables numerically with the included canonical data. New results go to ignored `build/reproduced/`; canonical inputs are never silently replaced. Eigenvector coordinates can depend on the numerical environment; use the pinned environment for the strict table comparison.

## Reproduce manuscript figures

```sh
make figures
```

All five main and four supplementary figures regenerate from included tables. Input manifests identify their data and generator hashes.

## Build manuscript

```sh
make paper
```

The build verifies data, regenerates figures/tables and compiles both documents, rejecting missing inputs, unresolved references and overfull boxes. `make clean` removes disposable builds, caches and LaTeX auxiliary files, retaining final figures and PDFs.

## Data

[DATA_MANIFEST.md](DATA_MANIFEST.md) lists every included dataset and checksum. [REPRODUCIBILITY.md](REPRODUCIBILITY.md) gives the clean-clone workflow. [Claim audit](docs/final_claims.md) maps quantitative paper statements to values, selections and run IDs. Canonical physical analysis: `20260919T021436Z-40be24e6`.

## Scientific scope and limitations

The calculation covers the anomalous-velocity channel in a homogeneous, zero-temperature, phenomenological four-valley model with reversible Born scattering and restricted disorder ensembles. It does not calculate all nonlinear Hall channels, predict quantitative SnTe transport, establish universal error thresholds, or provide a complete experimental device signal. Scalar/shape pieces are projection diagnostics, not separately measurable microscopic currents.
