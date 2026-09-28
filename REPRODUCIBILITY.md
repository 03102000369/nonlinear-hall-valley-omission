# Reproducibility

Clone the standalone repository, then work from its root. Tested with Python 3.12 and the exact versions in `requirements-lock.txt`; no access to the original research directory is required. Install a TeX distribution providing `pdflatex`, `bibtex`, Latin Modern, natbib, geometry, microtype and the standard mathematics/graphics packages.

```sh
git clone https://github.com/03102000369/nonlinear-hall-valley-omission.git
cd nonlinear-hall-valley-omission
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements-lock.txt
python -m pip install -e '.[test]'
make test
make validate
make reproduce-core
make figures
make paper
```

`make` uses `python3`; override with `make test PYTHON=path/to/python` when necessary. All targets set one BLAS/OpenMP thread. The plots use a local disposable Matplotlib cache under `build/`.

- **test**: independent Hamiltonian/Berry/symmetry/quadrature, conservation/PSD/zero modes, constant-lifetime, full/odd/Schur, frequency sign, texture covariance, mass, physical mode and checksum mutation tests.
- **validate**: exact included-data hashes; 46 independent band/transport checks; an N=256 finite-range M1–M4 regression; fixed-budget dependent-point sign/budget assertions.
- **reproduce-core**: directly constructs all 3,600 principal points, 300 two-kernel rows, five representatives, four spectra, 64 mass rows and four mechanism coefficient rows. Produces 13 central tables in `build/reproduced/processed/` and a numerical comparison report. No archived source snapshot or previous calculation output is a computational input. The canonical tables are used only for comparison after reconstruction.
- **figures**: regenerates all five main figures, four selected supplement figures, extended supplemental tables and the claim ledger from included canonical tables. Supplementary supporting controls are retained results, not rerun by this target.
- **paper**: runs figures, compiles main and supplement, and rejects missing inputs, overfull boxes and unresolved references/citations.
- **clean**: removes generated scratch outputs, caches and LaTeX auxiliaries. Canonical data and final figures/PDFs remain.

The comparison uses the full schema and keyed rows, absolute tolerance 1e-9 plus relative tolerance 1e-7 for floating-point cells, and exact categorical decisions. The production environment reproduces these comparisons; tiny residuals are not claimed as portable significant digits. Individual eigenvector signs/degenerate coordinates may change with LAPACK, so strict spectral-table reproduction should use the pinned environment. The physically meaningful signed sums, reconstruction identities and leakage powers are checked separately.

Canonical data have a fail-closed SHA-256 inventory in `data/manifest.json`. Deliberate scientific changes require reviewing the data, configuration, prose and manifest together; the pipeline never silently accepts new inputs. New calculations do not replace canonical files. To draw figures from a changed calculation, make an explicit reviewed data update rather than altering paths implicitly.

Only essential main calculations are rerun. The selected weak-series, convergence, certificate-cost, texture/gauge and open-circuit supporting tables carry original run and source hashes in `docs/source_provenance.json`; their figures regenerate directly from these included tables. Historical full scans and source snapshots remain in the original local scientific workspace and are unnecessary to use this release. No external archive is required.

See `docs/release_verification.md` for measured runtime, independent-copy execution, numerical comparison, page counts and visual inspection. Allow approximately 1 GB RAM and several minutes for warm numerical and document builds; dependency downloads and first font-cache construction take additional time.
