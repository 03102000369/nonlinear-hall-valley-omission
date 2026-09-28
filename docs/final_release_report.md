# Final standalone release report

Final state: **PUBLISHED_PUBLIC**. Following the user's explicit instruction to publish publicly on 2026-09-28, the separate repository [https://github.com/03102000369/nonlinear-hall-valley-omission](https://github.com/03102000369/nonlinear-hall-valley-omission) was created and the validated release was pushed to `main`.

1. **Original project:** `../nonlinear_hall_transport`, retained read-only. Paths here are relative to avoid committing personal machine paths.
2. **Clean release:** `nonlinear_hall_transport_release` (this repository), a separate sibling.
3. **BH/NS preservation:** no BH/NS repository was entered, written, committed, merged or used as a source.
4. **Separation:** no BH/NS files, remote or history is present. No `.git` directory was copied. Initial parent and original-project identity checks found no Git repository or remote.
5. **Canonical run:** `20260919T021436Z-40be24e6`; scientific parent `20260918T215035Z-4260f56e`; selected controls `20260917T094513Z-083977ff`. Selection and verified hashes are in `release_provenance.md` and `source_provenance.json`.
6. **Included:** source, compact tests, explicit configuration, portable scripts, 23 processed tables, one concise canonical validation summary, final figures/previews, LaTeX and two PDFs, claims and reproducibility documentation. Every payload filename is listed in `release_files.txt`.
7. **Excluded:** original history and `.git`, all old releases/runs, dense raw arrays, duplicate long-form ablation, obsolete plots and regression comparison, full historical ledgers, notebooks, caches, virtual environments, IDE files, archives and LaTeX auxiliaries. The original keeps these files; none is needed by the independently tested release.
8. **Payload files:** 106, excluding Git metadata and ignored local environments/builds.
9. **Payload size:** 17561611 bytes (17.562 decimal MB). This is the published content, excluding `.git` objects and ignored dependency installations.
10. **Largest files:** see table below. No file approaches the normal 100 MB GitHub limit; no external large archive is required.
11. **Tests:** 40 passed, including checksum-mutation protection.
12. **Validation:** 23 hashes; 46 fresh independent checks; N=256 representative regression; 22,133 recomputed projection gates; 786,294 comparison cells. Every numerical difference was zero in the pinned environment and categorical decisions agree exactly. Twelve central CSVs are byte-identical; the thirteenth differs only in row order.
13. **Main compile:** PASS; no overfull boxes or undefined references/citations; all pages visually inspected.
14. **Supplement compile:** PASS; same checks and full visual inspection.
15. **Title:** When Can Hall-Silent Valleys Be Neglected in Anomalous-Velocity Nonlinear Hall Transport?
16. **Main pages:** 11 (including references).
17. **Supplement pages:** 12.
18. **Final figures:** five main and four supplementary, listed below. All main PNGs are byte-identical to their canonical originals; all nine agree with the independent release build.
19. **Final processed datasets:** 23, listed below and individually hashed in `DATA_MANIFEST.md`.
20. **Human decisions:** authors/order/affiliations and declarations; target journal; archival deposition and DOI; confirmation of software license and author-approved manuscript/data/figure licenses; complete citation metadata; optional final-full-text comparison for Liu et al. beyond its verified published abstract. Author scientific sign-off remains required before submission.
21. **Git branch:** main. This is fresh, unrelated history.
22. **Published scientific commits:** listed below. A subsequent documentation commit records the public repository URL and publication status without changing scientific content.
23. **Remote status:** `origin` is `https://github.com/03102000369/nonlinear-hall-valley-omission.git`; visibility is public; default branch is `main`. No tag or formal GitHub Release was created.
24. **Access:** clone using the command below. Journal submission and archival deposition are not part of GitHub publication.
25. **State:** PUBLISHED_PUBLIC.

## Largest payload files

| Bytes | Path |
|---:|---|
| 7534449 | `data/processed/final_mode_ablation.csv` |
| 3193649 | `data/processed/final_tolerance_cancellation.csv` |
| 2106630 | `data/processed/final_invariant_subspace_residual.csv` |
| 822451 | `data/processed/final_leakage_spectrum.csv` |
| 583163 | `manuscript/main.pdf` |
| 516937 | `supplement/supplement.pdf` |
| 379444 | `data/processed/weak_coupling_remainder.csv` |
| 267417 | `figures/supplement/figure_S2_weak_control.png` |
| 176436 | `figures/main/figure_3_leakage.png` |
| 163755 | `docs/final_claims.json` |

## Final figures

- `figures/main/figure_1_physical_modes`
- `figures/main/figure_2_ablation`
- `figures/main/figure_3_leakage`
- `figures/main/figure_4_scalar_shape`
- `figures/main/figure_5_decision_map`
- `figures/supplement/figure_S1_spectral_weights`
- `figures/supplement/figure_S2_weak_control`
- `figures/supplement/figure_S3_mass`
- `figures/supplement/figure_S4_successive_convergence`

Each stem has `.pdf` and `.png`. Input and generator hashes are recorded in `figures/manifest.json`.

## Processed dataset list

- `data/processed/certificate_decision_utility.csv`
- `data/processed/final_basis_transformation.csv`
- `data/processed/final_cancellation_statistics.csv`
- `data/processed/final_convergence_summary.csv`
- `data/processed/final_invariant_subspace_residual.csv`
- `data/processed/final_leakage_error_correlations.csv`
- `data/processed/final_leakage_spectrum.csv`
- `data/processed/final_leakage_summary.csv`
- `data/processed/final_mass_scaling.csv`
- `data/processed/final_mode_ablation.csv`
- `data/processed/final_mode_characterization.csv`
- `data/processed/final_mode_decision_agreement.csv`
- `data/processed/final_physical_validation_summary.csv`
- `data/processed/final_reduced_mode_coefficients.csv`
- `data/processed/final_reduced_model_validation.csv`
- `data/processed/final_tolerance_cancellation.csv`
- `data/processed/gauge_invariance.csv`
- `data/processed/internal_structure_robustness.csv`
- `data/processed/measurement_level_omission.csv`
- `data/processed/reduction_cost_benchmark.csv`
- `data/processed/successive_resolution.csv`
- `data/processed/weak_coupling_intervals.csv`
- `data/processed/weak_coupling_remainder.csv`

## Published scientific commit sequence

```text
dc4f2b5 Add reproducible nonlinear Hall transport source and tests
d0c5047 Add final manuscript datasets and figures
df009db Add manuscript, supplement, and reproducibility documentation
```

Publication documentation commit subject: `Record public GitHub publication and clone instructions`.

## Access the public repository

```sh
git clone https://github.com/03102000369/nonlinear-hall-valley-omission.git
cd nonlinear-hall-valley-omission
```

The remote belongs only to this nonlinear-Hall release. No BH/NS repository, history or remote was used or modified. The scientific payload was preserved through publication.

No tag, formal GitHub Release, DOI or duplicate archive was created. The previously suggested tag remains `v1.0.0-manuscript`; it has not been applied. Verification evidence is in `release_verification.md`, `release_verification.json` and `validation_summary.md`. Scientific assessment is in `scientific_self_review.md`; journal and archival metadata remain in `HUMAN_TODO.md`.
