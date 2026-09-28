# Final dataset manifest

Only datasets physically included in this repository are listed. All paths are relative. Canonical tables are read-only inputs to the release pipeline; reconstruction outputs are disposable. Source-selection hashes and aggregation rules are in `docs/source_provenance.json`.

| Path | Rows | Purpose / figure | Source run | SHA-256 |
|---|---:|---|---|---|
| `data/processed/certificate_decision_utility.csv` | 24 | Supplement certificate table and coverage | `20260917T094513Z-083977ff` | `5c0375c30eed7c66dece89897acdb56e20bb210b602b29a695e41d8fe90223ac` |
| `data/processed/final_basis_transformation.csv` | 1824 | Raw-to-orthonormal mode transform; Sec. III and Supplement mode definitions | `20260919T021436Z-40be24e6` | `9e1b275dd905c6676ebc9c7edb48797db3fedf64de0b9635350e88cc212f3663` |
| `data/processed/final_cancellation_statistics.csv` | 30 | Main Table I; abstract and conclusion; all-grid Supplement table | `20260919T021436Z-40be24e6` | `a90da2b966ff0186bfd323160c3d57246520b9c953b982d61f7e119763c87bbd` |
| `data/processed/final_convergence_summary.csv` | 41 | Supplement convergence summary, grouped maxima without historical per-point ledger | `20260917T094513Z-083977ff` | `ec0a72b74c1ca5a3f12e397774600d1f6337f9657ceb434198deb115a6f6fc2e` |
| `data/processed/final_invariant_subspace_residual.csv` | 3905 | Main Fig. 3; closure and error diagnostics | `20260919T021436Z-40be24e6` | `983971cb2789609610c0d4cd026d675d371ca0d15d2104b72b63dee38e13a754` |
| `data/processed/final_leakage_error_correlations.csv` | 9 | Sec. V.B descriptive rank correlations | `20260919T021436Z-40be24e6` | `629d426a1ba5285ba38dec11c9fea2ab9841efa0f2c863a09199663437d85e46` |
| `data/processed/final_leakage_spectrum.csv` | 2048 | Main Fig. 3d; Fig. S1 signed collision-mode error and weights | `20260919T021436Z-40be24e6` | `71c248c5882ae6463306e7677dde34e9ba9828a212ce33d1336d933840e97fb7` |
| `data/processed/final_leakage_summary.csv` | 4 | Sec. V.B spectral and primal/dual diagnostics | `20260919T021436Z-40be24e6` | `92a4914e3088fb6b0c3c5281175b6ec49fc58d34270dce3052edd19c4245a450` |
| `data/processed/final_mass_scaling.csv` | 64 | Sec. V.E; Fig. S3 and mass table | `20260919T021436Z-40be24e6` | `b0ae701e5c570a35d0a28a73f7bdd0ad00337c8344906e00047563ff0778feea` |
| `data/processed/final_mode_ablation.csv` | 3905 | Main Fig. 2; Sec. V.A and extended ablation table | `20260919T021436Z-40be24e6` | `23d7534de0c539eb589c669f1419f3b1a86f3f5361075495a2ba603583036d7b` |
| `data/processed/final_mode_characterization.csv` | 456 | Main Fig. 1 / Sec. III; physical symmetry and overlap definitions | `20260919T021436Z-40be24e6` | `de8dbb99913279850a56bc8e9c45e7f6042bb4ff27e711d1db8fd280bf8bc2f3` |
| `data/processed/final_mode_decision_agreement.csv` | 12 | Sec. V.A projection-based omission agreement | `20260919T021436Z-40be24e6` | `27cdeb596c8d1eddb3e898989a960b8377b2a5189a34036f56e5ae93ad7b8f37` |
| `data/processed/final_reduced_mode_coefficients.csv` | 4 | Main Fig. 4; Sec. V.C; physical forcing and stiffness | `20260919T021436Z-40be24e6` | `aa748bd68ecd7e2ec1693981770a57362b144542ab60b8793708b4c9dbf7b72e` |
| `data/processed/final_reduced_model_validation.csv` | 72 | Sec. V.A; aggregate projection errors | `20260919T021436Z-40be24e6` | `dafc551c86068c3d9b92c02c58d8000b210c0b45e87951f4e56792788498f4a4` |
| `data/processed/final_tolerance_cancellation.csv` | 10800 | Main Fig. 5; fixed-budget decision map | `20260919T021436Z-40be24e6` | `8a1bce731f934ff299ec166bb0a13064904fe3b990b994041cede12b81f99782` |
| `data/processed/gauge_invariance.csv` | 75 | Supplement basis covariance | `20260917T094513Z-083977ff` | `9dfd9f793f6ae9c7d5e401c95e5f4e5d93d7971c9ee474ba0f0306ddd2cac534` |
| `data/processed/internal_structure_robustness.csv` | 35 | Supplement texture and kernel controls | `20260917T094513Z-083977ff` | `172023f8445fe7a8b184284bb2826410cd1a9093b3d01fb26175a9f3b6efd73d` |
| `data/processed/measurement_level_omission.csv` | 50 | Sec. V.F and Supplement open-circuit algebra; selected 50 rows | `20260917T094513Z-083977ff` | `e4eb6a1d8b97b4d0ca686c507c32974d9d09ea6b5879ce034c8c4d08d01e26db` |
| `data/processed/reduction_cost_benchmark.csv` | 20 | Supplement dense solve / certificate cost assessment | `20260917T094513Z-083977ff` | `c32bb936f9d5562bf39a65e9089d3ddd20d71447d4f5ca684f942541c40614bc` |
| `data/processed/successive_resolution.csv` | 150 | Fig. S4; successive angular-grid convergence | `20260917T094513Z-083977ff` | `60a3e8d66b51078916321d05041f8cfcdafc480481398988fec6f3b4d5bba602` |
| `data/processed/weak_coupling_intervals.csv` | 144 | Supplement sampled weak-coupling passing intervals | `20260917T094513Z-083977ff` | `99b48d226a16a89804a407091a3791dce5abf064ac7b60cc0f1dcd6dcda8f550` |
| `data/processed/weak_coupling_remainder.csv` | 1152 | Fig. S2; weak-coupling remainder controls | `20260917T094513Z-083977ff` | `b6bd33e65bdf93949d38bf143f9b283231a25f704000514bd8e6504bb99576b3` |
| `data/processed/final_physical_validation_summary.csv` | 7 | Sec. IV rate envelope; Supplement validity extrema | `20260917T094513Z-083977ff` | `a04ede177131074a56f0a28c37e72cb189e5a8a19c3e8ef87290cc2766c3004d` |

`validation/final_validation_summary.csv`: concise aggregation of all 25,733 canonical checks; source run `20260919T021436Z-40be24e6`; SHA-256 `82f91efcffbe7c54a4c205b5114f567cc6055769b8fae440265df89b63b38a4a`.

Mass and physical-coefficient tables originated in `20260918T215035Z-4260f56e` and were carried forward byte-for-byte into the canonical mode run. They are independently recomputed by this release. Supporting controls originate in `20260917T094513Z-083977ff`. No external large dataset is required.
