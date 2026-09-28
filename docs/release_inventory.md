# Release inventory and selection policy

This classification was written before copying any project content. The source project is retained read-only. No other research project is a source.

| Classification | Files / families | Decision |
|---|---|---|
| ESSENTIAL | core, referee_numerics, final_physics, mode_ablation, io; current physical figures; final mode/configuration tables; main and supplement TeX | Include one copy; refactor orchestration and paths only. |
| SUPPORTING | independent manuscript_numerics and validation solvers; five physics test modules; mass, texture, gauge, weak-series, certificate, convergence and open-circuit controls | Include only controls cited in the final paper, compact convergence summaries and selected measurement rows. |
| LOCAL/HISTORICAL | results/, old release/, notebooks/, prior revision plans/replies/audits; old source snapshots and old run namespaces | Keep exclusively in the original working project. Preserve relevant source hashes and run IDs in release provenance. |
| REDUNDANT | long-form ablation (derivable from wide table), processed mirrors, duplicate manuscript supplement PDF, numbered duplicate source files, superseded figures, unused TeX macros | Exclude. |
| LARGE/EXCLUDED | dense arrays, full per-point validation ledger, all historical maps and backups, environments and caches | Exclude; no large external dataset is required by the curated study. |

The descriptive regime map, earlier scalar-regression comparison and extra forcing plot are omitted with their manuscript discussion. The retained physical hierarchy directly answers the central question. Selected supplementary figures are spectrum, weak-series control, mass and successive convergence. Final concrete file and size inventories are in final_release_report.md.
