# Independent release verification

Executed on 2026-09-28. The release was copied to a separate temporary directory outside the scientific workspace, given its own clean virtual environment, installed from the pinned requirements and installed as an editable package. Its source, tests, scripts and configuration hashes match the release. No original-project path was supplied to any calculation or build. The tested Python version is 3.12.14; core numerical versions match the canonical environment.

| Check | Result |
|---|---|
| Package installation from the isolated copy | PASS |
| `make test` | 40 passed |
| `make validate` | 23 dataset hashes; 46 fresh independent physics checks; N=256 finite-range M1–M4 regression; fixed-budget sign/budget checks passed |
| `make reproduce-core` | PASS; 3,600 principal points, 3,510 positive; 22,133 projection/identity gates; 13 central tables recomputed |
| Numerical comparison | 786,294 cells checked; maximum numeric difference 0 in the pinned environment; categorical fields identical |
| Byte comparison | 12/13 reconstructed tables byte-identical; physical coefficient table differs only in row order, with all keyed values identical |
| Main figures versus original canonical figures | All five PNGs byte-identical |
| Independent figure build versus release | All nine PNGs byte-identical |
| `make figures` | Five main and four supplementary figures, tables and claim ledger regenerated |
| `make paper` | Main 11 pages; Supplement 12 pages; no overfull boxes or unresolved references/citations |
| `make clean`, followed by `make validate` | PASS; canonical inputs and final PDFs retained |
| Visual inspection | Every final main and supplemental page inspected; one orphan reference fixed; no clipped equations, overlapping labels or missing glyphs |

The isolated core reconstruction took **24.412 seconds** with one BLAS thread; dependency installation, font-cache preparation and TeX builds are excluded. Timing is descriptive, not a performance guarantee. `release_verification.json` contains the per-table comparison, exact figure hashes and environment.

The scientific source was read only. All 57 distinct selected source-file hashes were rechecked. The complete 10,676-file original path set and file sizes were unchanged. Read-triggered cloud hydration rounded four timestamps by less than one microsecond; their content hashes still match. No original source, dataset, plot or manuscript was written, deleted or reorganized. No BH/NS directory was entered or used, and no Git metadata was copied.

The retained input manifest is intentionally much smaller than the prior complete working package. Historical runs, full raw matrices, long-form duplicate ablation, legacy control scans, old figures, archives and development logs are excluded. Their absence was tested by the independent execution above. Supporting figures use included validated tables; the core target does not rerun every historical supporting scan.

The user authorized public GitHub publication on 2026-09-28. The separate repository `03102000369/nonlinear-hall-valley-omission` is public, uses `main` as its default branch, and received the validated scientific commits. Human authorship, target-journal requirements, archival-deposition metadata and nonsource licensing still require author decisions, as listed in `HUMAN_TODO.md`.

The Git staging audit covered all 106 intended release files and verified that each staged object preserved the working file bytes. No historical result directories, symlinks, oversized files, inherited history or remotes were present. The largest payload file is 7,534,449 bytes. Git attributes disable automatic line-ending conversion so checksum-verified inputs and generators retain their exact bytes on checkout. CSV CRLF and Markdown hard-break whitespace were intentionally preserved.
