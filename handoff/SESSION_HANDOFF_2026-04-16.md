# Session Handoff (2026-04-16)

## Scope completed
- Added standalone KL product proof-of-concept program: `kl_structure_poc.c`.
- Added bulk generator for all products: `kl_structure_bulk.c`.
- Added optimized bulk generator: `kl_structure_bulk_opt.c`.
- Added standalone validator: `kl_structure_validate.c`.
- Added/updated build targets in `makefile`.
- Generated structure-constant outputs for S4/S5/S6 in text and CSV forms.

## Main deliverables in workspace
- Source files:
  - `kl_structure_poc.c`
  - `kl_structure_bulk.c`
  - `kl_structure_bulk_opt.c`
  - `kl_structure_validate.c`
- Generated outputs (key):
  - `S4_structure_constants_bulk_opt.txt`
  - `S4_structure_constants_bulk_opt.csv`
  - `S5_structure_constants_bulk.txt`
  - `S5_structure_constants_bulk.csv`
  - `S6_structure_constants_bulk_opt.txt`
  - `S6_structure_constants_bulk_opt.csv`

## Validation status
- Identity check in validator: passes.
- Inversion-symmetry and associativity checks under naive labeling assumptions: fail frequently.
- Deep investigation indicates this is likely a convention mismatch, not random data corruption.

## Important convention note
Observed behavior in produced CSV labels:
- `e * y = y` (left identity in displayed labels)
- `x * e = x^{-1}` (right identity appears inverted in displayed labels)

This means the currently emitted labels are not in the standard two-sided unit convention for structure constants, so direct associativity checks on displayed labels can fail even when internal algebra computation is coherent.

## Current repo state notes
- `.gitignore` now includes `*.csv` in addition to `*.txt`, `*.exe`, and `*.o`.
- Existing generated outputs and binaries are present in the workspace.

## Recommended handoff decision
Choose one path before final publication of tables:
1. Keep current generator output convention and update validator semantics accordingly.
2. Normalize generator output to standard structure-constant labeling (recommended for downstream use and publication).
3. Support both conventions with an explicit switch and document both output schemas.

## Performance session update (same date)

### What was added
- Added optional benchmark mode in `kl_structure_bulk_opt.c` via `--bench`.
- Benchmark mode now reports wall time, CPU time, and per-phase timing for:
  - build Cx/Cy
  - multiply
  - decompose
  - write TXT
  - write CSV
- Added 1 MiB buffered output for TXT/CSV streams using `setvbuf`.

### Multiply-path optimizations implemented
- Kept specialized lower-term multiply path (`v^-1 - v`) in the simple-generator correction path.
- Reworked sparse product accumulation to use a dense stamped accumulator (reduces repeated sparse map churn).
- Added bounded left-multiply cache for reused `(rightId, leftIndex)` expansions.
- Wired these into `MultiplyHeckeSparse` and main compute loop allocation/free lifecycle.

### Results observed
- `--bench` runs on S5 show multiply as the dominant phase by far.
- Representative non-PGO result after serial optimizations:
  - wall: about 3.16 s (NUL output benchmark style run)
  - multiply phase: about 2.66 s
- Experimental PGO run (temporary) improved further to about 2.44 s wall, but training cost was significant.

### Final build workflow decision
- PGO targets were tested and then removed from `makefile` for a minimal workflow appropriate for one-off large runs.
- Current recommended path is the normal optimized build (`-Ofast -march=native`) with optional `--bench` for diagnostics.

### Current integrity checkpoint
- Optimization markers confirmed present in `kl_structure_bulk_opt.c`:
  - `DenseAccum`
  - `LeftMulCache`
  - `PolyMultiplyLowerTerm`
  - `MultiplyHeckeSparse` integration
  - `--bench` CLI support
- No PGO-related targets/flags remain in `makefile`.

## Update (2026-04-17)

### Additional scope completed
- Added parallel optimized generator: `kl_structure_bulk_opt_parallel.c`.
- Added weighted row scheduling support to mitigate heavier late x-rows:
  - `--row-schedule block|weighted` (default `weighted`).
- Added runtime sampling estimator (separate tool): `kl_runtime_sample_estimate.c`.
- Added build targets in `makefile` for:
  - `kl_structure_bulk_opt_parallel`
  - `kl_runtime_sample_estimate`

### Parallel implementation notes
- Parallelization is by x-row.
- Workspaces and caches are thread-local in the main compute loop.
- Output strategy uses per-thread temporary files and deterministic merge.
- Correctness was checked against sequential output:
  - CSV exact match
  - TXT body match (header banner differs)

### Performance observations
- S5:
  - weighted scheduling improves low-to-mid thread counts versus block.
  - at very high threads, block may be competitive or better on small n.
- S6:
  - reported full runtime in this workspace/settings: about `841s`.

### Runtime estimator notes
- `kl_runtime_sample_estimate` performs sampled x-row evaluation (full y-loop per sampled row), then extrapolates full runtime.
- Supports optional baseline fitting for exponent `a` in `T ~ (n!)^a` and predicts `S(n+1), S(n+2)`.
- Early check indicates small sample sizes can overestimate S6 runtime; calibration is pending.

### Recommended immediate next steps
1. Calibrate estimator on S6 using multiple sample sizes (`K=24,48,96,192`) and compute correction/confidence bands.
2. Use calibrated estimator to forecast S7/S8 runtime with uncertainty.
3. Keep `weighted` as default row schedule; add `auto` schedule selection if needed.

## Update (2026-04-17, later)

### Parallel generator changes
- `kl_structure_bulk_opt_parallel.c` now defaults to CSV output only.
- Added optional TXT output switch:
  - `--with-txt` (default txt path)
  - `--with-txt out.txt` (explicit txt path)
- Legacy positional mode is still accepted for compatibility:
  - `./kl_structure_bulk_opt_parallel n out.txt out.csv`

### Bench/reporting change
- Parallel `--bench` now reports per-phase timing (build Cx/Cy, multiply, decompose, write TXT, write CSV).
- Report uses aggregate thread-time in parallel mode; phase sums can exceed total wall time.

### Cleanup/refactor change
- Reduced repeated error cleanup loops by introducing a shared helper for removing per-thread temp files.

### Additional optimization attempt
- Updated left-multiply cache lookup from full linear scan to bounded hashed probe set.
- Intended effect: lower lookup overhead in hot multiply path while preserving bounded memory.

### Measured runtime snapshots (S5, weighted, 8 threads)
- Baseline before change (TXT+CSV): about `2.4s` wall.
- After changes, TXT+CSV (`--with-txt`): about `2.5s` wall.
- After changes, CSV-only default: about `2.1s` wall.

Interpretation:
- The new default output mode provides meaningful end-to-end speedup by avoiding TXT writes when not needed.
- Kernel/cache changes did not improve the TXT+CSV path in this quick S5 check; additional profiling is recommended before keeping or extending cache-policy changes.

## Update (2026-04-20)

### Multiply-path profiling and cache decision
- Added fine-grained multiply instrumentation in `kl_structure_bulk_opt_parallel.c` (cache lookup/build, polynomial multiply, accumulation, commit).
- Profiling on S5 (`n=5`, 4 threads, `--bench`) showed no practical cache reuse in the tested configuration.
- Cache hit/miss counters observed effectively zero hits.
- Removed cache lookup/build from `MultiplyHeckeSparse` hot path and switched to direct `LeftMultiplyByIndex(...)`.

### Timing interpretation update
- Confirmed that parallel CPU-time reports can significantly exceed wall time because they aggregate work across threads.
- Confirmed that heavy internal timers in `--bench` materially increase wall time.
- Non-bench run remains the reliable mode for wall-time comparison across implementation variants.

### Current runtime snapshot (S5)
- With bench instrumentation and cache removed: approximately `3.0s` wall at 4 threads.
- Without bench instrumentation: approximately `1.0s` wall at 4 threads.

### Workspace cleanup completed
- Removed temporary benchmark stderr artifacts (`*.err`) generated during profiling sweeps.
- No temporary `tmp_kl_parallel_*` directories remained after cleanup.

## Update (2026-04-20, validator follow-up)

### Validator diagnostics improvement
- Enhanced `kl_structure_validate.c` to print concrete failing cases (not only counts):
  - identity counterexamples (expected vs actual expansion)
  - inversion symmetry mismatches (specific mirrored pair and coefficient mismatch)
  - associativity counterexamples (specific `(x,y,z)` triple and first mismatching coefficient)

### Validator semantic fix
- Corrected validator product-key usage to match CSV semantics:
  - right identity now checks `x * e` directly
  - associativity expansion paths now use direct pair keys instead of inverse-indexed keys

### Validation results after fix
- `S4_structure_constants_bulk_opt_parallel.csv` (freshly generated) validates with:
  - identity: OK
  - inversion symmetry sampling: completed
  - associativity sampling: OK
- `S5_structure_constants_bulk_opt_parallel.csv` initially showed right-identity failures on an older file snapshot.
- After regenerating S5 with current binary, validator reports:
  - identity: OK
  - inversion symmetry sampling: completed
  - associativity sampling: OK

### Practical note for future runs
- After any generator code changes, regenerate target `S<n>_structure_constants_bulk_opt_parallel.csv` before validation to avoid stale-file false alarms.

## Update (2026-04-20, triple-check latency and spot checks)

### Validator single-triple path
- `kl_structure_validate.c` now supports a fast single-triple path for `--triple x y z`.
- Legacy indexed behavior remains available as `--triple-legacy x y z` for comparison/debug.

### Observed S6 runtime difference (same triple)
- Fast path (`--triple`): approximately `9.4s` wall.
- Legacy indexed path (`--triple-legacy`): approximately `114.0s` wall.

### Recheck of previously problematic S6 triples
Re-ran the three triples against `S6_structure_constants_bulk_opt_parallel.csv`:
1. `(694, 688, 324)`
2. `(477, 690, 664)`
3. `(591, 693, 691)`

Current result: all three report associativity holds.
