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
