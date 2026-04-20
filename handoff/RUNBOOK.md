# Runbook

## Build
In project root:

```powershell
make
```

## Programs

### 1) POC single-pair style exploration
```powershell
./kl_structure_poc 4
./kl_structure_poc 5
```

### 2) Bulk generation (baseline)
```powershell
./kl_structure_bulk 4
./kl_structure_bulk 5
```

### 3) Bulk generation (optimized)
```powershell
./kl_structure_bulk_opt 4
./kl_structure_bulk_opt 5
./kl_structure_bulk_opt 6
```

Default outputs from optimized run:
- `S<n>_structure_constants_bulk_opt.txt`
- `S<n>_structure_constants_bulk_opt.csv`

### 4) Bulk generation (optimized, parallel x-row)
```powershell
./kl_structure_bulk_opt_parallel 5
./kl_structure_bulk_opt_parallel 6 out.csv
./kl_structure_bulk_opt_parallel 6 out.csv --with-txt
./kl_structure_bulk_opt_parallel 6 out.csv --with-txt out.txt --bench
./kl_structure_bulk_opt_parallel 6 --threads 8 --row-schedule weighted
./kl_structure_bulk_opt_parallel 6 --threads 12 --row-schedule block
```

Parallel arguments:
- `--with-txt [path]` enables TXT output in addition to CSV (path optional)
- `--threads N` sets OpenMP thread count (default: all available)
- `--row-schedule block|weighted`
	- `weighted` (default) balances heavier late rows better
	- `block` can be better at high thread counts on smaller `n`
- `--bench` prints per-phase timing report in parallel mode
	- Note: phase wall times are aggregate thread-time and can exceed total wall time

Default outputs from parallel run:
- `S<n>_structure_constants_bulk_opt_parallel.csv`

Default TXT output when enabled:
- `S<n>_structure_constants_bulk_opt_parallel.txt`

### 5) Runtime sampling estimator (no output tables)
```powershell
./kl_runtime_sample_estimate 6 --sample-rows 24 --threads 8 --row-schedule weighted
./kl_runtime_sample_estimate 7 --sample-rows 96 --threads 8 --row-schedule weighted --baseline-seconds 841 --baseline-n 6
```

Estimator arguments:
- `--sample-rows K`: number of sampled x-rows
- `--threads N`: worker threads
- `--row-schedule block|weighted`: row sampling mode (`weighted` default)
- `--baseline-seconds T --baseline-n M`: optional fit of exponent `a` in `T ~ (n!)^a` and predictions for `S(n+1), S(n+2)`

### 6) Validation
```powershell
./kl_structure_validate 4 S4_structure_constants_bulk_opt.csv
./kl_structure_validate 5 S5_structure_constants_bulk.csv
./kl_structure_validate 6 S6_structure_constants_bulk_opt.csv 5000 0
./kl_structure_validate 6 S6_structure_constants_bulk_opt_parallel.csv --triple 694 688 324
./kl_structure_validate 6 S6_structure_constants_bulk_opt_parallel.csv --triple-legacy 694 688 324
```

Validator arguments:
- arg1: `n`
- arg2: csv path
- arg3 (optional): inversion sample count
- arg4 (optional): associativity sample count (`0` skips)
- `--triple x y z`: single-triple associativity check using the fast targeted-scan path
- `--triple-legacy x y z`: single-triple check using the full index-build path (debug/reference)

Validator note:
- For large files (for example S6 parallel CSV), `--triple` is much faster than the legacy indexed path.

## Typical performance notes observed this session
- Optimized bulk path gives substantial speedup vs baseline on S5.
- S6 optimized full run has already completed in this workspace and produced files.
- Parallel x-row path is confirmed correct against sequential outputs (CSV exact match; TXT body match ignoring header banner).
- For S5, weighted row schedule improves low-to-mid thread counts (not always best at very high threads).
- Observed full S6 runtime on this machine/settings: about 841 seconds.
- Sampling estimator exists but currently can overestimate for small sample sizes on S6; use larger sample sizes for forecasting.

## Cleanup
```powershell
make clean
# if a previous parallel run was interrupted, remove leftover temp directories:
Get-ChildItem -Directory -Filter "tmp_kl_parallel_*" | Remove-Item -Recurse -Force
```

Note: generated `*.txt` and `*.csv` outputs are ignored by `.gitignore`.
