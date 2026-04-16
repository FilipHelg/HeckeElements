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

### 4) Validation
```powershell
./kl_structure_validate 4 S4_structure_constants_bulk_opt.csv
./kl_structure_validate 5 S5_structure_constants_bulk.csv
./kl_structure_validate 6 S6_structure_constants_bulk_opt.csv 5000 0
```

Validator arguments:
- arg1: `n`
- arg2: csv path
- arg3 (optional): inversion sample count
- arg4 (optional): associativity sample count (`0` skips)

## Typical performance notes observed this session
- Optimized bulk path gives substantial speedup vs baseline on S5.
- S6 optimized full run has already completed in this workspace and produced files.

## Cleanup
```powershell
make clean
```

Note: generated `*.txt` and `*.csv` outputs are ignored by `.gitignore`.
