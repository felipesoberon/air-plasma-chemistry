# Results Folder Notes

All subfolders in `results/` except `results/rh_scan_01s/` were produced with
an older model version that had known bugs and did not include the Peng
electron-impact reactions (`R625-R673`).

Those legacy subfolders are kept only to preserve test setup context
(scripts/plot templates), not as trusted numerical results.

## Legacy Subfolder Purpose

- `results/dry/`
  - Early dry-air baseline run (`[H2O]=0`) using older executable (`airGM.x`),
    staged over increasing total times.
- `results/dry2/`
  - Alternate dry-air baseline run (`[H2O]=0`) using older executable
    (`airGM2`) and a similar staged-time approach.
- `results/humid/`
  - Early humid-air baseline run using older executable (`airGM.x`) at fixed
    low `Te` with staged total times.
- `results/dry_scanTe/`
  - Dry-air (`[H2O]=0`) scan over electron temperature (`Te`) with staged-time
    runs, plus summary plotting scripts.
- `results/humid_scanTe/`
  - Humid-air scan over electron temperature (`Te`) with staged-time runs,
    plus summary plotting scripts.
- `results/dry_scanNO2/`
  - Dry-air (`[H2O]=0`) scan over initial `NO2` values at fixed `Te`
    (`5.75 eV` and `6 eV` cases in subfolders).
- `results/humid_scanNO2/`
  - Humid-air scan over initial `NO2` values at fixed `Te`
    (`5.75 eV` and `6 eV` cases in subfolders).

## Current Trusted Scan

Use `results/rh_scan_01s/` for the current workflow and outputs.
