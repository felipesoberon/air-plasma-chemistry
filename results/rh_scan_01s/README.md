# RH Scan at 0.1 s (Te = 5.0 eV)

## Purpose

Evaluate the effect of relative humidity on the air plasma chemistry model by
running five RH conditions while keeping other inputs fixed.

## Configuration

- `Te = 5.0 eV`
- `totaltime = 0.1 s`
- `RH = 0, 25, 50, 75, 100 %`
- Model chemistry: Sakiyama (`R1-R624`) + Peng (`R625-R673`)

## Files

- `run_rh_scan_01s.sh`
  - Runs all 5 RH cases.
  - Saves each run to `outputs/run_RH000.csv`, `run_RH025.csv`, ..., `run_RH100.csv`.
- `plot_top10_ppb.py`
  - Reads all scan CSV files.
  - Computes final PPB at the end of each run using:
    - `PPB = species_density / (N2 + O2 + H2O) * 1e9`
  - Excludes `N2`, `O2`, `H2O`.
  - Selects top 10 species by highest **final** PPB across RH runs.
  - Produces one combined plot: RH on x-axis, final PPB on y-axis, with
    scatter markers and connecting lines for each species, plus legend.
- `outputs/`
  - Storage folder for run CSV outputs.

## Usage

From repository root:

```bash
bash results/rh_scan_01s/run_rh_scan_01s.sh
python results/rh_scan_01s/plot_top10_ppb.py
```

## Generated Figures

- `results/rh_scan_01s/top10_species_final_ppb_vs_rh.png`
- `results/rh_scan_01s/top10_species_final_ppb_vs_rh.pdf`

## Notes

- The solver resumes from `output.csv` if it exists. The run script removes
  root `output.csv` before each case so each RH run starts clean.
