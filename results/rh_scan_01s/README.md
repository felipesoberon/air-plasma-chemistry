# RH Scan at 0.1 s (Te = 4.8 eV)

## Purpose

Evaluate the effect of relative humidity on the air plasma chemistry model by
running 11 RH conditions while keeping other inputs fixed.

## Configuration

- `Te = 4.8 eV`
- `totaltime = 0.1 s`
- `RH = 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 %`
- Model chemistry: Sakiyama (`R1-R624`) + Peng (`R625-R673`)

## Files

- `run_rh_scan_01s.sh`
  - Runs all 11 RH cases.
  - Saves each run to `outputs/run_RH000.csv`, `run_RH010.csv`, ..., `run_RH100.csv`.
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

## Observations (Current Top-10 Figure)

Based on `top10_species_final_ppb_vs_rh` from the current scan:

- `O3` increases steadily with RH, from about `1.69 ppb` at `0% RH` to about
  `2.23 ppb` at `100% RH`.
- `OH` also increases strongly with RH, from essentially `0 ppb` in dry air to
  about `0.223 ppb` at `100% RH`.
- Relative to `O3` and `OH`, most other top species are comparatively flat
  across RH. Small changes are present, but the RH sensitivity is weaker.
- `NO` shows a modest upward trend (`~0.142 -> ~0.192 ppb`) with RH.
- Some water-related products (`HO2`, `H2O2`, `HNO2`) increase with RH, but
  remain at lower absolute PPB than `O3`, `OH`, and `NO`.
