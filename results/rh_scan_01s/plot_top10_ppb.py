#!/usr/bin/env python3
"""Plot final PPB vs RH for top 10 species across RH scan outputs.

Top-10 selection rule:
- Exclude background gases N2, O2, H2O.
- Rank by highest PPB at the final time point of each RH run.
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt


N2_IDX = 50  # 0-based in output rows
O2_IDX = 51
H2O_IDX = 52
EXCLUDED = {N2_IDX, O2_IDX, H2O_IDX}


def parse_rh_from_filename(path: Path) -> int:
    stem = path.stem  # e.g., run_RH025
    token = stem.split("RH")[-1]
    return int(token)


def read_output_csv(path: Path) -> tuple[list[str], list[float], list[list[float]]]:
    header_species: list[str] = []
    times: list[float] = []
    values: list[list[float]] = []

    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                cols = [c.strip() for c in line[1:].split(",") if c.strip()]
                # Header format: 53 species + Time(s) + StepNo
                header_species = cols[:53]
                continue

            cols = [float(x) for x in line.split(",") if x != ""]
            if len(cols) < 55:
                continue
            times.append(cols[53])
            values.append(cols[:53])

    if not header_species:
        raise ValueError(f"Missing species header in {path}")
    return header_species, times, values


def to_ppb(density: float, n2: float, o2: float, h2o: float) -> float:
    denom = n2 + o2 + h2o
    if denom <= 0.0:
        return math.nan
    return (density / denom) * 1.0e9


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    output_dir = script_dir / "outputs"
    csv_files = sorted(output_dir.glob("run_RH*.csv"))
    if not csv_files:
        raise SystemExit(f"No input files found in {output_dir}. Run run_rh_scan_01s.sh first.")

    per_rh_final_ppb: dict[int, list[float]] = {}
    species_names: list[str] | None = None
    final_max_ppb = [0.0] * 53

    for csv_path in csv_files:
        rh = parse_rh_from_filename(csv_path)
        names, times, rows = read_output_csv(csv_path)
        if species_names is None:
            species_names = names
        elif species_names != names:
            raise SystemExit(f"Species header mismatch in {csv_path}")

        if not rows:
            raise SystemExit(f"No data rows found in {csv_path}")
        last = rows[-1]
        n2_last, o2_last, h2o_last = last[N2_IDX], last[O2_IDX], last[H2O_IDX]
        final_ppb = [to_ppb(last[i], n2_last, o2_last, h2o_last) for i in range(53)]
        per_rh_final_ppb[rh] = final_ppb

        for i in range(53):
            if i in EXCLUDED:
                continue
            ppb_last = final_ppb[i]
            if not math.isnan(ppb_last) and ppb_last > final_max_ppb[i]:
                final_max_ppb[i] = ppb_last

    assert species_names is not None

    ranked = sorted(
        (i for i in range(53) if i not in EXCLUDED),
        key=lambda i: final_max_ppb[i],
        reverse=True,
    )
    top10 = ranked[:10]

    rhs = sorted(per_rh_final_ppb.keys())
    fig, ax = plt.subplots(figsize=(11, 7))
    for species_idx in top10:
        y = [per_rh_final_ppb[rh][species_idx] for rh in rhs]
        ax.plot(rhs, y, marker="o", linewidth=1.8, markersize=6, label=species_names[species_idx])

    ax.set_xlabel("Relative Humidity [%]")
    ax.set_ylabel("Final Concentration [PPB]")
    ax.set_title("Top 10 Species by Final PPB (excluding N2, O2, H2O) | Te=5.0 eV, totaltime=0.1 s")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    fig.tight_layout()

    out_png = script_dir / "top10_species_final_ppb_vs_rh.png"
    out_pdf = script_dir / "top10_species_final_ppb_vs_rh.pdf"
    fig.savefig(out_png, dpi=200)
    fig.savefig(out_pdf)
    print(f"Saved: {out_png}")
    print(f"Saved: {out_pdf}")
    print("Top 10 species:", ", ".join(species_names[i] for i in top10))


if __name__ == "__main__":
    main()
