"""
Parse BOLSIG+ results_Peng.dat and produce a clean CSV with columns
labelled by Peng reaction number (R625-R673).

Usage:
    python parse_bolsig.py

Input:  02_output/results_Peng.dat
Output: 02_output/rates_Peng.csv
"""

import re, csv, sys, os

# ── C-number to Peng reaction mapping ──
# From 02_output/C_to_Peng_mapping.txt
# Entries with None are extra processes not in Peng's table
C_TO_PENG = {
    1: "R625",  2: "R626",  3: None,    4: "R627",  5: "R628",
    6: "R629",  7: "R630",  8: None,    9: "R631", 10: "R633",
   11: "R634", 12: "R636", 13: "R638", 14: "R639", 15: "R640",
   16: "R641", 17: "R642", 18: "R653_C18", 19: "R655", 20: "R656",
   21: "R657", 22: "R643", 23: "R644", 24: None,   25: None,
   26: "R645", 27: "R646", 28: "R647", 29: "R648", 30: "R649",
   31: "R650", 32: "R651", 33: "R652", 34: "R653", 35: "R658",
   36: "R659", 37: "R660", 38: "R661", 39: "R663", 40: "R665",
   41: None,   42: None,   43: "R666", 44: "R667", 45: "R672",
   46: "R668", 47: "R669", 48: "R670", 49: "R671", 50: "R673",
}

# Note: C18 (11.0 eV) may correspond to R654 rather than R653.
# C34 (8.9 eV) is the true R653.  Label C18 as R653_C18 to flag this.

BASE = os.path.dirname(os.path.abspath(__file__))
INPUT  = os.path.join(BASE, "02_output", "results_Peng.dat")
OUTPUT = os.path.join(BASE, "02_output", "rates_Peng.csv")

# ── Parse rate coefficient blocks ──
with open(INPUT, "r") as f:
    lines = f.readlines()

# Find the start of the rate coefficient section (first C-number after
# cross section echo).  Rate blocks look like:
#   C1    H2O    Attachment
#   Energy (eV)\tRate coefficient (m3/s)
#   0.100005\t0.00000
#   ...  (20 rows)
pattern = re.compile(r"^C(\d+)\s+\S+\s+\S+")

blocks = {}   # c_number -> [(energy, rate), ...]
i = 0
in_rate_section = False

while i < len(lines):
    line = lines[i].rstrip()
    # Detect start of rate coefficient section: first "Energy (eV)" header
    if not in_rate_section:
        if "Rate coefficient (m3/s)" in line:
            in_rate_section = True
            # Back up to find the C-number header
            j = i - 1
            while j >= 0 and not pattern.match(lines[j].rstrip()):
                j -= 1
            m = pattern.match(lines[j].rstrip())
            c_num = int(m.group(1))
            # Read data rows
            data = []
            i += 1
            while i < len(lines) and lines[i].strip():
                parts = lines[i].strip().split()
                if len(parts) == 2:
                    energy = float(parts[0])
                    rate = float(parts[1])
                    data.append((energy, rate))
                i += 1
            blocks[c_num] = data
            continue
        i += 1
        continue

    # Already in rate section — look for next C-number header
    m = pattern.match(line)
    if m:
        c_num = int(m.group(1))
        i += 1  # skip "Energy (eV)..." line
        if i < len(lines) and "Rate coefficient" in lines[i]:
            i += 1
        data = []
        while i < len(lines) and lines[i].strip():
            parts = lines[i].strip().split()
            if len(parts) == 2:
                energy = float(parts[0])
                rate = float(parts[1])
                data.append((energy, rate))
            i += 1
        blocks[c_num] = data
        continue
    i += 1

# ── Build CSV ──
# Columns: mean_energy_eV, Te_K, then one column per Peng reaction
# Only include C-numbers that map to a Peng reaction (skip extras)
peng_cols = []  # (c_num, peng_label)
for c in sorted(blocks.keys()):
    label = C_TO_PENG.get(c)
    if label is not None:
        peng_cols.append((c, label))

# Use energy values from the first block (all blocks share the same grid)
energies = [row[0] for row in blocks[peng_cols[0][0]]]
n_points = len(energies)

with open(OUTPUT, "w", newline="") as f:
    writer = csv.writer(f)
    header = ["mean_energy_eV", "Te_eV"] + [label for _, label in peng_cols]
    writer.writerow(header)
    for idx in range(n_points):
        e = energies[idx]
        Te = (2.0 / 3.0) * e   # Te(eV) = (2/3) * mean_energy(eV)
        row = [f"{e:.6f}", f"{Te:.6f}"]
        for c, _ in peng_cols:
            rate = blocks[c][idx][1]
            row.append(f"{rate:.6E}")
        writer.writerow(row)

print(f"Wrote {OUTPUT}")
print(f"  {n_points} energy points, {len(peng_cols)} Peng reactions")
print(f"  Energy range: {energies[0]:.4f} - {energies[-1]:.4f} eV")
print(f"  Columns: {', '.join(label for _, label in peng_cols)}")
