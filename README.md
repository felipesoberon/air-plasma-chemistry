# airGM

A zero-dimensional plasma chemistry model for humid air discharge at atmospheric pressure. The model tracks 53 species and 673 reactions: 624 from Sakiyama et al. (2012) plus 49 electron-impact reactions (R625--R673) from Peng et al. (2022) whose rate coefficients are computed from cross-section data using the BOLSIG+ Boltzmann equation solver.

## Directory structure
- `src/` - C++ source and solver (`airGM2.1`)
- `python/` - Python wrapper and GUI (C++ backend)
- `bolsigplus072024-win/` - BOLSIG+ solver, input cross-section files, and output rate tables
- `documentation/` - detailed README in Org/HTML and related docs
- `results/` - example scripts, output data, and gnuplot files

## C++ version

### Build
Requires `make` and a C++11 compiler (tested with GCC g++).

From `src/`:

```bash
make
sudo make install    # optional
make clean           # remove build artifacts
```

### Run

```bash
./airGM2.1 [options]
```

Supported options:
- `-totaltime <s>` total simulation time
- `-Te <eV>` peak electron temperature
- `-dt <s>` fixed time step (if provided, adaptive timestep is disabled)
- `-RH <percent>` relative humidity in % (preferred humidity input)
- `-[H2O] <m^-3>` water concentration (alternative humidity input)
- `-plasmatime <s>` plasma pulse duration
- `-metricmin <value>` adaptive dt lower relative-change limit (default `0.05`)
- `-metricmax <value>` adaptive dt upper relative-change limit (default `0.10`)

Examples:

Adaptive dt (no `-dt`):

```bash
./airGM2.1 -Te 5.5 -totaltime 1 -RH 50 -metricmin 0.01 -metricmax 0.05
```

Fixed dt:

```bash
./airGM2.1 -Te 5.5 -totaltime 1 -RH 50 -dt 1e-12
```

Notes:
- If `output.csv` exists, the solver resumes from its last data row.
- Species density floor is `10 m^-3` in the current C++ model.

## Python interface

### Requirements
- Python 3.10+
- `matplotlib` (install via `requirements.txt`)
- WSL with C++ build tools (`g++`, `make`) available

### CLI wrapper
From repository root:

```bash
python -m python.main [options]
```

This forwards flags to the C++ solver in WSL (`src/airGM2.1`).

Example:

```bash
python -m python.main -Te 5.5 -totaltime 1 -RH 50
```

### GUI
Install dependencies:

```bash
pip install -r requirements.txt
```

Launch:

```bash
python -m python.gui
```

GUI behavior:
- Uses C++ solver backend via WSL and streams `output.csv` for plotting.
- Inputs shown:
- `Te [eV]`
- `Total time [s]`
- `RH [%]`
- `Metric min`
- `Metric max`
- Default GUI values:
- `Te = 5.5`
- `Total time = 1`
- `RH = 50`
- `Metric min = 0.01`
- `Metric max = 0.05`
- Controls: `Start`, `Stop`, `Reset`
- Plot:
- log scale on both axes
- zoom/pan via Matplotlib toolbar
- per-species visibility checkboxes
- extra checkbox curves: `SUM + ions`, `SUM - ions`
- Live indicators:
- simulation time `t`
- estimated `max rel dC/C` between saved rows
- estimated saved-step `dt` between rows

## Results
Simulation data are written to `output.csv`.

## Documentation
See `documentation/README.org` (or `documentation/README.html`) for detailed model notes.

## BOLSIG+ rate coefficients (Peng reactions)

Reactions R625--R673 use electron-impact rate coefficients that depend on electron temperature rather than analytical formulas. These rates were obtained as follows:

Important runtime scope: the Peng reaction set (R625--R673) is intended to run only during early simulation time, up to `1E-4 s`.

1. **Cross-section data** from the LXCat database (IST-Lisbon and Morgan sets) for e + H2O, e + N2, and e + O2 collisions were used as input to the BOLSIG+ Boltzmann equation solver (v07/2024), available at `https://www.bolsig.laplace.univ-tlse.fr/download.html`.
2. **BOLSIG+** solves the electron energy distribution function (EEDF) for a range of reduced electric fields and outputs rate coefficients as a function of mean electron energy.
3. The resulting rate coefficients (m3/s) vs electron temperature (eV) are stored in `bolsigplus072024-win/02_output/rates_Peng.csv` (40 data points, 44 reaction columns).
4. At runtime, the C++ model **loads this CSV once** and obtains rate coefficients by **linear interpolation** on Te (eV), clamping to the table bounds.

The Peng reactions include dissociative attachment, dissociation, ionization, electronic excitation to tracked excited states (N2(A), N2(B), O2(a1D), O(1D)), and energy-loss-only channels (elastic and vibrational) that consume electron energy but do not change species densities.

## References
- Yukinori Sakiyama, David B. Graves, Hung-Wen Chang, Tetsuji Shimizu, and Gregor E. Morfill. "Plasma chemistry model of surface microdischarge in humid air and dynamics of reactive neutral species," *Journal of Physics D: Applied Physics*, **45** (2012) 425201. doi:10.1088/0022-3727/45/42/425201.
- Peng et al. "Modeling analysis of humid air plasma chemistry in a streamer-like pulsed discharge," *Plasma Science and Technology*, **24** (2022) 055404.
- Hagelaar, G.J.M. and Pitchford, L.C. "Solving the Boltzmann equation to obtain electron transport coefficients and rate coefficients for fluid models," *Plasma Sources Science and Technology*, **14** (2005) 722--733.
