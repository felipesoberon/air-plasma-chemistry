# airGM

A zero-dimensional plasma chemistry model for humid air discharge at atmospheric pressure. The model tracks 53 species and 624 reactions derived from Sakiyama et al. (2012).

## Directory structure
- `src/` - C++ source (version 2.1)
- `python/` - Python port of the model
- `documentation/` - detailed README in Org/HTML and a TODO list
- `results/` - example scripts, output data and gnuplot plotting files

## C++ version

### Build

Requires `make` and a C++11 compiler (tested with GCC g++). From the `src` directory:

```bash
make
sudo make install    # optional
make clean           # remove build artifacts
```

### Run

Run the executable:

```bash
./airGM2.1 [options]
```

Supported options:

- `-totaltime <s>`    total simulation time
- `-Te <eV>`          peak electron temperature
- `-dt <s>`           time step
- `-[H2O] <m^-3>`     water concentration
- `-plasmatime <s>`   plasma pulse duration

Example:

```bash
./airGM2.1 -Te 2.45 -totaltime 1E-5
```

If `output.csv` exists, the simulation resumes from it; otherwise it begins with default values (1 ms, 2.6 eV, 50 ps, [H2O] = 1.2E24 m^-3).

## Python version

### Requirements

- Python 3.10+
- `matplotlib` (install via `requirements.txt` for GUI support)

### Run

From the repository root:

```bash
python -m python.main [options]
```

The Python CLI matches the C++ flags:

- `-totaltime <s>`    total simulation time
- `-Te <eV>`          peak electron temperature
- `-dt <s>`           time step
- `-[H2O] <m^-3>`     water concentration
- `-plasmatime <s>`   plasma pulse duration

Example:

```bash
python -m python.main -Te 2.45 -totaltime 1E-5
```

The Python version reads and appends `output.csv` with the same schema and column order as C++.

### GUI

Install dependencies:

```bash
pip install -r requirements.txt
```

Launch the GUI:

```bash
python -m python.gui
```

GUI features:
- Two-pane layout:
- Left pane: model inputs and controls
- Right pane: live species concentration plot
- Inputs exposed in GUI:
- `Te [eV]`
- `Total time [s]`
- `RH [%]` (H2O density is computed internally at fixed `T=298 K`, `P=101325 Pa`)
- Adaptive timestep limits:
- `Metric min` (default `0.05`)
- `Metric max` (default `0.5`)
- Runtime controls:
- `Start`, `Stop`, `Reset`
- Live indicators:
- simulation time `t`
- max relative concentration change (`max rel dC/C`)
- current timestep `dt`
- Plot behavior:
- log scale on both axes
- zoom/pan via Matplotlib toolbar
- species visibility via `Show all species` and per-species checkboxes

## Results

Simulation data are saved in `output.csv`. The scripts in `results/` demonstrate batch runs and plotting species densities with gnuplot.

## Documentation

See `documentation/README.org` (or `README.html`) for a comprehensive explanation of the model and full source listings.

## Reference

- Yukinori Sakiyama, David B. Graves, Hung-Wen Chang, Tetsuji Shimizu, and Gregor E. Morfill. "Plasma chemistry model of surface microdischarge in humid air and dynamics of reactive neutral species," *Journal of Physics D: Applied Physics*, **45** (2012) 425201. doi:10.1088/0022-3727/45/42/425201.
