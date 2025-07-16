# airGM

A zero-dimensional plasma chemistry model for humid air discharge at atmospheric pressure. The code tracks 53 species and 624 reactions derived from the work of Sakiyama et al. (2012).

## Directory structure
- `src/` – C++ source (version 2.1)
- `documentation/` – detailed README in Org/HTML and a TODO list
- `results/` – example scripts, output data and gnuplot plotting files

## Build

Requires `make` and a C++11 compiler (tested with GCC g++). From the `src` directory:

```bash
make
sudo make install    # optional
make clean           # remove build artifacts
```

## Usage

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

## Results

Simulation data are saved in `output.csv`. The scripts in `results/` demonstrate batch runs and plotting species densities with gnuplot.

## Documentation

See `documentation/README.org` (or `README.html`) for a comprehensive explanation of the model and full source listings.

## Reference

- Yukinori Sakiyama, David B. Graves, Hung-Wen Chang, Tetsuji Shimizu, and Gregor E. Morfill. "Plasma chemistry model of surface microdischarge in humid air and dynamics of reactive neutral species," *Journal of Physics D: Applied Physics*, **45** (2012) 425201. doi:10.1088/0022-3727/45/42/425201.
