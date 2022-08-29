# photocatalyst-TDDFT-benchmark

The scripts can be used to generate the database from the Gaussian or ORCA outputs located in the subfolders named after the molecules in the benchmark set. It extracts the excitations from the output files of excited state calculations and fits Gaussians using bandwidth and linear scaling factor parameters. The fitting is performed with respect to the reference spectra located in each subdirectory. The details about the process is reported in our paper (DOI: tbd).

## Requirements

Python3 (tested on 3.8+) with standard scientific libraries (numpy, pandas, seaborn, matplotlib, scipy) and palettable. An exported `environment.yml` file for conda is included.

## Usage

Run `run_optimizer.py` to obtain sub-databases then run `merge_all.py` to get the full database.

- By default, the optimizer script uses `max_workers = 37` (37 is the number of molecules) and runs fully parallelized if 37 cores are available.
- The fitting uses a brute-force algorithm to evaluate the fitting errors on a grid. Using default parameters the runtime should be around three hours. To reduce this time the grid resolution (default: 70x90) can be reduced or the number of error functions (default: `[RMSLE, MAE, MSE, r^2, R^2]`) can be decreased.

> **_NOTE:_**  Currently, to change the number of cores (`max_workers`), grid resolution (`SHIFTS`,`SIGMAS`) and the list of error functions (`error_functions`) the `run_optimizer.py` and `bench_opt.py` files need to be modified manually.

## Analysis

The analysis of the database can be done using the Jupyter Notebook located in the for_analysis/ subdirectory.
