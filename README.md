# photocatalyst-TDDFT-benchmark

The scripts can be used to generate the database from the Gaussian or ORCA outputs located in the subfolders named after the molecules in the benchmark set. The approach is reported in our paper (DOI: tbd, [ChemRxiv](https://doi.org/10.26434/chemrxiv-2022-pfvhc)) and it includes the following steps:

1. excitation energies oscillator strengths are extracted from the output files of excited state calculations 
2. the extracted stick spectrum is shifted and broadened using a sum of Gaussians with bandwidth and linear scaling factor parameters
3. the two parameters are optimized until the calculated lineshape offers the best possible fit (can be evaluated using any error metric) to the experimental reference spectrum located in a given subdirectory
4. the optimized parameters and lineshapes for every molecule/functional/error_metric combination are merged into a database file that can be analyzed using our Jupyter Notebook located in the for_analysis/ subdirectory


## Requirements

Python3 (tested on 3.8+) with standard scientific libraries (numpy, pandas, seaborn, matplotlib, scipy) and palettable. An exported `environment.yml` file for conda is included.

## Usage

Run `run_optimizer.py` to obtain sub-databases then run `merge_all.py` to get the full database.

- By default, the optimizer script uses `max_workers = 37` (37 is the number of molecules) and runs fully parallelized if 37 cores are available.
- The fitting uses a brute-force algorithm to evaluate the fitting errors on a grid. Using default parameters the runtime should be around three hours. To reduce this time the grid resolution (default: 70x90) can be reduced or the number of error functions (default: `[RMSLE, MAE, MSE, r^2, R^2]`) can be decreased.

> **_NOTE:_**  Currently, to change the number of cores (`max_workers`), grid resolution (`SHIFTS`,`SIGMAS`) and the list of error functions (`error_functions`) the `run_optimizer.py` and `bench_opt.py` files need to be modified manually.

## Analysis

The analysis of the database can be performed with the Jupyter Notebook located in the for_analysis/ subdirectory. We recommend using Google Colab for which we have set up default file paths that work in a cloned git repository.
