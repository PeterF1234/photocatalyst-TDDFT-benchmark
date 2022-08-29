#!/usr/bin/env python3

import os
import bench_opt as BO
import add_labels_optimized_defonly as labels
# import matplotlib.pyplot as plt
import numpy as np
import sys
import concurrent.futures
import time
from multiprocessing import cpu_count

def absoluteFilePaths(directory): # just in case we need it
    for dirpath,_,filenames in os.walk(directory):
        for f in filenames:
            yield os.path.abspath(os.path.join(dirpath, f))

def opt_spect_molecule(folder_path): # folder path example: '/home/peter/Munka/bechmark/v1.0_multicore/PDI_2'
    start = time.perf_counter() # start timing
    molname = os.path.basename(folder_path)
    files = os.listdir(folder_path)
    files = [os.path.join(folder_path, i) for i in files]
    ref_csv = molname + ".csv" # the program expects that every folder contains a reference spectrum file, e.g., PA_1.csv in the PA_1 folder
    ref_csv = os.path.join(folder_path, ref_csv)
    header = "reference,molecule,functional,basis,error_function,error_value,bandwidth,shift_factor,Wavelength [nm],Intensity [a.u.]\n"
    error_functions = [BO.rmsle, BO.mae, BO.mse, BO.r_square, BO.R_square] # these will be the error functions

    databases = [] # each list element will contain a data file for 15 "functionals" and 1 error metric
    for i in error_functions:
        errfn_name = i.__name__
        error_dir = os.path.join(folder_path,errfn_name)
        if not os.path.exists(error_dir):
            print(f'The program was unable to find the following folder: {error_dir}\n')
            sys.exit(4)
        dbname = ref_csv[:-4] + "_database_" + errfn_name + ".csv" # ref_csv has full path included
        databases.append(dbname)

    # Merge all results into a single file per molecule (i.e., 15 "functionals", 5 error metric, and the reference data)
    finaldb_name = "database_merged.csv"
    finaldb_name = os.path.join(folder_path, finaldb_name)
    labels.merge_db(databases,merged_name=finaldb_name)
    ref_filename = ref_csv[:-4] + "_integerx.csv" # BO.save_expt creates a file with this name
    labels.add_reference(ref_filename, finaldb_name)

    # finish timing
    finish = time.perf_counter()
    print(f'Finished with database building for molecule {molname} in *** {round(finish-start, 2)} second(s) ***!')


def main():

    molecules = labels.molecule # The list of molecules is defined in the add_labels_optimized_defonly.py file

    pth = os.getcwd() # Change this if running from another folder
    folders = []
    for i in os.listdir(pth):
        if i in molecules:
            folders.append(i)
    print('*** Spectrum optimizer program initialized! ***')
    print("Job started at " + time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + '\n')
    print(f'The following folders will be processed:\n{folders}\n')
    #sys.exit(1)
    
    # max_workers=cpu_count()
    # max_workers sould be set to the number of folders is enough CPU is available
    print("number of processors available: " + str(cpu_count()) + "\n")
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:
        e = executor.map(opt_spect_molecule, folders)
        for i in e:
            print(i)

    print("Job finished at " + time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + '\n')


if __name__ ==  '__main__':
    main()

