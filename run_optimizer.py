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
    start_molecule_timing = time.perf_counter()
    molname = os.path.basename(folder_path)
    print(f'Optimizing the spectra of molecule {molname}...')
    lambdarange = np.arange(200,700)
    SIGMA = 0.2
    lambda_shift = 1
    ref_pad = 50
    files = os.listdir(folder_path)
    files = [os.path.join(folder_path, i) for i in files]
    ref_csv = molname + ".csv" # the program expects that every folder contains a reference spectrum file, e.g., PA_1.csv in the PA_1 folder
    ref_csv = os.path.join(folder_path, ref_csv)
    BO.save_expt(ref_csv,ref_pad=50) # reformat the reference spectrum for the final database
    header = "reference,molecule,functional,basis,error_function,error_value,bandwidth,shift_factor,Wavelength [nm],Intensity [a.u.]\n"
    error_functions = [BO.rmsle, BO.mae, BO.mse, BO.r_square, BO.R2]

    databases = [] # each list element will be a data file for 15 "functionals" and 1 error metric (filename includes full path)
    for i in error_functions:
        start = time.perf_counter() # timing for each error metric calculaton
        errfn_name = i.__name__
        error_dir = os.path.join(folder_path,errfn_name)
        if not os.path.exists(error_dir):
            os.mkdir(error_dir)
        dumpfile_name = 'optimized_params_'+errfn_name+'.csv' # no need to add folder_path here, it is done by BO.print_result()
        # optimize the spectra
        print(f'process for {errfn_name} has started for molecule {molname}')
        for logfile in files:
            if logfile[-4:] == '.log' or logfile[-4:] == '.out':
                print(f'Working on file: {logfile}')
                print(f'Ref csv: {ref_csv}')
                print(f'Error directory: {error_dir}')
                BO.print_result(ref_csv,logfile,lambdarange,SIGMA,lambda_shift,ref_pad,dumpfile=dumpfile_name,writeopt=True,errorfunc=i,workdir=error_dir,logscale_hmap=True)
            # sys.exit(1)

        # build the database
        print(f"\nMerging results for {molname} and error function {errfn_name}")
        opt_csv_files = [] # one file for each functional (per molecule and error metric) 
        outfile = ref_csv[:-4] + "_database_" + errfn_name + ".csv" # merge all functionals into one file per error metric
        for i in os.listdir(error_dir):
            if i[-4:] == ".csv" and i != dumpfile_name:
                opt_csv_files.append(i)
            else:
                print(f'Skipping file {i} ...')
        labels.build_filedata(opt_csv_files,outfile,dumpfile_name,header,workdir=error_dir)
        finish = time.perf_counter() # finish timing
        print(f'Finished with {errfn_name} for molecule {molname} in *** {round(finish-start, 2)} second(s) ***!')
        databases.append(outfile)

    # Merge all results into a single file per molecule (i.e., 15 "functionals", 5 error metric, and the reference data)
    finaldb_name = "database_merged.csv"
    finaldb_name = os.path.join(folder_path, finaldb_name)
    labels.merge_db(databases,merged_name=finaldb_name)
    ref_filename = ref_csv[:-4] + "_integerx.csv" # BO.save_expt() creates a file with this name
    labels.add_reference(ref_filename, finaldb_name)

    # finish timing
    finish_molecule_timing = time.perf_counter()
    print(f'\nDatabase built for molecule {molname} in *** {round(finish_molecule_timing-start_molecule_timing, 2)} second(s) ***!\n')


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
    N_folders = len(folders )
    print(f'The number of folders to be processed: {N_folders}\n')
    # sys.exit(1)
    
    # max_workers=cpu_count()
    # max_workers sould be set to the number of folders is enough CPU is available
    print("number of processors available: " + str(cpu_count()) + "\n")
    with concurrent.futures.ProcessPoolExecutor(max_workers=37) as executor:
        e = executor.map(opt_spect_molecule, folders)
        for i in e:
            print(i)

    print("Job finished at " + time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + '\n')


if __name__ ==  '__main__':
    main()

