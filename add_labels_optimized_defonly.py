#!/usr/bin/env python3

import os

# These strings will be searched in the filenames
molecule = ["NCE_1","NCE_2","NCE_3","NCE_4","PTZ_1","PTZ_2","PTZ_3",
            "PA_1","PA_2","PA_3","PDI_1","PDI_2","PDI_3","POZ_1","POZ_2","POZ_3","POZ_4","POZ_5",
            "BOH-Acr_m","BOH-Acr_o","BOH-Acr_p","BF3-Acr_m","BF3-Acr_o","BF3-Acr_p",
            "Ph-Acr_1","Ph-Acr_2","Me2-Acr_1","Me2-Acr_2","Me2-Acr_3","Mes-Acr_1",
            "CA_1","CA_2","CA_3","Eos_1","Eos_2","Eos_3","Eos_Y","Rh_6G","Rh_B"]

def build_filedata(files,outfile,dumpfile,header,workdir='.'):

    functionals = ["DSD-BLYP","_B2PLYP","wB2PLYP","STEOM-DLPNO-CCSD",
                   "_B3LYP-D3","CAM-B3LYP-D3","M06_","M062X","M06L",
                   "PBE-D3","PBE0","wB97XD","TPSS","B97D3","TDHF"]
    basis = ["_def2TZVP","_TZVP","_def2SVP","_DEF2-TZVP"]
    # Build the database
    print("Building file data...")
    G = {}
    for i in files:
        # extract info from the filenames
        for m in molecule:
            if m in i:
                G[i] = [m]
                break
        try:
            len(G[i]) == 1 # can also test proper conditions an if statement
        except KeyError:
            print("Warning: molecule identification was not found in " + i)
        for f in functionals:
            if f in i:
                if f == "wB97XD": # change to proper functional name
                    G[i].append("$\omega$-B97XD")
                elif f == "wB2PLYP": # change to proper functional name
                    G[i].append("$\omega$-B2PLYP")
                elif f == "_B2PLYP": # avoid conflict with wB2PLYP
                    G[i].append("B2PLYP")
                elif f == "_B3LYP-D3": # avoid conflict with CAM-B3LYP
                    G[i].append("B3LYP-D3")
                elif f == "M06_": # avoid conflict with M062X and M06L
                    G[i].append("M06")
                else:
                    G[i].append(f)
        if len(G[i]) != 2:
            print("Warning: functional was not found in " + i)
        for b in basis:
            if b in i:
                if b == "_def2TZVP" or b == "_DEF2-TZVP":
                    G[i].append("def2-TZVP")
                if b == "_TZVP":
                    G[i].append("TZVP")
                if b == "_def2SVP":
                    G[i].append("def2-SVP")
        if len(G[i]) != 3:
            print("Warning: basis was not found in " + i)
        # extract the optimized parameters (error, bandwidth, shift factor, etc.)
        opt_params = os.path.join(workdir,dumpfile)
        with open(opt_params) as inp:
            for line in inp:
                line = line.strip()
                dat = line.split(',')
                if dat[0] == i:
                    G[i].append(dat[1:])
    print("Building file data done!")
    # G[i] will have this ordering: [molecule, functional, basis, error_function, error_value, bandwidth, shift_factor]
    # default header looks like this: "reference,molecule,functional,basis,error_function,error_value,bandwidth,shift_factor,Wavelength [nm],Intensity [a.u.]\n"
    k = open(outfile, "w")
    k.write(header)
    print("Merging the following files:")
    for filename, data in G.items(): # i = 0,1,2,3
        filename = os.path.join(workdir,filename)
        with open(filename) as inp:
            print(filename)
            for line in inp:
                newline = "no" + "," + data[0] + "," + data[1] + "," + data[2] + "," + data[3][0] + "," + data[3][1] + "," + data[3][2] + "," + data[3][3] + "," + line
                k.write(newline)
    print("Merge complete!\n")
    k.close()

def merge_db(databases,merged_name="database_merged.csv"):
    k = open(merged_name, "w")
    header = ""
    for db in databases:
        firstline = True
        with open(db) as inp:
            for line in inp:
                if firstline == True and header == "":
                    header = line
                    firstline = False
                    k.write(header)
                    continue
                elif firstline == True and line == header:
                    firstline = False
                    continue
                if firstline == True and line != header:
                    raise ValueError("The first lines (headers) in the input databases do not seem to match!")
                elif firstline == False:
                    k.write(line)
    k.close()

def add_reference(ref_file,dbfile):
    molname = "Molecule name was not found in " + ref_file
    for m in molecule:
        if m in ref_file:
            molname = m

    db = open(dbfile, "a")
    with open(ref_file, "r") as ref:
        for line in ref:
            newline = "yes" + "," + molname + "," + "expt." + "," + "expt." + "," + "-" + "," + "-" + "," + "-" + "," + "-" + "," + line
            db.write(newline)
    db.close()


