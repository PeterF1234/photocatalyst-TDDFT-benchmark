#!/usr/bin/env python3

import os


def merge_fulldb(db_names="database_merged.csv",pth='.'):
    molecules = ["NCE_1","NCE_2","NCE_3","NCE_4","PTZ_1","PTZ_2","PTZ_3",
                 "PA_1","PA_2","PA_3","PDI_1","PDI_2","PDI_3","POZ_1","POZ_2","POZ_3","POZ_4","POZ_5",
                 "BOH-Acr_m","BOH-Acr_o","BOH-Acr_p","BF3-Acr_m","BF3-Acr_o","BF3-Acr_p",
                 "Ph-Acr_1","Ph-Acr_2","Me2-Acr_1","Me2-Acr_2","Me2-Acr_3","Mes-Acr_1",
                 "CA_1","CA_2","CA_3","Eos_1","Eos_2","Eos_3","Eos_Y","Rh_6G","Rh_B"]
    fulldb_name = "OPC_database.csv"
    folders = next(os.walk(pth))[1]
    k = open(fulldb_name, "w")
    header = ""
    for i in folders:
        if i not in molecules:
            print("Skipping folder \"" + i + "\" because it does not contain a recognized molecule name.")
            continue
        else:
            db_fullpath = os.path.join(i,db_names) # all databases should have the same name
            firstline = True
            with open(db_fullpath) as inp:
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
    print("Database has been built!")

def main():
    merge_fulldb(db_names="database_merged.csv",pth='.')

if __name__ == "__main__":
    main()
