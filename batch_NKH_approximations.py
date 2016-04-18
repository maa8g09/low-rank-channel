#!/usr/bin/env python3
import argparse
import os
parser = argparse.ArgumentParser(description="Run nkh searches for W03 EQ")

parser.add_argument("-f",
                    "--File",
                    metavar='\b',
                    help="File to approximate",
                    required=True)


args = parser.parse_args()
name=str(args.File)
directory="/home/arslan/Documents/work/cfd-channelflow_solutions/w03_EQ/"+name+"/"
os.chdir(directory)
# Deconstruct
deconstruct = "eDeconstructField.py -d eq1_Details.txt -f "+name+".h5"
os.system(deconstruct)
# Construct
# Rank: Full
rankfull = "eConstructField.py -f "+name+"_deconstructed.h5 -r 1000"
os.system(rankfull)
fileName = name+"_rank_full"
rankfull_d = directory + fileName
os.chdir(rankfull_d)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 2
r=2
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 4
r=4
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 10
r=10
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 20
r=20
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)

