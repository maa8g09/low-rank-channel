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
print(directory)
#### Plot flow field
plot = "ePlotField.py -f "+name+".h5 -d eq1_Details.txt -i 0 -n 0"
print(plot)
os.system(plot)
plot = "ePlotField.py -f "+name+".h5 -d eq1_Details.txt -i 0 -n 2"
print(plot)
os.system(plot)
#### Deconstruct
deconstruct = "eDeconstructField.py -d eq1_Details.txt -f "+name+".h5"
print(deconstruct)
os.system(deconstruct)

#### Rank: Full
#### ---Construct
rankfull = "eConstructField.py -f "+name+"_deconstructed.h5 -r 1000"
os.system(rankfull)
fileName = name+"_rank_full"
rankfull_d = directory + fileName
os.chdir(rankfull_d)
#### ---Plot flow field
plot = "ePlotField.py -f "+fileName+".h5 -d ../eq1_Details.txt -i 0 -n 0"
print(plot)
os.system(plot)
plot = "ePlotField.py -f "+fileName+".h5 -d ../eq1_Details.txt -i 0 -n 2"
print(plot)
os.system(plot)
#### ---NKH Search
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
print("\n\n::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n\n")

#### Loop through ranks
ranks = [2,4,10,20]
for r in range(0, len(ranks)):
    rank = ranks[r]
    #### Change into main folder
    os.chdir(directory)
    #### Construct field into rank folder
    rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(rank)
    os.system(rank)
    #### Change into rank folder
    fileName = name+"_rank_" + str(rank)
    rank_dir = directory + fileName
    os.chdir(rank_dir)
    #### Plot flow field
    plot = "\nePlotField.py -f "+fileName+".h5 -d ../eq1_Details.txt -i 0 -n 0"
    print(plot)
    os.system(plot)
    plot = "\nePlotField.py -f "+fileName+".h5 -d ../eq1_Details.txt -i 0 -n 2"
    print(plot)
    os.system(plot)
    #### NKH Search
    nkh = "\nfindsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
    print(nkh)
    os.system(nkh)
    print("\n\n::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n\n")
