#!/usr/bin/env python3
import argparse
import os
import sys
sys.path.append("/".join(sys.path[0].split("/")[:-1]))
import Utils as ut
parser = argparse.ArgumentParser(description="Run nkh searches")

parser.add_argument("-f",
                    "--File",
                    metavar='\b',
                    help="File to approximate",
                    required=True)
parser.add_argument("-dir",
                    "--Directory",
                    metavar='\b',
                    help="Directory to perform searches in.",
                    required=True)
parser.add_argument("-d",
                    "--Details",
                    metavar='\b',
                    help="Details file.",
                    required=True)
args = parser.parse_args()
name=str(args.File)
directory = ut.format_Directory_Path(args.Directory)
directory+=name
directory = ut.format_Directory_Path(directory)
os.chdir(directory)
print(directory)
#### Plot flow field
plot = "e_Plot_Field.py -f "+name+".h5 -d "+args.Details+" -i 0 -n 0 -a"
print(plot)
os.system(plot)
##### ---NKH Search
T = 20
if name == "tw2":
    T = 13.9
elif name == "tw3":
    T = 11.9
nkh = "findsoln -eqb -xrel -T "+str(T)+" -sn -log nkh-"+name+".log "+name+".h5"
print(nkh)
os.system(nkh)
#### Deconstruct
deconstruct = "e_Deconstruct_Field.py -d "+args.Details+" -f "+name+".h5"
print(deconstruct)
os.system(deconstruct)

#### Rank: Full
#### ---Construct
print("\n======================================================================")
print("RANK: Full")
print("======================================================================")
rankfull = "\n\ne_Construct_Field.py -f "+name+"_deconstructed.h5 -r 1000"
print(rankfull)
os.system(rankfull)
fileName = name+"_rank_full"
rankfull_d = directory + fileName
os.chdir(rankfull_d)
#### ---Plot flow field
plot = "e_Plot_Field.py -f "+fileName+".h5 -d ../"+args.Details+" -i 0 -n 0 -a"
print(plot)
os.system(plot)
#### ---NKH Search
nkh = "findsoln -eqb -xrel -T "+str(T)+" -sn -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)

#### Loop through ranks
ranks = [2,4,10,20,40,60]
for r in range(0, len(ranks)):
    rank = ranks[r]
    print("\n======================================================================")
    print("RANK: " + str(rank))
    print("======================================================================")
    #### Change into main folder
    os.chdir(directory)
    #### Construct field into rank folder
    construct = "\n\ne_Construct_Field.py -f "+name+"_deconstructed.h5 -r " + str(rank)
    print(construct)
    os.system(construct)
    #### Change into rank folder
    fileName = name+"_rank_" + str(rank).zfill(2) 
#    fileName2= name+"_rank_" + str(rank)
    rank_dir = directory + fileName
    os.chdir(rank_dir)
    #### Plot flow field
    plot = "e_Plot_Field.py -f "+fileName+".h5 -d ../"+args.Details+" -i 0 -n 0 -a"
    print(plot)
    os.system(plot)
    #### NKH Search
    nkh = "\nfindsoln -eqb -xrel -T "+str(T)+" -sn -log nkh-"+fileName+".log "+fileName+".h5"
    print(nkh)
    os.system(nkh)

