#!/usr/bin/env python3
import argparse
import os
import sys
sys.path.append("/".join(sys.path[0].split("/")[:-1]))
import Utils as ut
parser = argparse.ArgumentParser(description="Run nkh searches and approx searches and post-process at the end.")
parser.add_argument("-dir",
                    "--Directory",
                    metavar='\b',
                    help="Directory to perform searches in.",
                    required=True)
args = parser.parse_args()
#===================================================================#
#### ---Change to the case directory                             ####
#===================================================================#
directory = ut.format_Directory_Path(args.Directory)
os.chdir(directory)
print(directory)
#===================================================================#
#### Get Reynolds number                                         ####
#===================================================================#
for root, sub_dirs, files in os.walk(args.Directory):
    for file in files:
        if file.endswith("tails.txt"):
            details_file = file
details = ut.read_Details(details_file)
reynolds = details['Re']
#===================================================================#
#### Get the name of the case                                    ####
#===================================================================#
name = directory.split("/")[-2]
#===================================================================#
#### ---Plot flow field                                          ####
#===================================================================#
plot = "e_Plot_Field.py -f "+name+".h5 -d "+details_file+" -i 0 -n 0 -a "
print(plot)
os.system(plot)
#===================================================================#
#### ---NKH Original                                             ####
#===================================================================#
nkh = "findsoln -eqb -R "+str(reynolds)+" -Nn 40 -sn -log nkh-"+name+".log "+name+".h5"
print(nkh)
os.system(nkh)
#===================================================================#
#### Deconstruct                                                 ####
#===================================================================#
deconstruct = "e_Deconstruct_Field.py -d "+details_file+" -f "+name+".h5"
print(deconstruct)
os.system(deconstruct)





#===================================================================#
#### Rank: Full                                                  ####
#===================================================================#
#### ---Construct                                                ####
#===================================================================#
print("\n======================================================================")
print("RANK: Full")
print("======================================================================")
rankfull = "\n\ne_Construct_Field.py -f "+name+"_deconstructed.h5 -r 1000"
print(rankfull)
os.system(rankfull)
fileName = name+"_rank_full"
rankfull_d = directory + fileName
os.chdir(rankfull_d)
#===================================================================#
#### ---Plot flow field                                          ####
#===================================================================#
plot = "e_Plot_Field.py -f "+fileName+".h5 -d ../"+details_file+" -i 0 -n 0 -a"
print(plot)
os.system(plot)
#===================================================================#
#### ---NKH Full rank                                            ####
#===================================================================#
nkh = "findsoln -eqb -R "+str(reynolds)+" -Nn 40 -sn -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)




#===================================================================#
#### Loop through ranks                                          ####
#===================================================================#
ranks = [2,4,10,20,40,60]
for r in range(0, len(ranks)):
    rank = ranks[r]
    print("\n======================================================================")
    print("RANK: " + str(rank))
    print("======================================================================")
    #### Change into main folder
    os.chdir(directory)
    #### Construct field into rank folder
    construct = "\ne_Construct_Field.py -f "+name+"_deconstructed.h5 -r " + str(rank)
    print(construct)
    os.system(construct)
    #### Change into rank folder
    fileName = name+"_rank_" + str(rank).zfill(2)
    rank_dir = directory + fileName
    os.chdir(rank_dir)
    #### Plot flow field
    plot = "e_Plot_Field.py -f "+fileName+".h5 -d ../"+details_file+" -i 0 -n 0 -a"
    print(plot)
    os.system(plot)
    #### NKH Search
    nkh = "\nfindsoln -eqb -R "+str(reynolds)+" -Nn 40 -sn -log nkh-"+fileName+".log "+fileName+".h5"
    print(nkh)
    os.system(nkh)



#===================================================================#
#### Post process                                                ####
#===================================================================#
#### Change into main folder
os.chdir(directory)
command = "e_Plot_Convergence_NKH.py -dir " + directory
print(command)
os.system(command)
