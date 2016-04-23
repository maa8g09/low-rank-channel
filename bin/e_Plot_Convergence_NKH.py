#!/usr/bin/env python3
import argparse
import os
from .. import Utils as ut
parser = argparse.ArgumentParser(description="Convergence plot of NKH data.")
parser.add_argument("-f",
                    "--File",
                    metavar='\b',
                    help="File to read.",
                    required=True)
parser.add_argument("-dir",
                    "--Directory",
                    metavar='\b',
                    help="Directory to perform searches in.")
args = parser.parse_args()

# If a convergence file given:
# plot it

# If a directory given:
# plot each convergence file in each subdirectory
# 
# check to see if a convergence file exists in the directory given, if there is... use it as a constant line.

if args.Directory: # if a directory is given...
    directory = ut.format_Directory_Path(args.Directory)
    os.chdir(directory)
    # the last directory name is the case name.
    fileName = args.Directory.split("/")
    fileName = fileName[-1]
    # Read convergence files in this directory and all sub directories...
    all_convergence_data = {}
    for root, sub_dirs, files in os.walk(directory):
        for file in files:
            if file.endswith("convergence.asc"):
                 print(os.path.join(root, file))
                 # read in convergence file
                 # initialise convergence dictionary
#                 all_convergence_data[fileName] = {}
                 # populate the dictionary
#                 all_convergence_data[fileName] = ut.read_NKH_convergence(file)
        

        for sub_dir in sub_dirs:
            os.chdir(sub_dir)
            for file in os.listdir(sub_dir):
                if file.endswith("convergence.asc"):
                    print(os.path.join(sub_dir, file))
                    # read in convergence file
                    # initialise convergence dictionary
#                    all_convergence_data[fileName] = {}
                    # populate the dictionary
#                    all_convergence_data[fileName] = ut.read_NKH_convergence(file)
           
            
            
#    os.chdir(directory)