#!/usr/bin/env python
import argparse
import os
parser = argparse.ArgumentParser(description="Rank approximate all files in given folder")
parser.add_argument("-r",
                    "--Rank",
                    metavar='\b',
                    help="Number of velocity modes to use to recreate flow field.",
                    required=True,
                    type=int)
args = parser.parse_args()
# Run from ff_files
pwd = os.getcwd()
if pwd[-1] != "/":
    pwd += "/"
# List all files ith format .ff
files = [fi for fi in os.listdir(pwd) if os.path.isfile(os.path.join(pwd,fi))]
files = sorted(files)
# Loop through and run the rank command.
for file in files:
    file = str(file)
    if file[-3:] == ".ff":
        command = "eRankApproximation.py"
        command+=" -d " + pwd
        command+=" -v turbdeviation.txt"
        command+=" -f " + file
        command+=" -r " + str(args.Rank)
        command+=" -s"
        print("")
        print(command)
        print("")
        os.system(command)
print("Done")
