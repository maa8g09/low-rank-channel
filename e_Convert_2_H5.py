#!/usr/bin/env python


import argparse
import os

parser = argparse.ArgumentParser(description="Convert all files in directory given from .ff to .h5")
parser.add_argument("-d",
                    "--Directory",
                    metavar='\b',
                    help="Directory where all .ff files are kept.",
                    required=True)

args = parser.parse_args()

directory = args.Directory
if directory[-1] != "/":
    directory += "/"
print("")
print(directory)
print("")
# List all files in the directory,
# convert each files with channelflow command convertfield
files = [fi for fi in os.listdir(directory) if os.path.isfile(os.path.join(directory,fi))]
files = sorted(files)
print('\nConverting...')
for file in files:
    file = str(file)
    if file[-3:] == ".ff":
        file_ff = file
        file_h5 = file[:-3] + ".h5"
        
        command = "\nfieldconvert " + file_ff + " " + file_h5
        print(command)
        os.system(command)

print('\nConverted all flow fields to HDF5 format.')
