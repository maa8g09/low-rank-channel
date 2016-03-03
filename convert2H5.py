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

# List all files in the directory,
# convert each files with channelflow command convertfield
files = [fi for fi in os.listdir(directory) if os.path.isfile(os.path.join(directory,fi))]
files = sorted(files)
print('\nConverting...')
for k in files:
    k = str(k)
    if k[0] == "u" and k[-3:] == ".ff":
        file_ff = k
        file_h5 = k[:-3]
        
        command = "\nfieldconvert " + file_ff + " " + file_h5
        print(command)
        os.system(command)

print('\nConverted all flow fields to HDF5 format.')
