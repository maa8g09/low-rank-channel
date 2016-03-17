#!/usr/bin/env python

import argparse
import os

parser = argparse.ArgumentParser(description="Convert all files in directory given from .ff to .h5")
parser.add_argument("-d",
                    "--Directory",
                    metavar='\b',
                    help="Directory where all rank folders are kept.",
                    required=True)

args = parser.parse_args()

pwd = args.Directory
if pwd[-1] != "/":
    pwd += "/"

#### List all directories in directory passed in...
dirs = []
for directory in os.listdir(pwd):
    tmp = os.path.join(pwd, directory)
    if os.path.isdir(tmp):
        print("")
        print(tmp)
        dirs.append(tmp)



#### Loop through all directories
for directory in dirs:
    #### Change into directory
    os.chdir(directory)
    print("\nChanged to"+directory)
    #### List files in directory
    files = [fi for fi in os.listdir(directory) if os.path.isfile(os.path.join(directory,fi))]
    print('\nConverting...')
    for file in files:
        file = str(file)
        #### If the file is a channelflow binary...
        if file[-3:] == ".ff":
            file_ff = file
            file_h5 = file[:-3] + ".h5"
            ### Convert it
            command = "\nfieldconvert " + file_ff + " ../" + file_h5
            print(command)
            os.system(command)

