#!/usr/bin/env python3
import numpy as np
import os
import Utils as ut
output_directory = os.getcwd()
output_directory = ut.format_Directory_Path(output_directory)
# Make a directory to store the files with integer names.
ints_directory = output_directory + 'ints'
# if directory exists, delete it:
if os.path.exists(ints_directory):
    print("\nRemoving:\t" + str(ints_directory))
    command = "rm -rf " + ints_directory
    os.system(command)
# if directory doesn't exist, make it:
if not os.path.exists(ints_directory):
    print("Making:\t\t" + str(ints_directory))
    os.mkdir(ints_directory)
ints_directory = output_directory + 'ints'
T0 = 0
T1 = 1400
steps = T1 - T0 + 1
TRange = np.linspace(T0, T1, steps)
#### Change directory down
# change into the ints folder and list files and rename the float files to integer files.
os.chdir(ints_directory)

files = [fi for fi in os.listdir(output_directory) if os.path.isfile(os.path.join(output_directory,fi))]
files = sorted(files)
for k in files:
    if str(k)[0] == 'u' and str(k).find("asc") == -1:
        oldFile = k
        k = k[1:-3]
        k = float(k)

        if k in TRange:
            k = int(k)
            newFile = 'u' + str(k) + '.h5'
            oldFile = '../' + oldFile
            command = 'fieldconvert ' + str(oldFile) + ' ' + str(newFile)
            print(command)
            os.system(command)

print('\nRenamed and stored flow fields with integers in names.')