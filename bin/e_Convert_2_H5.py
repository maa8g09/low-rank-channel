#!/usr/bin/env python3
"""
=====================================================================
Convert to H5
=====================================================================
The decomposed flow field is reconstructed using its resolvent modes,
singular values and amplitude coefficients.

At each Fourier mode pair:
u_hat[mx,:,mz] = resolvent_modes[mx,mz,:,rank] * singular_values[mx,mz,rank] *
        amplitude_coefficients[mx,mz,rank]

The decomposed field's components are saved in an HDF5 file
with the following structure:
    deconstructed_field/:
        resolvent_modes
        forcing_modes
        singular_values
        coefficients


Muhammad Arslan Ahmed
maa8g09@soton.ac.uk
Room 5069
Building 13
Aerodynamics and Flight Mechanics
University of Southampton

"""
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
