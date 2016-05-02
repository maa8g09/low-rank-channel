#!/usr/bin/env python3
"""
=====================================================================
Average fields
=====================================================================
Average given fields. Pretty obvious what this does.


Muhammad Arslan Ahmed
maa8g09@soton.ac.uk
Room 5069
Building 13
Aerodynamics and Flight Mechanics
University of Southampton

"""
import argparse
import os
import numpy as np
import sys
sys.path.append("/".join(sys.path[0].split("/")[:-1]))
import FlowField as ffClass
import Utils as ut

descriptionString = "Average given flows."
parser = argparse.ArgumentParser(description=descriptionString)
parser.add_argument("-f", "--Files", help="Files to approximate",
                    metavar='\b', nargs='+', required=True)
args = parser.parse_args()
pwd = ut.format_Directory_Path(os.getcwd())
tmp_f = args.Files[0]
tmp_f, tmp_o = ut.read_H5(tmp_f)
cumulative_sum = np.zeros((3, tmp_f['Nx'], tmp_f['Ny'], tmp_f['Nz']))
# Add files together to get cumulative sum
print("Averaging...")
for file in args.Files:
    file_info, original_attrs = ut.read_H5(file)
    cumulative_sum += file_info['ff']
    print("\t"+file)
average = np.zeros((3, tmp_f['Nx'], tmp_f['Ny'], tmp_f['Nz']))
average = cumulative_sum * (1.0 / len(args.Files))
# Initialise flow field class for average
ff_average = ffClass.FlowFieldChannelFlow( 3,
                                           tmp_f['Nx'],
                                           tmp_f['Ny'],
                                           tmp_f['Nz'],
                                           tmp_f['Lx'],
                                           tmp_f['Lz'],
                                           tmp_f['alpha'],
                                           tmp_f['beta'],
                                           0.0,
                                           "",
                                           0,
                                           average,
                                           "pp")
fileName = "averaged-"+args.Files[0][:-3]+"-"+args.Files[-1][:-3]
ut.write_H5(ff_average, tmp_o, fileName)
ut.print_EndMessage()