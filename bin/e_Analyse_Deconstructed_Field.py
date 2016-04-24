#!/usr/bin/env python3
import argparse
import os
import numpy as np
import sys
sys.path.append("/".join(sys.path[0].split("/")[:-1]))
import Utils as ut
import FlowField as ffClass
import Tests
parser = argparse.ArgumentParser(description="Analyse deconstructed flow field.")
parser.add_argument("-f", "--File", help="Deconstructed file to analyse",
                    metavar='\b', required=True)
parser.add_argument("-s", "--SingVals", help="Plot/Extract singular values",
                    action='store_true')
parser.add_argument("-a", "--AmplCoeffs", help="Plot/Extract amplitude coefficients",
                    action='store_true')
args = parser.parse_args()
#===================================================================#
#### Format current directory path                               ####
#===================================================================#
pwd = ut.format_Directory_Path(os.getcwd())
#===================================================================#
#### Read deconstructed field                                    ####
#===================================================================#
deconstructed_field = ut.read_H5_Deconstructed(args.File)
#===================================================================#
#### Chck what we are working with?                              ####
#===================================================================#
if args.SingVals:
    # Extract the norm singular values...
    deconstructed_field['singular_values']
    
    # Plot magnitude of singular values at each Fourier pair. 
    for 