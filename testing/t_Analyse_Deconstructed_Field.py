#!/usr/bin/env python3
import argparse
import os
import numpy as np
import sys
sys.path.append("/".join(sys.path[0].split("/")[:-1]))
import Utils as ut
import FlowField as ffClass
import Tests
def main(File, s):
    #===================================================================#
    #### Format current directory path                               ####
    #===================================================================#
    pwd = ut.format_Directory_Path(os.getcwd())
    #===================================================================#
    #### Read deconstructed field                                    ####
    #===================================================================#
    deconstructed_field = ut.read_H5_Deconstructed(File)
    #===================================================================#
    #### Chck what we are working with?                              ####
    #===================================================================#
    if s:
        # Extract the norm singular values...
        sing_vals = np.zeros((len(deconstructed_field["Mx"]),
                              len(deconstructed_field["Mz"])),
                              dtype=float)
        Mx_shifted = []
        if deconstructed_field["Nx"] % 2 == 0:
            # even Nx
            Mx_shifted = list(deconstructed_field["Mx"][deconstructed_field["Nx"]/2+1:])
            Mx_shifted.extend(list(deconstructed_field["Mx"][:deconstructed_field["Nx"]/2+1]))
        else:
            # odd Nx
            Mx_shifted = list(deconstructed_field["Mx"][(deconstructed_field["Nx"]+1)/2:])
            Mx_shifted.extend(list(deconstructed_field["Mx"][:(deconstructed_field["Nx"]-1)/2+1]))
    
        Mz_shifted = []
        if deconstructed_field["Nz"] % 2 == 0:
            # even Nz
            Mz_shifted = list(deconstructed_field["Mz"][deconstructed_field["Nz"]/2+1:])
            Mz_shifted.extend(list(deconstructed_field["Mz"][:deconstructed_field["Nz"]/2+1]))
        else:
            # odd Nz
            Mz_shifted = list(deconstructed_field["Mz"][(deconstructed_field["Nz"]+1)/2:])
            Mz_shifted.extend(list(deconstructed_field["Mz"][:(deconstructed_field["Nz"]-1)/2+1]))
        
        Mx_shifted = int(Mx_shifted)
        Mz_shifted = int(Mz_shifted)
        
        for mx in range(0,len(Mx_shifted)):
            for mz in range(0,len(Mz_shifted)):
                # 2-norm to get magnitude of singular values
                sing_vals[mx,mz] = np.linalg.norm(deconstructed_field['singular_values'][mx,mz,:], 2)
        
        # sing_vals is a 2D matrix of data
        # plot a contour now
        fileName = str(File[:-3]) + "_singular_values_2Norm"
        max_val = max(sing_vals)
        min_val = min(sing_vals)
        colorbar_label = "$|| \sigma ||$"
        levels = 16
        ut.plot_Contour_coeffs(pwd, fileName,
                               Mx, Mz, sing_vals,
                               "Mx", "Mz",
                               colorbar_label, max_val, min_val,
                               levels)
    
directory = "/Users/arslan/temp/w03_EQ/eq1"
os.chdir(directory)
file = "eq1_deconstructed.h5"
main(file, True)