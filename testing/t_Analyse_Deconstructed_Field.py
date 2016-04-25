#!/usr/bin/env python3
import argparse
import os
import numpy as np
import sys
sys.path.append("/".join(sys.path[0].split("/")[:-1]))
import Utils as ut
import FlowField as ffClass
import Tests
def main(File, s, a, contourPlot):
    print(File)
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
    #### Singular values
    if s:
        print("Printing singular values")
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
        
        Mx_shifted = np.array(Mx_shifted).astype(int)
        Mz_shifted = np.array(Mz_shifted).astype(int)
        
        for mx in range(0,len(Mx_shifted)):
            for mz in range(0,len(Mz_shifted)):
                oldIndex_mx = Mx_shifted[mx]
                oldIndex_mz = Mz_shifted[mz]
                # 2-norm to get magnitude of singular values
                sing_vals[mx,mz] = np.linalg.norm(deconstructed_field['singular_values'][oldIndex_mx,oldIndex_mz,:], 2)
                
                if oldIndex_mx == 0 and oldIndex_mz == 0:
                    sing_vals[mx,mz] = 1e-17
        
        # sing_vals is a 2D matrix of data
        # plot a contour now
        fileName = str(File[:-3]) + "_singular_values_2Norm"
        max_val = sing_vals.max()
        min_val = sing_vals.min()
        colorbar_label = "$|| \sigma ||$"
        levels = 16
        title = "$|| \sigma ||$ at each Fourier mode"
        ut.plot_Contour_coeffs(pwd, fileName,
                               Mx_shifted, Mz_shifted, sing_vals,
                               "Mx", "Mz",
                               colorbar_label, max_val, min_val,
                               levels, title, contourPlot)
    #### Amplitude coefficients
    if a:
        print("Printing amplitude coefficients values")
        # Extract the norm singular values...
        amp_coefs = np.zeros((len(deconstructed_field["Mx"]),
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
        
        Mx_shifted = np.array(Mx_shifted).astype(int)
        Mz_shifted = np.array(Mz_shifted).astype(int)
        
        for mx in range(0,len(Mx_shifted)):
            for mz in range(0,len(Mz_shifted)):
                oldIndex_mx = Mx_shifted[mx]
                oldIndex_mz = Mz_shifted[mz]
                # 2-norm to get magnitude of amplitude coefficients
                amp_coefs[mx,mz] = np.linalg.norm(deconstructed_field['coefficients'][oldIndex_mx,oldIndex_mz,:], 2)
                
                if oldIndex_mx == 0 and oldIndex_mz == 0:
                    amp_coefs[mx,mz] = 1e-17
        
        # sing_vals is a 2D matrix of data
        # plot a contour now
        fileName = str(File[:-3]) + "_amplitude_coeffs_2Norm"
        max_val = amp_coefs.max()
        min_val = amp_coefs.min()
        colorbar_label = "$|| \Xi ||$"
        levels = 16
        title = "$|| \Xi ||$ at each Fourier mode"
        ut.plot_Contour_coeffs(pwd, fileName,
                               Mx_shifted, Mz_shifted, amp_coefs,
                               "Mx", "Mz",
                               colorbar_label, max_val, min_val,
                               levels, title, contourPlot)
    #### Combined
    if a and s:
        print("Printing multiple of above values")
        # Extract the norm singular values...
        chi = np.zeros((len(deconstructed_field["Mx"]),
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
        
        Mx_shifted = np.array(Mx_shifted).astype(int)
        Mz_shifted = np.array(Mz_shifted).astype(int)
        
        for mx in range(0,len(Mx_shifted)):
            for mz in range(0,len(Mz_shifted)):
                oldIndex_mx = Mx_shifted[mx]
                oldIndex_mz = Mz_shifted[mz]
                # 2-norm to get magnitude of singular values
                s_val = np.asmatrix(deconstructed_field['singular_values'][oldIndex_mx,oldIndex_mz,:])
                a_val = np.asmatrix(deconstructed_field['coefficients'][oldIndex_mx,oldIndex_mz,:]).T
                tmp = s_val * a_val
                chi[mx,mz] = np.linalg.norm(tmp, 2)
                
                if oldIndex_mx == 0 and oldIndex_mz == 0:
                    chi[mx,mz] = 1e-10
        
        # sing_vals is a 2D matrix of data
        # plot a contour now
        fileName = str(File[:-3]) + "_chi_2Norm"
        max_val = chi.max()
        min_val = chi.min()
        colorbar_label = "$|| \chi ||$"
        levels = 16
        title = "$|| \chi ||$ at each Fourier mode"
        ut.plot_Contour_coeffs(pwd, fileName,
                               Mx_shifted, Mz_shifted, chi,
                               "Mx", "Mz",
                               colorbar_label, max_val, min_val,
                               levels, title, contourPlot)
    print("Finished")

name="eq11"
directory = "/home/arslan/Documents/work/cfd-channelflow_solutions/w03_EQ/" + name
os.chdir(directory)
file = name+"_deconstructed.h5"
main(file, True, True, False)