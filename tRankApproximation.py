#!/usr/bin/env python

import os

import Utils as ut
import FlowField as ffClass
import ChannelResolvent as cr

import numpy as np


def main(File, Rank, Directory, MeanProfile, sparse):
    #================================================================
    #### Create a temporary folder
    #================================================================
    parent_directory = os.getcwd()
    
    if parent_directory[-1] != "/":
        parent_directory += "/"
    
    temp_rank_folder = "rank-temp/"
    temp_rank_folder = parent_directory + temp_rank_folder
    
    #if a temporary directory exists, delete it:
    if os.path.exists(temp_rank_folder):
        command = "rm -rf " + temp_rank_folder
        os.system(command)
    
    #if a temporary directory doesn't exist, create one:
    if not os.path.exists(temp_rank_folder):
        os.mkdir(temp_rank_folder)
    
    # All work is done from the temporary directory
    os.chdir(temp_rank_folder)
    
    
    #================================================================
    #### Check file type
    #================================================================
    if File[-3:] == ".h5": # H5 file type
        #------------------------------------------------
        # Read H5 file
        #------------------------------------------------
        print("HDF5 file given.")
        #------------------------------------------------
        # Initialize instance of flow field class (Ati's class)
        #------------------------------------------------
    
        #------------------------------------------------
        # Fourier transform the velocity field
        #------------------------------------------------
    
        # (Ati has written a make_spectral method which you can use.)
    
    
    
    elif File[-3:] == ".ff": # channelflow binary file type
        print("\nA channelflow binary file given...")
        #------------------------------------------------
        #### Convert the binary file to ascii
        #------------------------------------------------
        command = "field2ascii -p ../" + str(File) + " " + str(File)[:-3]
        print(command)
        os.system(command)
    
        #------------------------------------------------
        #### Read physical ascii file
        #------------------------------------------------
        var = ut.read_ASC_PP(temp_rank_folder, str(File)[:-3])
        
        details_directory = Directory
        if details_directory[-1] != "/":
            details_directory += "/"
        
        var2 = ut.read_Details(details_directory, "u0")
        
        #------------------------------------------------
        #### Initialize an instance of FlowField class (Fluctuations)
        #------------------------------------------------
        ffcf_flucs = ffClass.FlowFieldChannelFlow2( var['Nd'],
                                                    var['Nx'],var['Ny'],var['Nz'],
                                                    var['Lx'],var['Lz'],
                                                    var['alpha'],var['beta'],
                                                    var2['c'],var2['bf'],var2['Re'],
                                                    var['ff'],
                                                    "pp")

    else: # No file type given.
        ut.error("Invalid file given.")
    

    #================================================================
    #### Check mean file
    #================================================================
    turb_mean_profile = []
    if MeanProfile:
        #------------------------------------------------
        #### Read velocity profile
        #------------------------------------------------
        turb_deviation_profile = ut.read_Vel_Profile(parent_directory, MeanProfile)
        if str(MeanProfile).find("mean") == -1:
            print("Turbulent deviation profile given.\nAdding Laminar profile to it.")
            # Laminar profile
            lam = 1.0 - ffcf_flucs.y**2.0
            # Add turbulent deviation profile to the parabolic laminar base flow profile
            turb_mean_profile = turb_deviation_profile + lam
    
        elif str(MeanProfile).find("mean") != -1:
            print("Turbulent mean profile given.")
            turb_mean_profile = turb_deviation_profile
    
        #------------------------------------------------
        #### Construct 4D array of turb profile
        #------------------------------------------------
        turb_mean = ut.make_mean_ff_pp(turb_mean_profile, ffcf_flucs.Nd, ffcf_flucs.Nx, ffcf_flucs.Nz)
    
        #------------------------------------------------
        #### Initialize mean instance FlowField class
        #------------------------------------------------
        ffmean = ffClass.FlowFieldChannelFlow2( var['Nd'],
                                                var['Nx'],var['Ny'],var['Nz'],
                                                var['Lx'],var['Lz'],
                                                var['alpha'],var['beta'],
                                                var2['c'],var2['bf'],var2['Re'],
                                                turb_mean,
                                                "pp")
        ffmean.set_ff(turb_mean, "pp")
        instantaneous_field = turb_mean + var['ff']
        ffcf = ffClass.FlowFieldChannelFlow2(var['Nd'],
                                            var['Nx'],var['Ny'],var['Nz'],
                                            var['Lx'],var['Lz'],
                                            var['alpha'],var['beta'],
                                            var2['c'],var2['bf'],var2['Re'],
                                            instantaneous_field,
                                            "pp")
    else:
        #------------------------------------------------
        #### Use original file as the mean
        #------------------------------------------------
        ffmean = ffClass.FlowFieldChannelFlow2( var['Nd'],
                                                var['Nx'],var['Ny'],var['Nz'],
                                                var['Lx'],var['Lz'],
                                                var['alpha'],var['beta'],
                                                var2['c'],var2['bf'],var2['Re'],
                                                var['ff'],
                                                "pp")
    
    
    #================================================================
    # Approximate the file w/regards to specified rank
    #================================================================
    print("\nStarting approximation...\n")
    approximated_inst_ff, alpha_beta_chi = cr.resolvent_approximation2(ffcf, Rank, turb_mean_profile, ffmean, sparse)
    # The flow field that is saved in the instance is (s,p)
    print("State of approximated field.")
    print(ffcf.state)
    
    
    if_u = instantaneous_field[0, :, :, :].real
    if_v = instantaneous_field[1, :, :, :].real
    if_w = instantaneous_field[2, :, :, :].real

    af_u = approximated_inst_ff[0, :, :, :].real
    af_v = approximated_inst_ff[1, :, :, :].real
    af_w = approximated_inst_ff[2, :, :, :].real

    # Original - New = fluctuations.
    diff = instantaneous_field.real - approximated_inst_ff.real

    diff_u = diff[0, :, :, :]
    diff_v = diff[1, :, :, :]
    diff_w = diff[2, :, :, :]

    ndiff_u = np.linalg.norm(diff_u)
    ndiff_v = np.linalg.norm(diff_v)
    ndiff_w = np.linalg.norm(diff_w)

    nodiff_u = np.linalg.norm(var['ff'][0, :, :, :].real)
    nodiff_v = np.linalg.norm(var['ff'][1, :, :, :].real)
    nodiff_w = np.linalg.norm(var['ff'][2, :, :, :].real)

    flucs = approximated_inst_ff.real - turb_mean.real

#    flucs_u = af_u - turb_mean[0, :, :, :].real
#    flucs_v = af_v - turb_mean[1, :, :, :].real
#    flucs_w = af_w - turb_mean[2, :, :, :].real
#
    orig_fu = if_u - turb_mean[0, :, :, :].real
    orig_fv = if_v - turb_mean[1, :, :, :].real
    orig_fw = if_w - turb_mean[2, :, :, :].real

    
    diff_na = np.linalg.norm(diff)
    diff_no = np.linalg.norm(var['ff'].real)
    flucs_diffs = var['ff'].real - diff
    flucs_diffsu = flucs_diffs[0,:,:,:]
    flucs_diffsv = flucs_diffs[1,:,:,:]
    flucs_diffsw = flucs_diffs[2,:,:,:]
    
    flucs_diffs_n= np.linalg.norm(flucs_diffs)
    
    

    approx_field = ffClass.FlowFieldChannelFlow2( var['Nd'],
                                            var['Nx'],var['Ny'],var['Nz'],
                                            var['Lx'],var['Lz'],
                                            var['alpha'],var['beta'],
                                            var2['c'],var2['bf'],var2['Re'],
                                            diff,
                                            "pp")


    #================================================================
    # Create a folder to store the approximated velocity field in
    #================================================================
    os.chdir(parent_directory) # Go up one directory
    rank_folder = File[:-3]+"_rank_" + str(ffcf.rank) + "/"
    rank_folder = parent_directory + rank_folder
    
    #if a rank directory does exist, delete it:
    if os.path.exists(rank_folder):
        command = "rm -rf " + rank_folder
        os.system(command)
    
    #if a rank directory doesn't exist, create one:
    if not os.path.exists(rank_folder):
        os.mkdir(rank_folder)
    
    # Change into the new rank directory
    os.chdir(rank_folder)
    
    
    #================================================================
    # Save flow field to file
    #================================================================
    # Check file type
    if File[-3:] == ".h5":
        #------------------------------------------------
        # Inverse Fourier transform the velocity field
        #------------------------------------------------
        print("Inverse Fourier transform the velocity field.")
        #------------------------------------------------
        # Write the file to disk in H5 format
        #------------------------------------------------
    
    
    elif File[-3:] == ".ff":
        
    #    #------------------------------------------------
    #    # Write spectral ascii file
    #    #------------------------------------------------
    #    sp_ASC_fileName = ut.write_approximated_ASC(ffcf, rank_folder, ffcf.rank)
    #    
    #    #------------------------------------------------
    #    # Convert ascii to binary flowfield
    #    #------------------------------------------------
    #    print(os.getcwd())
    #    approximationFileName = args.File[:-3]+"_rank_"+str(ffcf.rank)
    #    command = "ascii2field -p false -ge ../rank-temp/" + str(args.File)[:-3] + ".geom " + sp_ASC_fileName + " " + approximationFileName
    #    print(command)
    #    os.system(command)
    
        #------------------------------------------------
        # Write physical ascii file
        #------------------------------------------------
        fileName = File[:-3] + "_rnk_" + str(ffcf.rank)
        ut.write_ASC(approx_field, rank_folder, fileName)
        command = "ascii2field -p false -ge ../rank-temp/" + str(File)[:-3] + ".geom " + fileName + ".asc " + fileName
        print(command)
        os.system(command)
        
        #------------------------------------------------
        # Write amplitude coefficients for each Fourier mode combination
        #------------------------------------------------
        fileName = File[:-3] + "_coeffs"
        ut.write_amplitude_coefficients(approx_field, rank_folder, fileName, alpha_beta_chi)
    
        #------------------------------------------------
        # Remove ascii file and temporary folder
        #------------------------------------------------
    #    os.system("rm *.asc")
        os.chdir(parent_directory)
        command = "rm -rf " + temp_rank_folder
        os.system(command)
    
    print("Done!")




dirc = "/home/arslan/Desktop/ati-modes-copy/modes/B/kz2/ff_files"
#dirc = "/Users/arslan/Desktop/ati-modes-copy/modes/B/kz2/ff_files"
os.chdir(dirc)
main("mode-00.ff", 
     2, 
     dirc, 
     "turbdeviation.txt",
     True)
