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
    # Add slash at the end of the string if there isn't one already
    if parent_directory[-1] != "/":
        parent_directory += "/"

    temp_rank_folder = "rank-temp/"
    temp_rank_folder = parent_directory + temp_rank_folder

    #if a temporary directory exists, delete it.
    if os.path.exists(temp_rank_folder):
        command = "rm -rf " + temp_rank_folder
        os.system(command)

    #if a temporary directory doesn't exist, create one.
    if not os.path.exists(temp_rank_folder):
        os.mkdir(temp_rank_folder)

    # All work is done from the temporary directory.
    os.chdir(temp_rank_folder)

    #================================================================
    #### Check file type
    #================================================================
    if File[-3:] == ".h5": # H5 file type
        print("HDF5 file given.")
        #------------------------------------------------
        #### Convert it to binary flow field
        #------------------------------------------------
        command = "fieldconvert ../" + str(File) + " ../" + str(File)[:-3] + ".ff"
        print(command)
        os.system(command)

        #------------------------------------------------
        #### Convert binary to ASCII
        #------------------------------------------------
        command = "field2ascii -p -g ../" + str(File)[:-3] + ".ff " + str(File)[:-3]
        print(command)
        os.system(command)

        #------------------------------------------------
        #### Read physical ascii file
        #------------------------------------------------
        file_info = ut.read_ASC_channelflow(temp_rank_folder, str(File)[:-3])

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
        file_info = ut.read_ASC_PP(temp_rank_folder, str(File)[:-3])
        details = ut.read_Details(parent_directory, "u0_Details.txt")


    else: # No file type given.
        ut.error("Invalid file given.")
    

    #================================================================
    #### Initialise original flow field object
    #================================================================
    ff_original = ffClass.FlowFieldChannelFlow( file_info['Nd'],
                                                file_info['Nx'],
                                                file_info['Ny'],
                                                file_info['Nz'],
                                                file_info['Lx'],
                                                file_info['Lz'],
                                                file_info['alpha'],
                                                file_info['beta'],
                                                details['c'],
                                                details['bf'],
                                                details['Re'],
                                                file_info['ff'],
                                                "pp")

    test_u = file_info['ff'][0,:,:,:].real

    # Remove wall boundaries
    ff_original.remove_wall_boundaries()
    test_u_walls = ff_original.velocityField[0,:,:,:].real

    # FFT
    ff_original.make_xz_spectral()
    test_u_fft = ff_original.velocityField[0,:,:,:].real

    # Stack
    ff_original.stack_ff_in_y()
    test_u_stk = ff_original.velocityField[:,:,:].real

    # Unstack
    ff_original.unstack_ff()
    test_u_ustk = ff_original.velocityField[0,:,:,:].real

    # IFFT
    ff_original.make_xz_physical()
    test_u_ifft = ff_original.velocityField[0,:,:,:].real

    # Add wall boundaries
    ff_original.add_wall_boundaries()
    test_u_walls2 = ff_original.velocityField[0,:,:,:].real

    test_u3 = ff_original.velocityField[0,:,:,:].real
    d = test_u3.real - file_info['ff'][0,:,:,:].real
    dnorm = np.linalg.norm(d)


    #================================================================
    #### Check velocity profile
    #================================================================
    # Create empty 4D array to store mean flow field
    mean = np.zeros((file_info['Nd'], file_info['Nx'], file_info['Ny'], file_info['Nz']))
    mean_profile = []

    if MeanProfile: # Velocity profile given
        #------------------------------------------------
        #### Read velocity profile
        #------------------------------------------------
        vel_profile = ut.read_Vel_Profile(parent_directory, MeanProfile)
        # Check to see if it is a mean profile or a deviation profile.
        deviation = any(n < 0 for n in vel_profile)
        if deviation: # Deviation profile given
            # Add baseflow to deviation
            baseflow = []
            if details['bf'] == "lam": # Laminary base flow
                baseflow = 1.0 - ff_original.y**2.0
            elif details['bf'] == "cou": # Couette base flow
                baseflow = ff_original.y

            # Add baseflow to deviation
            mean_profile = vel_profile + np.asarray(baseflow)

        else: # Mean profile given
            mean_profile = vel_profile

        #------------------------------------------------
        #### Construct 4D array from mean_profile
        #------------------------------------------------
        mean = ut.make_ff_from_profile(mean_profile, 
                                       ff_original.Nd, 
                                       ff_original.Nx, 
                                       ff_original.Nz)
    else:
        #------------------------------------------------
        #### Use base flow only
        #------------------------------------------------
        # Add baseflow to deviation
        baseflow = []
        if details['bf'] == "lam": # Laminary base flow
            baseflow = 1.0 - ff_original.y**2.0
        elif details['bf'] == "cou": # Couette base flow
            baseflow = ff_original.y

        #------------------------------------------------
        #### Construct 4D array from mean_profile
        #------------------------------------------------
        mean = ut.make_ff_from_profile(np.asarray(baseflow), ff_original.Nd, ff_original.Nx, ff_original.Nz)


    #================================================================
    #### Initialize mean flow field object
    #================================================================
    ff_mean = ffClass.FlowFieldChannelFlow( file_info['Nd'],
                                            file_info['Nx'],
                                            file_info['Ny'],
                                            file_info['Nz'],
                                            file_info['Lx'],
                                            file_info['Lz'],
                                            file_info['alpha'],
                                            file_info['beta'],
                                            details['c'],
                                            details['bf'],
                                            details['Re'],
                                            mean,
                                            "pp")


    #================================================================
    #### Remove the wall boundaries
    #================================================================
    # Removing the xz-planes at y=1 and y=-1,
    # so that the chebyshev nodes can be used to construct 
    # the transfer function.
    ff_original.remove_wall_boundaries()
    ff_mean.remove_wall_boundaries()


    #================================================================
    #### Fourier transform original and mean velocity fields in xz directions
    #================================================================
    ff_original.make_xz_spectral()
    ff_mean.make_xz_spectral()


    #================================================================
    #### Stack velocity fields in the wall-normal direction
    #================================================================
    ff_original.stack_ff_in_y()
    ff_mean.stack_ff_in_y()


    #================================================================
    #### Create arrays of Fourier modes to use
    #================================================================
    # Modes multiplied with fundamental wavenumbers
    #(Modes: the physical modes, i.e. the grid points)
    kx_array = ff_original.Mx * ff_original.alpha
    kz_array = ff_original.Mz * ff_original.beta


    #================================================================
    #### Ensure valid rank is specified
    #================================================================
    rank = min(Rank, 3*ff_original.numModes)


    #================================================================
    #### Approximate the file w/regards to specified rank
    #================================================================
    deconstructed_field = cr.deconstruct_field(ff_original.velocityField,
                                              kx_array,
                                              kz_array,
                                              ff_original.numModes,
                                              ff_original.c,
                                              ff_original.Re,
                                              ff_original.baseflow,
                                              rank,
                                              mean_profile,
                                              sparse)


    approximated_ff_spectral = cr.construct_field(deconstructed_field['resolvent_modes'],
                                                  deconstructed_field['singular_values'],
                                                  deconstructed_field['coefficients'],
                                                  ff_mean.velocityField,
                                                  kx_array,
                                                  kz_array,
                                                  ff_original.numModes)


    #================================================================
    #### Initialize approximated flow field object
    #================================================================
    ff_approximated = ffClass.FlowFieldChannelFlow( file_info['Nd'],
                                                    file_info['Nx'],
                                                    ff_original.numModes, # the velocity field is missing wall boundaries
                                                    file_info['Nz'],
                                                    file_info['Lx'],
                                                    file_info['Lz'],
                                                    file_info['alpha'],
                                                    file_info['beta'],
                                                    details['c'],
                                                    details['bf'],
                                                    details['Re'],
                                                    approximated_ff_spectral,
                                                    "sp")


    #================================================================
    #### Unstack velocity fields in the wall-normal direction
    #================================================================
    if ff_approximated.is_stacked_in_y:
        ff_approximated.unstack_ff()

    if ff_mean.is_stacked_in_y:
        ff_mean.unstack_ff()

    if ff_original.is_stacked_in_y:
        ff_original.unstack_ff()


    #================================================================
    #### Inverse Fourier transform approximated and mean velocity fields in xz directions
    #================================================================
    ff_approximated.make_xz_physical()
    ff_mean.make_xz_physical()
    ff_original.make_xz_physical()


    #================================================================
    #### Add wall boundaries
    #================================================================
    ff_approximated.add_wall_boundaries()
    ff_mean.add_wall_boundaries()
    ff_original.add_wall_boundaries()

    test_u3 = ff_original.velocityField[0,:,:,:].real
    doriginal = test_u3.real - file_info['ff'][0,:,:,:].real
    dnorm = np.linalg.norm(doriginal)
#    if dnorm >= 1e-10:
#        message = "The original field has not been retrieved after FFT, stacking and removing wall boundaries and then reversing those operations... ||delta||" + str(dnorm)
#        ut.error(message)
        
    test_u3 = ff_mean.velocityField[0,:,:,:].real
    dmean = test_u3.real - mean[0,:,:,:].real
    dnorm = np.linalg.norm(dmean)
#    if dnorm >= 1e-10:
#        message = "The original mean has not been retrieved after FFT, stacking and removing wall boundaries and then reversing those operations... ||delta||" + str(dnorm)
#        ut.error(message)


    #================================================================
    #### Retrieve approximated fluctuations
    #================================================================  
    tmp = ff_approximated.velocityField[:,:,:,:] # - mean.real The mean is not set...
    app_u = tmp[0,:,:,:].real
    app_v = tmp[1,:,:,:].real
    app_w = tmp[2,:,:,:].real
    org_u = ff_original.velocityField[0,:,:,:].real
    org_v = ff_original.velocityField[1,:,:,:].real
    org_w = ff_original.velocityField[2,:,:,:].real
    ff_approximated.set_ff(tmp.real, "pp")


    #================================================================
    # Create a folder to store the approximated velocity field in
    #================================================================
    os.chdir(parent_directory) # Go up one directory
    rank_folder = File[:-3]+"_rank_" + str(rank) + "/"
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
        #------------------------------------------------
        # Write physical ascii file
        #------------------------------------------------
        fileName = File[:-3] + "_rnk_" + str(rank)
        ut.write_ASC(ff_approximated, rank_folder, fileName)

        #------------------------------------------------
        # Write binary file
        #------------------------------------------------
        command = "ascii2field -p false -ge ../rank-temp/" + str(File)[:-3] + ".geom " + fileName + ".asc " + fileName
        print(command)
        os.system(command)

        #------------------------------------------------
        # Write amplitude coefficients for each Fourier mode combination
        #------------------------------------------------
        fileName = File[:-3] + "_coeffs"
        ut.write_amplitude_coefficients(ff_approximated, rank_folder, fileName, deconstructed_field['coefficients'])

        #------------------------------------------------
        # Remove ascii file and temporary folder
        #------------------------------------------------
    #    os.system("rm *.asc")
        os.chdir(parent_directory)
        command = "rm -rf " + temp_rank_folder
        os.system(command)
    
    print("\nFinished")




dirc = "/home/arslan/Desktop/ati-modes-copy/modes/A/kz2/ff_files/"
os.chdir(dirc)
main("mode-00.ff",
     2,
     dirc, 
     "turbdeviation.txt",
     True)
