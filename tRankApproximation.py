#!/usr/bin/env python

import os
import sys
import h5py
import Utils as ut
import FlowField as ffClass
import ChannelResolvent as cr
import Tests

import numpy as np


def main(File, Rank, Directory, MeanProfile, Sparse, Testing):

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
    if File[-3:] == ".h5":
        print("HDF5 file given.")
        #------------------------------------------------
        #### Read the HDF5 and details file
        file_info = ut.read_H5(parent_directory, File)
        details = ut.read_Details(parent_directory, "u0_Details.txt")
    
    
    elif File[-3:] == ".ff":
        #------------------------------------------------
        #### Convert the binary file to ascii
        command = "field2ascii -p ../" + str(File) + " " + str(File)[:-3]
        print(command)
        os.system(command)
        #------------------------------------------------
        #### Read ASCII file and details file
        file_info = ut.read_ASC_channelflow(temp_rank_folder, str(File)[:-3])
        details = ut.read_Details(parent_directory, "u0_Details.txt")
    
    
    elif File[-3:] == "asc":
        #------------------------------------------------
        #### Read ASCII file and details file
        file_info = ut.read_ASC_PP(parent_directory, str(File)[:-7])
        details = ut.read_Details(parent_directory, "u0_Details.txt")
    
    
    else: # No file type given.
        ut.error("Invalid file given.")
    
    
    #================================================================
    #### Initialise flow field object for field (to approximate)
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
    
    Tests.fft_ifft(ff_original)
    
    
    #================================================================
    #### Check velocity profile
    #================================================================
    # Create empty 4D array to store mean flow field
    mean = np.zeros((file_info['Nd'], file_info['Nx'], file_info['Ny'], file_info['Nz']))
    mean_profile = []
    
    if MeanProfile: # Velocity profile given
        #------------------------------------------------
        #### Read velocity profile
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
    
            # Add baseflow to deviation to get turbulent mean profile
            mean_profile = vel_profile + np.asarray(baseflow)
    
        else: # Turbulent mean profile given
            mean_profile = vel_profile
    
        #------------------------------------------------
        #### Construct 4D array from mean_profile
        #------------------------------------------------
        mean = ut.make_ff_from_profile(vel_profile, 
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
    #### ---- Remove the wall boundaries
    #================================================================
    # Removing the xz-planes at y=1 and y=-1,
    # so that the chebyshev nodes can be used to construct 
    # the transfer function.
    ff_original.remove_wall_boundaries()
    ff_mean.remove_wall_boundaries()
    
    
    #================================================================
    #### ---- Fourier transform in xz directions
    #================================================================
    ff_original.make_xz_spectral()
    ff_mean.make_xz_spectral()
    
    
    #================================================================
    #### ---- Stack velocity fields in the wall-normal direction
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
    
        
#    Tests.checkHermitianSymmetry(ff_original.velocityField[:,0:ff_original.numModes,:],
#                                 ff_original.Nx, 
#                                 ff_original.Nz)
#    Tests.checkHermitianSymmetry(ff_original.velocityField[:,ff_original.numModes:ff_original.numModes*2,:],
#                                 ff_original.Nx, 
#                                 ff_original.Nz)
#    Tests.checkHermitianSymmetry(ff_original.velocityField[:,2*ff_original.numModes:ff_original.numModes*3,:],
#                                 ff_original.Nx, 
#                                 ff_original.Nz)

    
    #================================================================
    #### Ensure valid rank is specified
    #================================================================
    rank = min(Rank, 3*ff_original.numModes)
    
    
    
    #### -!-!- TESTING -!-!-:   Zeroth mode differences
    meanFF = ff_mean.velocityField
    origFF = ff_original.velocityField
    diffs = meanFF[0, :, 0] - origFF[0, :, 0]
    
    diffsn = np.linalg.norm(diffs)
    
    
    
    
    #### -!-!- TESTING -!-!-:   Synthesizing a Fourier domain flow field
    fake_field = ff_original.velocityField
#    fake_field *= 0.0+0.0j
#    fake_field[1,:,1] += 1.0+2.0j
#    fake_field[1,:,-1] += 1.0-2.0j
    
    
    
    #================================================================
    #### Deconstruct original flow field
    #================================================================
    deconstructed_field = cr.deconstruct_field_testing(fake_field,
                                                      kx_array,
                                                      kz_array,
                                                      ff_original.numModes,
                                                      ff_original.c,
                                                      ff_original.Re,
                                                      ff_original.baseflow,
                                                      mean_profile,
                                                      Sparse,
                                                      False)


    #================================================================
    #### Reconstruct approximated flow field
    #================================================================
    approximated_ff_spectral = cr.construct_field_testing(deconstructed_field['resolvent_modes'],
                                                          deconstructed_field['singular_values'],
                                                          deconstructed_field['coefficients'],
        #                                                  ff_mean.velocityField, # Use mean field at zeroth modes
                                                          fake_field, # Use original field at zeroth modes
                                                          kx_array,
                                                          kz_array,
                                                          ff_original.numModes,
                                                          rank,
                                                          fake_field)






    #### -!-!- TESTING -!-!-:   Synthesizing a Fourier domain flow field
    # The retrieved field should be the same as the fak_field...
    retrieved_difference = approximated_ff_spectral - fake_field
    retrieved_difference_n = np.linalg.norm(retrieved_difference)





    #### -!-!- TESTING -!-!-:   Full-Rank decomposition and recomposition (Fourier Domain)
    meanFF = ff_mean.velocityField
    origFF = ff_original.velocityField
    difference = np.linalg.norm(ff_original.velocityField[1:,:,1:] - approximated_ff_spectral[1:,:,1:])
    
    
    #### -!-!- TESTING -!-!-:   Zeroth mode differences
    difference2 = np.linalg.norm( approximated_ff_spectral[0,:,0] - ff_original.velocityField[0,:,0])
    difference3 = np.linalg.norm( approximated_ff_spectral[0,:,0] - meanFF[0, :, 0])
#    diffsn and difference2 should be the same...
    print("")
    print("The norm of the difference is " + str(difference))
    print("")


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
    

    
    
    # If not symmetric: you need to filter the negative frequencies in the approximated result.
    # This will introduce hermitian symmetry.


    #================================================================
    #### ---- Unstack velocity fields in the wall-normal direction
    #================================================================
    ff_approximated.unstack_ff()
    ff_mean.unstack_ff()
    ff_original.unstack_ff()
    
    
    #================================================================
    #### ---- Inverse Fourier transform approximated and mean velocity fields in xz directions
    #================================================================
    ff_approximated.make_xz_physical()
    ff_mean.make_xz_physical()
    ff_original.make_xz_physical()

    
    #================================================================
    #### ---- Add wall boundaries
    #================================================================
    ff_approximated.add_wall_boundaries()
    ff_mean.add_wall_boundaries()
    ff_original.add_wall_boundaries()
    

    #================================================================
    #### Remove mean flow from approximated field (if mean used to reconstruct)
    #================================================================
#    approximated_field = ff_approximated.velocityField.real - ff_mean.velocityField.real
#    ff_approximated.set_ff(approximated_field, "pp")
#
#    difference = np.linalg.norm(approximated_field - ff_approximated.velocityField)



    #### -!-!- TESTING -!-!-:   Full-Rank decomposition and recomposition: (Real Domain)
    u_a = ff_approximated.velocityField[0,:,:,:].real
    u_o = ff_original.velocityField[0,:,:,:].real
    difference = np.linalg.norm(ff_original.velocityField.real - ff_approximated.velocityField.real)
    print("")
    print("The norm of the difference is " + str(difference))
    print("")
    
    #### -!-!- TESTING -!-!-:   Projections
    if Testing:
        #### Remove boundaries
        ff_approximated.remove_wall_boundaries()
        #### FFT
        ff_approximated.make_xz_spectral()
        #### stack in y
        ff_approximated.stack_ff_in_y()
        #### Deconstruct the approximated flow field.
        deconstructed_approximated = cr.deconstruct_field(ff_approximated.velocityField,
                                                          kx_array,
                                                          kz_array,
                                                          ff_approximated.numModes,
                                                          ff_approximated.c,
                                                          ff_approximated.Re,
                                                          ff_approximated.baseflow,
                                                          rank,
                                                          mean_profile,
                                                          Sparse,
                                                          False)
        #### Reconstruct projected approximated flow field.
        doubly_approximated_ff_spectral = cr.construct_field( deconstructed_approximated['resolvent_modes'],
                                                              deconstructed_approximated['singular_values'],
                                                              deconstructed_approximated['coefficients'],
                #                                             ff_mean.velocityField,
                                                              kx_array,
                                                              kz_array,
                                                              ff_approximated.numModes)

        #### Initialize projected approximated flow field object.
        ff_doubly_approximated = ffClass.FlowFieldChannelFlow(  file_info['Nd'],
                                                                file_info['Nx'],
                                                                ff_approximated.numModes, # the velocity field is missing wall boundaries
                                                                file_info['Nz'],
                                                                file_info['Lx'],
                                                                file_info['Lz'],
                                                                file_info['alpha'],
                                                                file_info['beta'],
                                                                details['c'],
                                                                details['bf'],
                                                                details['Re'],
                                                                doubly_approximated_ff_spectral,
                                                                "sp")

        #### unstack
        ff_doubly_approximated.unstack_ff()
        ff_approximated.unstack_ff()
        #### IFFT
        ff_doubly_approximated.make_xz_physical()
        ff_approximated.make_xz_physical()
        #### add wall boundaries
        ff_doubly_approximated.add_wall_boundaries()
        ff_approximated.add_wall_boundaries()


        #### check differences between modes/sing vals/coeffs
        del_sigma = deconstructed_approximated['singular_values'] - deconstructed_field['singular_values']
        del_sigmaR = del_sigma.real
        del_sigmaI = del_sigma.imag
        del_sigma_n = np.linalg.norm(del_sigma)
        print("\n\n||sigma_tilde - sigma|| = " + str(del_sigma_n))
    
        del_psi = deconstructed_approximated['resolvent_modes'] - deconstructed_field['resolvent_modes']
        del_psiR = del_psi.real
        del_psiI = del_psi.imag
        del_psi_n = np.linalg.norm(del_psi)
        print("\n\n||psi_tilde - psi|| = " + str(del_psi_n))
    
        del_xi = deconstructed_approximated['coefficients'] - deconstructed_field['coefficients']
        del_xiR = del_xi.real
        del_xiI = del_xi.imag
        del_xi_n = np.linalg.norm(del_xi)
        print("\n\n||xi_tilde - xi|| = " + str(del_xi_n))
        
        
        del_ff = ff_doubly_approximated.velocityField.real - ff_approximated.velocityField.real
        dim_a = ff_approximated.velocityField.shape
        print("\nDimension of approximated ff: " + str(dim_a))

        dim_da = ff_doubly_approximated.velocityField.shape
        print("\nDimension of doubly approximated ff: " + str(dim_da))

        del_ff_n = np.linalg.norm(del_ff)
        print("\n\n||u_tilde_tilde - u_tilde|| = " + str(del_ff_n))
    
    
    
    
    
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
    
    
    
    
    elif File[-3:] == "asc":
        #------------------------------------------------
        # Write physical ascii file
        #------------------------------------------------
        fileName = File[:-3] + "_rnk_" + str(rank)
        ut.write_ASC_Py(ff_approximated, rank_folder, fileName)
    
    
    
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




dirc="/home/arslan/Documents/work/cfd-symmetry_scans/s_tw1_sigma_z_tau_x/Re600.0/KB/2016_02_18/002_theta_-0.5000/data-alt"

os.chdir(dirc)
fileName = "u975.000.h5"
vel_profile_file = "turbulent_deviation950-999.txt"
testing=False
sparse=True
main(fileName,
     20,
     dirc,
     vel_profile_file,
     sparse,
     testing)
