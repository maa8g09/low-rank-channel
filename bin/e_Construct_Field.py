#!/usr/bin/env python3
"""
=====================================================================
Construct Field
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
import numpy as np
import sys
sys.path.append("/".join(sys.path[0].split("/")[:-1]))
import ChannelResolvent as cr
import FlowField as ffClass
import Utils as ut
parser = argparse.ArgumentParser(description="Construct an (approximate)field at a given rank with a given mean profile (if specified).")
parser.add_argument("-f", "--File", help="File to re-construct.",
                    metavar='\b', required=True)
parser.add_argument("-r", "--Rank", help="Number of velocity modes to use to recreate flow field.",
                    metavar='\b', required=True, type=int)
parser.add_argument("-v", "--MeanProfile", help="Turbulent mean velocity profile as .txt file. Keep in same directory as file to approximate.",
                    metavar='\b')
args = parser.parse_args()
#===================================================================#
#### Read deconstructed field                                    ####
#===================================================================#
deconstructed_field = ut.read_H5_Deconstructed(args.File)
#===================================================================#
#### Create arrays of Fourier modes to use                       ####
#===================================================================#
# Modes multiplied with fundamental wavenumbers
#(Modes: the physical modes, i.e. the grid points)
kx_array = deconstructed_field['Mx'] * deconstructed_field['alpha']
kz_array = deconstructed_field['Mz'] * deconstructed_field['beta']
#===================================================================#
#### Ensure valid rank is specified                              ####
#===================================================================#
numModes = deconstructed_field['Ny'] - 2
rank = min(args.Rank, 3*numModes)
if rank == 3*numModes:
    rankStr = "full"
else:
    rankStr = str(rank).zfill(2) 
#===================================================================#
#### Check velocity profile                                      ####
#===================================================================#
# Create empty 4D array to store mean flow field
mean = np.zeros((deconstructed_field['Nd'], deconstructed_field['Nx'], 
                 deconstructed_field['Ny'], deconstructed_field['Nz']))
spectral_deviation_profile = []
deviation_profile = []
if args.MeanProfile: # Velocity profile given
    #----------------------------------------------------------------
    #### Read velocity profile                                      -
    #----------------------------------------------------------------
    vel_profile = ut.read_Vel_Profile(args.MeanProfile)
    # Check to see if it is a mean profile or a deviation profile.
    deviation = any(n < 0 for n in vel_profile)
    if deviation: # Deviation profile given
        deviation_profile = vel_profile
    else: # Turbulent mean profile given
        # Remove baseflow from turbulent mean profile
        baseflow_profile = []
        if deconstructed_field['bf'] == "lam": # Laminary base flow
            baseflow_profile = 1.0 - deconstructed_field['y']**2.0
        elif deconstructed_field['bf'] == "cou": # Couette base flow
            baseflow_profile = deconstructed_field['y']
        # Remove baseflow from mean to get turbulent deviation profile
        deviation_profile = vel_profile - np.asarray(baseflow_profile)    
    deviation_profile_sp = np.fft.fft(deviation_profile)
    spectral_deviation_profile = np.zeros((3*numModes),dtype=complex)
    spectral_deviation_profile[1:numModes] = deviation_profile_sp[1:numModes]
#===================================================================#
#### Reconstruct approximated flow field                         ####
#===================================================================#
approximated_ff_spectral = cr.construct_field(deconstructed_field['resolvent_modes'],
                                              deconstructed_field['singular_values'],
                                              deconstructed_field['coefficients'],
                                              kx_array,
                                              kz_array,
                                              numModes,
                                              rank,
                                              spectral_deviation_profile)
#===================================================================#
#### Initialize approximated flow field object                   ####
#===================================================================#
ff_approximated = ffClass.FlowFieldChannelFlow(deconstructed_field['Nd'],
                                               deconstructed_field['Nx'],
                                               numModes, # the velocity field is missing wall boundaries
                                               deconstructed_field['Nz'],
                                               deconstructed_field['Lx'],
                                               deconstructed_field['Lz'],
                                               deconstructed_field['alpha'],
                                               deconstructed_field['beta'],
                                               deconstructed_field['c'],
                                               deconstructed_field['bf'],
                                               deconstructed_field['Re'],
                                               approximated_ff_spectral,
                                               "sp")
#===================================================================#
#### ---- Unstack velocity field in the wall-normal direction    ####
#===================================================================#
ff_approximated.unstack_ff()
#===================================================================#
#### ---- Inverse Fourier transform approximated field in xz directions
#===================================================================#
ff_approximated.make_xz_physical()
#===================================================================#
#### ---- Add wall boundaries                                    ####
#===================================================================#
ff_approximated.add_wall_boundaries()
#===================================================================#
#### Write approximated HDF5 file                                ####
#===================================================================#
prefix = args.File.split("deconstructed")[0]
fileName = prefix + "rank_" + rankStr
# Make a folder to put it in
folder = ut.make_Folder(os.getcwd(), fileName, True)
os.chdir(folder)
ut.write_H5(ff_approximated, deconstructed_field['original_attrs'],fileName)
ut.print_EndMessage()