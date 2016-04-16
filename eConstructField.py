#!/usr/bin/env python3
import argparse
import os
import numpy as np
import ChannelResolvent as cr
import FlowField as ffClass
import Tests
import Utils as ut
parser = argparse.ArgumentParser(description="Construct an (approximate)field at a given rank with a given mean profile (if specified).")
parser.add_argument("-f",
                    "--File",
                    metavar='\b',
                    help="File to re-construct.",
                    required=True)
parser.add_argument("-r",
                    "--Rank",
                    metavar='\b',
                    help="Number of velocity modes to use to recreate flow field.",
                    required=True,
                    type=int)
parser.add_argument("-v",
                    "--MeanProfile",
                    metavar='\b',
                    help="Turbulent mean velocity profile as .txt file. Keep in same directory as file to approximate.")
args = parser.parse_args()
#================================================================
#### Read deconstructed field
#================================================================
deconstructed_field = ut.read_H5_Deconstructed(args.File)
#================================================================
#### Create arrays of Fourier modes to use
#================================================================
# Modes multiplied with fundamental wavenumbers
#(Modes: the physical modes, i.e. the grid points)
kx_array = deconstructed_field['Mx'] * deconstructed_field['alpha']
kz_array = deconstructed_field['Mz'] * deconstructed_field['beta']
#================================================================
#### Ensure valid rank is specified
#================================================================
numModes = deconstructed_field['Ny'] - 2
rank = min(args.Rank, 3*numModes)
#================================================================
#### Check velocity profile
#================================================================
# Create empty 4D array to store mean flow field
mean = np.zeros((deconstructed_field['Nd'], deconstructed_field['Nx'], 
                 deconstructed_field['Ny'], deconstructed_field['Nz']))
mean_profile = []
if args.MeanProfile: # Velocity profile given
    #------------------------------------------------------------
    #### Read velocity profile
    #------------------------------------------------------------
    vel_profile = ut.read_Vel_Profile(args.MeanProfile)
    # Check to see if it is a mean profile or a deviation profile.
    deviation = any(n < 0 for n in vel_profile)
    if deviation: # Deviation profile given
        # Add baseflow to deviation
        baseflow = []
        if deconstructed_field['bf'] == "lam": # Laminary base flow
            baseflow = 1.0 - deconstructed_field["y"].y**2.0
        elif deconstructed_field['bf'] == "cou": # Couette base flow
            baseflow = deconstructed_field["y"]
        # Add baseflow to deviation to get turbulent mean profile
        mean_profile = vel_profile + np.asarray(baseflow)
    else: # Turbulent mean profile given
        mean_profile = vel_profile
#================================================================
#### Reconstruct approximated flow field
#================================================================
approximated_ff_spectral = cr.construct_field(deconstructed_field['resolvent_modes'],
                                              deconstructed_field['singular_values'],
                                              deconstructed_field['coefficients'],
                                              kx_array,
                                              kz_array,
                                              numModes,
                                              rank,
                                              mean_profile,
                                              deconstructed_field['bf'],
                                              deconstructed_field["y"])
#================================================================
#### Initialize approximated flow field object
#================================================================
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
#================================================================
#### ---- Unstack velocity field in the wall-normal direction
#================================================================
ff_approximated.unstack_ff()
#================================================================
#### ---- Inverse Fourier transform approximated field in xz directions
#================================================================
ff_approximated.make_xz_physical()
#================================================================
#### ---- Add wall boundaries
#================================================================
ff_approximated.add_wall_boundaries()
#================================================================
#### Write approximated HDF5 file
#================================================================
fileName = args.File[:-3] + "_rank_" + str(rank)
ut.write_H5(ff_approximated, fileName) # need to write attributes to file...
print("Finished\n")