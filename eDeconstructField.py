#!/usr/bin/env python3
import argparse
import ChannelResolvent as cr
import FlowField as ffClass
import Tests
import Utils as ut
import numpy as np
parser = argparse.ArgumentParser(description="Project modes onto a velocity field to get a rank approximation.")
parser.add_argument("-f",
                    "--File",
                    metavar='\b',
                    help="File to approximate",
                    required=True)
parser.add_argument("-d",
                    "--Details",
                    metavar='\b',
                    help="Details file with Re, c and bf info as .txt file.",
                    required=True)
parser.add_argument("-v",
                    "--MeanProfile",
                    metavar='\b',
                    help="Turbulent mean velocity profile as .txt file. Keep in same directory as file to approximate.")
parser.add_argument("-s",
                    "--Sparse",
                    help="Use sparse SVD algorithm.",
                    action='store_true')
args = parser.parse_args()
#================================================================
#### Format the output directory
#================================================================
if args.File[-3:] == ".h5":
#    print("HDF5 file given.")
    #------------------------------------------------------------
    #### Read the HDF5 and details file
    #------------------------------------------------------------
    file_info, original_attrs = ut.read_H5(args.File)
    details = ut.read_Details(args.Details)
elif args.File[-3:] == ".ff":
#    print("Channelflow binary file given.")
    #------------------------------------------------------------
    #### Convert from .ff to .h5
    #------------------------------------------------------------
    command = "fieldconvert " + args.File + " " + args.File[:-3] + ".h5"
    #------------------------------------------------------------
    #### Read the HDF5 and details file
    #------------------------------------------------------------
    file_info = ut.read_H5(args.File)
    details = ut.read_Details(args.Details)
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
        if details['bf'] == "lam": # Laminary base flow
            baseflow = 1.0 - ff_original.y**2.0
        elif details['bf'] == "cou": # Couette base flow
            baseflow = ff_original.y
        # Add baseflow to deviation to get turbulent mean profile
        mean_profile = vel_profile + np.asarray(baseflow)
    else: # Turbulent mean profile given
        mean_profile = vel_profile
#================================================================
#### ---- Remove the wall boundaries
#================================================================
# Removing the xz-planes at y=1 and y=-1,
# so that the chebyshev nodes can be used to construct 
# the transfer function.
ff_original.remove_wall_boundaries()
#================================================================
#### ---- Fourier transform in xz directions
#================================================================
ff_original.make_xz_spectral()
#================================================================
#### ---- Stack velocity fields in the wall-normal direction
#================================================================
ff_original.stack_ff_in_y()
#================================================================
#### Create arrays of Fourier modes to use
#================================================================
# Modes multiplied with fundamental wavenumbers
#(Modes: the physical modes, i.e. the grid points)
kx_array = ff_original.Mx * ff_original.alpha
kz_array = ff_original.Mz * ff_original.beta
#================================================================
#### Deconstruct original flow field
#================================================================
deconstructed_field = cr.deconstruct_field(ff_original.velocityField,
                                           kx_array,
                                           kz_array,
                                           ff_original.numModes,
                                           ff_original.y,
                                           ff_original.c,
                                           ff_original.Re,
                                           ff_original.baseflow,
                                           mean_profile,
                                           args.Sparse)
#================================================================
#### Write decomposed flow field to disk as an HDF5 file
#================================================================
ut.write_H5_Deconstructed(deconstructed_field, original_attrs, ff_original, args.File[:-3])
#print("\nFinished\n")