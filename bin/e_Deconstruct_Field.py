#!/usr/bin/env python3
"""
=====================================================================
Deconstruct Field
=====================================================================
The resolvent formulation is used to decompose a given velocity
field into its resolvent modes, forcing modes, singular values and
amplitude coefficients.

The decomposed field's components are saved in an HDF5 file
with the following structure:
    deconstructed_field/:
        resolvent_modes
        forcing_modes
        singular_values
        coefficients

The decomposition of the flow field takes place in the Fourier
domain. The resolvent formulation calculates a transfer function,
which is singular value decomposed into resolvent modes, singular
values and forcing modes.

These vectors are then projected onto the Fourier transformed
velocity field, yielding the amplitude coefficients. 

The resolvent modes, singular values and amplitude coefficients
can be truncated and used to construct an approximation of the 
original velocity field, see e_Construct_Field.py


Muhammad Arslan Ahmed
maa8g09@soton.ac.uk
Room 5069
Building 13
Aerodynamics and Flight Mechanics
University of Southampton

"""
import argparse
import numpy as np
import sys
sys.path.append("/".join(sys.path[0].split("/")[:-1]))
import ChannelResolvent as cr
import FlowField as ffClass
import Tests
import Utils as ut
parser = argparse.ArgumentParser(description="Project modes onto a velocity field to get a rank approximation.")
parser.add_argument("-f", "--File", help="File to approximate",
                    metavar='\b', required=True)
parser.add_argument("-d", "--Details", help="Details file with Re, c and bf info as .txt file.",
                    metavar='\b', required=True)
parser.add_argument("-v", "--MeanProfile", help="Turbulent mean velocity profile as .txt file. Keep in same directory as file to approximate.",
                    metavar='\b')
parser.add_argument("-s", "--Sparse", help="Use sparse SVD algorithm.",
                    action='store_true')
args = parser.parse_args()
#===================================================================#
#### Format the output directory                                 ####
#===================================================================#
if args.File[-3:] == ".h5":
#    print("HDF5 file given.")
    #----------------------------------------------------------------
    #### Read the HDF5 and details file                             -
    #----------------------------------------------------------------
    file_info, original_attrs = ut.read_H5(args.File)
    details = ut.read_Details(args.Details)
elif args.File[-3:] == ".ff":
#    print("Channelflow binary file given.")
    #----------------------------------------------------------------
    #### Convert from .ff to .h5                                    -
    #----------------------------------------------------------------
    command = "fieldconvert " + args.File + " " + args.File[:-3] + ".h5"
    #----------------------------------------------------------------
    #### Read the HDF5 and details file                             -
    #----------------------------------------------------------------
    file_info = ut.read_H5(args.File)
    details = ut.read_Details(args.Details)
else: # No file type given.
    ut.error("Invalid file given.")
#===================================================================#
#### Initialise flow field object for field (to approximate)     ####
#===================================================================#
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
#===================================================================#
#### Check velocity profile                                      ####
#===================================================================#
# Create empty 4D array to store mean flow field
mean = np.zeros((file_info['Nd'], file_info['Nx'], file_info['Ny'], file_info['Nz']))
mean_profile = []
if args.MeanProfile: # Velocity profile given
    #----------------------------------------------------------------
    #### Read velocity profile                                      -
    #----------------------------------------------------------------
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
#===================================================================#
#### ---- Remove the wall boundaries                             ####
#===================================================================#
# Removing the xz-planes at y=1 and y=-1,
# so that the chebyshev nodes can be used to construct 
# the transfer function.
ff_original.remove_wall_boundaries()
#===================================================================#
#### ---- Fourier transform in xz directions                     ####
#===================================================================#
ff_original.make_xz_spectral()
#===================================================================#
#### ---- Stack velocity fields in the wall-normal direction     ####
#===================================================================#
ff_original.stack_ff_in_y()
#===================================================================#
#### Create arrays of Fourier modes to use                       ####
#===================================================================#
# Modes multiplied with fundamental wavenumbers
#(Modes: the physical modes, i.e. the grid points)
kx_array = ff_original.Mx * ff_original.alpha
kz_array = ff_original.Mz * ff_original.beta
#===================================================================#
#### Deconstruct original flow field                             ####
#===================================================================#
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
#===================================================================#
#### Write decomposed flow field to disk as an HDF5 file         ####
#===================================================================#
ut.write_H5_Deconstructed(deconstructed_field, original_attrs, ff_original, args.File[:-3])
ut.print_EndMessage()