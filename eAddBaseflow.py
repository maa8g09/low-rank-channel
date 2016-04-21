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
                    help="File to add to.",
                    required=True)
parser.add_argument("-bf",
                    "--Baseflow",
#                    metavar='\b',
                    help="The baseflow to add.",
                    choices=['lam', 'cou'])
parser.add_argument("-p",
                    "--Profile",
                    metavar='\b',
                    help="The mean/deviation/baseflow profile to add")
args = parser.parse_args()
#================================================================
#### Read the HDF5 and details file
#================================================================
file_info, original_attrs = ut.read_H5(args.File)
#================================================================
#### Declare flow field object
#================================================================
ff_original = ffClass.FlowFieldChannelFlow( file_info['Nd'],
                                            file_info['Nx'],
                                            file_info['Ny'],
                                            file_info['Nz'],
                                            file_info['Lx'],
                                            file_info['Lz'],
                                            file_info['alpha'],
                                            file_info['beta'],
                                            0.0,
                                            args.Baseflow,
                                            0.0,
                                            file_info['ff'],
                                            "pp")
#================================================================
#### Calculate mean/deviation/baseflow profile
#================================================================
suffix = ""
profile = []
if args.Baseflow: # given baseflow flag
    if args.Baseflow == "lam": # Laminary base flow
        profile = 1.0 - ff_original.y**2.0
        suffix = "bf_lam"
    elif args.Baseflow == "cou": # Couette base flow
        profile = ff_original.y
        suffix = "bf_cou"
elif args.Profile: # given velocity profile file.
    profile = ut.read_Vel_Profile(args.Profile)
    deviation = any(n < 0 for n in profile)
    if deviation:
        suffix = "turb_deviation"
    else:
        suffix = "turb_mean"
profile_ff = ut.make_ff_from_profile(profile,
                                     ff_original.Nd,
                                     ff_original.Nx,
                                     ff_original.Nz)
total_ff = profile_ff + ff_original.velocityField
ff_original.set_ff(total_ff, "pp")
fileName = args.File[:-3] + "_with_" + suffix
ut.write_H5(ff_original, original_attrs, fileName)
ut.print_EndMessage()