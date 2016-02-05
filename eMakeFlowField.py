#!/usr/bin/env python

# This is an executable file that should be run from the command line.
# The command line arguments determine the type of flow field created.


## Need to add all the usual things about yourself. Name, Institute etc.


import os
import argparse
import time
import FlowField as ffClass
import ChannelResolvent as cr
import Utils as ut
date = time.strftime("%Y_%m_%d")


####################################################################################################
# Parse the command line arguments (flag parameters)
#ut.print_ResolventHeader()
#ut.print_ResolventSubHeader()
parser = argparse.ArgumentParser(description="Create a flow field using the resolvent formulation.")
parser.add_argument("-wp",
                    "--Wavepacket",
#                    metavar='\b',
                    help="The wavepacket you want to create.",
                    required=True,
                    choices=['KA', 'KB', 'KC', 'KD', 'KE'])
#                    nargs=1
parser.add_argument("-nd",
                    "--Nd",
                    metavar='\b',
                    help="The dimensions of the flow field you want to create.",
                    required=True,
                    type=int)
parser.add_argument("-nx",
                    "--Nx",
                    metavar='\b',
                    help="Number of grid point in the streamwise direction.",
                    required=True,
                    type=int)
parser.add_argument("-ny",
                    "--Ny",
                    metavar='\b',
                    help="Number of grid point in the wall-normal direction.",
                    required=True,
                    type=int)
parser.add_argument("-nz",
                    "--Nz",
                    metavar='\b',
                    help="Number of grid point in the spanwise direction.",
                    required=True,
                    type=int)
parser.add_argument("-re",
                    "--Reynolds",
                    metavar='\b',
                    help="Reynolds number to generate the flow field.",
                    required=True,
                    type=float)
parser.add_argument("-c",
                    "--Wavespeed",
                    metavar='\b',
                    help="Wave speed of the flow field.",
                    required=True,
                    type=float)
parser.add_argument("-th",
                    "--theta",
                    metavar='\b',
                    help="Scaling of the amplitude coefficients where the scaled amplitude coefficient defined as chi_tilde = 10**(theta) * chi.",
                    required=True,
                    type=float)
parser.add_argument("-bf",
                    "--Baseflow",
#                    metavar='\b',
                    help="The baseflow to use.",
                    required=True,
                    choices=['lam', 'cou'])
parser.add_argument("-d",
                    "--Directory",
                    metavar='\b',
                    help="Output directory.",
                    required=True)
parser.add_argument("-i",
                    "--Iteration",
                    metavar='\b',
                    help="Iteration.",
                    type=int)
parser.add_argument("-dat",
                    "--SaveDAT",
                    help="Save flow field in DAT format.",
                    action='store_true')
parser.add_argument("-asc",
                    "--SaveASC",
                    help="Save flow field in ASCII format with an accompanying geometry file.",
                    action='store_true')
parser.add_argument("-h5",
                    "--SaveHDF5",
                    help="Save flow field in HDF5 format.",
                    action='store_true')
parser.add_argument("-ff",
                    "--SaveFF",
                    help="Save flow field in binary '.ff' format.",
                    action='store_true')
args = parser.parse_args()

####################################################################################################
# Pass these variables to the FlowField class to construct an instance of the geometry
flowFieldGeometry = ffClass.FlowFieldGeometry(args.Baseflow,
                                              args.Wavepacket,
                                              args.Nd,
                                              args.Nx,
                                              args.Ny,
                                              args.Nz,
                                              args.Reynolds,
                                              args.Wavespeed,
                                              args.theta)

####################################################################################################
# Make the correct directory to store the flow field in
if args.Iteration:
    print(args.Iteration)
    ut.print_ResolventSubHeader()
    output_directory = ut.make_FlowField_output_directory_wIteration(args.Directory, flowFieldGeometry, date, args.Iteration)
else:
    ut.print_ResolventHeader()
    ut.print_ResolventSubHeader()
    output_directory = ut.make_FlowField_output_directory(args.Directory, flowFieldGeometry, date)


####################################################################################################
# Pass the FlowField class to main_resolvent (this is where the flow field is generated)
ff = cr.resolvent_formulation(flowFieldGeometry)
ff = ffClass.FlowField(flowFieldGeometry, ff, "pp")
####################################################################################################
# Save flow field to output_directory
if args.SaveDAT:
    ut.write_DAT(ff, output_directory, "u0")

if args.SaveASC:
    ut.write_ASC(ff, output_directory, "u0")
    ut.write_GEOM(ff, output_directory, "u0")
    
    if args.SaveFF:
        ut.write_FF(output_directory, "u0")

#if args.saveHDF5:
#    ut.writeHDF5(ff, flowFieldGeometry, output_directory)

# Save details of the flow field in a file called u0_details.txt
ut.write_Details(ff, output_directory, "u0")

print("\nThe flow field has been stored at")
print(output_directory)
ut.print_EndMessage()
