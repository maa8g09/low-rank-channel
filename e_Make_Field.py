#!/usr/bin/env python3
import argparse
import ChannelResolvent as cr
import FlowField as ffClass
import numpy as np
import time
import Utils as ut
parser = argparse.ArgumentParser(description="Create a flow field using the resolvent formulation.")
parser.add_argument("-nx", "--Nx", help="Number of grid point in the streamwise direction.",
                    metavar='\b', required=True, type=int)
                    
parser.add_argument("-ny", "--Ny", help="Number of grid point in the wall-normal direction.",
                    metavar='\b', required=True, type=int)
                    
parser.add_argument("-nz", "--Nz", help="Number of grid point in the spanwise direction.",
                    metavar='\b', required=True, type=int)

parser.add_argument("-lx", "--Lx", help="Streamwise (x) box length",
                    metavar='\b', required=True, type=float)

parser.add_argument("-lz", "--Lz", help="Spanwise (z) box length",
                    metavar='\b', required=True, type=float)
                    
parser.add_argument("-re", "--Reynolds", help="Reynolds number.",
                    metavar='\b', required=True, type=float)
                    
parser.add_argument("-c", "--Wavespeed", help="Wave speed of flow field.",
                    metavar='\b', required=True, type=float)
                    
parser.add_argument("-th", "--Theta", help="Scaling of the amplitude coefficients (chi),\ndefined as chi_tilde = 10**(theta) * chi.",
                    metavar='\b', required=True, type=float)
                    
parser.add_argument("-wp", "--Wavepacket", help="The wavepacket to create.",
                                  required=True, choices=['KA', 'KB', 'KC', 'KD', 'KE'])
                    
parser.add_argument("-bf", "--Baseflow", help="The baseflow to use.",
                                  required=True, choices=['lam', 'cou'])
                                  
parser.add_argument("-d", "--Directory", help="Output directory.",
                    metavar='\b', required=True)
                    
parser.add_argument("-i", "--Iteration", help="Iteration.",
                    metavar='\b', type=int)
                    
parser.add_argument("-dat", "--SaveDAT",  help="Save flow field in DAT format.",
                    action='store_true')
                    
parser.add_argument("-asc", "--SaveASC", help="Save flow field in ASCII format with an accompanying geometry file.",
                    action='store_true')
                    
parser.add_argument("-h5", "--SaveHDF5", help="Save flow field in HDF5 format.",
                    action='store_true')
                    
parser.add_argument("-ff", "--SaveFF", help="Save flow field in binary '.ff' format.",
                    action='store_true')
                    
args = parser.parse_args()
date = time.strftime("%Y_%m_%d")
#====================================================================
#### Construct an instance of the FlowFieldGeometry class        ####
#====================================================================
tmp = np.zeros((3, args.Nx, args.Ny, args.Nz),dtype=float)
ff = ffClass.FlowField(args.Nx,
                        args.Ny,
                        args.Nz,
                        args.Lx,
                        args.Lz,
                        args.Wavepacket,
                        args.Wavespeed,
                        args.Theta,
                        args.Reynolds,
                        args.Baseflow,
                        tmp,
                        "pp")
#====================================================================
#### Generate flow field using resolvent formulation             ####
#====================================================================
velocity_field = cr.resolvent_formulation(ff)
ff.set_ff(velocity_field, "pp")
#====================================================================
#### ---- Unstack velocity field in the wall-normal direction    ####
#====================================================================
ff.unstack_ff(ff.velocityField)
#====================================================================
#### ---- Add wall boundaries                                    ####
#====================================================================
ff.add_wall_boundaries()
#====================================================================
#### Make the correct directory to store the flow field in       ####
#====================================================================
if args.Iteration:
    output_directory = ut.make_FlowField_output_directory_wIteration(args.Directory, ff, date, args.Iteration)
else:
    output_directory = ut.make_FlowField_output_directory(args.Directory, ff, date)
#====================================================================
#### Save flow field to output_directory                         ####
#====================================================================
if args.SaveASC:
    ut.write_ASC(ff, output_directory, "u0")
    ut.write_GEOM(ff, output_directory, "u0")
    if args.SaveFF:
        ut.ascii2field(output_directory, "u0", "ff")
    if args.SaveHDF5:
        ut.ascii2field(output_directory, "u0", "h5")
if args.SaveDAT:
    ut.write_DAT(ff, output_directory, "u0")
# Save details of the flow field in a file called u0_details.txt
ut.write_Details(ff, output_directory, "u0")
print("\nThe flow field has been stored in")
print(output_directory)
#====================================================================
# End of file                                                    ####
#====================================================================
ut.print_EndMessage()