#!/usr/bin/env python3
import argparse
import sys
sys.path.append("/".join(sys.path[0].split("/")[:-1]))
import Utils as ut
import FlowField as ffClass
parser = argparse.ArgumentParser(description="Calculate difference between two fields and output field in HDF5 format")
parser.add_argument("-f1", "--File1", help="File 1",
                    metavar='\b', required=True)
parser.add_argument("-f2", "--File2", help="File 1",
                    metavar='\b', required=True)
args = parser.parse_args()
#===================================================================#
#### Read the HDF5 and details file                              ####
#===================================================================#
file1_info, original_attrs1 = ut.read_H5(args.File1)
file2_info, original_attrs2 = ut.read_H5(args.File2)
#===================================================================#
#### Calculate difference                                        ####
#===================================================================#
difference = ut.calculate_Difference(file1_info['ff'], file2_info['ff'])
#===================================================================#
#### Declare flow field object                                   ####
#===================================================================#
ff_diff = ffClass.FlowFieldChannelFlow( file1_info['Nd'],
                                        file1_info['Nx'],
                                        file1_info['Ny'],
                                        file1_info['Nz'],
                                        file1_info['Lx'],
                                        file1_info['Lz'],
                                        file1_info['alpha'],
                                        file1_info['beta'],
                                        0.0,
                                        "lam",
                                        0.0,
                                        difference,
                                        "pp")
fileName = args.File1[:-3] + "_-_" + args.File2[:-3]
ut.write_H5(ff_diff, original_attrs1, fileName)
ut.print_EndMessage()