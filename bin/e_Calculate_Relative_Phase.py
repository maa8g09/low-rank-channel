#!/usr/bin/env python3
"""
.. - - - - - - 
"""
import argparse
from numpy import angle
from numpy import exp
from numpy import mean
from numpy.linalg import norm
import sys
sys.path.append("/".join(sys.path[0].split("/")[:-1]))
import FlowField as ffClass
import Tests
import Utils as ut
parser = argparse.ArgumentParser(description="Calculate difference between two fields and output field in HDF5 format")
parser.add_argument("-f1", "--File1", help="File 1",
                    metavar='\b', required=True)
parser.add_argument("-f2", "--File2", help="File 2",
                    metavar='\b', required=True)
args = parser.parse_args()
#===================================================================#
#### Read the HDF5 and details file                              ####
#===================================================================#
file1_info, original_attrs1 = ut.read_H5(args.File1)
file2_info, original_attrs2 = ut.read_H5(args.File2)
#===================================================================#
#### Declare flow field objects                                  ####
#===================================================================#
ff1 = ffClass.FlowFieldChannelFlow( file1_info['Nd'],
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
                                    file1_info['ff'],
                                    "pp")
ff2 = ffClass.FlowFieldChannelFlow( file2_info['Nd'],
                                    file2_info['Nx'],
                                    file2_info['Ny'],
                                    file2_info['Nz'],
                                    file2_info['Lx'],
                                    file2_info['Lz'],
                                    file2_info['alpha'],
                                    file2_info['beta'],
                                    0.0,
                                    "lam",
                                    0.0,
                                    file2_info['ff'],
                                    "pp")
#===================================================================#
#### ---- Remove the wall boundaries                             ####
#===================================================================#
# Removing the xz-planes at y=1 and y=-1,
# so that the chebyshev nodes can be used to construct 
# the transfer function.
ff1.remove_wall_boundaries()
ff2.remove_wall_boundaries()
#===================================================================#
#### ---- Fourier transform in xz directions                     ####
#===================================================================#
ff1.make_xz_spectral()
ff2.make_xz_spectral()
#===================================================================#
#### ---- Stack velocity fields in the wall-normal direction     ####
#===================================================================#
ff1.stack_ff_in_y()
ff2.stack_ff_in_y()
#===================================================================#
#### Calculate amplitude/phase difference                        ####
#===================================================================#
tolerance = 1e-9
print("==========================================")
print("Differences between " + args.File1[:-3] + " and " + args.File2[:-3])
print("==========================================")
message = "Amplitude difference between " + args.File1[:-3] + " and " + args.File2[:-3]
diff_A = Tests.difference(abs(ff1.velocityField), abs(ff2.velocityField), tolerance, message)
diff_P = mean(angle(ff1.velocityField,deg=True) - angle(ff2.velocityField,deg=True))
print("Amplitude diff \t\t%.2E" % diff_A)
print("")
print("Field 1 (deg) \t\t%.6f (avg)" % mean(angle(ff1.velocityField,deg=True)))
print("Field 2 (deg) \t\t%.6f (avg)" % mean(angle(ff2.velocityField,deg=True)))
print("Phase diff (deg) \t%.6f (avg)" % diff_P)
print("")
#===================================================================#
#### Euler form of Fourier transformed field                     ####
#===================================================================#
print("==========================================")
print("Reconstructed using Euler formula")
print("==========================================")
test_ff1 = abs(ff1.velocityField) * exp(1j * angle(ff1.velocityField))
test_ff2 = abs(ff2.velocityField) * exp(1j * angle(ff2.velocityField))
message = "Difference between spectral " + args.File1[:-3] + " and " + args.File1[:-3] + "_test"
diff_ff1 = Tests.difference(ff1.velocityField, test_ff1, tolerance, message)
message = "Difference between spectral " + args.File2[:-3] + " and " + args.File2[:-3] + "_test"
diff_ff2 = Tests.difference(ff2.velocityField, test_ff2, tolerance, message)
print("Flow field 1 diff \t%.2E" % diff_ff1)
print("Flow field 2 diff \t%.2E" % diff_ff2)
print("")
#===================================================================#
#### Phase shift ff2 -> ff1                                      ####
#===================================================================#
print("==========================================")
print("Differences between " + args.File1[:-3] + " and " + args.File2[:-3] + "_shifted")
print("==========================================")
ff2_shifted = abs(ff2.velocityField) * exp(1j * angle(ff1.velocityField))
message = "Amplitude difference between " + args.File1[:-3] + " and " + args.File2[:-3] + "_shifted"
diff_A_shifted = Tests.difference(abs(ff1.velocityField), abs(ff2_shifted), tolerance, message)
diff_P_shifted = mean(angle(ff1.velocityField,deg=True) - angle(ff2_shifted,deg=True))   # should be the same(ish)
print("Amplitude diff \t\t%.2E" % diff_A_shifted)
print("")
print("Field 1 (deg) \t\t%.6f (avg)" % mean(angle(ff1.velocityField,deg=True)))
print("Field 2 (deg) \t\t%.6f (avg)" % mean(angle(ff2_shifted,deg=True)))
print("Phase diff (deg) \t%.6f (avg)" % diff_P_shifted)
print("")
#===================================================================#
#### Initialize shifted flow field object                        ####
#===================================================================#
ff2_shifted = ffClass.FlowFieldChannelFlow( file2_info['Nd'],
                                            file2_info['Nx'],
                                            file2_info['Ny']-2,
                                            file2_info['Nz'],
                                            file2_info['Lx'],
                                            file2_info['Lz'],
                                            file2_info['alpha'],
                                            file2_info['beta'],
                                            0.0,
                                            "lam",
                                            0.0,
                                            ff2_shifted,
                                            "sp")
#===================================================================#
#### ---- Unstack velocity field in the wall-normal direction    ####
#===================================================================#
ff2_shifted.unstack_ff()
#===================================================================#
#### ---- Inverse Fourier transform approximated field in xz directions
#===================================================================#
ff2_shifted.make_xz_physical()
#===================================================================#
#### ---- Add wall boundaries                                    ####
#===================================================================#
ff2_shifted.add_wall_boundaries()
#===================================================================#
#### Write shifted HDF5 file                                     ####
#===================================================================#
fileName = args.File2[:-3] + "_shifted_to_" + args.File1[:-3]
ut.write_H5(ff2_shifted, original_attrs2, fileName)
ut.print_EndMessage()