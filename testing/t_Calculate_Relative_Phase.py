#!/usr/bin/env python3
"""
=====================================================================
Add profile to field
=====================================================================
Add given velocity profile to the streamwise component of the 
given flow field.

#### To do:
####    - Interpolate velocity profile if len(profile) != field.Ny

Muhammad Arslan Ahmed
maa8g09@soton.ac.uk
Room 5069
Building 13
Aerodynamics and Flight Mechanics
University of Southampton

"""
import os
from numpy import angle
from numpy import exp
from numpy.linalg import norm
import sys
sys.path.append("/".join(sys.path[0].split("/")[:-1]))
import FlowField as ffClass
import Utils as ut

directory = "/home/arslan/Documents/work/cfd-channelflow_solutions/w03_EQ_nsb1/testing_diffs"
os.chdir(directory)
File1 = "eq15-1.h5"
File2 = "eq15-2.h5"
#===================================================================#
#### Read the HDF5 and details file                              ####
#===================================================================#
file1_info, original_attrs1 = ut.read_H5(File1)
file2_info, original_attrs2 = ut.read_H5(File2)
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
amplitude_ff1 = abs(ff1.velocityField)
amplitude_ff2 = abs(ff2.velocityField)
diff_A = norm(amplitude_ff1 - amplitude_ff2)
phase_ff1 = angle(ff1.velocityField)
phase_ff2 = angle(ff2.velocityField)
diff_P = norm(phase_ff1 - phase_ff2)
phase_difference = phase_ff1 - phase_ff2
#===================================================================#
#### Euler form of Fourier transformed field                     ####
#===================================================================#
test_ff1 = amplitude_ff1 * exp(1j * phase_ff1)
test_ff2 = amplitude_ff2 * exp(1j * phase_ff2)
diff_ff1 = norm(ff1.velocityField - test_ff1)
diff_ff2 = norm(ff2.velocityField - test_ff2)

test_ff1_phase = angle(test_ff1)
diff_P_test = norm(phase_ff1 - test_ff1_phase)
#===================================================================#
#### Phase shift ff2 -> ff1                                      ####
#===================================================================#
ff2_shifted = amplitude_ff2 * exp(1j * phase_ff1)
diff_A_shifted = norm(amplitude_ff1 - abs(ff2_shifted)) # should be the same
diff_P_shifted = norm(phase_ff1 - angle(ff2_shifted))   # should be the same
#===================================================================#
#### Initialize shifted flow field object                        ####
#===================================================================#
ff2_shifted = ffClass.FlowFieldChannelFlow( file2_info['Nd'],
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
fileName = File2[:-3] + "_shifted_to_" + File1[:-3]
ut.write_H5(ff2_shifted, original_attrs2, fileName)
ut.print_EndMessage()