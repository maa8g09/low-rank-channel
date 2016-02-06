#!/usr/bin/env python
import argparse
import os
import numpy as np
import Utils as ut
import FlowField as ffClass

parser = argparse.ArgumentParser(description="Calculate mean velocity field of converged data.")

parser.add_argument("-T0",
                    "--t_start",
                    metavar='\b',
                    help="Plotting start time unit.",
                    required=True,
                    type=int)
parser.add_argument("-T1",
                    "--t_end",
                    metavar='\b',
                    help="Plotting end time unit.",
                    required=True,
                    type=int)
parser.add_argument("-d",
                    "--Directory",
                    metavar='\b',
                    help="Directory where the DNS data is stored. (Give full directory)",
                    required=True)

ut.print_Start_Bar()
args = parser.parse_args()

# The directory should be the data-X/ folder.
dns_data_directory = args.Directory
if dns_data_directory[-1] != '/':
    dns_data_directory += '/'


#### Check directory
# Now check to see if the directory has data after the penultimate '/'
present_directory = dns_data_directory.split("/")[-2]
if present_directory.find("data") != -1:
    print("Correct directory given")

else:
    message = "Invalid directory flag. Please add data-X/ folder at the end. Present directory: " + present_directory
    ut.error()


#### Make Directory
# Make a directory to store the files with integer names.
tmp_directory = dns_data_directory + 'tmp/'
# if directory exists, delete it:
if os.path.exists(tmp_directory):
    print("\nRemoving:\t" + str(tmp_directory))
    command = "rm -rf " + tmp_directory
    os.system(command)

# if directory doesn't exist, make it:
if not os.path.exists(tmp_directory):
    print("Making:\t\t" + str(tmp_directory))
    os.mkdir(tmp_directory)


# construct a flowField class for the mean (read the u0_Details file)
var = ut.calculate_Mean(dns_data_directory, tmp_directory, args.t_start, args.t_end)
var2= ut.read_Details(dns_data_directory[:dns_data_directory.find("data")], "u0")
ffg = ffClass.FlowFieldGeometry(var2['bf'],
                                var2['wp'],
                                var2['Nd'],
                                var2['Nx'],
                                var2['Ny'],
                                var2['Nz'],
                                var2['Re'],
                                var2['c'],
                                var2['theta'])

ff = ffClass.FlowField(ffg, var['ff'], "pp")


# write mean file to disc as asc in the following folder:
# ..../data-X/mean/uMean.asc
# if directory doesn't exist, make it:
mean_directory = dns_data_directory + "mean/"
if not os.path.exists(mean_directory):
    print("Making:\t\t" + str(mean_directory))
    os.mkdir(mean_directory)

fileName = "uMean"
ut.write_ASC(ff, mean_directory, fileName)
ut.write_GEOM(ff, mean_directory, fileName)
ut.write_Details(ff, mean_directory, fileName)
ut.write_FF(mean_directory, fileName)








#### Calculate velocity profile
# Calculate the mean profile
cumulative = np.zeros((ff.Nx, ff.Ny))
vel_profile = np.zeros((ff.Ny))

z_avgd_mean = cumulative


# Average in spanwise direction
for nz in range(0, ff.Nz):
    cumulative[:, :] += ff.ff[args.VelComponent, :, :, nz]

z_avgd_mean[:, :] = cumulative[:, :] * (1.0/ff.Nz)
    
cumulative = vel_profile
# Average in streamwise direction
for nx in range(0, ff.Nx):
    cumulative[:] += z_avgd_mean[nx, :]
    
vel_profile[:] = cumulative[:] * (1.0/ff.Nx)


# Write the mean profile in wall-normal direction
ut.write_Vel_Profile(vel_profile, mean_directory, "turb_mean")


# maybe plot contour plots


command = "rm -rf " + tmp_directory
os.system(command)

ut.print_EndMessage()