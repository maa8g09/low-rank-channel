#!/usr/bin/env python
import argparse
import os
import numpy as np
import Utils as ut
import FlowField as ffClass

parser = argparse.ArgumentParser(description="Calculate mean velocity field of converged data.\n\nThis function should be executed from the directory in which the DNS data is stored, i.e. .../data-X ")

parser.add_argument("-T0",
                    "--T0",
                    metavar='\b',
                    help="Calculate mean from time unit.",
                    required=True,
                    type=int)
parser.add_argument("-T1",
                    "--T1",
                    metavar='\b',
                    help="Calculate mean till time unit.",
                    required=True,
                    type=int)

ut.print_Start_Bar()
args = parser.parse_args()

#================================================================
# Add a '/' to the end of the directory string.
#================================================================
dns_data_directory = os.getcwd()
if dns_data_directory[-1] != '/':
    dns_data_directory += '/'



#================================================================
#### Ensure we are in the DNS directory
#================================================================
present_directory = dns_data_directory.split("/")[-2]
if present_directory.find("data") != -1:
    print("Correct directory given")

else:
    message = "Execute from the DNS directory. Present directory is \n\n" + present_directory
    print(message)
    ut.error()


#================================================================
#### Make a directory to store the files with integer names.
#================================================================
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


#================================================================
#### Construct a FlowField object for the temporal mean (read the u0_Details file)
#================================================================
var = ut.calculate_Temporal_Mean(dns_data_directory, tmp_directory, args.T0, args.T1)
var2= ut.read_Details(dns_data_directory[:dns_data_directory.find("data")], "u0_Details.txt")
ffg = ffClass.FlowFieldGeometry(var2['bf'],
                                var2['wp'],
                                var2['Nd'],
                                var2['Nx'],
                                var2['Ny'],
                                var2['Nz'],
                                var2['Re'],
                                var2['c'],
                                var2['theta'])
                                
ff_t = ffClass.FlowField(ffg, var['ff'], "pp")

#================================================================
#### Calculate velocity profiles U, u, v, w
#================================================================
vel_profiles = ut.calculate_Vel_Profiles(ff_t)

spatial_mean = np.zeros((ffg.Nd, ffg.Nx, ffg.Ny, ffg.Nz))

tmp = np.zeros((ffg.Ny))

for nx in range(0, ffg.Nx):
    for nz in range(0, ffg.Nz):
        spatial_mean[0, nx, :, nz] = vel_profiles['u']
        spatial_mean[1, nx, :, nz] = vel_profiles['v']
        spatial_mean[2, nx, :, nz] = vel_profiles['w']

ff_ts = ffClass.FlowField(ffg, spatial_mean, "pp")

# Add the laminar base flow to the velocity profiles to get a turbulent mean. Otherwise you are outputting deviation.

#================================================================
#### Write mean file to disc as ASC 
# in the following folder:
# ..../data-X/mean/uMean.asc
#================================================================
# if directory doesn't exist, make it:
mean_directory = dns_data_directory + "mean/"
if not os.path.exists(mean_directory):
    print("Making:\t\t" + str(mean_directory))
    os.mkdir(mean_directory)

fileName = "uMean_" + str(args.T0) + "-" + str(args.T1)
ut.write_ASC(ff_ts, mean_directory, fileName)
ut.write_GEOM(ff_ts, mean_directory, fileName)
ut.write_Details(ff_ts, mean_directory, fileName)
ut.write_FF(mean_directory, fileName)


#================================================================
#### Write the mean profile in wall-normal direction
#================================================================
fileName = "turbulent_deviation" + str(args.T0) + "-" + str(args.T1)
ut.write_Vel_Profile(vel_profiles['u'], mean_directory, fileName)


#================================================================
#### Plot the velocity profiles
#================================================================
fileName = "uMean_" + str(args.T0) + "-" + str(args.T1) + "-u"
ut.plot_Vel_Profile(mean_directory, fileName, vel_profiles['u'], vel_profiles['y'], "u", "y")

fileName = "uMean_" + str(args.T0) + "-" + str(args.T1) + "-v"
ut.plot_Vel_Profile(mean_directory, fileName, vel_profiles['v'], vel_profiles['y'], "v", "y")

fileName = "uMean_" + str(args.T0) + "-" + str(args.T1) + "-w"
ut.plot_Vel_Profile(mean_directory, fileName, vel_profiles['w'], vel_profiles['y'], "w", "y")

#fileName = "uMean_" + str(args.T0) + "-" + str(args.T1) + "-U"
#ut.plot_Vel_Profile(mean_directory, fileName, vel_profiles['U'], vel_profiles['y'], "U", "y")


#================================================================
#### Remove the temporary directory
#================================================================
command = "rm -rf " + tmp_directory
os.system(command)


ut.print_EndMessage()
