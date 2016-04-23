#!/usr/bin/env python

# This is an executable file that should be run from the command line.
# The command line arguments determine the type of flow field created.


## Need to add all the usual things about yourself. Name, Institute etc.

## This file should be run from the directory where the 

import argparse
import time
import os
import numpy as np

import Utils as ut
import FlowField as ffClass
import Tests

# Project flow field. 
# this should work by taking a flow field and a rank and retrning a rank-approximated flowfield 
# for you to use. 


date = time.strftime("%Y_%m_%d")

####################################################################################################
# Parse the command line arguments (flag parameters)
#ut.print_ResolventHeader()
#ut.print_ResolventSubHeader()
parser = argparse.ArgumentParser(description="Plot flow field slice at time given. \nThis should be executed from the directory where the flow fields are saved.")
parser.add_argument("-t",
                    "--Time",
                    metavar='\b',
                    help="File to plot at time t.",
                    required=True,
                    type=int)
parser.add_argument("-i",
                    "--VelComponent",
                    metavar='\b',
                    help="Velocity component to plot (integers i.e. 0 = u, ...)",
                    required=True,
                    type=int)
parser.add_argument("-n",
                    "--SpatialNorm",
                    metavar='\b',
                    help="Spatial direction to plot in (integers i.e. 0 = x, ...),\nE.g. selecting 0 will plot yz planes.",
                    required=True,
                    type=int)
parser.add_argument("-d",
                    "--Directory",
                    metavar='\b',
                    help="Directory where u0_Details.txt is kept.",
                    required=True)
parser.add_argument("-o",
                    "--OutputDirectory",
                    metavar='\b',
                    help="Directory where slices are output."
                    )
args = parser.parse_args()





# The directory should be the data-X/ folder.
output_directory = os.getcwd()
if output_directory[-1] != '/':
    output_directory += '/'



#================================================================
#### Ensure we are in the DNS directory
#================================================================
present_directory = output_directory.split("/")[-2]
if present_directory.find("data") != -1:
    print("Correct directory given")

else:
    message = "Execute from the DNS directory. Present directory is \n\n" + present_directory
    print(message)
    ut.error()


#================================================================
# Make a temporary folder:
#================================================================
tmp_directory = output_directory + "tmp/"
if os.path.exists(tmp_directory):
    print("\nRemoving:\t" + str(tmp_directory))
    command = "rm -rf " + tmp_directory
    os.system(command)
# if directory doesn't exist, make it:
if not os.path.exists(tmp_directory):
    print("Making:\t\t" + str(tmp_directory))
    os.mkdir(tmp_directory)


#================================================================
# Change into temp folder and convert ff to ascii.
#================================================================
os.chdir(tmp_directory)
fileName = "u" + str(args.Time) + ".000"
command = "field2ascii -p ../" + fileName + ".ff " + fileName
os.system(command)


#================================================================
# Read the ascii and construct a class to use for plotting the slice.
#================================================================
var = ut.read_ASC_PP(tmp_directory, fileName)
var2= ut.read_Details(output_directory[:output_directory.find("data")], "u0")
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


#================================================================
# This is where the image gets saved.
#================================================================
slice_directory = ""
if str(args.OutputDirectory):
    slice_directory = args.OutputDirectory

else:
    slice_directory = output_directory + "slices/"

if os.path.exists(slice_directory):
    print("\nThis directory already exists:\t" + str(slice_directory))
    print("\nThis slce will be saved there.")

if not os.path.exists(slice_directory):
    print("Making:\t\t" + str(slice_directory))
    os.mkdir(slice_directory)



i = args.VelComponent
n = args.SpatialNorm
m = n+1
p = n-1

if m==3:
    m=0
if p==-1:
    p=2

velName = ""
if i==0:
    velName = "u"
elif i==1:
    velName = "v"
elif i==2:
    velName = "w"


ffdata = np.zeros((ff.Nx, ff.Ny, ff.Nz))
ffdata = ff.velocityField[i, :, :, :]

n = args.SpatialNorm

#### Spatial component check 
fileName = "u" + str(args.Time).zfill(5)
if n == 0:

    print("\nPlotting contourplots in x direction at averaged in the streamwise direction.")

    # Average the flow field in teh streamwise direction.
    x_averaged_ff = np.zeros((ff.Nd, ff.Ny, ff.Nz))
    for nd in range(0, ff.Nd):
        for nx in range(0, ff.Nx):
            x_averaged_ff[nd, :, :] += ff.velocityField[nd, nx, :, :]
    
    x_averaged_ff *= (1.0 / ff.Nx)
    
    vl_max = np.amax(x_averaged_ff[i, :, :])
    vl_min = np.amin(x_averaged_ff[i, :, :])

    ut.plot_Contour(slice_directory, fileName, 
                    ff.z, ff.y, x_averaged_ff[i, :, :], 
                    "avg", ":", ":", 
                    ff.Re, 
                    "x-avg",  "-", "-",
                    x_averaged_ff[m, :, :], 
                    x_averaged_ff[p, :, :], 
                    "z", "y", velName,
                    vl_max, vl_min, True)


elif n == 1:

    print("\nPlotting contourplots in y direction at averaged in the wall-normal direction.")

    # Average the flow field in teh streamwise direction.
    y_averaged_ff = np.zeros((ff.Nd, ff.Nx, ff.Nz))
    for nd in range(0, ff.Nd):
        for ny in range(0, ff.Ny):
            y_averaged_ff[nd, :, :] += ff.velocityField[nd, :, ny, :]

    y_averaged_ff *= (1.0 / ff.Ny)

    vl_max = np.amax(y_averaged_ff[i, :, :])
    vl_min = np.amin(y_averaged_ff[i, :, :])
    
    ut.plot_Contour(slice_directory, fileName, 
                    ff.z, ff.x, y_averaged_ff[i, :, :], 
                    ":", "avg", ":", 
                    ff.Re, 
                    "-", "y-avg", "-", 
                    y_averaged_ff[m, :, :], 
                    y_averaged_ff[p, :, :], 
                    "z", "x", velName,
                    vl_max, vl_min, True)

elif n == 2:

    print("\nPlotting contourplots in z direction at averaged in the spanwise direction.")

    # Average the flow field in teh streamwise direction.
    z_averaged_ff = np.zeros((ff.Nd, ff.Nx, ff.Ny))
    for nd in range(0, ff.Nd):
        for nz in range(0, ff.Nz):
            z_averaged_ff[nd, :, :] += ff.velocityField[nd, :, :, nz]

    z_averaged_ff *= (1.0 / ff.Nz)

    vl_max = np.amax(z_averaged_ff[i, :, :])
    vl_min = np.amin(z_averaged_ff[i, :, :])
    
    ut.plot_Contour(slice_directory, fileName, 
                    ff.x, ff.y, z_averaged_ff[i, :, :].T, 
                    ":", ":", "avg", 
                    ff.Re, 
                    "-", "-", "z-avg", 
                    z_averaged_ff[m, :, :].T, 
                    z_averaged_ff[p, :, :].T, 
                    "x", "y", velName,
                    vl_max, vl_min, False)



# Remove the temporary directory
os.chdir(output_directory)
command = "rm -rf " + tmp_directory
os.system(command)

ut.print_EndMessage()
