#!/usr/bin/env python3

# This is an executable file that should be run from the command line.
# The command line arguments determine the type of flow field created.


## Need to add all the usual things about yourself. Name, Institute etc.

## This file should be run from the directory where the 

import argparse
import time
import os
import numpy as np
import sys
sys.path.append("/".join(sys.path[0].split("/")[:-1]))
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
parser.add_argument("-coord",
                    "--Coordinate",
                    metavar='\b',
                    help="Spatial co-ordinate of normal plane.",
                    required=True,
                    type=float)
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
if str(args.OutputDirectory) != "None":
    print(str(args.OutputDirectory))
    slice_directory = args.OutputDirectory

else:
    slice_directory = output_directory + "slices/"
    print("\n\n Hello \n\n")

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

    print("\nPlotting some awesome contourplots in x direction at co-ordinate: " + str(args.Coordinate))
    
    coords = Tests.indices(ff.x, lambda m: m > float(args.Coordinate))
    x_coord = coords[0]

    vl_max = np.amax(ff.velocityField[i, :, :, :])
    vl_min = np.amin(ff.velocityField[i, :, :, :])

    ut.plot_Contour(slice_directory, fileName, 
                    ff.z, ff.y, ff.velocityField[i, x_coord, :, :], 
                    format(ff.x[x_coord], '.2f'), ":", ":", 
                    ff.Re, 
                    str(x_coord),  "-", "-",
                    ff.velocityField[m, x_coord, :, :], 
                    ff.velocityField[p, x_coord, :, :], 
                    "z", "y", velName,
                    vl_max, vl_min, True)


elif n == 1:

    print("\nPlotting some awesome contourplots in y direction at co-ordinate: " + str(args.Coordinate))
    
    coords = Tests.indices(ff.y, lambda m: m > float(args.Coordinate))
    y_coord = coords[-1]

    vl_max = np.amax(ff.velocityField[i, :, :, :])
    vl_min = np.amin(ff.velocityField[i, :, :, :])
    
    ut.plot_Contour(slice_directory, fileName, 
                    ff.z, ff.x, ff.velocityField[i, :, y_coord, :], 
                    ":", format(ff.y[y_coord], '.2f'), ":", 
                    ff.Re, 
                    "-", str(y_coord), "-", 
                    ff.velocityField[m, :, y_coord, :], 
                    ff.velocityField[p, :, y_coord, :], 
                    "z", "x", velName,
                    vl_max, vl_min, True)

elif n == 2:

    print("\nPlotting some awesome contourplots in z direction at co-ordinate: " + str(args.Coordinate))
    
    coords = Tests.indices(ff.z, lambda m: m > float(args.Coordinate))
    z_coord = coords[0]

    vl_max = np.amax(ff.velocityField[i, :, :, :].T)
    vl_min = np.amin(ff.velocityField[i, :, :, :].T)

    ut.plot_Contour(slice_directory, fileName, 
                    ff.x, ff.y, ff.velocityField[i, :, :, z_coord].T, 
                    ":", ":", format(ff.z[z_coord], '.2f'), 
                    ff.Re, 
                    "-", "-", str(z_coord), 
                    ff.velocityField[m, :, :, z_coord].T, 
                    ff.velocityField[p, :, :, z_coord].T, 
                    "x", "y", velName,
                    vl_max, vl_min, False)



# Remove the temporary directory
os.chdir(output_directory)
command = "rm -rf " + tmp_directory
os.system(command)

ut.print_EndMessage()
