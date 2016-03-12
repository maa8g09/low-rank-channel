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

# Project flow field. 
# this should work by taking a flow field and a rank and retrning a rank-approximated flowfield 
# for you to use. 


date = time.strftime("%Y_%m_%d")

####################################################################################################
# Parse the command line arguments (flag parameters)
#ut.print_ResolventHeader()
#ut.print_ResolventSubHeader()
parser = argparse.ArgumentParser(description="Plot flow field.")
parser.add_argument("-f",
                    "--File",
                    metavar='\b',
                    help="File to plot",
                    required=True)
parser.add_argument("-i",
                    "--VelComponent",
                    metavar='\b',
                    help="Velocity component to plot (integers i.e. 0 = u, ...)",
                    required=True,
                    type=int)
parser.add_argument("-n",
                    "--SpatialComponent",
                    metavar='\b',
                    help="Spatial direction to plot in (integers i.e. 0 = x, ...),\nE.g. selecting 0 will plot yz planes.",
                    required=True,
                    type=int)
parser.add_argument("-d",
                    "--Directory",
                    metavar='\b',
                    help="Directory where u0_Details.txt is kept.",
                    required=True)
args = parser.parse_args()

# The directory should be the data-X/ folder.
output_directory = os.getcwd()
if output_directory[-1] != '/':
    output_directory += '/'

#### Make output directory:
# Make a direcotry to store all the information
images_directory = output_directory + "images_" + str(args.File)[:-3] + "/"
x_directory = images_directory + "x/"
y_directory = images_directory + "y/"
z_directory = images_directory + "z/"

# Check to see if the images directory exists.
if os.path.exists(images_directory):
    print("\nThis directory already exists:\t" + str(images_directory))
#    ut.error("Images directory for this flow field already exists, check it before running this again!")

# if directory doesn't exist, make it:
if not os.path.exists(images_directory):
    print("Making:\t\t" + str(images_directory))
    os.mkdir(images_directory)


#### Change down
# Change into the images directory
os.chdir(images_directory)


#### Convert to ascii
# Convert the binary flow field parsed in, to ascii
command = "field2ascii -p ../" + str(args.File) + " " + str(args.File)[:-3]
os.system(command)


#### Read into class
var = ut.read_ASC_PP(images_directory, str(args.File)[:-3])

details_directory = args.Directory
if details_directory[-1] != "/":
    details_directory += "/"

var2= ut.read_Details(details_directory, str(args.File[:-3]))

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


i = args.VelComponent
n = args.SpatialComponent
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

n = args.SpatialComponent

#### Spatial component check 
if n == 0:
    
    if os.path.exists(x_directory):
        print("\nThis directory already exists:\t" + str(x_directory))
        command = "rm -rf " + x_directory
    if not os.path.exists(x_directory):
        print("Making:\t\t" + str(x_directory))
        os.mkdir(x_directory)


    print("\nPlotting some awesome contourplots in x direction")
    
    vl_max = np.amax(ff.velocityField[i, :, :, :])
    vl_min = np.amin(ff.velocityField[i, :, :, :])
    
    for j in range(0, ff.Nx):
        ut.plot_Contour(x_directory, str(args.File[:-3]), 
                        ff.z, ff.y, ff.velocityField[i, j, :, :], 
                        format(ff.x[j], '.2f'), ":", ":", 
                        ff.Re, 
                        str(j),  "-", "-",
                        ff.velocityField[m, j, :, :], 
                        ff.velocityField[p, j, :, :], 
                        "z", "y", velName,
                        vl_max, vl_min, True)
        print(j)


elif n == 1:
    
    if os.path.exists(y_directory):
        print("\nThis directory already exists:\t" + str(y_directory))
        command = "rm -rf " + y_directory
    if not os.path.exists(y_directory):
        print("Making:\t\t" + str(y_directory))
        os.mkdir(y_directory)


    print("\nPlotting some awesome contourplots in y direction")
    
    vl_max = np.amax(ff.velocityField[i, :, :, :])
    vl_min = np.amin(ff.velocityField[i, :, :, :])
    
    for j in range(0, ff.Ny):
        ut.plot_Contour(y_directory, str(args.File[:-3]), 
                        ff.z, ff.x, ff.velocityField[i, :, j, :], 
                        ":", format(ff.y[j], '.2f'), ":", 
                        ff.Re, 
                        "-", str(j), "-", 
                        ff.velocityField[m, :, j, :], 
                        ff.velocityField[p, :, j, :], 
                        "z", "x", velName,
                        vl_max, vl_min, True)
        print(j)

elif n == 2:
    
    if os.path.exists(z_directory):
        print("\nThis directory already exists:\t" + str(z_directory))
        command = "rm -rf " + z_directory
    if not os.path.exists(z_directory):
        print("Making:\t\t" + str(z_directory))
        os.mkdir(z_directory)


    print("\nPlotting some awesome contourplots in z direction")
    
    vl_max = np.amax(ff.velocityField[i, :, :, :].T)
    vl_min = np.amin(ff.velocityField[i, :, :, :].T)
    
    for j in range(0, ff.Nz):
        ut.plot_Contour(z_directory, str(args.File[:-3]), 
                        ff.x, ff.y, ff.velocityField[i, :, :, j].T, 
                        ":", ":", format(ff.z[j], '.2f'), 
                        ff.Re, 
                        "-", "-", str(j), 
                        ff.velocityField[m, :, :, j].T, 
                        ff.velocityField[p, :, :, j].T, 
                        "x", "y", velName,
                        vl_max, vl_min, True)
        print(j)



# Delete the asc and geom files.
os.chdir(images_directory)
command = "rm *.asc *.geom"
os.system(command)


ut.print_EndMessage()
