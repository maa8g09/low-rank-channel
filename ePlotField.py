#!/usr/bin/env python3
import argparse
import time
import os
import numpy as np
import Utils as ut
import FlowField as ffClass
date = time.strftime("%Y_%m_%d")
parser = argparse.ArgumentParser(description="Plot flow field.")
parser.add_argument("-f",
                    "--File",
                    metavar='\b',
                    help="File to plot",
                    required=True)
parser.add_argument("-d",
                    "--Details",
                    metavar='\b',
                    help="Details file.",
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
parser.add_argument("-q",
                    "--Quiver",
                    help="Plot vector arrows as well.",
                    action='store_true')
args = parser.parse_args()
pwd = ut.format_Directory_Path(os.getcwd())
file_info, original_attrs = ut.read_H5(args.File)
details = ut.read_Details(args.Details)
images_directory = ut.make_Folder(pwd, "images_" + str(args.File)[:-3], False)
os.chdir(images_directory)
ff = ffClass.FlowFieldChannelFlow( file_info['Nd'],
                                            file_info['Nx'],
                                            file_info['Ny'],
                                            file_info['Nz'],
                                            file_info['Lx'],
                                            file_info['Lz'],
                                            file_info['alpha'],
                                            file_info['beta'],
                                            details['c'],
                                            details['bf'],
                                            details['Re'],
                                            file_info['ff'],
                                            "pp")
i = args.VelComponent
n = args.SpatialComponent
m=n+1; p=n-1
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
print("Plotting...")
vl_max = np.amax(ff.velocityField[i, :, :, :])
vl_min = np.amin(ff.velocityField[i, :, :, :])
if n == 0:
    x_directory = ut.make_Folder(images_directory, "x", False)
    for j in range(0, ff.Nx):
        ut.plot_Contour(x_directory, str(args.File[:-3]), 
                        ff.z, ff.y, ff.velocityField[i, j, :, :], 
                        format(ff.x[j], '.2f'), ":", ":", 
                        ff.Re, 
                        str(j),  "-", "-",
                        ff.velocityField[m, j, :, :], 
                        ff.velocityField[p, j, :, :], 
                        "z", "y", velName,
                        vl_max, vl_min, args.Quiver)


elif n == 1:
    y_directory = ut.make_Folder(images_directory, "y", False)
    for j in range(0, ff.Ny):
        ut.plot_Contour(y_directory, str(args.File[:-3]), 
                        ff.z, ff.x, ff.velocityField[i, :, j, :], 
                        ":", format(ff.y[j], '.2f'), ":", 
                        ff.Re, 
                        "-", str(j), "-", 
                        ff.velocityField[m, :, j, :], 
                        ff.velocityField[p, :, j, :], 
                        "z", "x", velName,
                        vl_max, vl_min, args.Quiver)
elif n == 2:
    z_directory = ut.make_Folder(images_directory, "z", False)
    for j in range(0, ff.Nz):
        ut.plot_Contour(z_directory, str(args.File[:-3]), 
                        ff.x, ff.y, ff.velocityField[i, :, :, j].T, 
                        ":", ":", format(ff.z[j], '.2f'), 
                        ff.Re, 
                        "-", "-", str(j), 
                        ff.velocityField[m, :, :, j].T, 
                        ff.velocityField[p, :, :, j].T, 
                        "x", "y", velName,
                        vl_max, vl_min, args.Quiver)
ut.print_EndMessage()
