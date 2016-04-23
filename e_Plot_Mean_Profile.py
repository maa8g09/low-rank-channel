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
#import pylab 
from matplotlib import pyplot as plt
#from matplotlib import animation
from matplotlib import cm as cm

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
                    help="File to approximate",
                    required=True)
parser.add_argument("-i",
                    "--VelComponent",
                    metavar='\b',
                    help="Velocity component to plot.",
                    required=True,
                    type=int)
parser.add_argument("-re",
                    "--Reynolds",
                    metavar='\b',
                    help="Reynolds number of flow field.",
                    required=True,
                    type=float)
parser.add_argument("-d", "--delta", help="plot delta",
                    action="store_true")
args = parser.parse_args()

####################################################################################################
# Create a temporary folder in which to do all the magic in
parent_directory = os.getcwd()

if parent_directory[-1] != "/":
    parent_directory += "/"

temp_plot_folder = "tmp_plots/"
temp_plot_folder = parent_directory + temp_plot_folder

#if a temporary directory exists, delete it:
if os.path.exists(temp_plot_folder):
    command = "rm -rf " + temp_plot_folder
    os.system(command)

#if a temporary directory doesn't exist, create one:
if not os.path.exists(temp_plot_folder):
    os.mkdir(temp_plot_folder)

# All work is done from the temporary directory
os.chdir(temp_plot_folder)

####################################################################################################
# Convert the binary flow field parsed in, to ascii
command = "field2ascii -p ../" + str(args.File) + " " + str(args.File)[:-3]
os.system(command)

####################################################################################################
# Read the physical ascii and grab a slice to plot. Boom!
var = ut.read_ASC_PP(temp_plot_folder, str(args.File)[:-3])
ffcf = ffClass.FlowFieldChannelFlow(var['Nd'],
                                    var['Nx'],
                                    var['Ny'],
                                    var['Nz'],
                                    var['Lx'],
                                    var['Lz'],
                                    var['alpha'],
                                    var['beta'],
                                    0.0,
                                    "lam",
                                    0.0,
                                    var['ff'],
                                    "pp")

####################################################################################################
# Calculate the mean profile
#### CASE
cumulative = np.zeros((ffcf.Nx, ffcf.Ny))
vel_stable = np.zeros((ffcf.Ny))

z_avgd_mean = cumulative


# Average in spanwise direction
for nz in range(0, ffcf.Nz):
    cumulative[:, :] += ffcf.ff[args.VelComponent, :, :, nz]

z_avgd_mean[:, :] = cumulative[:, :] * (1.0/ffcf.Nz)
    
cumulative = vel_stable
# Average in streamwise direction
for nx in range(0, ffcf.Nx):
    cumulative[:] += z_avgd_mean[nx, :]
    
vel_stable[:] = cumulative[:] * (1.0/ffcf.Nx)


#### LAMINAR
x = np.linspace(0.0, ffcf.Lx, ffcf.Nx)
y = np.linspace(-1.0, 1.0, ffcf.Ny)

if args.VelComponent == 0:
    
u_lam = np.linspace(-1.0, 1.0, ffcf.Ny)
for ny in range(0, ffcf.Ny):                             # laminar profile
    y[ny] = np.cos(ny*np.pi/(ffcf.Ny-1))
    u_lam[ny] = 1.0 - y[ny]**2.0


#### Plot 1
u_lam = u_lam[ffcf.Ny/2.0:-1].T
vel = vel_stable[ffcf.Ny/2.0:-1].T
half_channel = y[ffcf.Ny/2.0:-1]

plt.figure()
title="Re " + str(args.Reynolds)
plt.plot(vel, half_channel, 'r-', u_lam, half_channel, 'g-')
plt.xlabel('u')
plt.ylabel('y')
plt.ylim([-1.0, 0.0])
plt.grid(True)
legnd = [str(args.File)[:-3], 'Laminar']
plt.legend(legnd, loc='best') 
plt.title(title)
plt.savefig(temp_plot_folder + 'mean_profile_'+str(args.File)[:-3]+'.png')


#### Plot 2
delta = vel - u_lam
if args.delta:
    plt.figure()
    title="Re " + str(args.Reynolds)
    plt.plot(delta, half_channel, 'r-')
    plt.xlabel('u')
    plt.ylabel('y')
    plt.ylim([-1.0, 0.0])
    plt.grid(True)
    legnd = ["(unewt - lam)"]
    plt.legend(legnd, loc='best') 
    plt.title(title)
    plt.savefig(temp_plot_folder + 'mean_profile_delta_'+str(args.File)[:-3]+'.png')


