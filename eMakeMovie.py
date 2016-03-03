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
parser = argparse.ArgumentParser(description="Make Movie using slices specified, between time range specified.\nThis program should be run from the directory where the flow fields are saved.")


parser.add_argument("-i",
                    "--VelComponent",
                    metavar='\b',
                    help="Velocity component to plot (integers i.e. 0 = u, ...)",
                    required=True,
                    type=int)
parser.add_argument("-n",
                    "--SpatialComponent",
                    metavar='\b',
                    help="Spatial direction of normal (integers i.e. 0 = x, ...),\nE.g. selecting 0 will plot yz planes.",
                    required=True,
                    type=int)
parser.add_argument("-d",
                    "--Directory",
                    metavar='\b',
                    help="Directory where u0_Details.txt is kept.",
                    required=True)
parser.add_argument("-T0",
                    "--T0",
                    metavar='\b',
                    help="Plotting start time unit.",
                    required=True,
                    type=int)
parser.add_argument("-T1",
                    "--T1",
                    metavar='\b',
                    help="Plotting end time unit.",
                    required=False,
                    type=int)

args = parser.parse_args()


output_directory = os.getcwd()
if output_directory[-1] != '/':
    output_directory += '/'


T0 = str(args.T0)
if args.T1:
    T1 = ""
else:
    T1 = str(args.T1)
    
movie_directory = output_directory + "movie_" + T0 + "-" + T1 + "/"


# List all files in this directory,
# Then we make sure that they are within the range specified.
# If yes, then plot and save
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

