#!/usr/bin/env python

import argparse
import os
import numpy as np

import Utils as ut

from images2gif import writeGif
from PIL import Image
from time import time 


####################################################################################################
# Parse the command line arguments (flag parameters)
#ut.print_ResolventHeader()
#ut.print_ResolventSubHeader()
parser = argparse.ArgumentParser(description="Plot flow field slice at time given. \nThis should be executed from the directory where the flow fields are saved.")
parser.add_argument("-T0",
                    "--T0",
                    metavar='\b',
                    help="Start movie at T0.",
                    required=True,
                    type=int)
parser.add_argument("-T1",
                    "--T1",
                    metavar='\b',
                    help="End movie at T1.",
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
parser.add_argument("-a",
                    "--SpatiallyAveraged",
                    help="Spatially averaged?.",
                    action='store_true')

args = parser.parse_args()

time_range = int(args.T1) - int(args.T0) + 1
movie_time = np.linspace(int(args.T0), int(args.T1), time_range)

# The directory should be the data-X/ folder.
pwd = os.getcwd()
if pwd[-1] != '/':
    pwd += '/'


pwd2 = pwd[:pwd.find("data")]
if pwd2[-1] != '/':
    pwd2 += '/'

normal = ""
if args.SpatialNorm == 0:
    normal = "x"
elif args.SpatialNorm == 1:
    normal = "y"
elif args.SpatialNorm == 2:
    normal = "z"
vel = ""
if args.VelComponent == 0:
    vel = "u"
elif args.VelComponent == 1:
    vel = "v"
elif args.VelComponent == 2:
    vel = "w"

movie_directory = pwd2 + "movie_" + str(args.T0) + "-" + str(args.T1) + "_" + normal + str(args.Coordinate) + "_" + vel + "/"

for t in range(0, time_range):
    time_unit = int(movie_time[t])
    # If spatially averaged?
    # change the command.
    command = ""
    if args.SpatiallyAveraged:
        command+= "ePlotSpatialAverage.py"
    else:
        command+= "ePlotSlice.py"
        command+= " -coord " + str(args.Coordinate)

    command+= " -d " + str(args.Directory)
    command+= " -o " + movie_directory
    command+= " -n " + str(args.SpatialNorm)
    command+= " -i " + str(args.VelComponent)
    command+= " -t " + str(time_unit)
    print("")
    print(command)
    print("")
    os.system(command)

if args.SpatiallyAveraged:
    command = "mv " + movie_directory + " " + pwd2 + "movie_" + str(args.T0) + "-" + str(args.T1) + "_" + normal + str(args.Coordinate) + "_" + vel + "_avgd/"
    os.system(command)

os.chdir(movie_directory)
#slices = sorted((fn for fn in os.listdir(movie_directory) if fn.endswith('.png')))
#
#images = [Image.open(str(movie_directory) + fn) for fn in slices]
#
#runningtime = 10.0
#runningtime = 5.0
#fileName = movie_directory + "movie_" + str(args.T0) + "-" + str(args.T1) + ".gif"
#writeGif(fileName, images, duration=runningtime, dither=1, nq = 1)
fileName = movie_directory + "movie_" + str(args.T0) + "-" + str(args.T1) + ".gif"
command = "convert -delay 50 *.png " + fileName
os.system(command)

