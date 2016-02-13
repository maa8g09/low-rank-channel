#!/usr/bin/env python


import argparse
import time
import os

import Utils as ut
import FlowField as ffClass
import ChannelResolvent as cr


parser = argparse.ArgumentParser(description="Project modes onto a flow field in order to yield a rank approximation of the velocity field.")
parser.add_argument("-f",
                    "--File",
                    metavar='\b',
                    help="File to convert",
                    required=True)
parser.add_argument("-type",
                    "--Type",
#                    metavar='\b',
                    help="Convert to type: ",
                    required=True,
                    choices=['dat', 'asc', 'h5'])
#                    nargs=1
args = parser.parse_args()




# Make a temporary directory:
pwd = os.getcwd()

if pwd[-1] != "/":
    pwd += "/"

temp_folder = "rank-temp/"
temp_folder = pwd + temp_folder

#if a temporary directory exists, delete it:
if os.path.exists(temp_folder):
    command = "rm -rf " + temp_folder
    os.system(command)

#if a temporary directory doesn't exist, create one:
if not os.path.exists(temp_folder):
    os.mkdir(temp_folder)

# All work is done from the temporary directory
os.chdir(temp_folder)


####################################################################################################
# Convert the binary flow field parsed in, to ascii
command = "field2ascii -p ../" + str(args.File) + " " + str(args.File)[:-3]
os.system(command)


####################################################################################################
# Read and store the spectral file.
var = ut.read_ASC_PP(temp_folder, str(args.File)[:-3])
var2 = ut.read_Details(pwd, str(args.File)[:-3])

command = "rm -rf " + temp_folder
os.system(command)

ffcf = ffClass.FlowFieldChannelFlow(var['Nd'],
                                    var['Nx'],
                                    var['Ny'],
                                    var['Nz'],
                                    var['Lx'],
                                    var['Lz'],
                                    var['alpha'],
                                    var['beta'],
                                    var2['c'],
                                    var2['bf'],
                                    var2['Re'],
                                    var['ff'],
                                    "pp")

if args.Type == "dat":
    ut.write_DAT(ffcf, pwd, str(args.File)[:-3])
