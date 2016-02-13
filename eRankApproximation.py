#!/usr/bin/env python

# This is an executable file that should be run from the command line.
# The command line arguments determine the type of flow field created.


## Need to add all the usual things about yourself. Name, Institute etc.

## This file should be run from the directory where the 

import argparse
import time
import os

import Utils as ut
import FlowField as ffClass
import ChannelResolvent as cr


# Project flow field. 
# this should work by taking a flow field and a rank and retrning a rank-approximated flowfield 
# for you to use. 

date = time.strftime("%Y_%m_%d")

####################################################################################################
# Parse the command line arguments (flag parameters)
#ut.print_ResolventHeader()
#ut.print_ResolventSubHeader()
parser = argparse.ArgumentParser(description="Project modes onto a flow field in order to yield a rank approximation of the velocity field.")
parser.add_argument("-f",
                    "--File",
                    metavar='\b',
                    help="File to approximate",
                    required=True)
parser.add_argument("-r",
                    "--Rank",
                    metavar='\b',
                    help="Number of velocity modes to use to recreate flow field.",
                    required=True,
                    type=int)
parser.add_argument("-d",
                    "--Directory",
                    metavar='\b',
                    help="Directory where u0_Details.txt is kept.",
                    required=True)
parser.add_argument("-v",
                    "--TurbMeanProfile",
                    metavar='\b',
                    help="Turbulent mean velocity profile. (Prefix with full directory)",
                    required=True)
parser.add_argument("-m",
                    "--MeanFile",
                    metavar='\b',
                    help="(S,P) Mean ascii file. (Prefix with full directory)",
                    required=True)

args = parser.parse_args()


####################################################################################################
# Create a temporary folder in which to do all the magic in
parent_directory = os.getcwd()

if parent_directory[-1] != "/":
    parent_directory += "/"

temp_rank_folder = "rank-temp/"
temp_rank_folder = parent_directory + temp_rank_folder

#if a temporary directory exists, delete it:
if os.path.exists(temp_rank_folder):
    command = "rm -rf " + temp_rank_folder
    os.system(command)

#if a temporary directory doesn't exist, create one:
if not os.path.exists(temp_rank_folder):
    os.mkdir(temp_rank_folder)

# All work is done from the temporary directory
os.chdir(temp_rank_folder)


####################################################################################################
# Convert the binary flow field parsed in, to ascii
command = "field2ascii -p ../" + str(args.File) + " " + str(args.File)[:-3]
os.system(command)


####################################################################################################
# Read and store the spectral file.
var = ut.read_ASC_SP(temp_rank_folder, str(args.File)[:-3])

details_directory = args.Directory
if details_directory[-1] != "/":
    details_directory += "/"

var2 = ut.read_Details(details_directory, str(args.File)[:-3])

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
                                    "sp")


# Read velocity profile
turb_mean = ut.read_Vel_Profile(args.TurbMeanProfile)

# Read mean asc file
tmp = args.MeanFile.find("uMean")
meanDir = args.MeanFile[:tmp]
var = ut.read_ASC_SP(meanDir, "uMean")
ffmean = ffClass.FlowFieldGeometry(var2['bf'],
                                   var2['wp'],
                                   var2['Nd'],
                                   var2['Nx'],
                                   var2['Ny'],
                                   var2['Nz'],
                                   var2['Re'],
                                   var2['c'],
                                   var2['theta'])

ffmean = ffClass.FlowField(ffmean, var['ff'], "sp")


####################################################################################################
# Constrct an approximation.
ffcf = cr.resolvent_approximation(ffcf, args.Rank, turb_mean, ffmean)


####################################################################################################
# Write new velocity field to disk as ascii.

os.chdir(parent_directory)      # Go up one directory
rank_folder = args.File[:-3]+"_rank_" + str(ffcf.rank) + "/"
rank_folder = parent_directory + rank_folder

#if a rank directory does exist, delete it:
if os.path.exists(rank_folder):
    command = "rm -rf " + rank_folder
    os.system(command)

#if a rank directory doesn't exist, create one:
if not os.path.exists(rank_folder):
    os.mkdir(rank_folder)

# Change into the new rank directory
os.chdir(rank_folder)


####################################################################################################
# Save the ascii file here
sp_ASC_fileName = ut.write_approximated_ASC(ffcf, rank_folder, ffcf.rank)


####################################################################################################
# Convert ascii to binary flowfield.
print(os.getcwd())
approximationFileName = args.File[:-3]+"_rank_"+str(ffcf.rank)
command = "ascii2field -p false -ge ../rank-temp/" + str(args.File)[:-3] + ".geom " + sp_ASC_fileName + " " + approximationFileName
print(command)
os.system(command)


####################################################################################################
# Delete the ascii files
os.system("rm *.asc")
os.chdir(parent_directory)
command = "rm -rf " + temp_rank_folder
os.system(command)

print("Done!")
