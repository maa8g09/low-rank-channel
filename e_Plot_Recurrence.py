#!/usr/bin/env python3
import argparse
import os
import numpy as np
import Utils as ut
from matplotlib import cm as cm
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser(description="Plot recurrence file.")
parser.add_argument("-f",
                    "--File",
                    metavar='\b',
                    help="File to plot",
                    required=True)
args = parser.parse_args()
seriesFile = open(args.File, 'r')
T0=-1;T1=-1;tmax=-1
data=[]
perLine = 0
# Populate data and perLine
for j, line in enumerate(seriesFile):
    if j==0:
        # Get details
        items = line.split()
        for k in range(0,len(items)):
            if str(items[k]) == "-T0":
                T0 = int(items[k+1])
            elif items[k] == "-T1":
                T1 = int(items[k+1])
            elif items[k] == "-tmax":
                tmax = int(items[k+1])
    else:
        values = line.split()
        vl_len=len(values)
        perLine = vl_len
        
        for i in range(0,vl_len):
            data.append(values[i])

data = np.asarray(data, dtype = np.float128)

lines = int(len(data) / perLine)
data = data.reshape((perLine, lines))
data = np.flipud(data)
y = np.arange(data.shape[0])
x = np.arange(data.shape[1]) + T0
print("Plotting...")
ut.plot_Recurrence(x, y, data, T0, T1, tmax, args.File[:-4])
ut.print_EndMessage()