#!/usr/bin/env python
import os
import numpy as np
import Utils as ut
from matplotlib import cm as cm
from matplotlib import pyplot as plt

def main(File):
    #================================================================
    #### Read output
    #================================================================
    seriesFile = open(File, 'r')
    T0=-1;T1=-1;tmax=-1
    data=[]
    perLine = 0
    # Populate data and perLine
    for j, line in enumerate(seriesFile):
        if j==0:
            # Get details
            items = line.split()
            T0 = int(items[3])
            T1 = int(items[5])
            tmax = int(items[9])
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
    
    ut.plot_Recurrence(x, y, data, T0, T1, tmax, File[:-4])

    
    
    
directory = "/home/arslan/temp"
os.chdir(directory)
file = "dist.asc"
main(file)