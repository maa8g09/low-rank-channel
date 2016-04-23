#!/usr/bin/env python
import argparse
import os
import numpy as np
import Utils as ut
from matplotlib import cm as cm
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser(description="Generate a recurrence plot, over time range with time shift (both specified).")
parser.add_argument("-d",
                    "--Directory",
                    metavar='\b',
                    help="Output directory. The directory should be the data-X/ folder.",
                    required=True)
parser.add_argument("-T0",
                    "--T0",
                    metavar='\b',
                    help="Plot starting time unit.",
                    required=True,
                    type=int)
parser.add_argument("-T1",
                    "--T1",
                    metavar='\b',
                    help="Plot ending time unit.",
                    required=True,
                    type=int)
parser.add_argument("-t",
                    "--tmax",
                    metavar='\b',
                    help="Plot shift unit.",
                    required=True,
                    type=int)
parser.add_argument("-kx",
                    "--kxmax",
                    metavar='\b',
                    help="L2Norm for |kx|<kxmax, 0=>all",
                    required=True,
                    type=int)
parser.add_argument("-kz",
                    "--kzmax",
                    metavar='\b',
                    help="L2Norm for |kz|<kzmax, 0=>all",
                    required=True,
                    type=int)
args = parser.parse_args()


# The directory should be the data-X/ folder.
output_directory = args.Directory
if output_directory[-1] != '/':
    output_directory += '/'


#### Check directory
# Now check to see if the directory has data after the penultimate '/'
present_directory = output_directory.split("/")[-2]
if present_directory.find("data") != -1:
    print("In the correct directory")
    
else:
    message = "Invalid directory flag. Please add data-X/ folder at the end. Present directory: " + present_directory
    ut.error()


#### Set variables
tmax = args.tmax
T0 = args.T0
T1 = args.T1
steps = T1 - T0 + 1
TRange = np.linspace(T0, T1, steps)


#### Make Directory
# Make a directory to store the files with integer names.
ints_directory = output_directory + 'ints'
# if directory exists, delete it:
if os.path.exists(ints_directory):
    print("\nRemoving:\t" + str(ints_directory))
    command = "rm -rf " + ints_directory
    os.system(command)
# if directory doesn't exist, make it:
if not os.path.exists(ints_directory):
    print("Making:\t\t" + str(ints_directory))
    os.mkdir(ints_directory)


#### Change directory down
# change into the ints folder and list files and rename the float files to integer files.
os.chdir(ints_directory)

files = [fi for fi in os.listdir(output_directory) if os.path.isfile(os.path.join(output_directory,fi))]

for k in files:
    if str(k)[0] == 'u' and str(k).find("asc") == -1:
        oldFile = k
        k = k[1:-3]
        k = float(k)

        if k in TRange:
            k = int(k)
            newFile = 'u' + str(k) + '.ff'
            oldFile = '../' + oldFile
            command = 'cp ' + str(oldFile) + ' ' + str(newFile)
            os.system(command)

print('\nRenamed and stored flow fields with integers in names.')


#### Change directory up
# change into the data-X directory, i.e. go up one directory from the ints directory.
print("\nCurrently in:\t" + str(ints_directory))
os.chdir(os.path.dirname(os.getcwd()))
print("Changed to:\t" + str(os.getcwd()))

#### Run Seriesdist
seriesFile = output_directory + "seriesdist-" + str(T0) + "-" + str(T1) + "-" + str(args.kxmax) + "-" + str(args.kzmax)
command = "seriesdist "
command+= " -T0 " + str(T0)             # start time
command+= " -T1 " + str(T1 - tmax)             # end time
command+=" -dT 1"                       # save interval
command+=" -tmax " + str(tmax)          # maximum t' for L2Dist(u(t), t(t+t'))
command+=" -da ints/"                   # flowfield series directory A
command+=" -db ints/"                   # flowfield series directory B
command+=" -xs"                         # translate ub by Lx/2
command+=" -o " + seriesFile
command+=" -dg 16"                      # # digits in output
command+=" -kx " + str(args.kxmax)      # L2Norm for |kx|<kxmax, 0=>all
command+=" -kz " + str(args.kzmax)      # L2Norm for |kz|<kzmax, 0=>all
print("\nRunning:")
print(command)
print("")
os.system(command)


#### Read output
print('\nNow generating recurrence plot.')
seriesFile += ".asc"
seriesFile = open(seriesFile, 'r')

data=[]
perLine = 0
# Populate data and perLine
for j, line in enumerate(seriesFile):
    if j!=0:
        values = line.split()
        vl_len=len(values)
        perLine = vl_len
        
        for i in range(0,vl_len):
            data.append(values[i])

data = np.asarray(data, dtype = np.float128)

lines = len(data) / float(perLine)
lines = int(lines)
data = data.reshape((perLine, lines))
y = np.arange(data.shape[0])
x = np.arange(data.shape[1])
x+= T0


#### Plot output
x, y = np.meshgrid(x, y)
origin = 'lower'
v = 7
#v = np.linspace(0.0, 0.4, 7, endpoint=True)
CS = plt.contourf(x, y, data, v, cmap=cm.jet, origin=origin)
cbar = fig.colorbar(CS,ticks=ticks_at,format='%1.2g')
plt.xlabel('$t$ (Time units)')
plt.ylabel('$T$ (shift)')
plt.title('$||u(t)  -  u(t+T)||$')
plt.grid(True)
plt.xlim([T0, T1-tmax])
plt.ylim([0, tmax])
#### Change up
os.chdir(os.path.dirname(os.getcwd()))
output_directory = os.getcwd() + "/"
plotName = output_directory + "plot_recurrence_" + str(T0) + "-" + str(T1) + "-" + str(args.kxmax) + "-" + str(args.kzmax) + ".png"
plt.savefig(plotName)

# Now you have to delete the "ints" folder.
# if directory exists, delete it:
if os.path.exists(ints_directory):
    print("Deleting:\t" + str(ints_directory))
    command = "rm -rf " + ints_directory
    os.system(command)

ut.print_EndMessage()