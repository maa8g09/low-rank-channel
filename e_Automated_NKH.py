#!/usr/bin/env python3

import argparse
from datetime import datetime
import os
import numpy as np
import Utils as ut

descriptionString = "Run a NKH search on given files."
descriptionString+= "This program should be run from the directory where the DNS data is stored."
parser = argparse.ArgumentParser(description=descriptionString)
parser.add_argument("-f",
                    "--Files",
                    metavar='\b',
                    nargs='+',
                    help="Files to approximate",
                    required=True)
parser.add_argument("-sym",
                    "--SymmetryFile",
                    metavar='\b',
                    help="Symmetry file",
                    required=True)
args = parser.parse_args()


parent_directory = os.getcwd()    
# Add slash at the end of the string if there isn't one already
if parent_directory[-1] != "/":
    parent_directory += "/"


#================================================================
#### Extract Reynolds number from path
#================================================================
tmp = parent_directory[parent_directory.find("Re"):]
Re = float(tmp[2:tmp.find("/")])


#================================================================
#### Loop over files given
#================================================================
for file in args.Files:
    #------------------------------------------------
    #### Start timing
    #------------------------------------------------
    startTime = datetime.now()
    print("\nSearching " + file)
    #------------------------------------------------
    #### Remove filetype from string
    #------------------------------------------------
    fileName = file[:-3]

    #------------------------------------------------
    #### Create an NKH directory to do the searches in
    #------------------------------------------------
    nkh_directory = parent_directory + "nkh-" + fileName + "/"

    # if a nkh directory exists, stop loop and start again.
    if os.path.exists(nkh_directory):
        continue

    # if an nkh directory doesn't exist, create one.
    if not os.path.exists(nkh_directory):
        os.mkdir(nkh_directory)

    #------------------------------------------------
    #### Change into NKH directory
    #------------------------------------------------
    os.chdir(nkh_directory)
    print("Changed directory into:\n\n" + nkh_directory)
    
    #------------------------------------------------
    #### Run NKH search
    #------------------------------------------------
    TW = True
    
    # NKH Settings
    command = "findsoln"
    if TW:
        command += " -eqb "
    else:
        command += " -orb "
    
    
    command += " -xrel"         # search over x phase shift for relative orbit or eqb
    command += " -R " + str(Re) # Extracted from parent directory
    command += " -sigma ../" + args.SymmetryFile  # file containing symmetry of relative solution (default == identity)
                                            # this sigma file is a copy of the symmetry file when searching for equilibriums
    # If you leave this flag empty, the symmetry of the final solution is 
    # shifted by half domain in the streamwise direction.
    # Run and double check for yourself...

    # The sigma file only has one generator, make sure it is the proper one...

    command += " -es 1e-12"      # stop search if L2Norm(s f^T(u) - u) < epsEQB
    command += " -symms ../" + args.SymmetryFile
    command += " -is smrk2 "     # timestepping algorithm for initializing multistep algorithms [cnfe1 cnrk2 smrk2 sbdf1]
    command += " -ts sbdf3 "     # timestepping algorithm [cnfe1 cnab2 smrk2 sbdf1 sbdf2 sbdf3 sbdf4]
    command += " -nl alt "       # method of calculating nonlinearity [rot conv div skew alt]
    command += " -mc bulkv "     # hold this fixed
    command += " -Ubulk " + str(2.0/3.0)
    command += " -Uwall 0"       # static walls
    command += " -vdt"           # variable dt
    command += " -CFLmin 0.1"    # minimum  CFL
    command += " -CFLmax 1.0"    # maximum CFL
    command += " -sn"            # save field at each Newton step
    command += " -log nkh-search.log"
    command += " ../" + file
    print("")
    print(command)
    os.system(command)
    print("Done with searching "+file)

    #------------------------------------------------
    #### Stop timing
    #------------------------------------------------
    endTime = datetime.now() - startTime
    
    #------------------------------------------------
    #### Output time taken to a time logfile
    #------------------------------------------------
    timeFileName = nkh_directory + "searchTime.log"
    timeFile = open(timeFileName, "w")
    timeFile.write("\n\nNKH time for " + nkh_directory + "\n\n")
    timeFile.write(str(endTime) + "\n\n")
    timeFile.close()

print("\nDone!")
