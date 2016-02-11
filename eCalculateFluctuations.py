#!/usr/bin/env python
import argparse
import os
import Utils as ut
import FlowField as ffClass

parser = argparse.ArgumentParser(description="Calculate mean velocity field of converged data.\n(Both files must be in the same directory)")

parser.add_argument("-m",
                    "--MeanFile",
                    metavar='\b',
                    help="Mean velocity field.",
                    required=True)
parser.add_argument("-f",
                    "--File",
                    metavar='\b',
                    help="Total velocity field.",
                    required=True)
parser.add_argument("-d",
                    "--Directory",
                    metavar='\b',
                    help="Directory where the details file is.",
                    required=True)

ut.print_Start_Bar()
args = parser.parse_args()


pwd = os.getcwd()
if pwd[-1] != '/':
    pwd += '/'
d = args.Directory
if d[-1] != '/':
    d += '/'

#================================================================
#### Read total (instantaneous) vel field
#================================================================
var = ut.read_FF(pwd, args.File[:-3])

var2 = ut.read_Details(d, "u0")


#================================================================
#### Construct instance of ffClass
#================================================================
ff = ffClass.FlowFieldGeometry(var2['bf'],
                               var2['wp'],
                               var2['Nd'],
                               var2['Nx'],
                               var2['Ny'],
                               var2['Nz'],
                               var2['Re'],
                               var2['c'],
                               var2['theta'])

ff = ffClass.FlowField(ff, var['ff'], "pp")


#================================================================
#### Read mean file
#================================================================
var = ut.read_FF(pwd, args.MeanFile[:-3])


#================================================================
#### Construct instance of ffClass
#================================================================
meanff = ffClass.FlowField(ff, var['ff'], "pp")


#================================================================
#### Calculate the difference
#================================================================

delta = ut.calculate_Difference(ff.velocityField, meanff.velocityField) # ff1 - ff2


#================================================================
# Construct instance of ffClass
#================================================================
delta = ffClass.FlowField(ff, delta, "pp")


#================================================================
#### Save ASC, GEOM, FF
#================================================================
fileName = args.File[:-3] + "_flucs"
ut.write_ASC(delta, pwd, fileName)
ut.write_GEOM(delta, pwd, fileName)
ut.write_FF(pwd, fileName)


#================================================================
#### Remove excess files
#================================================================
command = "rm " + fileName + ".asc " + fileName + ".geom "
os.system(command)


ut.print_EndMessage()
