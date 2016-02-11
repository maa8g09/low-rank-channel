#!/usr/bin/env python

# This is an executable file that should be run from the command line.
# The command line arguments determine the type of flow field created.


## Need to add all the usual things about yourself. Name, Institute etc.


import time
import FlowField as ffClass
import ChannelResolvent as cr
import Utils as ut

def main(wp, Nd, Nx, Ny, Nz, Re, c, th, bf, d, asc, dat, hdf5):

    date = time.strftime("%Y_%m_%d")

    ####################################################################################################
    # Pass these variables to the FlowField class to construct an instance of the geometry
    flowFieldGeometry = ffClass.FlowFieldGeometry(bf,wp, Nd, Nx, Ny, Nz, Re, c, th)

    ####################################################################################################
    # Make the correct directory to store the flow field in
    output_directory = ut.make_FlowField_output_directory(d, flowFieldGeometry, date)

    ####################################################################################################
    # Pass the FlowField class to main_resolvent (this is where the flow field is generated)
    ff = cr.resolvent_formulation(flowFieldGeometry)
    ff = ffClass.FlowField(flowFieldGeometry, ff, "pp")
    ####################################################################################################
    # Save flow field to output_directory
    if dat:
        ut.write_DAT(ff, output_directory, "u0")
    
    if asc:
        ut.write_ASC(ff, output_directory, "u0")
        ut.write_GEOM(ff, output_directory, "u0")
    
    #if hdf5:
    #    ut.writeHDF5(ff, flowFieldGeometry, output_directory)
    
    # Save details of the flow field in a file called u0_details.txt
    ut.write_Details(ff, output_directory, "u0")
    
    print("\nThe flow field has been stored at")
    print(output_directory)
    ut.print_EndMessage()
    return 0

asc = True
dat = True
hdf = False
main("KB", 3, 60, 30, 40, 1000.0, 0.3, -1.0, "lam", "/home/arslan/Desktop/test", asc, dat, hdf)