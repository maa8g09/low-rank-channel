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
    ffg = ffClass.FlowFieldGeometry(bf,wp, Nd, Nx, Ny, Nz, Re, c, th)

    ####################################################################################################
    # Make the correct directory to store the flow field in
    output_directory = ut.make_FlowField_output_directory(d, ffg, date)

    ####################################################################################################
    # Pass the FlowField class to main_resolvent (this is where the flow field is generated)
    vff, y = cr.resolvent_formulation(ffg)
    ff = ffClass.FlowField(ffg, vff, "pp")
    ff.set_y(y)
    ####################################################################################################
    # Save flow field to output_directory
    if dat:
        ut.write_DAT(ff, output_directory, "u0")
    
    if asc:
        ut.write_ASC(ff, output_directory, "u0")
        ut.write_ASC_Py(ff, output_directory, "u0")
        ut.write_GEOM(ff, output_directory, "u0")
    
    #if hdf5:
    #    ut.writeHDF5(ff, ffg, output_directory)
    
    # Save details of the flow field in a file called u0_details.txt
    ut.write_Details(ff, output_directory, "u0")
    
    print("\nThe flow field has been stored at")
    print(output_directory)
    ut.print_EndMessage()
    return 0

asc = True
dat = True
hdf = False
main("KB", 3, 35, 35, 35, 1200.0, 0.666666666666667, 0.0, "lam", "/home/arslan/Documents/work/channelflow-related/symmetry_scans/s_eq", asc, dat, hdf)
# NX MUST EQUAL NZ