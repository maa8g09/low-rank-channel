#!/usr/bin/env python

import argparse
import time
import os
from datetime import datetime
import numpy as np

import Utils as ut
date = time.strftime("%Y_%m_%d")

ut.print_ResolventHeader()
parser = argparse.ArgumentParser(description="Run DNS with initial flow field generated by the resolvent formulation.")
parser.add_argument("-wp",
                    "--Wavepacket",
#                    metavar='\b',
                    help="The wavepacket you want to create.",
                    required=True,
                    choices=['KA', 'KB', 'KC', 'KD', 'KE'])
#                    nargs=1
parser.add_argument("-nd",
                    "--Nd",
                    metavar='\b',
                    help="The dimensions of the flow field you want to create.",
                    required=True,
                    type=int)
parser.add_argument("-nx",
                    "--Nx",
                    metavar='\b',
                    help="Number of grid point in the streamwise direction.",
                    required=True,
                    type=int)
parser.add_argument("-ny",
                    "--Ny",
                    metavar='\b',
                    help="Number of grid point in the wall-normal direction.",
                    required=True,
                    type=int)
parser.add_argument("-nz",
                    "--Nz",
                    metavar='\b',
                    help="Number of grid point in the spanwise direction.",
                    required=True,
                    type=int)
parser.add_argument("-re",
                    "--Reynolds",
                    metavar='\b',
                    help="Reynolds number to generate the flow field.",
                    required=True,
                    type=float)
parser.add_argument("-c",
                    "--Wavespeed",
                    metavar='\b',
                    help="Wave speed of the flow field.",
                    required=True,
                    type=float)
parser.add_argument("-bf",
                    "--Baseflow",
#                    metavar='\b',
                    help="The baseflow to use.",
                    required=True,
                    choices=['lam', 'cou'])
parser.add_argument("-d",
                    "--Directory",
                    metavar='\b',
                    help="Output directory.",
                    required=True)
parser.add_argument("-th_min",
                    "--theta_min",
                    metavar='\b',
                    help="Minimum theta for scaling amplitude. ( chi_tilde = 10.0**(theta) * chi )",
                    required=True,
                    type=float)
parser.add_argument("-th_max",
                    "--theta_max",
                    metavar='\b',
                    help="Maximum theta for scaling amplitude.",
                    required=True,
                    type=float)
parser.add_argument("-th_steps",
                    "--theta_steps",
                    metavar='\b',
                    help="Number of cases to generate. (If 1 sample to generate, the maximum theta will be used for scaling.)\n(theta_max, theta_min, theta_steps)",
                    required=True,
                    type=int)
parser.add_argument("-T0",
                    "--t_start",
                    metavar='\b',
                    help="DNS starting time unit.",
                    required=True,
                    type=float)
parser.add_argument("-T1",
                    "--t_end",
                    metavar='\b',
                    help="DNS ending time unit.",
                    required=True,
                    type=float)

args = parser.parse_args()


#================================================================
# Format the output directory correctly.
#================================================================
output_directory = args.Directory
if output_directory[-1] != '/':
    output_directory += '/'

print('\nThe process will be carried out in:')
print(output_directory, '\n')


#================================================================
# Define the indices array
#================================================================
powers = np.linspace(args.theta_max, args.theta_min, args.theta_steps)


#================================================================
# Start timing: initialising
#================================================================
startTime_1 = datetime.now()


#================================================================
# Initialise initial velocity fields in their respective directories
#================================================================
directories = []

for i in range(0, args.theta_steps):
    command  = "eMakeFlowField.py"
    command += " -wp " + str(args.Wavepacket)
    command += " -nd " + str(args.Nd)
    command += " -nx " + str(args.Nx)
    command += " -ny " + str(args.Ny)
    command += " -nz " + str(args.Nz)
    command += " -re " + str(args.Reynolds)
    command += " -c  " + str(args.Wavespeed)
    command += " -bf " + str(args.Baseflow)
    command += " -d  " + str(output_directory)
    command += " -dat"
    command += " -asc"
    command += " -ff"
    command += " -th " + format(powers[i], '.16f')
    #------------------------------------------------
    # Iteration
    #------------------------------------------------
    if args.theta_steps > 1:
        command += " -i " + str(i+1)
        os.system(command)
        d = str(i+1).zfill(3) +'_theta_' + format(powers[i], '.4f') + '/'
        directories.append(d)
    #------------------------------------------------
    #  No Iteration
    #------------------------------------------------
    else:
        os.system(command)
        d = 'theta_' + format(powers[0], '.4f') + '/'
        directories.append(d)


#================================================================
# End timing: initialising
#================================================================
endTime_1 = datetime.now() - startTime_1


#================================================================
# Output directory where all the cases are saved
#================================================================
parent_directory = output_directory + 'Re' + str(args.Reynolds) + '/KB/' + date + "/"



#================================================================
# Setup each directory to run DNS
#================================================================
directories = sorted(directories)

for i in range(0, len(directories)):

    ut.print_DNSHeader()
    ut.print_DNSSubHeader()

    #------------------------------------------------
    # Start timing: DNS
    #------------------------------------------------
    startTime_2 = datetime.now()


    #------------------------------------------------
    # Format the case directory
    #------------------------------------------------
    case_directory = parent_directory + directories[i]
    print("\n" + case_directory)
    if case_directory[-1] != '/':
        case_directory += '/'


    #------------------------------------------------
    # Define symmetries file
    #------------------------------------------------
    # Filename is the name of parent directory (symmetry name)
    symmsFileName = output_directory.split("/")[-2]
    symmsFile = symmsFileName + ".asc"

    if symmsFileName == "s_eq":
        ut.write_Symms_File(case_directory, symmsFile, 4, ['1 1 1 1 0.0 0.0',   # e
                                                           '1 1 1 -1 0.5 0.0',  # sigma_z tau_x
                                                           '1 1 -1 1 0.0 0.5',  # sigma_y tau_z
                                                           '1 1 -1 -1 0.5 0.5']) # sigma_yz tau_xz

    elif symmsFileName == "s_tw1_sigma_z_tau_x":
        ut.write_Symms_File(case_directory, symmsFile, 2, ['1 1 1 1 0.0 0.0',   # e
                                                           '1 1 1 -1 0.5 0.0']) # sigma_z tau_x

    elif symmsFileName == "s_tw2_sigma_z_tau_xz":
        ut.write_Symms_File(case_directory, symmsFile, 2, ['1 1 1 1 0.0 0.0',   # e
                                                           '1 1 -1 1 0.5 0.5']) # sigma_z tau_xz

    elif symmsFileName == "s_po":
        ut.write_Symms_File(case_directory, symmsFile, 4, ['1 1 1 1 0.0 0.0',   # e
                                                           '1 1 1 1 0.5 0.0',   # tau_x
                                                           '1 1 1 1 0.0 0.5',   # tau_z
                                                           '1 1 1 1 0.5 0.5'])  # tau_xz

    else:
        ut.write_Symms_File(case_directory, symmsFile, 1, ['1 1 1 1 0.0 0.0'])  # e
#        ut.write_Symms_File(case_directory, symmsFile, 1, ['1 1 1 1 0.5 0.0'])  # tau_x for testing...

    symmsFile = case_directory + symmsFile


    #------------------------------------------------
    # Define DNS flags
    #------------------------------------------------
    # I use an aliased transform.
    nonlinearity = 'alt'
    # Method of calculating nonlinearity, one of [rot conv div skew alt]
    # 
    # Zang recommends using the skew-symmetric or alternating forms with aliased transforms
    # or
    # the rotational form with dealiased transforms.
    # padded is dealiased (i.e. 2/3-style dealiasing)
    #
    # T.A. Zang.   On the rotation and skew-symmetric forms for incompressible flow simulations.
    # Appl.Numer. Math., 7:27–40, 1991.
    # 
    # Preyet, 2002, Spectral Methods for Incompressible Viscous Flow, Springer.


    timestepping = 'sbdf3'
    # SBDF2, SBDF3, SBDF4:  
    # 2nd,  3rd,  and  4th-order  Semi-implicit  Backward  Differentiation Formulae,
    # requiring 1, 2 and 3 initialization steps.
    #
    # Gibson found the SBDF schemes to be the best-behaved timestepping schemes,
    # from: Crank-Nicolson, Foward-Euler 1,
    #       Crank-Nicolson, Adams-Bashforth 2,
    #       Crank-Nicolson, Runge-Kutta,
    #       Semi-Implicit Runge-Kutta,
    #       Semi-implicit  Backward  Differentiation Formulae
    #
    # i.e. from: [cnfe1 cnab2 smrk2 sbdf1 sbdf2 sbdf3 sbdf4]
    #
    #
    # When solving u^n+1 and p^n+1, SBDF schemes enforce divergence and 
    # momentum equations at t^n+1.
    #
    # This strongly implicit formulation poduces strong damping for
    # high-frequency modes and results in pressure field as accurate as the velocity field.
    # SBDF3 is particularly good: it has the strongest asympotitc decay of all 
    # 3rd-order implicit-explicit linear multistep schemes.
    #
    # For these reasons,SBDF3 is the default value of timestepping.
    # Peyret terms these algorithms AB/BDEk (kth-order Adams-Bashforth Backward-Differentiation).
    #
    # Preyet, 2002, Spectral Methods for Incompressible Viscous Flow, Springer.


    initstepping = 'cnrk2'
    # Some of the time-stepping algorithms listed above (SBDF in particular)
    # require data from N previous time steps.  
    # Supplying these past values to the DNS constructor would
    # entail a number of tedious practical problems,
    # so the DNS class instead takes its first N steps with an
    # initialization timestepping algorithm that requires no previous data. 
    # This initialization algorithm is specified by "initstepping".
    # Valid values are CNFE, CNRK2 and SMRK2.
    #
    # The default is CNRK2.
    #
    # In my version of channelflow, this option is unavailable. *wamp wamp waaaaamp*

    
    #------------------------------------------------
    # Define Output directory and data file name.
    #------------------------------------------------
    output_dir  = case_directory + 'data-' + nonlinearity + '/'
    output_file = case_directory + 'data_' + nonlinearity + '.txt'

    #------------------------------------------------
    # Define initial flow field name and start time
    #------------------------------------------------
    t_start = int(args.t_start)
    if t_start == 0:
        initialFF = case_directory + 'u0.ff'
    else:
        initialFF = output_dir + 'u'+str(t_start) + '.000.ff'

    #------------------------------------------------
    # Define time unit variables
    #------------------------------------------------
    t_end     = args.t_end
    t_dt      = 0.001
    # Sufficiently small delta_t
    # Preyet uses 0.01 and 0.005
    # Preyet, 2002, Spectral Methods for Incompressible Viscous Flow, Springer, pg 286-289.
    t_dtmin   = 0.0001
    t_dtmax   = 0.01
    t_dtsave  = 0.5

    #------------------------------------------------
    # Define CFL limits
    #------------------------------------------------
    cfl_min   = 0.001
    cfl_max   = 1.0
    # The Courant–Friedrichs–Lewy (CFL) condition: CFL < 1.

    #------------------------------------------------
    # Define Reynolds number
    #------------------------------------------------
    reynolds = args.Reynolds
    # This Reynolds number is based on the centreline velocity and
    # the channel half-height.
    # Re = (U_c*h)/nu

    #------------------------------------------------
    # Define Reynolds number
    #------------------------------------------------
    Ubulk = 2.0/3.0
    # This is the average bulk velocity
    # 0.5 * integral^{1}_{-1} 1 - y^2 dy
    #
    # This is based off the calculation used by Gibson as shown by 
    # Tuckerman et al., 2014, Turbulent-laminar patterns in plane Poiseuille flow.

    #------------------------------------------------
    # Define DNS command for channelflow
    #------------------------------------------------
    print('\nExecuting the following command:')
    command  = "couette --channel"
    command += " -T0 "      + str(t_start)
    command += " -T1 "      + str(t_end)
    command += " -dt "      + str(t_dt)
    command += " -dtmin "   + str(t_dtmin)
    command += " -dtmax "   + str(t_dtmax)
    command += " -dT "      + str(t_dtsave)
    command += " -CFLmin "  + str(cfl_min)
    command += " -CFLmax "  + str(cfl_max)
    command += " -nl "      + nonlinearity
#    command += " -ts "      + timestepping
    command += " -symms "   + symmsFile
    command += " --outdir " + output_dir
    command += " -R "        + str(reynolds)
    command += " -b -U "    + str(Ubulk)
    command += " -cfl -l2 -D -I -dv -u -Up -p" # Print these parameters
    command += " " + initialFF
    command += " >> " + output_file
    print(command + "\n")
    
    #------------------------------------------------
    # Issue Command
    #------------------------------------------------
    os.system(command)

    #------------------------------------------------
    # Stop timing
    #------------------------------------------------
    endTime_2 = datetime.now() - startTime_2

    #------------------------------------------------
    # Output time taken to a time logfile
    #------------------------------------------------
    timeFileName = case_directory + "time_" + nonlinearity + ".log"
    timeFile = open(timeFileName, "w")
    timeFile.write("\n\nDNS time for " + case_directory + "\n\n")
    timeFile.write(str(endTime_2) + "\n\n")
    timeFile.close()


    #------------------------------------------------
    # Plot Convergence  
    #------------------------------------------------
    command = "\nePlotConvergenceDNS.py"
    command+= " -T0 " + str(t_start)
    command+= " -T1 " + str(t_end)
    command+= " -f "  + str(output_file)
    print(command + "\n")
    os.system(command)


    #------------------------------------------------
    # Plot Recurrence  
    #------------------------------------------------
    rec_dir = output_dir
    rec_T0 = int(0.7*t_end)
    rec_T1 = int(0.9*t_end)
    rec_Tmax = int(0.1*t_end)
    rec_kx = int(args.Nx/2) # Look at all available modes in x
    rec_kz = int(args.Nz/2) # Look at all available modes in z

    command = "\nePlotRecurrence.py"
    command+= " -d "  + rec_dir
    command+= " -T0 " + str(rec_T0)
    command+= " -T1 " + str(rec_T1)
    command+= " -t "  + str(rec_Tmax)
    command+= " -kx " + str(rec_kx)
    command+= " -kz " + str(rec_kz)
    print(command + "\n")
    os.system(command)


    #------------------------------------------------
    # End of DNS
    #------------------------------------------------
    ut.print_EndMessage()


#================================================================
# End of file
#================================================================
ut.print_EndMessage()
