#!/usr/bin/env python
import argparse
import Utils as ut

ut.print_Convergence_DNS()

parser = argparse.ArgumentParser(description="Convergence plot of DNS data.")

parser.add_argument("-T0",
                    "--t_start",
                    metavar='\b',
                    help="Plotting start time unit.",
                    required=True,
                    type=int)
parser.add_argument("-T1",
                    "--t_end",
                    metavar='\b',
                    help="Plotting end time unit.",
                    required=True,
                    type=int)
parser.add_argument("-f",
                    "--File",
                    metavar='\b',
                    help="File to read.",
                    required=True)
args = parser.parse_args()

# Make sure that the numbers are integers for start and end time
T0 = int(args.t_start)
T1 = int(args.t_end)

data = ut.read_Output_DNS(args.File, T0, T1)
ut.plot_Convergence_DNS(data, T0, T1)
ut.plot_Convergence_DNS_log(data, T0, T1)
ut.print_EndMessage()