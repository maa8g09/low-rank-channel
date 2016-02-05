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
                    type=float)
parser.add_argument("-T1",
                    "--t_end",
                    metavar='\b',
                    help="Plotting end time unit.",
                    required=True,
                    type=float)
parser.add_argument("-f",
                    "--File",
                    metavar='\b',
                    help="File to read.",
                    required=True)
ut.print_Start_Bar()
args = parser.parse_args()
data = ut.read_Output_DNS(args.File, args.t_start, args.t_end)
ut.plot_Convergence_DNS(data, args.t_start, args.t_end)
ut.print_EndMessage()