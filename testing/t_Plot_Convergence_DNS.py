#!/usr/bin/env python
import os
import sys
sys.path.append("/".join(sys.path[0].split("/")[:-1]))
import Utils as ut

def main(T0, T1, File):
    # Make sure that the numbers are integers for start and end time
    T0 = int(T0)
    T1 = int(T1)
    
    data = ut.read_Output_DNS(File, T0, T1)
    ut.plot_Convergence_DNS(data, T0, T1)
    ut.plot_Convergence_DNS_log(data, T0, T1)
    ut.print_EndMessage()

directory = "/home/arslan/Documents/work/cfd-symmetry_scans/s_tw1_sigma_z_tau_x/Re1800.0/KB/2016_02_22/005_theta_-2.0000"
os.chdir(directory)
file = "data_alt.txt"
t0=0
t1=1400
main(t0,t1,file)