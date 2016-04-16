#!/usr/bin/env python3

import argparse
import time
import os
import numpy as np
import Utils as ut
import FlowField as ffClass
import ChannelResolvent as cr
import h5py
import Tests


parser = argparse.ArgumentParser(description="Construct an (approximate)field at a given rank with a given mean profile (if specified).")
parser.add_argument("-f",
                    "--File",
                    metavar='\b',
                    help="File to re-assemble",
                    required=True)
parser.add_argument("-r",
                    "--Rank",
                    metavar='\b',
                    help="Number of velocity modes to use to recreate flow field.",
                    required=True,
                    type=int)
parser.add_argument("-d",
                    "--Details",
                    metavar='\b',
                    help="Details file with Re, c and bf info as .txt file. Keep in same directory as file to approximate.",
                    required=True)
parser.add_argument("-v",
                    "--MeanProfile",
                    metavar='\b',
                    help="Turbulent mean velocity profile as .txt file. Keep in same directory as file to approximate.")
parser.add_argument("-s",
                    "--Sparse",
                    help="Use sparse SVD algorithm.",
                    action='store_true')
args = parser.parse_args()
