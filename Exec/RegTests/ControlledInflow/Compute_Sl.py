#!/usr/bin/env python3

import os
import sys
import shutil
import argparse
import fnmatch
import numpy as np
from itertools import islice
import re
import matplotlib.pyplot as plt

USAGE = """
    Template post-processing computing flame speed from fuel consumption rate integral and compare
    the result to the PMF file and the AC_history.
"""

def calc_sl(args):

    # User data
    fuel="H2"
    flow_direction=0    #0=X, 1=Y, 2=Z 
    Sl_ref=0.22804

    # Get a local copy of post-processing executable
    run_dir = os.getcwd()
    if ( not os.path.isfile(os.path.basename(args.fextract_exe)) ):
        shutil.copy(args.fextract_exe, run_dir)
    
    # Get the sorted list plt files
    pltfilelist=[]
    for f in os.listdir(run_dir):
        if ( not fnmatch.fnmatch(f, '*old*')):
            if (fnmatch.fnmatch(f,"*plt*")):
                pltfilelist.append(f)
    pltfilelist.sort()

    # Depending on plt mode use last or all
    pltfileiter=[]
    if ( args.plt_mode == "last" ):
        pltfileiter.append(pltfilelist[-1])
    elif ( args.plt_mode == "all" ):
        pltfileiter = pltfilelist
        TimeSl = 'Transient_Sl.txt'
        text_file = open(TimeSl, "w")
    else:
        print("Unknown plt_file argument. Should be either 'all' or 'last'")
        sys.exit()

    plot_time = []
    plot_Sl = []

    for filename in pltfileiter:
    
        # Run fextract to get an ASCII    
        tempfile = "tempAscii.txt"
        fextract_args = "-s {} -d {} -v \'density {}_ConsumptionRate Y({}) temp \' {}".format(tempfile,flow_direction,fuel,fuel,filename)
        os.system("./{} {}".format(os.path.basename(args.fextract_exe),fextract_args))

        with open(tempfile) as fin:
            line = list(islice(fin, 1, 2))
            time = line[0].split()

        x, rho, omega, Y_fuel, Temp = np.loadtxt(tempfile, unpack=True, comments='#')
    
        # Compute the value of Sl from the integral of fuel consumption rate    
        Sl = np.trapz(omega,x)/(rho[0]*Y_fuel[0])
    
        plot_time.append(float(time[3]))
        plot_Sl.append(float(Sl))
    
        if ( args.plt_mode == "all" ): text_file.write('%0.14e \t %f \t %f \n'%(float(time[3]),Sl,-np.trapz(omega,x)))
        os.system("rm {}".format(tempfile))
    
    # Plot data if plt_mode == all
    if ( args.plt_mode == "all" ):
        plt.plot(plot_time,plot_Sl) 
        plt.savefig("Flame_speed.png")
        plt.show()

    # Get the data from AC_history if required    
    if ( args.compare_AC and os.path.isfile("AC_History") ):
        with open("AC_History") as fAC:
            for line in fAC:
                pass
            last = line
            Sl_AC = float(last.split()[2])
        diffSl = abs(Sl-Sl_ref) 
        diffSlAC = abs(Sl_AC-Sl_ref) 
        os.system("echo Ref. AC Consumption > Flame_speeds.dat")
        os.system("echo {} {} {} >> Flame_speeds.dat".format(Sl_ref, Sl_AC, Sl))
        os.system("echo 1.0 {} {} >> Flame_speeds.dat".format(diffSlAC,diffSl))

def parse_args(arg_string=None):
    parser = argparse.ArgumentParser(description=USAGE)

    parser.add_argument("--plt_mode", type=str, default="last", metavar="plt_input_mode",
                        help="Can either treat all plt* in current folder if 'all' or just the last one if 'last'. Default = last")

    parser.add_argument("--fextract_exe", type=str, default="None", metavar="fextract.exe",
                        required=True,help="path to the executable required for the analysis.")

    parser.add_argument("--compare_AC", type=bool, default=True, metavar="compare_AC",
                        help="Perform a comparison of computed Sl with AC_History data. Default = True if AC_History file found.")

    if not arg_string is None:
        args, unknown = parser.parse_known_args(arg_string)
    else:
        args, unknown = parser.parse_known_args()

    return args   

if __name__ == "__main__":
    args = parse_args(arg_string=sys.argv[1:])
    calc_sl(args)
