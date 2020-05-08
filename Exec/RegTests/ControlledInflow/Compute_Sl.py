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

def dump_data(direction, Sl_ref, Sl_AC, Sl, diffSlAC, diffSl):
    os.system("echo Ref. AC Consumption > FlameSpeedsCompa_{}.dat".format(direction))
    os.system("echo {} {} {} >> FlameSpeedsCompa_{}.dat".format(Sl_ref, Sl_AC, Sl,direction))
    os.system("echo 1.0 {} {} >> FlameSpeedsCompa_{}.dat".format(diffSlAC,diffSl,direction))

def calc_sl(args, direction):

    # User data
    fuel="H2"
    #0=X, 1=Y, 2=Z
    if direction == "X":
        flow_direction=0
    elif direction == "Y":
        flow_direction=1
    else:
        assert direction == "Z"
        flow_direction=2

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
    max_step = 0
    cur_step = 0
    for filename in pltfileiter:
        # Run fextract to get an ASCII
        tempfile = f"tempAscii_{filename}.txt"
        cur_step = int(re.sub("plt","", filename))
        max_step = max(max_step, cur_step)
        if args.PeleC_nomenclature:
            #TODO: Choose homogeneous names of variables and output for PeleC and PeleLM
            fextract_args = "-s {} -d {} -v \'density rho_omega_{} Y({}) Temp \' {}".format(tempfile,flow_direction,fuel,fuel,filename)
        else:
            fextract_args = "-s {} -d {} -v \'density {}_ConsumptionRate Y({}) temp \' {}".format(tempfile,flow_direction,fuel,fuel,filename)
        ##TODO: check for output, ex: invalid direction"
        os.system("./{} {}".format(os.path.basename(args.fextract_exe),fextract_args))

        with open(tempfile) as fin:
            line = list(islice(fin, 1, 2))
            time = line[0].split()

        x, rho, H2_consumptionRate, Y_fuel, Temp = np.loadtxt(tempfile, unpack=True, comments='#')
        # Compute the value of Sl from the integral of fuel consumption rate
        Sl = np.trapz(H2_consumptionRate,x)/(rho[0]*Y_fuel[0])

        plot_time.append(float(time[3]))
        plot_Sl.append(float(Sl))

        if ( args.plt_mode == "all" ): text_file.write('%0.14e \t %f \t %f \n'%(float(time[3]),Sl,-np.trapz(H2_consumptionRate,x)))
        #os.system("rm {}".format(tempfile))
    # Plot data if plt_mode == all
    if ( args.plt_mode == "all" ):
        plt.plot(plot_time,plot_Sl)
        print(plot_time)
        print(plot_Sl)
        plt.xlabel("time (s)")
        plt.ylabel("Flame speed (integral of fuel (H2) consumption rate)")
        plt.title("    Flame speed from PMF regtest: max_step = {}".format(max_step))
        plt.savefig("Flame_speed_from{}script.png".format(direction))
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
        dump_data(direction, Sl_ref, Sl_AC, Sl, diffSlAC,diffSl)
        plotCompa_Sl = [Sl_ref,Sl_AC,Sl]
        fig, ax = plt.subplots()
        ind = np.arange(3)
        rects1 = ax.bar(ind-0.15, plotCompa_Sl)
        ax.set_ylabel("Flame speeds [m/s]")
        ax.set_ylim(0.98*Sl_ref,1.02*Sl_ref)
        ax.set_xticks(ind)
        ax.set_xticklabels(("Reference","Vin_AC","Sl_cons"))
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels)
        plt.savefig("FlameSpeedCompa_{}.png".format(direction))

def parse_args(arg_string=None):
    parser = argparse.ArgumentParser(description=USAGE)

    parser.add_argument("--plt_mode", type=str, default="last", metavar="plt_input_mode",
                        help="Can either treat all plt* in current folder if 'all' or just the last one if 'last'. Default = last")

    parser.add_argument("--fextract_exe", type=str, default="None", metavar="fextract.exe",
                        required=True,help="path to the executable required for the analysis.")

    parser.add_argument("--compare_AC", type=bool, default=True, metavar="compare_AC",
                        help="Perform a comparison of computed Sl with AC_History data. Default = True if AC_History file found.")

    parser.add_argument("--PeleC_nomenclature", type=bool, default=False, metavar="PeleC_nomenclature",
                        help="To apply the script on plot directory generated with PeleC, the nomenclature must be adapted.")

    if not arg_string is None:
        args, unknown = parser.parse_known_args(arg_string)
    else:
        args, unknown = parser.parse_known_args()

    return args
