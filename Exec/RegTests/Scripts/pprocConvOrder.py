#!/usr/bin/env python3

# Post-processing script for convergence analysis
# Must be used after multirun.py script
# Usage:
#   ./pprocConvOrder.py fcompare_exec


import sys
import os
import shutil
import numpy as np

import matplotlib
import matplotlib.pyplot as plt

def pproc(fcompare_exe):

    # User data
    vars=["x_velocity", "y_velocity", "avg_pressure"]
    resolution = [32,64,128]        

    # Get a local copy of fcompare
    run_dir = os.getcwd()
#    shutil.copy(fcompare_exe, run_dir)
    test_name = run_dir.split("/")[-1]

    # Run fcompare on each resolution
    errors = np.empty([len(resolution),len(vars)+1])
    for res in range(len(resolution)):
        case = resolution[res]
        errors[res,0] = case
        outfile = "error_{}.analysis.out".format(case)
        args_file = "fcompare_params_{}".format(case)
        f = open(args_file, "r")
        fcomp_args = f.read()
        f.close()
#        os.system("rm {}".format(args_file)) 
        os.system("./{} {} > {}".format(os.path.basename(fcompare_exe), fcomp_args.split("\n")[0], outfile))

        # Extract errors on each variable
        with open(outfile) as fp:
            for i, line in enumerate(fp):
                if (i >= 5):
                    var = line.split()[0]
                    for v in range(len(vars)):
                        if ( var == vars[v] ):
                            errors[res,v+1] = line.split()[1]
        # Dump data in file                                
        os.system("rm {}".format(outfile))

    # Plot data
    plotdata(errors, test_name, vars)
    writetex(errors, test_name, vars)

def plotdata(data, test_name, vars):
    # Evaluate 2nd order slope
    snd_order = data[:,1]*1.05
    for i in range(1,len(data[:,1])):
        snd_order[i] = snd_order[i-1]/np.exp(2.0*np.log(2.0))
    for i in range(0, len(vars)):    
        plt.plot(data[:,0], data[:,i+1], label="{}".format(vars[i]))
    plt.plot(data[:,0], snd_order[:],linestyle='--',color='k', label='2nd-order')
    plt.xlabel("Resolution")
    plt.ylabel("Error L2norm")
    plt.xscale("log")
    plt.yscale("log")
    plt.grid(which='both',color='k', linestyle=':', linewidth=1)
    plt.legend(bbox_to_anchor=(0.9, 0.9), loc=1, borderaxespad=0.)
    plt.savefig("CoVo_Convergence_{}.png".format(test_name))

def writetex(data, test_name, vars):
    # Evaluate order
    conv_order = np.empty([len(data[:,0])-1,len(vars)])
    for v in range(len(vars)):
        for i in range(len(conv_order[:,0])):
            conv_order[i,v] = np.log(data[i,v+1]/data[i+1,v+1])/np.log(2.0)
    fout = open("ConvTable.tex", "w")            
    fout.write("\\begin{table}[ht!]\n")
    fout.write("\centering\n")
    fout.write("\\begin{tabular}{l|")
    for i in range(len(conv_order[:,0])):
        fout.write("c ")
    fout.write("}\n")    
    fout.write("\hline\n")
    fout.write("Variable ")
    for i in range(len(conv_order[:,0])):
        fout.write("&  {}/{} ".format(data[i+1,0],data[i,0]))
    fout.write("\\\\\n\hline\hline\n")
    for v in range(len(vars)):
        fout.write("{} ".format(vars[v].replace("_","\_")))
        for i in range(len(conv_order[:,0])):
            fout.write("&  {:.3f} ".format(conv_order[i,v]))
        fout.write("\\\\\n")
    fout.write("\end{tabular}\n")
    fout.write("\caption{PeleLM convergence order}\n")
    fout.write("\label{table:conv}\n")
    fout.write("\end{table}\n")
    fout.close()


if __name__ == "__main__":
    fcompare_exe = str(sys.argv[1])
    pproc(fcompare_exe)
