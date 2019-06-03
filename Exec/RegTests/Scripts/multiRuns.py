#!/usr/bin/env python3

# Script to lunch several times PeleLM at different resolution
# in order to evaluate the convergence order

# ./multiRuns.py testName PeleLMInputFile

import sys
import os
import shutil
import numpy as np

def multiRun(test_name,inputFile):
    
    print(" Scripted multiple runs for convergence study ")

    # Get the PeleLM exec
    run_dir = os.getcwd()
    for f in os.listdir(run_dir):
        if ( f.startswith("PeleLM") and f.endswith(".ex")):
               executable = f

    # Loop on /= resolutions, run & prepare fcompare input
    resolution = [32,64,128]        
    for case in resolution:
        outfile = "{}_{}.run.out".format(test_name,case)
        runtime_params = "amr.n_cell={} {} {} ".format(case,case,case)
        runtime_params += "amr.plot_file={}_plt_{}_".format(test_name,case)
#        os.system("mpiexec -n 1 ./{} {} {} > {}".format(executable, inputFile,runtime_params, outfile))
        pltfile=[]
        for f in os.listdir(run_dir):
            if (f.startswith("{}_plt_{}_".format(test_name,case))):
                pltfile.append(f)
        pltfile.sort()        
        f = open("fcompare_params_{}".format(case),"w")
        f.write(" -n 2 {} {}".format(pltfile[0], pltfile[-1]))
        f.close()


if __name__ == "__main__":
    test_name = str(sys.argv[1])
    inputFile = str(sys.argv[2])
    multiRun(test_name,inputFile)
