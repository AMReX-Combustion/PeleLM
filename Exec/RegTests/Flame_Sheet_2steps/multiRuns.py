#!/usr/bin/env python3

# A template script to lunch several times PeleLM at different resolution
# in order to evaluate the convergence order
# Used in the marex regression_testing framework

# Usage:
#   ./multiRuns.py testName PeleLMInputFile

# Input:
#   * testName: a TESTNAME that will be prependded to the plt files names
#   * PeleLMInputFile: the PeleLM input file

# "Internal" user input
#   * resolution : a list of the resolutions to run

# Head's up : 
#   * The PeleLM executable is searched for in the current directory.

import sys
import os
import shutil
import numpy as np

def multiRun(test_name,inputFile):
    
    print(" Scripted multiple runs for convergence study ")
    # User data
    resolution = [64,128,256,512,1024]        

    # Get the PeleLM exec
    run_dir = os.getcwd()
    for f in os.listdir(run_dir):
        if ( f.startswith("PeleLM") and f.endswith(".ex")):
               executable = f

    # Loop on /= resolutions, run 
    for case in resolution:
        width = case//16
        print(" Running {}x{} case".format(width,case))
        outfile = "{}_{}.run.out".format(test_name,case)
        runtime_params = "amr.n_cell={} {} {} ".format(width,case,case)
        runtime_params += "amr.plot_file={}_plt_{}_".format(test_name,case)
        os.system("mpiexec -n 1 ./{} {} {} > {}".format(executable, inputFile,runtime_params, outfile))


if __name__ == "__main__":
    test_name = str(sys.argv[1])
    inputFile = str(sys.argv[2])
    multiRun(test_name,inputFile)
