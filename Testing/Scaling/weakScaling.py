#!/usr/bin/env python3

import sys
import os
import shutil
import argparse
import socket
import random
import numpy as np

# Typical usage for a FlameSheet case:
#./weakScaling.py -n ScalingLM -e PeleLM3d.gnu.TPROF.MPI.CUDA.ex -i inputs.3d-regt_GPU -b batch.MACHINE --extra_files drm19_pmf.dat

# A set of known host is listed in getMachines
# This info is then used to generate the batch_script in setBatchScript
# and in launchRun.

USAGE = """    
    A template script to setup weak scaling of PeleLM.
    Given an initial problem provided by the input file, an increasing number of nodes
    is tested (2^n nodes) while increasing the problem size in x and y
"""

def parse_args(arg_string=None):
    parser = argparse.ArgumentParser(description=USAGE)

    parser.add_argument("-n","--test_name", type=str, default="Case1", metavar="test-name",
                        help="Name of the test. Default = Case1")

    parser.add_argument("-e","--exec", type=str, default="None", metavar="Pele-exec",
                        help="Pele executable name. Default = first Pele*.ex in current directory")

    parser.add_argument("-i","--input_file", type=str, default="None", metavar="Pele-input",
                        help="Pele input file name. Default = first inputs.* in current directory")

    parser.add_argument("-b","--batch_script", type=str, default="None", metavar="Cluster-batch",
                        help="The batch script name.", required=True)

    parser.add_argument("--extra_files", nargs='+', type=str, default="None", metavar="Extras",
                        help="Additional files reqs. for the run. Default = None")

    if not arg_string is None:
        args, unknown = parser.parse_known_args(arg_string)
    else:
        args, unknown = parser.parse_known_args()

    # Checks
    if (not args.exec== "None" ):
        if (not os.path.exists(args.exec)):
            sys.exit("Provided executable '{}' does not exists!".format(args.exec))

    if (not args.input_file== "None" ):
        if (not os.path.exists(args.input_file)):
            sys.exit("Provided inputs file '{}' does not exists!".format(args.input_file))

    if (not os.path.exists(args.batch_script)):
        sys.exit("Provided batch_script '{}' does not exists!".format(args.batch_script))

    if (not args.extra_files== "None" ):
        for files in args.extra_files:
            if (not os.path.exists(files)):
                sys.exit("Provided extra_files '{}' does not exists!".format(files))

    return args 

def getMachines():
    host = socket.getfqdn()

    site = "unknown"
    machine = "unknown"

    # handle NERSC
    nersc_host = os.environ.get('NERSC_HOST')
    if nersc_host is not None:
        site = "NERSC"
        if "cori" in nersc_host:
            machine = "Cori"
        if "perlmutter" in nersc_host:
            machine = "Perlmutter"
    
    if "olcf" in host:
        site = "OLCF"
        if "ascent" in host:
            machine = "Ascent"
        if "summit" in host:
            machine = "Summit"
    
    return [site,machine]

def getInitialSize(fileIn):
    ofile=open(fileIn,'r')
    sizeLine = "" 
    for line in ofile:
        if 'geometry.prob_hi' in line:
            sizeLine=line

    sizesStr=sizeLine.split("=")[1]
    sizesStr=sizesStr.split("#")[0]
    sizes=[float(x) for x in sizesStr.strip().split(" ")]
    
    ofile.close()
    return sizes

def getInitialNcell(fileIn):
    ofile=open(fileIn,'r')
    sizeLine = "" 
    for line in ofile:
        if 'amr.n_cell' in line:
            sizeLine=line

    ncellStr=sizeLine.split("=")[1]
    ncellStr=ncellStr.split("#")[0]
    ncells=[int(x) for x in ncellStr.strip().split(" ")]
    
    ofile.close()
    return ncells

def setInputFile(args,caseFolder,sizes,ncells):
    outfile="{}/{}".format(caseFolder,args.input_file)
    with open(args.input_file) as fin, open(outfile, 'w') as fout:
        for line in fin:
            lineout = line
            if 'amr.n_cell' in line:
                lineout = "amr.n_cell   = {} {} {}   # Level 0 number of cells \n".format(ncells[0],ncells[1],ncells[2])
            if 'geometry.prob_hi' in line:
                lineout = "geometry.prob_hi   = {} {} {} # x_hi y_hi (z_hi)\n".format(sizes[0],sizes[1],sizes[2])
            fout.write(lineout)

def setBatchScript(args,caseFolder,case,host):
    outfile="{}/{}".format(caseFolder,args.batch_script)
    with open(args.batch_script) as fin, open(outfile, 'w') as fout:
        for line in fin:
            lineout = line
            # Ascent or Summit
            if any(h in host[1] for h in ("Ascent","Summit")):
                if "-nnodes" in line:
                    lineout = "#BSUB -nnodes {} \n".format(case)
                if "-J" in line:
                    lineout = "#BSUB -J {}_{:04d} \n".format(args.test_name,case)
            # Perlmutter: 4 CPU/GPU couples / nodes
            if "Perlmutter" in host[1]:
                if "SBATCH -n" in line:
                    lineout = "#SBATCH -n {} \n".format(case*4)
                if "SBATCH -J" in line:
                    lineout = "#SBATCH -J {}_{:04d} \n".format(args.test_name,case)
            fout.write(lineout)

def launchRun(args,host):
    # Ascent or Summit
    if any(h in host[1] for h in ("Ascent","Summit")):
        os.system("bsub {}".format(args.batch_script))
    if any(h in host[1] for h in ("Cori","Perlmutter")):
        os.system("sbatch {}".format(args.batch_script))

if __name__ == "__main__":
    args = parse_args(arg_string=sys.argv[1:])

    # Setup 
    # number of system nodes
#    nodeCounts = [1,2,4,8,16,32,64,128,256,512]
    nodeCounts = [1,2,4,8,16,32,64,128,256,512,1024]

    # Get local info
    host = getMachines()
    print("Performing weak scaling on {} of {}".format(host[1],host[0]))

    # Get the Pele exec
    run_dir = os.getcwd()
    if (args.exec == "None"):
        for f in os.listdir(run_dir):
            if ( f.startswith("Pele") and f.endswith(".ex")):
                executable = f
                break
        print("Using {} executable".format(executable))
    else:
        executable = args.exec

    # Check for the input file: first inputs.* is default
    if ( args.input_file == "None" ):
        for f in os.listdir(run_dir):
            if ( f.startswith("inputs") ):
                args.input_file = f
                print("Using {} input file".format(args.input_file))
                break

    # Create test folder
    if os.path.exists(args.test_name):
        oldCase = "{}.old.{}".format(args.test_name,random.randint(10000, 99999))
        shutil.move(args.test_name,oldCase)
        print("Case folder already exists. Old folder moved to {}".format(oldCase))
    os.mkdir(args.test_name)

    # Get initial problem size. Weak scaling increases alternatively the x and y direction ncells/size.
    initNcells = getInitialNcell(args.input_file)
    initSizes  = getInitialSize(args.input_file)

    # Create case tree
    sizes  = initSizes
    ncells = initNcells
    for nCount in range(len(nodeCounts)):
        case = nodeCounts[nCount]       
        if case != 1:
            if (nCount % 2) == 0: 
                sizes[0] *= 2
                ncells[0] *= 2
            else:
                sizes[1] *= 2
                ncells[1] *= 2
        caseFolder = "{}/{:04d}".format(args.test_name,case)
        os.mkdir(caseFolder)
        shutil.copy(executable,"{}/{}".format(caseFolder,executable))
        if (not args.extra_files== "None" ):
            for files in args.extra_files:
                shutil.copyfile(files,"{}/{}".format(caseFolder,files))
        setInputFile(args,caseFolder,sizes,ncells)
        setBatchScript(args,caseFolder,case,host)

    # Launch the cases
    for case in nodeCounts:
        caseFolder = "{}/{:04d}".format(args.test_name,case)
        os.chdir(caseFolder)
        print("Launching {} nodes case".format(case))
        launchRun(args,host)
        os.chdir(run_dir)
