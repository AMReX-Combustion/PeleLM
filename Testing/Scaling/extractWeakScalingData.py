#!/usr/bin/env python3

import sys
import os
import shutil
import argparse
import socket
import numpy as np
import subprocess
import fnmatch

# Script used to extract weak scaling data from PeleLM log files
# Simply ./extractWeakScalingData.py in the folder created using weakScaling.py
# Adapt the "logprefix" below to your case

if __name__ == "__main__":

    # log file pattern
    logprefix = "slurm*"

    # get list of cases from folders
    cases = [ os.path.relpath(f.path,"./") for f in os.scandir(".") if f.is_dir() ]
    cases.sort()

    # data holder
    runTimes = []
    diffusionTimes = []
    reactionTimes = []
    MacProjTimes = []
    NodalProjTimes = []
    VelAdvTimes = []
    ScalAdvTimes = []
    SyncTimes = []

    for case in cases:
        print(case)
        folder = "./{}".format(case)

        logfile = ""
        for root, dirs, files in os.walk(folder):
            for name in files:
                if fnmatch.fnmatch(name, logprefix):
                    logfile = os.path.join(root, name)
                    break

        if (logfile == ""):
            print("WARNING ! Could not find logfile in {}".format(case))
            continue

        # Get runTime
        cmd = "cat {}".format(logfile)+" | grep 'Run time =' | awk -F= '{print $2}'"
        procRunTime = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        currRunTime = procRunTime.communicate()[0].decode("utf-8").strip()
        runTimes.append(currRunTime)
        # Get component times: reaction, diffusion, MacProj, NodalProj, ScalAdv, VelAdv, Sync
        cmd = "cat {}".format(logfile)+" | grep 'PeleLM::advance::reaction' | awk 'NR%2==0' | awk '{print $4}'"
        procRunTime = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        currReactTime = procRunTime.communicate()[0].decode("utf-8").strip()
        reactionTimes.append(currReactTime)
        cmd = "cat {}".format(logfile)+" | grep 'PeleLM::advance::diffusion' | awk 'NR%2==0' | awk '{print $4}'"
        procRunTime = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        currDiffTime = procRunTime.communicate()[0].decode("utf-8").strip()
        diffusionTimes.append(currDiffTime)
        cmd = "cat {}".format(logfile)+" | grep 'PeleLM::advance::mac' | awk 'NR%2==0' | awk '{print $4}'"
        procRunTime = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        currMACTime = procRunTime.communicate()[0].decode("utf-8").strip()
        MacProjTimes.append(currMACTime)
        cmd = "cat {}".format(logfile)+" | grep 'PeleLM::advance::project' | awk 'NR%2==0' | awk '{print $4}'"
        procRunTime = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        currProjTime = procRunTime.communicate()[0].decode("utf-8").strip()
        NodalProjTimes.append(currProjTime)
        cmd = "cat {}".format(logfile)+" | grep 'PeleLM::advance::scalars_adv' | awk 'NR%2==0' | awk '{print $4}'"
        procRunTime = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        currScalAdvTime = procRunTime.communicate()[0].decode("utf-8").strip()
        ScalAdvTimes.append(currScalAdvTime)
        cmd = "cat {}".format(logfile)+" | grep 'PeleLM::advance::velocity_adv' | awk 'NR%2==0' | awk '{print $4}'"
        procRunTime = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        currVelAdvTime = procRunTime.communicate()[0].decode("utf-8").strip()
        VelAdvTimes.append(currVelAdvTime)
        cmd = "cat {}".format(logfile)+" | grep 'PLM::mac_sync()' | awk 'NR%2==0' | awk '{print $4}'"
        procRunTime = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        currSyncTime = procRunTime.communicate()[0].decode("utf-8").strip()
        SyncTimes.append(currSyncTime)

    fout = open("ScalingData.dat", "w")
    fout.write(" NodeCount RunTime Reaction Diffusion MacProj NodalProj ScalAdv VelAdv Sync \n")
    for n in range(len(cases)):
        fout.write("{} {} {} {} {} {} {} {} {} \n".format(cases[n], runTimes[n], reactionTimes[n], diffusionTimes[n], MacProjTimes[n], NodalProjTimes[n], ScalAdvTimes[n], VelAdvTimes[n], SyncTimes[n]))
    fout.close()
