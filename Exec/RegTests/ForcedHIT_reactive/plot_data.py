import yt 
import fnmatch
import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
import math
from yt.units import mm

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.size'] = 16

yt.enable_parallelism()

kernel_height   = 1.01e-01
kernel_y_loc    = 1.92
cathode_tip = kernel_y_loc-kernel_height/2.0
anode_tip   = kernel_y_loc+kernel_height/2.0 
# Select the time 
istart = 0 # Start index
iskip  = 1 # Skip index
path = os.getcwd()

list_output=fnmatch.filter(os.listdir(path), 'plt_*')
list_output.sort()    
print("Found output   : {}".format(list_output))

list_process=list_output[istart:len(list_output):iskip]
print("To be processed: {} \n".format(list_process))

fig1, ax1 = plt.subplots()

# list_process = ['results_turbulence_forcing/plt_02000']
list_process = ['plt_00225']
scalars = ['temp','density','pressure', 'magvel','x_velocity','y_velocity', 'magvort','vfrac','Y(O2)','Y(O)','Y(CH2O)','Y(CO)','Y(OH)']
scalars = ['x_velocity','y_velocity','z_velocity','mag_vort','temp','Y(CH2O)','HeatRelease','IC8_ConsumptionRate']
for i,name in enumerate(list_process):
        ds = yt.load("{}/{}".format(path,name))
        ad = ds.all_data()
        print('In {} the minimum temperature is: {}'.format(name,ad['temp'].min()))
        print('In {} the maximum temperature is: {}'.format(name,ad['temp'].max()))
        for scalar in scalars:
            
            slc = yt.SlicePlot(ds, 'x', '{}'.format(scalar), center=[0.012/2., 1.0e-2, 0.01], width=(1.2e-2))
            
            if(scalar == 'density'):
                slc.set_cmap('{}'.format(scalar), 'RdBu')
#            elif(scalar == 'HeatRelease'):
#                ds.surface(slc, 'temp', 1500.0)
            else:
                slc.set_cmap('{}'.format(scalar), 'RdBu_r')
            slc.annotate_timestamp()
 #           slc.set_zlim('temp',300.0,2500.0)
            slc.set_log('{}'.format(scalar),log=False)
            slc.save()


 #           slc = yt.SlicePlot(ds, 'y', '{}'.format(scalar), center=[0.6e-2, 0.6e-2, 0.01], width=(1.2e-2))
 #           slc.set_log('{}'.format(scalar),log=False)
 #           slc.set_cmap('{}'.format(scalar), 'RdBu_r')
 #           slc.annotate_timestamp()
 #           slc.save()

