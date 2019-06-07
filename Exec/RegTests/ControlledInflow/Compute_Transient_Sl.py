import os
import subprocess
import numpy as np
import sys
import glob
from itertools import islice
import re
import matplotlib.pyplot as plt


numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

listing = sorted(glob.glob('plt*'),key=numericalSort)

TimeSl = 'Transient_Sl.txt'

text_file = open(TimeSl, "w")

plot_time = []
plot_Sl = []

for filename in listing:

    print(filename)

    exeName = 'fextract.Linux.gfortran.exe '

    cmdArgs = '-p ' + filename + ' -s ' + 'tempAscii.txt' +' -d 2 -v \'density H2_ConsumptionRate Y(H2) y_velocity temp \''

    callArgs = exeName+cmdArgs

    subprocess.call(callArgs, shell=True)
    
    with open('tempAscii.txt') as fin:
        line = list(islice(fin, 1, 2))
        time = line[0].split()

    x,rho,omega,Y_fuel,xVel,Temp = np.loadtxt('tempAscii.txt', unpack=True, comments='#')

    Sl = 100*-np.trapz(omega,x)/(rho[0]*Y_fuel[0])

    plot_time.append(float(time[3]))
    plot_Sl.append(float(Sl))

    text_file.write('%0.14e \t %f \t %f \n'%(float(time[3]),Sl,-np.trapz(omega,x)))


plt.plot(plot_time,plot_Sl) 
plt.show()   
