
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from collections import OrderedDict
import seaborn
seaborn.set(style="white", context="notebook", font_scale=1.5,
            rc={"axes.grid": True, "legend.frameon": False,
                "lines.markeredgewidth": 1.4, "lines.markersize": 10})
#import generateUfile
#import generateControlDict
import os
from subprocess import call
import sys
# the mock-0.3.1 dir contains testcase.py, testutils.py & mock.py


F = lambda U, Cd, rho, S, N2gfCoeffcient=0.00980665: .5*rho*U**2*Cd*S/N2gfCoeffcient
c2 = lambda rho1, rho2, h1, h2, g=9.81: g*(rho2 - rho1)/rho2*h1*h2/(h1 + h2) #Maximume wave speed squared
c2Hester = lambda rho1, rho2, h1, h2, g=9.81: (rho2 - rho1)*g/(rho2/h2 + rho1/h1)
c2Grue = lambda rho1, rho2, h1, h2, g=9.81: (rho2 - rho1)*g/rho2 * h1*h2/(h1 + h2)
c2Emsmailpour = lambda rho1, rho2, h1, g=9.81: g*h1*(rho2 - rho1)/rho1
Fr = lambda U, c2: U/c2**.5

def load_force_coeffs(filename): 
    """Load force coefficients into a DataFrame for the given time directory."""
    data = np.loadtxt(filename, skiprows=9)
    df = pd.DataFrame()
    df["time"] = data[:, 0]
    df["cl"] = data[:, 3]
    df["cd"] = data[:, 2]
    df["cm"] = data[:, 1]
    wantedIndexStart  = df.time[df.time > 1].index[0]
    wantedIndexStop  = df.time[df.time == 40].index[0]
    df2 = df[wantedIndexStart:wantedIndexStop] #wantedIndexStart
    return df2.time, df2.cd

#%% 
us = [0.14, 0.22]
rho1 = 997 # densisty fresh
rho2 = 1024 # salt
h1s = np.asarray([0.1, 0.2]) #fresh water height
 

interfaces = ['01', '02']
velocities = ['014', '022']
meshes = ['coarse', 'medium', 'fine']

resultsDict = OrderedDict()


for interface in interfaces:
    resultsDict[interface] = OrderedDict()
    for velocity in velocities:
        resultsDict[interface][velocity] = OrderedDict()
        for mesh in meshes:
            filename  = '/home/karl/OpenFOAM/karl-5.0/run/MT/MT_PK/twoPhaseIF2/freshSaltMixingMeshSensitivity/%s/interface%s/forceCoeffs/forceCoeffVelocity%s' %(mesh, interface, velocity)
            resultsDict[interface][velocity][mesh] = OrderedDict()
            resultsDict[interface][velocity][mesh]['time'] = load_force_coeffs(filename)[0]
            resultsDict[interface][velocity][mesh]['cd'] = load_force_coeffs(filename)[1]


#%% Plotting
fontSize = 23
legends = ['Coarse', 'Medium', 'Fine']
titles = ['0.1', '0.2']
velocityCounter = 0
for velocity in velocities:
    counterTitles = 0
    for interface in interfaces:
        fig, ax = plt.subplots()
        ax.tick_params(axis='both', which='major', labelsize=fontSize)
        #ax.tick_params(axis='both', which='minor', labelsize=10)
        ax.set_ylim(0.05, 0.14)
        counter = 0
        for mesh in meshes:
            ax.plot(resultsDict[interface][velocity][mesh]['time'], \
                 resultsDict[interface][velocity][mesh]['cd'], linewidth=4.0,\
                 label=legends[counter])
            counter += 1
        plt.legend(fontsize = fontSize)
        ax.set_xlabel('Time', fontsize = fontSize)
        ax.set_ylabel(r'$C_d$', fontsize = fontSize*1.25, rotation = 90)
        c2twoPhaseIf = c2Emsmailpour(rho1, rho2, h1s[counterTitles])
        FrS = Fr(us[velocityCounter], c2twoPhaseIf)
        ax.set_title('Pycnocline %sm \n Fr %.2f' %(titles[counterTitles], FrS),\
                     fontsize=fontSize)
        #fig.tight_layout()
        fig.savefig('/home/karl/OpenFOAM/karl-5.0/run/MT/MT_PK/plots/meshSensitivity/U%sInterface%s' %(velocity, interface),\
                    bbox_inches = 'tight')
        counterTitles += 1
    velocityCounter += 1

