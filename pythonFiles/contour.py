#%%

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from collections import OrderedDict
import seaborn
seaborn.set(style="white", context="notebook", font_scale=1.5,
            rc={"axes.grid": True, "legend.frameon": False,
                "lines.markeredgewidth": 1.4, "lines.markersize": 4})
#import generateUfile
#import 
import os
from subprocess import call
import sys
import matplotlib.ticker as plticker
import matplotlib as mpl

# the mock-0.3.1 dir contains testcase.py, testutils.py & mock.py

#%% Functions

F = lambda U, Cd, rho, S, N2gfCoeffcient=0.00980665: .5*rho*U**2*Cd*S/N2gfCoeffcient
c2 = lambda rho1, rho2, h1, h2, g=9.81: g*(rho2 - rho1)/rho2*h1*h2/(h1 + h2) #Maximume wave speed squared
c2Hester = lambda rho1, rho2, h1, h2, g=9.81: (rho2 - rho1)*g/(rho2/h2 + rho1/h1)
c2Grue = lambda rho1, rho2, h1, h2, g=9.81: (rho2 - rho1)*g/rho2 * h1*h2/(h1 + h2)
c2Emsmailpour = lambda rho1, rho2, h1, g=9.81: g*h1*(rho2 - rho1)/rho1
Fr = lambda U, c2: U/c2**.5
FrSingle = lambda U: U/np.sqrt(9.81*0.6)


def loadPoints(filename, y): 
    tol = 1e-4
    df = pd.read_csv(filename)
    wantedIndexes = df.Points1[np.abs(df.Points1 - y) < tol]
    wantedIndexes2 = wantedIndexes.index[0:-1]
    df2 = pd.DataFrame()
    df2['x'] = df.Points0[wantedIndexes2]
    df2['z'] = df.Points2[wantedIndexes2]
    return df2



def loadVelocityBox(filename, xValue, xtol, ytol): 
    df = pd.read_csv(filename)
    tolFactor = 6.25
    numberOfPoints = 15
    wantedZStart = -0.085
    wantedZStop = -0.225
    wantedIndexes = df.Points2[df.Points2 < wantedZStart] 
    wantedIndexes = wantedIndexes.index[0:-1]
    
    wantedIndexes2 = df.Points2[df.Points2 > wantedZStop]
    wantedIndexes2 = wantedIndexes2.index[0:-1]
    wantedIndexesX = [2]
    wantedIndexesY = [2]
    while len(wantedIndexesX) < numberOfPoints or len(wantedIndexesY) < numberOfPoints:
        if len(wantedIndexesX) < numberOfPoints:
            wantedIndexesX = df.Points0[abs(df.Points0 - xValue) < xtol] 
            xtol *= tolFactor
        if len(wantedIndexesY) < numberOfPoints:
            wantedIndexesY = df.Points1[abs(df.Points1) < ytol] 
            ytol *= tolFactor
    print 'x ', xValue, 'len x ', len(wantedIndexesX), 'len y', len(wantedIndexesY), '\n'

    wantedIndexesY = wantedIndexesY.index[0:-1]
    wantedIndexesX = wantedIndexesX.index[0:-1]
        
    wantedIndexes = wantedIndexes.intersection(wantedIndexes2)
    wantedIndexes = wantedIndexes.intersection(wantedIndexesY)
    wantedIndexes = wantedIndexes.intersection(wantedIndexesX)
    df2 = pd.DataFrame()
    df2['x'] = df.Points0[wantedIndexes]
    df2['y'] = df.Points1[wantedIndexes]
    df2['z'] = df.Points2[wantedIndexes]
    df2['Ux'] = df.U0[wantedIndexes]
    df2 = df2.sort_values('z')
    return df2

def loadPressure(filename, xValue, xtol, ytol): 
    df = pd.read_csv(filename0b)
    tolFactor = 6.25
    numberOfPoints = 15
    
    wantedZStart = 0.0
    wantedZStop = -0.125
    
    wantedXStart = 0.7
    wantedXStop = -0.1
    
    wantedIndexesZ = df.Points2[df.Points2 < wantedZStart] 
    wantedIndexesZ = wantedIndexesZ.index[0:-1]
    
    wantedIndexesZ2 = df.Points2[df.Points2 > wantedZStop]
    wantedIndexesZ2 = wantedIndexesZ2.index[0:-1]
    
    wantedIndexesX = df.Points0[df.Points0 < wantedXStart] 
    wantedIndexesX = wantedIndexesX.index[0:-1]
    
    wantedIndexesX2 = df.Points2[df.Points0 > wantedXStop]
    wantedIndexesX2 = wantedIndexesX2.index[0:-1]
    
    wantedIndexesY = [2]
    while len(wantedIndexesY) < numberOfPoints:
        wantedIndexesY = df.Points1[abs(df.Points1) < ytol] 
        ytol *= tolFactor
    
    wantedIndexesY = wantedIndexesY.index[0:-1]
    
    wantedIndexes = wantedIndexesZ   
    wantedIndexes = wantedIndexes.intersection(wantedIndexesZ2)
    wantedIndexes = wantedIndexes.intersection(wantedIndexesY)
    wantedIndexes = wantedIndexes.intersection(wantedIndexesX)
    wantedIndexes = wantedIndexes.intersection(wantedIndexesX2)
    df2 = pd.DataFrame()
    df2['x'] = df.Points0[wantedIndexes]
    df2['y'] = df.Points1[wantedIndexes]
    df2['z'] = df.Points2[wantedIndexes]
    df2['p_rgh'] = df.p_rgh[wantedIndexes]
    return df2



fontSize = 23

h2 = 3.8# 4.#16 #salt "
rho1 = 997 # densisty fresh
rho2 = 1024 # salt
g = 9.81

#%%  Eta Interface 02
yWidth= 0
h1s = 0.2 #fresh water height
Us01 = np.asarray([0.1, 0.22])
Us02 = np.asarray([0.08, 0.16, 0.22])
Us = Us02
U2inCm2 = [(Us[i]*100)**2 for i in xrange(len(Us))] # U**2, in cm
c2twoPhaseIf = c2Emsmailpour(rho1, rho2, h1s)
FrS = Fr(Us, c2twoPhaseIf) 

filename0 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/zg/interface02U008.csv'
filename1 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/zg/interface02U016.csv'
#filename1 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/interface02U022.csv'
filename2 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/zg/interface02U022.csv'
filename3 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/interface02SingleU016.csv'


#filename1 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/interface01U01.csv'
#filename2 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/interface01U022.csv'
    
filenames = [filename0, filename1, filename2]#, filename3]
x1 = np.arange(0, 0.6, .1)
y1 = np.ones_like(x1)*(-.1)
y2 = np.arange(-0.1, 0+.01, .01)
x2 = np.ones_like(y2)*.6
x3 = x1[::-1]
y3 = np.zeros_like(y1)
y4 = y2[::-1]
x4 = np.zeros_like(y4)
x = np.concatenate((x1, x2, x3, x4))
y = np.concatenate((y1, y2, y3, y4))

xForLine = np.array((0, 0))
yForLine = np.array((-.1, -.25))


fig, ax = plt.subplots()
ax.plot(x, y, 'k', label = 'Barge')
ax.plot(xForLine, yForLine,  'k--')
markers = ['o', 'x', "^", 'P']
colors = ['b', 'g', 'r', 'c']
counter = 0
for filename, markerType, colorType in zip(filenames, markers, colors):
    points= loadPoints(filename, yWidth)
    if filename != filename3:
        ax.plot(points.x, points.z, markerType,  label='$Fr_h$ %.2f' %FrS[counter], ms=5)
    else:
        ax.plot(points.x, points.z, 'mP', label='Non-stratified')
    counter += 1
#plt.axis('equal')
ax.text(-0.8, -0.095, 'Stern', fontsize=23)
ax.set_xlim(-1, 1.5)
ax.set_ylim(-0.25, -0.075)
ax.set_ylabel(r'$z$ [m]', fontsize = fontSize*1.25)
ax.set_xlabel('x-position [m]', fontsize = fontSize*1.25)
#legends = ['Barge', 'Fr %.2f' %FrS[0], 'Fr %.2f' %FrS[1] ]
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 1.0, box.height*.8])
ax.legend(loc='center left', bbox_to_anchor=(0.0, -0.45), \
           fontsize = fontSize*1.25, ncol=2)
#ax.axis('equal')
#fig5.tight_layout()
ax.tick_params(axis='both', which='major', labelsize=fontSize*1.25)
xtickValues = np.arange(-2.0,1.0+.1, 1.0)
ytickValues = np.arange(-.2,-0.125+0.05, 0.05)
ax.set_xticks(xtickValues)
ax.set_yticks(ytickValues)
#ax5.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
ax.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
ax.grid(b=True, which='minor', color='grey', linewidth=0.2)
ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.25))
ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.025))
fig.savefig('/home/peterek/OpenFOAM/PETER-5.0/MT_PK/plots/3Dtop/Surface'\
                , bbox_inches = 'tight')


    
#%% Velocity profiles interface 02
filename0b = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/data/zg/interface02U008box.csv'
filename1b = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/data/zg/interface02U016box.csv'
filename2b = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/data/zg/interface02U022box.csv'
filename3b = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/data/interface02SingleU016box.csv'
filenames2 = [filename0b, filename1b, filename2b, filename3b]
U0 = [0.08, 0.16, .22, 0.16]
xForLine = np.array((-1, -1))
yForLine = np.array((-.225, -.075))

colors = ['b', 'g', 'r', 'm']
eps = 0
xValues = np.asarray((0.0, 0.1, 0.3, 0.4, 0.5, 0.6))+eps
#xValue = 0.1
counter = 0

xtol = 1e-3
ytol = 1e-3

for xValue in xValues:
    uCounter = 0
    fig, ax = plt.subplots()
    ax.plot(xForLine, yForLine,  'k--')
    for filename, colorType in zip(filenames2, colors):
        zAndU = loadVelocityBox(filename, xValue, xtol, ytol)#loadVelocity(filename)
        #zAndU.sort('z')
        if filename != filename3b:
            ax.plot(zAndU.Ux/U0[uCounter], zAndU.z, colorType, \
                    label='$Fr_h$ %.2f' %FrS[uCounter], \
                    linewidth=4.0)
        else:
            ax.plot(zAndU.Ux/U0[uCounter], zAndU.z, colorType, \
                label='Non-stratified', \
                linewidth=4.0)
        uCounter += 1
    ax.set_title('x position %.2f [m]' %xValue, fontsize = fontSize*1.25)
    ax.legend(fontsize = fontSize)
    ax.set_xlabel(r'$U_x/U_0$ ', fontsize = fontSize*1.25)
    ax.set_ylabel('z-position [m]', fontsize = fontSize*1.25)
    ax.tick_params(axis='both', which='major', labelsize=fontSize*1.25)
    xtickValues = np.arange(-1.2,0+.1, .4)
    ytickValues = np.arange(-.25,-0.1+0.05, 0.05)
    ax.set_xticks(xtickValues)
    ax.set_yticks(ytickValues)
    #ax5.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
    ax.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
    ax.grid(b=True, which='minor', color='grey', linewidth=0.2)
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.2))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.025))
    fig.savefig('/home/peterek/OpenFOAM/PETER-5.0/MT_PK/plots/3Dtop/VelocityProfileX%d'\
                    %counter, bbox_inches = 'tight')
    counter += 1
    
#%%  Eta Interface 015
yWidth = 0.0
h1s = 0.15 #fresh water height
Us01 = np.asarray([0.1, 0.22])
Us02 = np.asarray([0.08, 0.12, 0.22])
Us = Us02
U2inCm2 = [(Us[i]*100)**2 for i in xrange(len(Us))] # U**2, in cm
c2twoPhaseIf = c2Emsmailpour(rho1, rho2, h1s)
FrS = Fr(Us, c2twoPhaseIf) 
    
    
filename0 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/zg/interface015U008.csv'
filename1 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/zg/interface015U012.csv'
#filename1 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/interface02U022.csv'
filename2 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/zg/interface015U022.csv'
filename3 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/interface02SingleU016.csv'


#filename1 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/interface01U01.csv'
#filename2 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/interface01U022.csv'
    
filenames = [filename0, filename1, filename2]#, filename3]
x1 = np.arange(0, 0.6, .1)
y1 = np.ones_like(x1)*(-.1)
y2 = np.arange(-0.1, 0+.01, .01)
x2 = np.ones_like(y2)*.6
x3 = x1[::-1]
y3 = np.zeros_like(y1)
y4 = y2[::-1]
x4 = np.zeros_like(y4)
x = np.concatenate((x1, x2, x3, x4))
y = np.concatenate((y1, y2, y3, y4))

xForLine = np.array((0, 0))
yForLine = np.array((-.1, -.25))

fig, ax = plt.subplots()
ax.plot(x, y, 'k', label = 'Barge')
ax.plot(xForLine, yForLine,  'k--')
markers = ['o', 'x', "^", 'P']
colors = ['b', 'g', 'r', 'c']
counter = 0
for filename, markerType, colorType in zip(filenames, markers, colors):
    points= loadPoints(filename, yWidth)
    if filename != filename3:
        ax.plot(points.x, points.z, markerType,  label='$Fr_h$ %.2f' %FrS[counter], ms=5)
    else:
        ax.plot(points.x, points.z, 'mP', label='Non-stratified')
    counter += 1
#plt.axis('equal')
ax.text(-0.8, -0.095, 'Stern', fontsize=20)
ax.set_xlim(-1, 1.5)
ax.set_ylim(-0.25, -0.075)
ax.set_ylabel(r'$z$ [m]', fontsize = fontSize*1.25)
ax.set_xlabel('x-position [m]', fontsize = fontSize*1.25)
#legends = ['Barge', 'Fr %.2f' %FrS[0], 'Fr %.2f' %FrS[1] ]
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 1.0, box.height*.8])
ax.legend(loc='center left', bbox_to_anchor=(0.0, -0.45), \
           fontsize = fontSize*1.25, ncol=2)
#ax.axis('equal')
#fig5.tight_layout()
ax.tick_params(axis='both', which='major', labelsize=fontSize*1.25)
xtickValues = np.arange(-2.0,1.0+.1, 1.0)
ytickValues = np.arange(-.2,-0.125+0.05, 0.05)
ax.set_xticks(xtickValues)
ax.set_yticks(ytickValues)
#ax5.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
ax.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
ax.grid(b=True, which='minor', color='grey', linewidth=0.2)
ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.25))
ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.025))
fig.savefig('/home/peterek/OpenFOAM/PETER-5.0/MT_PK/plots/3Dtop/Surface15'\
                , bbox_inches = 'tight')

#%% Velocity profiles interface 015
filename0b = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/data/zg/interface015U008box.csv'
filename1b = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/data/zg/interface015U012box.csv'
filename2b = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/data/zg/interface015U022box.csv'
filename3b = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/data/interface02SingleU016box.csv'
filenames2 = [filename0b, filename1b, filename2b, filename3b]
U0 = [0.08, 0.12, .22, 0.16]
xForLine = np.array((-1, -1))
yForLine = np.array((-.225, -.075))

colors = ['b', 'g', 'r', 'm']
xValues = [0.0, 0.1, 0.3, 0.4, 0.5, 0.6]
#xValue = 0.1
counter = 0

xtol = 1e-3
ytol = 1e-3

for xValue in xValues:
    uCounter = 0
    fig, ax = plt.subplots()
    ax.plot(xForLine, yForLine,  'k--')
    for filename, colorType in zip(filenames2, colors):
        zAndU = loadVelocityBox(filename, xValue, xtol, ytol)#loadVelocity(filename)
        #zAndU.sort('z')
        if filename != filename3b:
            ax.plot(zAndU.Ux/U0[uCounter], zAndU.z, colorType, \
                    label='$Fr_h$ %.2f' %FrS[uCounter], \
                    linewidth=4.0)
        else:
            ax.plot(zAndU.Ux/U0[uCounter], zAndU.z, colorType, \
                label='Non-stratified', \
                linewidth=4.0)
        uCounter += 1
    ax.set_title('x position %.2f [m]' %xValue, fontsize = fontSize*1.25)
    ax.legend(fontsize = fontSize)
    ax.set_xlabel(r'$U_x/U_0$ ', fontsize = fontSize*1.25)
    ax.set_ylabel('z-position [m]', fontsize = fontSize*1.25)
    ax.tick_params(axis='both', which='major', labelsize=fontSize*1.25)
    xtickValues = np.arange(-1.2,0+.1, .4)
    ytickValues = np.arange(-.25,-0.1+0.05, 0.05)
    ax.set_xticks(xtickValues)
    ax.set_yticks(ytickValues)
    #ax5.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
    ax.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
    ax.grid(b=True, which='minor', color='grey', linewidth=0.2)
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.2))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.025))
    fig.savefig('/home/peterek/OpenFOAM/PETER-5.0/MT_PK/plots/3Dtop/VelocityProfileX%dIF15'\
                    %counter, bbox_inches = 'tight')
    counter += 1

#%%  Eta Interface 01
yWidth = 0.0
h1s = 0.15 #fresh water height
Us01 = np.asarray([0.1, 0.22])
Us02 = np.asarray([0.08, 0.1, 0.22])
Us = Us02
U2inCm2 = [(Us[i]*100)**2 for i in xrange(len(Us))] # U**2, in cm
c2twoPhaseIf = c2Emsmailpour(rho1, rho2, h1s)
FrS = Fr(Us, c2twoPhaseIf) 
    
    
filename0 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/zg/interface01U008.csv'
filename1 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/zg/interface01U01.csv'
#filename1 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/interface02U022.csv'
filename2 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/zg/interface01U022.csv'
filename3 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/interface02SingleU016.csv'


#filename1 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/interface01U01.csv'
#filename2 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/interface01U022.csv'
    
filenames = [filename0, filename1, filename2]#, filename3]
x1 = np.arange(0, 0.6, .1)
y1 = np.ones_like(x1)*(-.1)
y2 = np.arange(-0.1, 0+.01, .01)
x2 = np.ones_like(y2)*.6
x3 = x1[::-1]
y3 = np.zeros_like(y1)
y4 = y2[::-1]
x4 = np.zeros_like(y4)
x = np.concatenate((x1, x2, x3, x4))
y = np.concatenate((y1, y2, y3, y4))

xForLine = np.array((0, 0))
yForLine = np.array((-.1, -.25))

fig, ax = plt.subplots()
ax.plot(x, y, 'k', label = 'Barge')
ax.plot(xForLine, yForLine,  'k--')
markers = ['o', 'x', "^", 'P']
colors = ['b', 'g', 'r', 'c']
counter = 0
for filename, markerType, colorType in zip(filenames, markers, colors):
    points= loadPoints(filename, yWidth)
    if filename != filename3:
        ax.plot(points.x, points.z, markerType,  label='Fr %.2f' %FrS[counter])
    else:
        ax.plot(points.x, points.z, 'mP', label='Non-stratified')
    counter += 1
#plt.axis('equal')
ax.text(-0.9, -0.075, 'Stern', fontsize=20)
ax.set_xlim(-2, 1)
ax.set_ylim(-0.2, 0)
ax.set_ylabel(r'$z$ [m]', fontsize = fontSize*1.25)
ax.set_xlabel('x-position [m]', fontsize = fontSize*1.25)
#legends = ['Barge', 'Fr %.2f' %FrS[0], 'Fr %.2f' %FrS[1] ]
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 1.0, box.height*.8])
ax.legend(loc='center left', bbox_to_anchor=(0.0, -0.45), \
           fontsize = fontSize*1.25, ncol=2)
#ax.axis('equal')
#fig5.tight_layout()
ax.tick_params(axis='both', which='major', labelsize=fontSize*1.25)
xtickValues = np.arange(-2.0,1.0+.1, 1.0)
ytickValues = np.arange(-.2,-0.025+0.05, 0.05)
ax.set_xticks(xtickValues)
ax.set_yticks(ytickValues)
#ax5.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
ax.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
ax.grid(b=True, which='minor', color='grey', linewidth=0.2)
ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.25))
ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.025))
fig.savefig('/home/peterek/OpenFOAM/PETER-5.0/MT_PK/plots/3Dtop/Surface1'\
                , bbox_inches = 'tight')

#%% Velocity profiles interface 01
filename0b = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/data/zg/interface01U008box.csv'
filename1b = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/data/zg/interface01U01box.csv'
filename2b = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/data/zg/interface01U022box.csv'
filename3b = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/data/interface02SingleU016box.csv'
filenames2 = [filename0b, filename1b, filename2b, filename3b]
U0 = [0.08, 0.1, .22, 0.16]
xForLine = np.array((-1, -1))
yForLine = np.array((-.225, -.075))

colors = ['b', 'g', 'r', 'm']
xValues = [0.0, 0.1, 0.3, 0.4, 0.5, 0.6]
#xValue = 0.1
counter = 0

xtol = 1e-3
ytol = 1e-3

for xValue in xValues:
    uCounter = 0
    fig, ax = plt.subplots()
    ax.plot(xForLine, yForLine,  'k--')
    for filename, colorType in zip(filenames2, colors):
        zAndU = loadVelocityBox(filename, xValue, xtol, ytol)#loadVelocity(filename)
        #zAndU.sort('z')
        if filename != filename3b:
            ax.plot(zAndU.Ux/U0[uCounter], zAndU.z, colorType, \
                    label='Fr %.2f' %FrS[uCounter], \
                    linewidth=4.0)
        else:
            ax.plot(zAndU.Ux/U0[uCounter], zAndU.z, colorType, \
                label='Non-stratified', \
                linewidth=4.0)
        uCounter += 1
    ax.set_title('x position %.2f [m]' %xValue)
    ax.legend()
    ax.set_xlabel(r'$U_x/U_0$ ', fontsize = fontSize*1.25)
    ax.set_ylabel('z-position [m]', fontsize = fontSize)
    ax.tick_params(axis='both', which='major', labelsize=fontSize)
    xtickValues = np.arange(-1.2,0+.1, .4)
    ytickValues = np.arange(-.25,-0.1+0.05, 0.05)
    ax.set_xticks(xtickValues)
    ax.set_yticks(ytickValues)
    #ax5.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
    ax.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
    ax.grid(b=True, which='minor', color='grey', linewidth=0.2)
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.2))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.025))
    fig.savefig('/home/peterek/OpenFOAM/PETER-5.0/MT_PK/plots/3Dtop/VelocityProfileX%dIF1'\
                    %counter, bbox_inches = 'tight')
    counter += 1

#%% Pressure
'''
filename0b = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/data/interface02U008box.csv'
filename1b = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/data/interface02U016box.csv'
filename2b = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/data/interface02U022box.csv'
xValue = 0.0
xtol = 1e-3
ytol = 1e-3
data = loadPressure(filename0b, xValue, xtol, ytol)

x2,z2 = np.meshgrid(data.x,data.z, indexing="ij") # Grid for x- og y-verdiene (km)
p_rgh = np.zeros_like(x2)
xDim, yDim = np.shape(p_rgh)
for i in xrange(xDim):
    for j in xrange(yDim):
        xVal = x2[i,j]
        zVal = z2[i,j]  
        xIndexes = data.x[data.x == xVal]
        xIndexes = xIndexes.index
        zIndexes = data.z[data.z == zVal]
        zIndexes = zIndexes.index
        wantedIndexes = xIndexes
        wantedIndexes = wantedIndexes.intersection(zIndexes)
        if len(wantedIndexes) == 0:
            p_rgh[i,j] = np.nan
        else:
            p_rgh[i,j] = np.asscalar(data.p_rgh[wantedIndexes].values)
#%%
B = np.ma.masked_where(np.isnan(p_rgh),p_rgh)
cmap = plt.cm.get_cmap("winter")
cmap.set_under("magenta")
cmap.set_over("yellow")
cmap.set_bad("green")
#plt.pcolor(C)
C=plt.contourf(x2,z2, B, cmap=cmap)

plt.colorbar(C)
plt.show()
#p = p0 - dp/(1+(x**2+y**2)/R**2)     # Beregn trykket p (hPa)
#C=plt.pcolor(x,y,p)
#plt.colorbar(C)
#plt.axis("equal")
#plt.xlabel("x") # Sett aksenavn
#plt.ylabel("y")
#plt.show()



#%%

R = 50    # Utstrekningen av lavtrykket (km)
p0 = 1000 # Lufttrykket langt borte fra sentrum (hPa)
dp = 40   # Trykkfallet inn mot sentrum (hPa)
tx = np.linspace(-150,1,151)
ty = np.linspace(0,150,151)
x,y = np.meshgrid(tx,ty, indexing="ij") # Grid for x- og y-verdiene (km)
p = p0 - dp/(1+(x**2+y**2)/R**2)     # Beregn trykket p (hPa)
C=plt.pcolor(x,y,p)
plt.colorbar(C)
plt.axis("equal")
plt.xlabel("x") # Sett aksenavn
plt.ylabel("y")
plt.show()

#%%
array = np.random.rand(4,10)
array[array<0.5]=np.nan
m = np.ma.masked_where(np.isnan(array),array)
plt.pcolor(m)
plt.colorbar(orientation='horizontal')       
plt.show()
'''

#%%  Eta Interface 02 Width side
yWidth= 0.225
h1s = 0.2 #fresh water height
Us01 = np.asarray([0.1, 0.22])
Us02 = np.asarray([0.08, 0.16, 0.22])
Us = Us02
U2inCm2 = [(Us[i]*100)**2 for i in xrange(len(Us))] # U**2, in cm
c2twoPhaseIf = c2Emsmailpour(rho1, rho2, h1s)
FrS = Fr(Us, c2twoPhaseIf) 

filename0 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/zg/interface02U008.csv'
filename1 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/zg/interface02U016.csv'
#filename1 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/interface02U022.csv'
filename2 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/zg/interface02U022.csv'
filename3 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/interface02SingleU016.csv'


#filename1 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/interface01U01.csv'
#filename2 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/interface01U022.csv'
    
filenames = [filename0, filename1, filename2]#, filename3]
x1 = np.arange(0, 0.6, .1)
y1 = np.ones_like(x1)*(-.1)
y2 = np.arange(-0.1, 0+.01, .01)
x2 = np.ones_like(y2)*.6
x3 = x1[::-1]
y3 = np.zeros_like(y1)
y4 = y2[::-1]
x4 = np.zeros_like(y4)
x = np.concatenate((x1, x2, x3, x4))
y = np.concatenate((y1, y2, y3, y4))

xForLine = np.array((0, 0))
yForLine = np.array((-.1, -.25))


fig, ax = plt.subplots()
ax.plot(x, y, 'k', label = 'Barge')
ax.plot(xForLine, yForLine,  'k--')
markers = ['o', 'x', "^", 'P']
colors = ['b', 'g', 'r', 'c']
counter = 0
for filename, markerType, colorType in zip(filenames, markers, colors):
    points= loadPoints(filename, yWidth)
    if filename != filename3:
        ax.plot(points.x, points.z, markerType,  label='Fr %.2f' %FrS[counter])
    else:
        ax.plot(points.x, points.z, 'mP', label='Non-stratified')
    counter += 1
#plt.axis('equal')
ax.text(-0.8, -0.095, 'Stern', fontsize=20)
ax.set_xlim(-1, 1.5)
ax.set_ylim(-0.25, -0.075)
ax.set_ylabel(r'$z$ [m]', fontsize = fontSize*1.25)
ax.set_xlabel('x-position [m]', fontsize = fontSize)
#legends = ['Barge', 'Fr %.2f' %FrS[0], 'Fr %.2f' %FrS[1] ]
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 1.0, box.height*.8])
ax.legend(loc='center left', bbox_to_anchor=(0.0, -0.45), \
           fontsize = fontSize, ncol=2)
#ax.axis('equal')
#fig5.tight_layout()
ax.tick_params(axis='both', which='major', labelsize=fontSize)
xtickValues = np.arange(-2.0,1.0+.1, 1.0)
ytickValues = np.arange(-.2,-0.125+0.05, 0.05)
ax.set_xticks(xtickValues)
ax.set_yticks(ytickValues)
#ax5.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
ax.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
ax.grid(b=True, which='minor', color='grey', linewidth=0.2)
ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.25))
ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.025))
fig.savefig('/home/peterek/OpenFOAM/PETER-5.0/MT_PK/plots/3Dtop/SurfaceWidth'\
                , bbox_inches = 'tight')


    
#%% Velocity profiles interface 02
filename0b = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/data/zg/interface02U008box.csv'
filename1b = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/data/zg/interface02U016box.csv'
filename2b = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/data/zg/interface02U022box.csv'
filename3b = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/data/interface02SingleU016box.csv'
filenames2 = [filename0b, filename1b, filename2b, filename3b]
U0 = [0.08, 0.16, .22, 0.16]
xForLine = np.array((-1, -1))
yForLine = np.array((-.225, -.075))

colors = ['b', 'g', 'r', 'm']
xValues = [0.0, 0.1, 0.3, 0.4, 0.5, 0.6]
#xValue = 0.1
counter = 0

xtol = 1e-3
ytol = 1e-3

for xValue in xValues:
    uCounter = 0
    fig, ax = plt.subplots()
    ax.plot(xForLine, yForLine,  'k--')
    for filename, colorType in zip(filenames2, colors):
        zAndU = loadVelocityBox(filename, xValue, xtol, ytol)#loadVelocity(filename)
        #zAndU.sort('z')
        if filename != filename3b:
            ax.plot(zAndU.Ux/U0[uCounter], zAndU.z, colorType, \
                    label='Fr %.2f' %FrS[uCounter], \
                    linewidth=4.0)
        else:
            ax.plot(zAndU.Ux/U0[uCounter], zAndU.z, colorType, \
                label='Non-stratified', \
                linewidth=4.0)
        uCounter += 1
    ax.set_title('x position %.2f [m]' %xValue)
    ax.legend()
    ax.set_xlabel(r'$U_x/U_0$ ', fontsize = fontSize*1.25)
    ax.set_ylabel('z-position [m]', fontsize = fontSize)
    ax.tick_params(axis='both', which='major', labelsize=fontSize)
    xtickValues = np.arange(-1.2,0+.1, .4)
    ytickValues = np.arange(-.25,-0.1+0.05, 0.05)
    ax.set_xticks(xtickValues)
    ax.set_yticks(ytickValues)
    #ax5.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
    ax.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
    ax.grid(b=True, which='minor', color='grey', linewidth=0.2)
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.2))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.025))
    fig.savefig('/home/peterek/OpenFOAM/PETER-5.0/MT_PK/plots/3Dtop/VelocityProfileX%dWidth'\
                    %counter, bbox_inches = 'tight')
    counter += 1
    
#%%  Eta Interface 015 Width side
yWidth = 0.225
h1s = 0.15 #fresh water height
Us01 = np.asarray([0.1, 0.22])
Us02 = np.asarray([0.08, 0.12, 0.22])
Us = Us02
U2inCm2 = [(Us[i]*100)**2 for i in xrange(len(Us))] # U**2, in cm
c2twoPhaseIf = c2Emsmailpour(rho1, rho2, h1s)
FrS = Fr(Us, c2twoPhaseIf) 
    
    
filename0 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/zg/interface015U008.csv'
filename1 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/zg/interface015U012.csv'
#filename1 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/interface02U022.csv'
filename2 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/zg/interface015U022.csv'
filename3 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/interface02SingleU016.csv'


#filename1 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/interface01U01.csv'
#filename2 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/interface01U022.csv'
    
filenames = [filename0, filename1, filename2]#, filename3]
x1 = np.arange(0, 0.6, .1)
y1 = np.ones_like(x1)*(-.1)
y2 = np.arange(-0.1, 0+.01, .01)
x2 = np.ones_like(y2)*.6
x3 = x1[::-1]
y3 = np.zeros_like(y1)
y4 = y2[::-1]
x4 = np.zeros_like(y4)
x = np.concatenate((x1, x2, x3, x4))
y = np.concatenate((y1, y2, y3, y4))

xForLine = np.array((0, 0))
yForLine = np.array((-.1, -.25))

fig, ax = plt.subplots()
ax.plot(x, y, 'k', label = 'Barge')
ax.plot(xForLine, yForLine,  'k--')
markers = ['o', 'x', "^", 'P']
colors = ['b', 'g', 'r', 'c']
counter = 0
for filename, markerType, colorType in zip(filenames, markers, colors):
    points= loadPoints(filename, yWidth)
    if filename != filename3:
        ax.plot(points.x, points.z, markerType,  label='Fr %.2f' %FrS[counter])
    else:
        ax.plot(points.x, points.z, 'mP', label='Non-stratified')
    counter += 1
#plt.axis('equal')
ax.text(-0.8, -0.095, 'Stern', fontsize=20)
ax.set_xlim(-1, 1.5)
ax.set_ylim(-0.25, -0.075)
ax.set_ylabel(r'$z$ [m]', fontsize = fontSize*1.25)
ax.set_xlabel('x-position [m]', fontsize = fontSize)
#legends = ['Barge', 'Fr %.2f' %FrS[0], 'Fr %.2f' %FrS[1] ]
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 1.0, box.height*.8])
ax.legend(loc='center left', bbox_to_anchor=(0.0, -0.45), \
           fontsize = fontSize, ncol=2)
#ax.axis('equal')
#fig5.tight_layout()
ax.tick_params(axis='both', which='major', labelsize=fontSize)
xtickValues = np.arange(-2.0,1.0+.1, 1.0)
ytickValues = np.arange(-.2,-0.125+0.05, 0.05)
ax.set_xticks(xtickValues)
ax.set_yticks(ytickValues)
#ax5.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
ax.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
ax.grid(b=True, which='minor', color='grey', linewidth=0.2)
ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.25))
ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.025))
fig.savefig('/home/peterek/OpenFOAM/PETER-5.0/MT_PK/plots/3Dtop/Surface15Width'\
                , bbox_inches = 'tight')

#%% Velocity profiles interface 015 Width side
filename0b = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/data/zg/interface015U008box.csv'
filename1b = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/data/zg/interface015U012box.csv'
filename2b = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/data/zg/interface015U022box.csv'
filename3b = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/data/interface02SingleU016box.csv'
filenames2 = [filename0b, filename1b, filename2b, filename3b]
U0 = [0.08, 0.12, .22, 0.16]
xForLine = np.array((-1, -1))
yForLine = np.array((-.225, -.075))

colors = ['b', 'g', 'r', 'm']
xValues = [0.0, 0.1, 0.3, 0.4, 0.5, 0.6]
#xValue = 0.1
counter = 0

xtol = 1e-3
ytol = 1e-3

for xValue in xValues:
    uCounter = 0
    fig, ax = plt.subplots()
    ax.plot(xForLine, yForLine,  'k--')
    for filename, colorType in zip(filenames2, colors):
        zAndU = loadVelocityBox(filename, xValue, xtol, ytol)#loadVelocity(filename)
        #zAndU.sort('z')
        if filename != filename3b:
            ax.plot(zAndU.Ux/U0[uCounter], zAndU.z, colorType, \
                    label='Fr %.2f' %FrS[uCounter], \
                    linewidth=4.0)
        else:
            ax.plot(zAndU.Ux/U0[uCounter], zAndU.z, colorType, \
                label='Non-stratified', \
                linewidth=4.0)
        uCounter += 1
    ax.set_title('x position %.2f [m]' %xValue)
    ax.legend()
    ax.set_xlabel(r'$U_x/U_0$ ', fontsize = fontSize*1.25)
    ax.set_ylabel('z-position [m]', fontsize = fontSize)
    ax.tick_params(axis='both', which='major', labelsize=fontSize)
    xtickValues = np.arange(-1.2,0+.1, .4)
    ytickValues = np.arange(-.25,-0.1+0.05, 0.05)
    ax.set_xticks(xtickValues)
    ax.set_yticks(ytickValues)
    #ax5.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
    ax.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
    ax.grid(b=True, which='minor', color='grey', linewidth=0.2)
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.2))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.025))
    fig.savefig('/home/peterek/OpenFOAM/PETER-5.0/MT_PK/plots/3Dtop/VelocityProfileX%dIF15Width'\
                    %counter, bbox_inches = 'tight')
    counter += 1

#%%  Eta Interface 01 Width side
eps = 0.03
yWidth =  .225/2
h1s = 0.15 #fresh water height
Us01 = np.asarray([0.1, 0.22])
Us02 = np.asarray([0.08, 0.1, 0.22])
Us = Us02
U2inCm2 = [(Us[i]*100)**2 for i in xrange(len(Us))] # U**2, in cm
c2twoPhaseIf = c2Emsmailpour(rho1, rho2, h1s)
FrS = Fr(Us, c2twoPhaseIf) 
    
    
filename0 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/zg/interface01U008.csv'
filename1 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/zg/interface01U01.csv'
#filename1 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/interface02U022.csv'
filename2 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/zg/interface01U022.csv'
filename3 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/interface02SingleU016.csv'


#filename1 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/interface01U01.csv'
#filename2 = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/eta/interface01U022.csv'
    
filenames = [filename0, filename1, filename2]#, filename3]
x1 = np.arange(0, 0.6, .1)
y1 = np.ones_like(x1)*(-.1)
y2 = np.arange(-0.1, 0+.01, .01)
x2 = np.ones_like(y2)*.6
x3 = x1[::-1]
y3 = np.zeros_like(y1)
y4 = y2[::-1]
x4 = np.zeros_like(y4)
x = np.concatenate((x1, x2, x3, x4))
y = np.concatenate((y1, y2, y3, y4))

xForLine = np.array((0, 0))
yForLine = np.array((-.1, -.25))

fig, ax = plt.subplots()
ax.plot(x, y, 'k', label = 'Barge')
ax.plot(xForLine, yForLine,  'k--')
markers = ['o', 'x', "^", 'P']
colors = ['b', 'g', 'r', 'c']
counter = 0
for filename, markerType, colorType in zip(filenames, markers, colors):
    points= loadPoints(filename, yWidth)
    if filename != filename3:
        ax.plot(points.x, points.z, markerType,  label='Fr %.2f' %FrS[counter])
    else:
        ax.plot(points.x, points.z, 'mP', label='Non-stratified')
    counter += 1
#plt.axis('equal')
ax.text(-0.9, -0.075, 'Stern', fontsize=20)
ax.set_xlim(-2, 1)
ax.set_ylim(-0.2, 0)
ax.set_ylabel(r'$z$ [m]', fontsize = fontSize*1.25)
ax.set_xlabel('x-position [m]', fontsize = fontSize)
#legends = ['Barge', 'Fr %.2f' %FrS[0], 'Fr %.2f' %FrS[1] ]
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 1.0, box.height*.8])
ax.legend(loc='center left', bbox_to_anchor=(0.0, -0.45), \
           fontsize = fontSize, ncol=2)
#ax.axis('equal')
#fig5.tight_layout()
ax.tick_params(axis='both', which='major', labelsize=fontSize)
xtickValues = np.arange(-2.0,1.0+.1, 1.0)
ytickValues = np.arange(-.2,-0.025+0.05, 0.05)
ax.set_xticks(xtickValues)
ax.set_yticks(ytickValues)
#ax5.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
ax.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
ax.grid(b=True, which='minor', color='grey', linewidth=0.2)
ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.25))
ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.025))
fig.savefig('/home/peterek/OpenFOAM/PETER-5.0/MT_PK/plots/3Dtop/Surface1Width'\
                , bbox_inches = 'tight')

#%% Velocity profiles interface 01 Width side
filename0b = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/data/zg/interface01U008box.csv'
filename1b = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/data/zg/interface01U01box.csv'
filename2b = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/data/zg/interface01U022box.csv'
filename3b = '/home/peterek/OpenFOAM/PETER-5.0/MT_PK/python/data/interface02SingleU016box.csv'
filenames2 = [filename0b, filename1b, filename2b, filename3b]
U0 = [0.08, 0.1, .22, 0.16]
xForLine = np.array((-1, -1))
yForLine = np.array((-.225, -.075))

colors = ['b', 'g', 'r', 'm']
xValues = [0.0, 0.1, 0.3, 0.4, 0.5, 0.6]
#xValue = 0.1
counter = 0

xtol = 1e-3
ytol = 1e-3

for xValue in xValues:
    uCounter = 0
    fig, ax = plt.subplots()
    ax.plot(xForLine, yForLine,  'k--')
    for filename, colorType in zip(filenames2, colors):
        zAndU = loadVelocityBox(filename, xValue, xtol, ytol)#loadVelocity(filename)
        #zAndU.sort('z')
        if filename != filename3b:
            ax.plot(zAndU.Ux/U0[uCounter], zAndU.z, colorType, \
                    label='Fr %.2f' %FrS[uCounter], \
                    linewidth=4.0)
        else:
            ax.plot(zAndU.Ux/U0[uCounter], zAndU.z, colorType, \
                label='Non-stratified', \
                linewidth=4.0)
        uCounter += 1
    ax.set_title('x position %.2f [m]' %xValue)
    ax.legend()
    ax.set_xlabel(r'$U_x/U_0$ ', fontsize = fontSize*1.25)
    ax.set_ylabel('z-position [m]', fontsize = fontSize)
    ax.tick_params(axis='both', which='major', labelsize=fontSize)
    xtickValues = np.arange(-1.2,0+.1, .4)
    ytickValues = np.arange(-.25,-0.1+0.05, 0.05)
    ax.set_xticks(xtickValues)
    ax.set_yticks(ytickValues)
    #ax5.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
    ax.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
    ax.grid(b=True, which='minor', color='grey', linewidth=0.2)
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.2))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.025))
    fig.savefig('/home/peterek/OpenFOAM/PETER-5.0/MT_PK/plots/3Dtop/VelocityProfileX%dIF1Width'\
                    %counter, bbox_inches = 'tight')
    counter += 1