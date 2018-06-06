
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
import matplotlib.ticker as plticker
import matplotlib as mpl

# the mock-0.3.1 dir contains testcase.py, testutils.py & mock.py

#%% Input for files
nutType = 'SpaldingRefType2kLowRe'#'LowYplus'
topOrNot = '3Dtop'

#%% Functions
def runMultipleCases(velocity, barge, phase):
    if phase == 'freshSaltKOmegaSST': 
        call(["./Allclean"])
        sys.path.append(os.getcwd()+'/'+ phase)
        import generateUfile
        import generateControlDict
        generateUfile.writeUfile(velocity, barge, phase)
        generateControlDict.writeControlDictFile(velocity)
        call(["./Allrun", str(velocity).replace(".", "")])

    elif phase == 'singleKOmegaSST':
        call(["./AllcleanSingle"])
        sys.path.append(os.getcwd()+'/'+ phase)
        import generateUfile
        import generateControlDict
        generateUfile.writeUfile(velocity, barge, phase)
        generateControlDict.writeControlDictFile(velocity)
        call(["./Allrun", str(velocity).replace(".", "")])


F = lambda U, Cd, rho, S, N2gfCoeffcient=0.00980665: .5*rho*U**2*Cd*S/N2gfCoeffcient
c2 = lambda rho1, rho2, h1, h2, g=9.81: g*(rho2 - rho1)/rho2*h1*h2/(h1 + h2) #Maximume wave speed squared
c2Hester = lambda rho1, rho2, h1, h2, g=9.81: (rho2 - rho1)*g/(rho2/h2 + rho1/h1)
c2Grue = lambda rho1, rho2, h1, h2, g=9.81: (rho2 - rho1)*g/rho2 * h1*h2/(h1 + h2)
c2Emsmailpour = lambda rho1, rho2, h1, g=9.81: g*h1*(rho2 - rho1)/rho1
Fr = lambda U, c2: U/c2**.5
FrSingle = lambda U: U/np.sqrt(9.81*0.6)

def load_force_coeffs(filename): 
    """Load force coefficients into a DataFrame for the given time directory."""
    data = np.loadtxt(filename, skiprows=9)
    df = pd.DataFrame()
    df["time"] = data[:, 0]
    df["cl"] = data[:, 3]
    df["cd"] = data[:, 2]
    df["cm"] = data[:, 1]
    return df.tail(1).cd

def load_force_coeffsAverage(filename): 
    """Load force coefficients into a DataFrame for the given time directory."""
    data = np.loadtxt(filename, skiprows=9)
    df = pd.DataFrame()
    df["time"] = data[:, 0]
    df["cl"] = data[:, 3]
    df["cd"] = data[:, 2]
    df["cm"] = data[:, 1]

    #t0 = df.time.max() - 1000
    #df = df[df.time >= t0]
    #df["cd"] = np.average(df.cd.values)    
    wantedIndex  = df.time[df.time == 40].index[0]
    #wantedCd = df["cd"][wantedIndex]
    #return wantedCd
    df2 = df[wantedIndex-1:wantedIndex+1]
    return df2.tail(1).cd #df.mean() 





#%% Gou et al esimtated
filenumbersGou = ["001", "002", "003", "004", "005", "006", "007", "008", "009"]
numberOfFilesGou = len(filenumbersGou)
UvaluesGou = [0.24, .2, .18, .16, .14, .12, .10, .08, .06]
U2inCm2Gou = [(UvaluesGou[i]*100)**2 for i in xrange(len(UvaluesGou))] # U**2, in cm
Gou = [.4952 + .2501*U2inCm2Gou[i] for i in xrange(numberOfFilesGou)]

from scipy.optimize import curve_fit

def func(x, a, b):
    return a + b*x

xdata = np.array((U2inCm2Gou[0], U2inCm2Gou[1], U2inCm2Gou[4], U2inCm2Gou[5],
              U2inCm2Gou[6], U2inCm2Gou[7], U2inCm2Gou[8]))
ydata = np.array((forces[0], forces[1], forces[2], forces[3], forces[4],
              forces[5], forcesNoSlipWidth[3]))
ydata = np.squeeze(ydata)

popt, pcov = curve_fit(func, xdata, ydata)

#%% Gou Testing Force calculations in Gou et al, Figures 5 and 6
F = lambda U, Cd, rho, S, N2gfCoeffcient=0.00980665: .5*rho*U**2*Cd*S/N2gfCoeffcient
CdGou = np.array((.025, .05, .085, .08, .075, .07, .06, .04))
FGou = F(Us[:-2], CdGou, rho1, S2)
plt.plot(Us[:-2]*100, FGou)
plt.ylim(0,200)


#%% Gou Plot, same as Gou et al figure 2
fig, ax = plt.subplots()
ax.plot(U2inCm2, forces, 'o')
ax.plot(U2inCm2NoSlipWidth, forcesNoSlipWidth, 'x')
ax.plot(U2inCm2Gou, Gou, 'b-')
ax.plot(xdata, func(xdata, *popt), 'r-')#label='fit: a=%5.3f, b=%5.3f' % tuple(popt))
ax.set_xlabel(r'$U^2 (\frac{cm^2}{s^2})$')
ax.set_ylabel(r'$F_d(gf)$')
ax.set_title('Non-stratified')
ax.legend(['CFD slip width', 'CFD no-slip width', 'Gou et al', 'CFD Curvefit'], loc=0) 
plt.show()
fig.tight_layout()
outfilename = "singplePhase.png"
fig.savefig(outfilename)
#plt.close()

#%% Running
d = 0.1
h1s = [0.1, 0.15, 0.2] #fresh water height
dtoh1squared = [(d/h1)**2 for h1 in h1s]
h2 = 3.8# 4.#16 #salt "
#us1 = np.linspace(.06, 0.24, 10)
us1 = np.linspace(.08, 0.22, 8)
#us2 = np.linspace(0.3, .6, 4)
Us = us1 #np.concatenate((us1, us2))
numberOfUs = len(Us)
rho1 = 997 # densisty fresh
rho2 = 1024 # salt
g = 9.81

U2inCm2 = [(Us[i]*100)**2 for i in xrange(len(Us))] # U**2, in cm
#c2twoPhaseIf = c2(rho1, rho2, h1, h2)
c2twoPhaseIf = c2Emsmailpour(rho1, rho2, h1)
FrS = Fr(Us, c2twoPhaseIf) 
rho = 997




barges = [0.10, 0.15, 0.20]
phases = ['freshSaltKOmegaSST','singleKOmegaSST' ]
numberOfBarges = len(barges)



resultsDict = OrderedDict()
resultsDictDiff = OrderedDict()

us2 = ['008', '01', '012', '014', '016', '018','02', '022']#, \
#us2 = ['006', '008', '01', '012', '014','016', '018','02', '022','024']#, \
       #'03','04', '05','06']

CdArrayBarge01Single = []
FArrayBarge01Single = []
CdArrayBarge01Two = []
FArrayBarge01Two = []

CdArrayBarge015Single = []
FArrayBarge015Single = []
CdArrayBarge015Two = []
FArrayBarge015Two = []

CdArrayBarge02Single = []
FArrayBarge02Single = []
CdArrayBarge02Two = []
FArrayBarge02Two = []


CdDiffBarge01 = []
FDiffBarge01 = []
CdDiffBarge015 = []
FDiffBarge015 = []
CdDiffBarge02 = []
FDiffBarge02 = []

for barge in barges:
    S = 0.6 * 0.225 + 0.6*barge + 2 *  0.225 * barge# Halv sufrace (Used in openFoam)
    #S2 = 0.6 * 0.45 + 2*0.6*barge + 2*0.45*barge# Full surface, wetted
    S2 = 0.6 * 0.45 + 2*0.6*barges[0] + 2*0.45*barges[0]# Full surface, wetted
    resultsDict[barge] = OrderedDict()
    resultsDictDiff[barge] = OrderedDict()
    counter = 0
    for U in Us:
        resultsDict[barge][U] = OrderedDict()
        resultsDictDiff[barge][U] = OrderedDict()
        for phase in phases:
            resultsDict[barge][U][phase] = OrderedDict()
            resultsDictDiff[barge][U] = OrderedDict()
        
            #runMultipleCases(U, barge, phase)
            if phase == 'singleKOmegaSST':
                filename = '/home/karl/OpenFOAM/karl-5.0/run/MT/MT_PK/twoPhaseIF2/freshSaltMixingSingle/forceCoeffs/forceCoeffVelocity%s' %(us2[counter])
            else:    
                if str(barge) == '0.1':
                    filename = '/home/karl/OpenFOAM/karl-5.0/run/MT/MT_PK/twoPhaseIF2/freshSaltMixing/forceCoeffsInterface01/forceCoeffVelocity%s' %us2[counter]
                elif str(barge) == '0.15':
                    if phase == 'freshSaltKOmegaSST':
                        filename = '/home/karl/OpenFOAM/karl-5.0/run/MT/MT_PK/twoPhaseIF2/freshSaltMixing/forceCoeffsInterface015/forceCoeffVelocity%s' %us2[counter]
                    #elif phase == 'singleKOmegaSST':
                    #    filename = '/home/karl/OpenFOAM/karl-5.0/run/MT/MT_PK/twoPhaseIF2/%sAltBc3smallMesh/forceCoeffsSpaldingInterface010/forceCoeffVelocity%s' %(phase, us2[counter])
                elif str(barge) == '0.2':
                    if phase == 'freshSaltKOmegaSST':
                        filename = '/home/karl/OpenFOAM/karl-5.0/run/MT/MT_PK/twoPhaseIF2/freshSaltMixing/forceCoeffsInterface02/forceCoeffVelocity%s' %us2[counter]
                    #elif phase == 'singleKOmegaSST':
                     #   filename = '/home/karl/OpenFOAM/karl-5.0/run/MT/MT_PK/twoPhaseIF2/%sAltBc3smallMesh/forceCoeffsSpaldingInterface010/forceCoeffVelocity%s' %(phase, us2[counter])
                
            
            #print filename
            if phase == 'freshSaltKOmegaSST':
                resultsDict[barge][U][phase]['Cd'] = load_force_coeffsAverage(filename)#load_force_coeffs(filename)
            else: 
                resultsDict[barge][U][phase]['Cd'] = load_force_coeffs(filename)#load_force_coeffs(filename)
            if str(barge) == '0.1' and phase == 'singleKOmegaSST':
                CdArrayBarge01Single.append(np.asscalar(resultsDict[barge][U][phase]['Cd'].values))
                FArrayBarge01Single.append(F(U, resultsDict[barge][U][phase]['Cd'].values, rho=rho1, S=S2))
            elif str(barge) == '0.1' and phase == 'freshSaltKOmegaSST':
                CdArrayBarge01Two.append(resultsDict[barge][U][phase]['Cd'].values)
                FArrayBarge01Two.append(F(U, resultsDict[barge][U][phase]['Cd'].values, rho=rho1, S=S2))
            elif str(barge) == '0.15' and phase == 'singleKOmegaSST':
                CdArrayBarge015Single.append(resultsDict[barge][U][phase]['Cd'].values)
                FArrayBarge015Single.append(F(U, resultsDict[barge][U][phase]['Cd'].values, rho=rho1, S=S2))
                #CdArrayBarge02Single.append(CdArrayBarge01Single[-1])
                #FArrayBarge02Single.append(FArrayBarge01Single[-1])
            elif str(barge) == '0.15' and phase == 'freshSaltKOmegaSST':
                CdArrayBarge015Two.append(resultsDict[barge][U][phase]['Cd'].values)
                FArrayBarge015Two.append(F(U, resultsDict[barge][U][phase]['Cd'].values, rho=rho1, S=S2))
            elif str(barge) == '0.2' and phase == 'singleKOmegaSST':
                CdArrayBarge02Single.append(resultsDict[barge][U][phase]['Cd'].values)
                FArrayBarge02Single.append(F(U, resultsDict[barge][U][phase]['Cd'].values, rho=rho1, S=S2))
                #CdArrayBarge02Single.append(CdArrayBarge01Single[-1])
                #FArrayBarge02Single.append(FArrayBarge01Single[-1])
            elif str(barge) == '0.2' and phase == 'freshSaltKOmegaSST':
                CdArrayBarge02Two.append(resultsDict[barge][U][phase]['Cd'].values)
                FArrayBarge02Two.append(F(U, resultsDict[barge][U][phase]['Cd'].values, rho=rho1, S=S2))
            
            resultsDict[barge][U][phase]['F'] = F(U, resultsDict[barge][U][phase]['Cd'], rho=rho1, S=S2)
        resultsDictDiff[barge][U]['Cd'] =  resultsDict[barge][U]['freshSaltKOmegaSST']['Cd'].values - resultsDict[barge][U]['singleKOmegaSST']['Cd'].values
        resultsDictDiff[barge][U]['F'] =  resultsDict[barge][U]['freshSaltKOmegaSST']['F'].values - resultsDict[barge][U]['singleKOmegaSST']['F'].values
        if str(barge) == '0.1':
            CdDiffBarge01.append(resultsDictDiff[barge][U]['Cd'])
            FDiffBarge01.append(resultsDictDiff[barge][U]['F'])
        elif str(barge) == '0.15':
            CdDiffBarge015.append(resultsDictDiff[barge][U]['Cd'])
            FDiffBarge015.append(resultsDictDiff[barge][U]['F'])
        elif str(barge) == '0.2':
            CdDiffBarge02.append(resultsDictDiff[barge][U]['Cd'])
            FDiffBarge02.append(resultsDictDiff[barge][U]['F'])
        
        counter += 1

CdArrayBarge01Single = np.asarray(CdArrayBarge01Single)
FArrayBarge01Single = np.asarray(FArrayBarge01Single)
CdArrayBarge01Two = np.asarray(CdArrayBarge01Two)
FArrayBarge01Two = np.asarray(FArrayBarge01Two)

CdArrayBarge015Single = np.asarray(CdArrayBarge015Single)
FArrayBarge015Single = np.asarray(FArrayBarge015Single)
CdArrayBarge015Two = np.asarray(CdArrayBarge015Two)
FArrayBarge015Two = np.asarray(FArrayBarge015Two)

CdArrayBarge02Single = np.asarray(CdArrayBarge02Single)
FArrayBarge02Single = np.asarray(FArrayBarge02Single)
CdArrayBarge02Two = np.asarray(CdArrayBarge02Two)
FArrayBarge02Two = np.asarray(FArrayBarge02Two)


CdDiffBarge01 = np.asarray(CdDiffBarge01)
FDiffBarge01 = np.asarray(FDiffBarge01)
CdDiffBarge015 = np.asarray(CdDiffBarge015)
FDiffBarge015 = np.asarray(FDiffBarge015)
CdDiffBarge02 = np.asarray(CdDiffBarge02)
FDiffBarge02 = np.asarray(FDiffBarge02)



#%% Teasting different autors formulae
c2GouV = c2(rho1, rho2, h1, h2)
c2HesterV = c2Hester(rho1, rho2, h1, h2)
c2GrueV = c2Grue(rho1, rho2, h1, h2) 
c2EmsmailpourV = c2Emsmailpour(rho1, rho2, h1)
c2s = [c2GouV, c2HesterV, c2GrueV, c2EmsmailpourV]

#%% Gou: Check if Gou et al get constant Cd for singlePhase, as Esmailpour
U2inCm2 = [(Us[i]*100)**2 for i in xrange(len(Us))] # U**2, in cm
GouEsimator = np.asarray([.4952 + .2501*U2inCm2[i] for i in xrange(numberOfUs)])
N2gfCoeffcient=0.00980665
GouEsimator *= N2gfCoeffcient
GouCdSingle = GouEsimator/(.5*rho1*Us**2*S2)
GouCdAdd = np.asarray([.02, .04, .075, .07, .065, .055, .045, .03, .02, .02])
gouCdTwo = GouCdSingle + GouCdAdd

fig0, ax0 = plt.subplots()
ax0.plot(FrS, GouCdSingle, 'bo', label='Single')
ax0.plot(FrS, gouCdTwo, 'x', label='Two')
ax0.set_ylabel(r'$C_d$')
ax0.set_xlabel(r'$Fr $')
ax0.set_ylim(0, 0.3)
ax0.legend()    
ax0.set_title('Cd single Gou')
fig0.tight_layout()
#fig4.savefig('/home/karl/OpenFOAM/karl-5.0/run/MT/MT_PK/plots/Gou2KO')

  

#%% PLOTTING
fontSize = 20



# Comparison with Gou et al's estimation results, Figure 2
U2inCm2 = [(Us[i]*100)**2 for i in xrange(len(Us))] # U**2, in cm
GouEsimator = [.4952 + .2501*U2inCm2[i] for i in xrange(numberOfUs)]

fig4, ax4 = plt.subplots()
ax4.plot(U2inCm2, GouEsimator, '-r', label='Gou et al.')
ax4.plot(U2inCm2, FArrayBarge02Single, 'bo', label='Current study')

ax4.set_ylabel(r'$F [gf]$', fontsize = fontSize, rotation=90)
ax4.set_xlabel(r'$u^2 [cm^2/s^2]$', fontsize = fontSize)
#ax4.set_ylim(0, 160)
ax4.legend(fontsize = fontSize)    
#ax4.set_title('Gou Figure 2. %s. %s' %(nutType, topOrNot))
#fig4.tight_layout()
ax4.tick_params(axis='both', which='major', labelsize=fontSize)
xtickValues = np.arange(0,500+.1, 100)
ytickValues = np.arange(0,120+.1, 20)
ax4.set_xticks(xtickValues)
ax4.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
ax4.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
ax4.grid(b=True, which='minor', color='grey', linewidth=0.2)
ax4.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(50))
ax4.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(10))
fig4.savefig('/home/karl/OpenFOAM/karl-5.0/run/MT/MT_PK/plots/%s/Gou2KO%s' \
             %(topOrNot, nutType), bbox_inches = 'tight')
    
#Plots as in Gou et al figures 5
fig, ax = plt.subplots()
ax.plot(Us, FArrayBarge02Two, 'o', Us, FArrayBarge02Single, 'x' ,  Us, FDiffBarge02, '^')
ax.set_ylabel(r'$F [gf]$', fontsize = fontSize, rotation=90)
ax.set_xlabel('U', fontsize = fontSize)
#ax.set_ylim(0, 200)
ax.set_title('Pycnocline 0.2m', \
             fontsize = fontSize)
legends = (['Stratified', 'Non-stratified',  'Difference'])
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(legends, loc='center left', bbox_to_anchor=(1, 0.5),\
          fontsize = fontSize)
#fig.tight_layout()
ax.tick_params(axis='both', which='major', labelsize=fontSize)
fig.savefig('/home/karl/OpenFOAM/karl-5.0/run/MT/MT_PK/plots/%s/Gou5KO%s'\
            %(topOrNot, nutType), bbox_inches = 'tight')

#Plots as in Gou et al figures 5b
fig, ax = plt.subplots()
ax.plot(Us, FArrayBarge015Single, 'x' , Us, FArrayBarge015Two, 'o', Us, FDiffBarge015, '^')
ax.set_ylabel(r'$F [gf]$', fontsize = fontSize, rotation=90)
ax.set_xlabel('U', fontsize = fontSize)
#ax.set_ylim(0, 200)
ax.set_title('Pycnocline 0.15m', \
             fontsize = fontSize)
legends = (['Stratified', 'Non-stratified',  'Difference'])
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(legends, loc='center left', bbox_to_anchor=(1, 0.5),\
          fontsize = fontSize)
#fig.tight_layout()
ax.tick_params(axis='both', which='major', labelsize=fontSize)
fig.savefig('/home/karl/OpenFOAM/karl-5.0/run/MT/MT_PK/plots/%s/Gou5bKO%s'\
            %(topOrNot, nutType), bbox_inches = 'tight')

#Plots as in Gou et al figures 5b Interface 0.10
fig, ax = plt.subplots()
ax.plot(Us, FArrayBarge01Single, 'x' , Us, FArrayBarge01Two, 'o', Us, FDiffBarge01, '^')
ax.set_ylabel(r'$F [gf]$', fontsize = fontSize, rotation=90)
ax.set_xlabel('U', fontsize = fontSize)
#ax.set_ylim(0, 200)
ax.set_title('Pycnocline 0.1m', \
             fontsize = fontSize)
legends = (['Stratified', 'Non-stratified',  'Difference'])
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(legends, loc='center left', bbox_to_anchor=(1, 0.5),\
          fontsize = fontSize)
#fig.tight_layout()
ax.tick_params(axis='both', which='major', labelsize=fontSize)
fig.savefig('/home/karl/OpenFOAM/karl-5.0/run/MT/MT_PK/plots/%s/Gou5cKO%s'\
            %(topOrNot, nutType), bbox_inches = 'tight')



#Plots as in Gou et al figures 6a
fig2, ax2 = plt.subplots()
#seaborn.set(rc={'xticks': np.linspace(0,1.5, 10)})
c2twoPhaseIf = c2Emsmailpour(rho1, rho2, np.asarray(h1s))
FrS = Fr(Us, c2twoPhaseIf[0]) 
ax2.plot(FrS, CdDiffBarge01, 'bo')
FrS = Fr(Us, c2twoPhaseIf[1]) 
ax2.plot(FrS, CdDiffBarge015, 'g^')
FrS = Fr(Us, c2twoPhaseIf[2]) 
ax2.plot(FrS, CdDiffBarge02, 'rx')
ax2.set_ylabel(r'$C_{d_{add}}$', fontsize = fontSize*1.25)
ax2.set_xlabel('Fr', fontsize = fontSize)
legends = (['Pycnocline 0.1', 'Pycnocline 0.15', 'Pycnocline 0.20'])
box = ax2.get_position()
ax2.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax2.legend(legends, loc='center left', bbox_to_anchor=(1, 0.5), \
           fontsize = fontSize)
ax2.tick_params(axis='both', which='major', labelsize=fontSize)
#xtickValues = np.concatenate((np.arange(0,0.5,0.2), np.arange(0.5, 0.8, 0.1), np.arange(0.8, 1.5, 0.2)))
xtickValues = np.arange(0,1.5+.1, 0.3)
ax2.set_xticks(xtickValues)
ax2.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
ax2.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
ax2.grid(b=True, which='minor', color='grey', linewidth=0.2)
ax2.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1))
ax2.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.00125))
fig2.savefig('/home/karl/OpenFOAM/karl-5.0/run/MT/MT_PK/plots/%s/Gou6aKO%s'\
             %(topOrNot, nutType), bbox_inches = 'tight')

#Plots as in Gou et al figures 6b
fig3, ax3 = plt.subplots()
FrS = Fr(Us, c2twoPhaseIf[0]) 
ax3.plot(FrS, CdDiffBarge01/dtoh1squared[0], 'bo')
FrS = Fr(Us, c2twoPhaseIf[1]) 
ax3.plot(FrS, CdDiffBarge015/dtoh1squared[1], 'g^')
FrS = Fr(Us, c2twoPhaseIf[2]) 
ax3.plot(FrS, CdDiffBarge02/dtoh1squared[2], 'rx')
ax3.set_ylabel(r'$\frac{C_{d_{add}}}{(d/h_1)^2}$',fontsize = fontSize*1.25)
ax3.set_xlabel('Fr',fontsize = fontSize)
#ax3.set_ylim(0, .014)
legends = (['Pycnocline 0.1', 'Pycnocline 0.15', 'Pycnocline 0.2'])
box = ax3.get_position()
ax3.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax3.legend(legends, loc='center left', bbox_to_anchor=(1, 0.5), \
           fontsize = fontSize)
#ax3.set_title('Gou Figure 6b. %s. %s. \n CHECK DENOMINATOR' %(nutType, topOrNot))
#fig3.tight_layout()
ax3.tick_params(axis='both', which='major', labelsize=fontSize)
xtickValues = np.arange(0,1.5+.1, 0.3)
ax3.set_xticks(xtickValues)
ax3.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
ax3.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
ax3.grid(b=True, which='minor', color='grey', linewidth=0.2)
ax3.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1))
ax3.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.0025))

fig3.savefig('/home/karl/OpenFOAM/karl-5.0/run/MT/MT_PK/plots/%s/Gou6bKO%s'\
             %(topOrNot, nutType), bbox_inches = 'tight')

# Esmaeilpour Fig 2
fig5, ax5 = plt.subplots()
FrS = Fr(Us, c2twoPhaseIf[0]) 
ax5.plot(FrS, CdArrayBarge01Two, 'bo')
FrS = Fr(Us, c2twoPhaseIf[1]) 
ax5.plot(FrS, CdArrayBarge015Two, 'g^')
FrS = Fr(Us, c2twoPhaseIf[2]) 
ax5.plot(FrS, CdArrayBarge02Two, 'rx')
FrS = Fr(Us, c2twoPhaseIf[1]) 
ax5.plot(FrS, CdArrayBarge01Single, 'sienna')
#ax5.set_title('Esmaeilpour fig 2')
ax5.set_ylabel(r'$C_D$', fontsize = fontSize*1.25)
ax5.set_xlabel('$Fr$', fontsize = fontSize)
#ax5.set_title('Esmaeilpour figure 2. %s. %s' %(nutType, topOrNot))
legends = ['$D/P = 1$', '$D/P = 2/3$',  '$D/P = 1/2$', \
           'Non-stratified']
box = ax5.get_position()
ax5.set_position([box.x0, box.y0, box.width * 1.0, box.height*.8])
ax5.legend(legends, loc='center left', bbox_to_anchor=(-0.2, -0.55), \
           fontsize = fontSize, ncol=2)
#fig5.tight_layout()
ax5.tick_params(axis='both', which='major', labelsize=fontSize)
xtickValues = np.arange(0.3,1.4+.1, 0.3)
ytickValues = np.arange(0.0825,.11, 0.005)
ax5.set_xticks(xtickValues)
ax5.set_yticks(ytickValues)
#ax5.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
ax5.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
ax5.grid(b=True, which='minor', color='grey', linewidth=0.2)
ax5.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1))
ax5.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.0025))
fig5.savefig('/home/karl/OpenFOAM/karl-5.0/run/MT/MT_PK/plots/%s/esmaeilpour%s'\
             %(topOrNot, nutType), bbox_inches = 'tight')

# Esmaeilpour Fig 2, percentage difference from single
fig6, ax6= plt.subplots()
FrS = Fr(Us, c2twoPhaseIf[0]) 
ax6.plot(FrS, (np.squeeze(CdArrayBarge01Two)/CdArrayBarge01Single-1)*100, 'bo')
FrS = Fr(Us, c2twoPhaseIf[1]) 
ax6.plot(FrS, (np.squeeze(CdArrayBarge015Two)/CdArrayBarge01Single-1)*100, 'g^')
FrS = Fr(Us, c2twoPhaseIf[2]) 
ax6.plot(FrS, (np.squeeze(CdArrayBarge02Two)/CdArrayBarge01Single-1)*100, 'rx')
ax6.set_ylabel(r'$(\frac{C_d}{C_{d, single}}-1) \cdot 100$', fontsize = fontSize*1)
ax6.set_xlabel('$Fr$', fontsize = fontSize)
#ax6.set_title('Esmaeilpour figure 2. %s. %s \n Cd. Percentage difference from single' %(nutType, topOrNot))
legends = ['$D/P = 1$', '$D/P = 2/3$',  '$D/P = 1/2$']
box = ax6.get_position()
ax6.set_position([box.x0, box.y0, box.width * 1.0, box.height*.8])
ax6.legend(legends, loc='center left', bbox_to_anchor=(-0.2, -0.55), \
           fontsize = fontSize, ncol=2)
ax6.tick_params(axis='both', which='major', labelsize=fontSize)
#fig6.tight_layout()
xtickValues = np.arange(0.3,1.4+.1, 0.3)
ytickValues = np.arange(0,22+1, 5)
ax6.set_xticks(xtickValues)
ax6.set_yticks(ytickValues)
ax6.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
ax6.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
ax6.grid(b=True, which='minor', color='grey', linewidth=0.2)
ax6.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1))
fig6.savefig('/home/karl/OpenFOAM/karl-5.0/run/MT/MT_PK/plots/%s/esmaeilpour%sPercentage' \
             %(topOrNot, nutType), bbox_inches = 'tight')
