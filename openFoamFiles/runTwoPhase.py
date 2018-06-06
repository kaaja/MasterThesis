#from generateUfile import writeUfile
import generateUfile
import generateKfile
import generateOmegafile
import generateNutfile
import generateControlDict
import os
from subprocess import call


def runMultipleCases():
    calculateK = lambda U, I = 0.01: 3./2*(U*I)**2
    calculateOmega = lambda k, l: k**.5/(0.09**.25*l)
    calculateEpsilon = lambda k, l, c_mu=0.09: c_mu**.75*k**(3./2)/l
    calculateNutKeps = lambda k, epsilon: 0.09*k**2/epsilon
    calculateNutKomega = lambda k, omega: k/omega#from https://turbmodels.larc.nasa.gov/wilcox.html
    #velocities = [0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24]#, 0.3, .4]
    #velocities = [0.10]#, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24, 0.3, .4, .5, .6]
    velocities = [0.08]#, 0.10]
    turbulentIntensity = .05
    shipLength = 0.6
    cylinderRadiusToShipLength = 1.
    cylinderRadius = cylinderRadiusToShipLength*shipLength
    turbulentLength = cylinderRadius/100.
    for velocity in velocities:
        call(["./Allclean"])
        k = calculateK(velocity, turbulentIntensity)
        omega = calculateOmega(k, turbulentLength)
        epsilon = calculateEpsilon(k, turbulentLength)
        nutKeps = calculateNutKeps(k, epsilon)
        nutKomega = calculateNutKomega(k, omega)
        generateUfile.writeUfile(velocity)
        generateKfile.writeKfile(k)
        generateNutfile.writeNutfile(nutKomega)
        generateOmegafile.writeOmegafile(omega)
        generateControlDict.writeControlDictFile(velocity)
        call(["./Allrun", str(velocity).replace(".", "")])
runMultipleCases()
