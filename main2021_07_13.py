import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
import scipy
import prepState
import timeEvo
from scipy.optimize import curve_fit


sigmaX = np.array([[0.,1.],[1.,0.]])
sigmaY = np.array([[0.,-1j],[1j,0.]])
sigmaZ = np.array([[1.,0.],[0.,-1.]])
sigmaPlus = .5 * (sigmaX + 1j*sigmaY)
sigmaMinus = .5 * (sigmaX - 1j*sigmaY)
kb = 8.617333262e-05

###
#
#prepare states and write it into a dictionary
#
###


N = 4 #chain length. counting starts at 0. arbitrary natural number.
#DO NOT do N > 10. takes ages to run.
states = {}
for i in range(N+1): #prepare every x state
	states["sigmaX", i] = prepState.stateX(N,i)
	states["sigmaY", i] = prepState.stateY(N,i)
	states["sigmaZ", i] = prepState.stateZ(N,i)

basis = {}
for i in range(N+1):
	basis["zero",i] = prepState.state0(N,i)
	basis["one",i] = prepState.state1(N,i)

###
#
#Hamiltonians of the system and partition function
#
###

Hint = 0
for i in range(N): # interaction hamiltonian, XY-Heisenberg
	Hint = Hint + np.matmul(states["sigmaX", i], states["sigmaX", i+1])
	Hint = Hint + np.matmul(states["sigmaY", i], states["sigmaY", i+1])

H0 = 0
for i in range(1,N+1):
		H0 = H0 - states["sigmaZ", i]

beta = 0.7 # 1/(kb*300) ----> choose a value for beta then calculate the temperature from that, avoids overflow error (kb is small)
rho = np.exp(-beta * np.diag(H0))
Z = np.sum(rho)
rho = np.diag(rho / Z)
first = prepState.state0(N,0) # get |0><0| state for first spin
rho = np.matmul(rho, first) # thermal state of the first spin


###
#
# Time evolution, very much unfinished.
#
###
t = np.linspace(0,100,100)
energyX = []
energyY = []
energyZ = []

"""

energies = {} #
for i in range(N+1):
	energyX = [] #
	energyY = [] # these are helper lists 
	energyZ = [] #
	for dt in t:
		dt = dt/(2*np.pi)
		rhoT = timeEvo.timeEvo(dt, rho, Hint)
		EoftX = timeEvo.energy(rhoT, states["sigmaX", i])
		EoftY = timeEvo.energy(rhoT, states["sigmaY", i])
		EoftZ = timeEvo.energy(rhoT, states["sigmaZ", i])
		energyX.append(EoftX)
		energyY.append(EoftY)
		energyZ.append(EoftZ)
	energies["energyX", i] = np.real(energyX) # somehow i get imaginary results here
	energies["energyY", i] = np.real(energyY) # or it gives me some 0.00000001j result.
	energies["energyZ", i] = energyZ # kinda useless, 


for i in range(N+1): #I chose hbar = 1 to avoid overflow errors. thats why it is time with units s/(Js) = 1/J.
	plt.plot(t/(2*np.pi), energies["energyX", i], label="x-energy")
	plt.xlabel("time [1/J]") 
	plt.ylabel("Energy [J]")
	plt.show()
	plt.plot(t/(2*np.pi), energies["energyY", i], label="y-energy")
	plt.xlabel("time [1/J]") 
	plt.ylabel("Energy [J]")
	plt.show()
	plt.plot(t/(2*np.pi), np.real(energies["energyZ", i]), label="z-energy, real part")
	plt.plot(t/(2*np.pi), np.imag(energies["energyZ", i]), label="z-energy, imag part")
	plt.xlabel("time [1/J]")
	plt.ylabel("Energy [J]")
	plt.show()
"""
# fidelity

fidelityList = []
rho_0 = partial_trace(rho) # partial trace to get state of i = 1 without i= 2,...,N
for dt in t:
	rhoF = timeEvo(dt, rho, Hint)
	rho_N = partial_trace(rhoF)
	fidelityList.append(timeEvo.fidelity(rho_0, rho_N))
plt.plot(t, fidelityList)
plt.show()
