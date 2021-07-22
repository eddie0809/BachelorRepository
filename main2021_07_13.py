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
for i in range(N): #
	Hint = Hint + np.matmul(states["sigmaX", i], states["sigmaX", i+1])
	Hint = Hint + np.matmul(states["sigmaY", i], states["sigmaY", i+1])

H0 = 0
for i in range(1,N+1):
		H0 = H0 - states["sigmaZ", i]

beta = 0.7 # 1/(kb*300) ----> choose a value for beta then calculate the temperature from that
rho = np.exp(-beta * np.diag(H0))
Z = np.sum(rho)
rho = np.diag(rho / Z)
first = prepState.state0(N,0)
rho = np.matmul(rho, first)


###
#
# Time evolution
#
###

t = np.linspace(0,100,100)
energyX = []
energyY = []
energyZ = []

#### X states
energies = {} # all the energies are "stored" in this dictionary.
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
	energies["energyX", i] = np.real(energyX) # this np.real is only to remove the +0.j from the numbers. no numbers are with imaginary part here
	energies["energyY", i] = np.real(energyY) # but they still have the +0.j.
	energies["energyZ", i] = energyZ # there are some energies with imaginary part (?) in this one (what does it mean? i don't know). i therefore didnt remove it here

#poptX, pcovX = curve_fit(timeEvo.sinFit, t/(2*np.pi), np.real(energyX), p0=[.7, 2, 0, 7]) <-- the function in X and Y should be sinusoidal.
#print(poptX)
#poptY, pcovY = curve_fit(timeEvo.sinFit, t/(2*np.pi), np.real(energyY), bounds=([.35, 2.2, 0, -.1],[.4, 2.6, 2, .1]))
#print(poptY)

for i in range(N+1):
	plt.plot(t/(2*np.pi), energies["energyX", i], label="x-energy")
	plt.xlabel("time [1/J]") #in timeEvo.py i chose hbar = 1 to avoid overflow errors. thats why it is time with units s/(Js) = 1/J.
	plt.ylabel("Energy [J]")
	plt.show()
	plt.plot(t/(2*np.pi), energies["energyY", i], label="y-energy")
	plt.xlabel("time [1/J]") #in timeEvo.py i chose hbar = 1 to avoid overflow errors. thats why it is time with units s/(Js) = 1/J.
	plt.ylabel("Energy [J]")
	plt.show()
	plt.plot(t/(2*np.pi), np.real(energies["energyZ", i]), label="z-energy, real part")
	plt.plot(t/(2*np.pi), np.imag(energies["energyZ", i]), label="z-energy, imag part")
	plt.xlabel("time [1/J]") #in timeEvo.py i chose hbar = 1 to avoid overflow errors. thats why it is time with units s/(Js) = 1/J.
	plt.ylabel("Energy [J]")
	plt.show()

# fidelity

fidelityList = []
rho_0 = partial_trace(rho)
for dt in t:
	rhoF = timeEvo(dt, rho, Hint)
	rho_N = partial_trace(rhoF)
	fidelityList.append(timeEvo.fidelity(rho_0, rho_N))
plt.plot(t, fidelityList)
plt.show()
