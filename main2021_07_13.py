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
"""

#some examples to see if it works
ax = prepState.prepStateX(3, 3) #this should be the sigmaX matrix along the main diagonal, 16x16 matrix 
ay = prepState.prepStateY(3, 0) 
az = prepState.prepStateZ(3, 3)
x0 = prepState.prepStateX(3, 0) #
bx = prepState.prepStateX(5, 4) #a larger chain length to see if it scales

print(ax)
print(ay)
#print(az)
print(bx)


#this following part is to visualise it in a plot and does not add anything
plt.matshow(bx)
plt.matshow(x0)
plt.matshow(ax)
#plt.matshow(ay) <- can't use matshow on a matrix with complex elements
plt.matshow(az)
plt.show() #plots 
"""
###
#
#prepare states and write it into a dictionary
#
###


N = 4
states = {}
for i in range(N+1): #prepare every x state
	states["sigmaX", i] = prepState.stateX(N,i)
	states["sigmaY", i] = prepState.stateY(N,i)
	states["sigmaZ", i] = prepState.stateZ(N,i)
	
#print(states)


###
#
#Hamiltonian of the system
#
###

Hint = 0
for i in range(N): #
	Hint = Hint + np.matmul(states["sigmaX", i], states["sigmaX", i+1])
	Hint = Hint + np.matmul(states["sigmaY", i], states["sigmaY", i+1])

H0 = 0
for i in range(1,N+1):
		H0 = H0 - states["sigmaZ", i]
#print(H)
# plt.matshow(np.real(Hint))
# plt.show()


###
#
# TEST!
#
###

#EigVa, EigVe = LA.eig(Hint)
#print(np.real(EigVa))
#plt.matshow(np.real(EigVe))
#plt.show()

###
#
# partition function ----> overflow error because kb is too small for python
#
###

beta = 0.7 # 1/(kb*300) ----> choose a value for beta then calculate the temperature from that
rho = np.exp(-beta * np.diag(H0))
Z = np.sum(rho)
rho = np.diag(rho / Z)
first = prepState.state0(N,0)
rho = np.matmul(rho, first)
#plt.matshow(rho)
#plt.show()

# print(Z)


###
#
# Time evolution
#
###


def sinFit(x, a, b, phi, c):
	return a*np.sin(b*x + phi) + c 

t = np.linspace(0,100,100)
energyX = []
energyY = []
energyZ = []

#### X states
"""
for dt in t:
	dt = dt/(2*np.pi)
	rhoT = timeEvo.timeEvo(dt, rho, Hint)
	EoftX = timeEvo.energy(rhoT, states["sigmaX", 1])
	EoftZ = timeEvo.energy(rhoT, states["sigmaZ", 1])
	energyX.append(EoftX)
	energyZ.append(EoftZ)
"""
energies = {}
for i in range(N+1):
	energyX = []
	energyY = []
	energyZ = []
	for dt in t:
		dt = dt/(2*np.pi)
		rhoT = timeEvo.timeEvo(dt, rho, Hint)
		EoftX = timeEvo.energy(rhoT, states["sigmaX", i])
		EoftY = timeEvo.energy(rhoT, states["sigmaY", i])
		EoftZ = timeEvo.energy(rhoT, states["sigmaZ", i])
		energyX.append(EoftX)
		energyY.append(EoftY)
		energyZ.append(EoftZ)
	energies["energyX", i] = np.real(energyX)
	energies["energyY", i] = np.real(energyY)
	energies["energyZ", i] = energyZ

poptX, pcovX = curve_fit(sinFit, t/(2*np.pi), np.real(energyX), p0=[.7, 2, 0, 7])
#print(poptX)
#poptY, pcovY = curve_fit(sinFit, t/(2*np.pi), np.real(energyY), bounds=([.35, 2.2, 0, -.1],[.4, 2.6, 2, .1]))
#print(poptY)
for i in range(N+1):
	plt.plot(t/(2*np.pi), energies["energyX", i], label="x-energy")
	plt.show()
	plt.plot(t/(2*np.pi), energies["energyY", i], label="y-energy")
	plt.show()
	plt.plot(t/(2*np.pi), np.real(energies["energyZ", i]), label="z-energy, real part")
	plt.plot(t/(2*np.pi), np.imag(energies["energyZ", i]), label="z-energy, imag part")
	plt.show()

print(sigmaPlus)
"""
plt.xlabel("time [1/J]")
plt.ylabel("Energy [J]")
plt.plot(t/(2*np.pi), sinFit(t/(2*np.pi), *poptX))
plt.plot(t/(2*np.pi), np.real(energyX))
#plt.plot(t/(2*np.pi), np.imag(energyX))
plt.show()
plt.xlabel("time [1/J]")
plt.ylabel("Energy [J]")
#plt.plot(t/(2*np.pi), sinFit(t/(2*np.pi), *poptY))
plt.plot(t/(2*np.pi), np.real(energyY))
#plt.plot(t/(2*np.pi), np.imag(energyY))
plt.show()
plt.xlabel("time [1/J]")
plt.ylabel("Energy [J]")
plt.plot(t/(2*np.pi), np.real(energyZ))
plt.plot(t/(2*np.pi), np.imag(energyZ))
plt.show()
"""
