### this file is just that the main file isnt filled with 2 line functions.
import numpy as np
from scipy import linalg as la


def timeEvo(dt, rho, Hint): #time evolution of an operator rho
	# The more effective way to make expontial of a matrix is using 'scipy.linalg.expm'
	U = -1.j * Hint * dt
	U = U.expm() 
	#U = la.expm(-1.j * Hint * dt)
	# In this case it's faster to define the unitary dagger than apply 'matrix.getH'(= dagger)
	Ud = 1.j * Hint * dt
	Ud = Ud.expm()
	#Ud = la.expm(1.j * Hint * dt)
	#return np.matmul(Ud, np.matmul(rho, U))
	rho = rho * Ud
	return U * rho

def energy(rho, sigma):
	return np.trace(np.matmul(sigma, rho))

"""
def fidelity(rho, sigma):
	sqrt_rho = la.sqrtm(rho)
	argument = np.matmul(sqrt_rho, np.matmul(sigma, sqrt_rho))
	return (np.trace(la.sqrtm(argument)))**2
	# fidelity is usually defined for pure states, this definition uses density matrices
	# so that mixed states can be parsed as well.
	# https://doi.org/10.1016/0034-4877(76)90060-4
"""


def qubitfidelity(rho, sigma):# only for 2x2 matrices
	tr = np.trace(np.matmul(rho,sigma)) + 2 * np.sqrt(la.det(rho) * la.det(sigma))


def stateTransfer(dt, Hint, first, last):
	U = la.expm(-1.j * Hint * dt)
	Fdings = np.matmul(U,first)
	F = np.matmul(last, Fdings)
	return F
	
