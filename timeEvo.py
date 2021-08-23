### this file is just that the main file isnt filled with 2 line functions.
import numpy as np
from scipy import linalg as la


def timeEvo(dt, rho, Hint): #time evolution of an operator rho
	# The more effective way to make expontial of a matrix is using `scipy.linalg.expm`
	U = la.expm(-1.j * Hint * dt)
	# In this case it's faster to define the unitary dagger than apply `matrix.getH`.
	Ud = la.expm(1.j * Hint * dt)
	return np.matmul(Ud, np.matmul(rho, U))


def energy(rho, sigma):
	return np.trace(np.matmul(sigma, rho))


def sinFit(x, a, b, phi, c): # defining fit function
	return a*np.sin(b*x + phi) + c

def fidelity(rho, sigma):
	"""Definition of fidelity from Wikipedia.
	"""
	sqrt_rho = la.sqrtm(rho)
	argument = np.matmul(sqrt_rho, np.matmul(sigma, sqrt_rho))
	return (np.trace(la.sqrtm(argument)))**2 
