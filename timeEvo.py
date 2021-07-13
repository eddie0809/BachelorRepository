### this file is just that the main file isnt filled with 2 line functions. 
import numpy as np
import math


def timeEvo(dt, rho, Hint): #time evolution of an operator rho
	U = np.exp(-1.j * Hint * dt)
	return np.matmul(np.matrix.getH(U), np.matmul(rho, U)) #np.matrix.getH creates the hermitian conjugate


def energy(rho, sigma):
	return np.trace(np.matmul(sigma, rho))


def sinFit(x, a, b, phi, c): # defining fit function
	return a*np.sin(b*x + phi) + c 
