### this file is just that the main file isnt filled with 2 line functions. one could argue that this is less efficient
import numpy as np
import math


def timeEvo(dt, rho, Hint): #time evolution of an operator rho
	U = np.exp(-1.j * Hint * dt)
	return np.matmul(np.matrix.getH(U), np.matmul(rho, U))


def energy(rho, sigma):
	return np.trace(np.matmul(sigma, rho))
