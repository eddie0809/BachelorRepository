import numpy as np
import math


def timeEvo(dt, rho, Hint): #time evolution of an operator rho
	U = np.exp(-1.j * Hint * dt)
	return np.matmul(np.matrix.getH(U), np.matmul(rho, U))


def energy(rho, sigma):
	return np.trace(np.matmul(sigma, rho))
