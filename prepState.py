import math
import numpy as np

def state0(N, i): #preparing the X state of the i-th element of a chain with length N
	if i == 0:
		sigma0 = np.array([[1.,0.],[0.,0.]])
	else:
		sigma0 = np.eye(2)
	if i != 0:
		for n in range(0,i-1):
			sigma0 = np.kron(sigma0, np.eye(2))
		sigma0 = np.kron(sigma0, np.array([[1.,0.],[0.,0.]]))
		for n in range(i,N):
			sigma0 = np.kron(sigma0, np.eye(2))
		return sigma0
	else:
		for n in range(0,N):
			sigma0 = np.kron(sigma0,np.eye(2))
		return sigma0

def stateX(N, i): #preparing the X state of the i-th element of a chain with length N
	if i == 0:
		sigma0 = np.array([[0.,1.],[1.,0.]])
	else:
		sigma0 = np.eye(2)
	if i != 0:
		for n in range(0,i-1):
			sigma0 = np.kron(sigma0, np.eye(2))
		sigma0 = np.kron(sigma0, np.array([[0.,1.],[1.,0.]]))
		for n in range(i,N):
			sigma0 = np.kron(sigma0, np.eye(2))
		return sigma0
	else:
		for n in range(0,N):
			sigma0 = np.kron(sigma0,np.eye(2))
		return sigma0


def stateY(N, i): #preparing the Y state of the i-th element of a chain with length N
	if i == 0:
		sigma0 = np.array([[0.,-1j],[1j,0.]])
	else:
		sigma0 = np.eye(2)
	if i != 0:
		for n in range(0,i-1):
			sigma0 = np.kron(sigma0, np.eye(2))
		sigma0 = np.kron(sigma0, np.array([[0.,-1j],[1j,0.]]))
		for n in range(i,N):
			sigma0 = np.kron(sigma0, np.eye(2))
		return sigma0
	else:
		for n in range(0,N):
			sigma0 = np.kron(sigma0,np.eye(2))
		return sigma0


def stateZ(N, i): #preparing the Y state of the i-th element of a chain with length N
	if i == 0:
		sigma0 = np.array([[1.,0.],[0.,-1.]])
	else:
		sigma0 = np.eye(2)
	if i != 0:
		for n in range(0,i-1):
			sigma0 = np.kron(sigma0, np.eye(2))
		sigma0 = np.kron(sigma0, np.array([[1.,0.],[0.,-1.]]))
		for n in range(i,N):
			sigma0 = np.kron(sigma0, np.eye(2))
		return sigma0
	else:
		for n in range(0,N):
			sigma0 = np.kron(sigma0,np.eye(2))
		return sigma0
