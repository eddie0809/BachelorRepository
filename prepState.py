### this file is just so that this massive block of just defining functions is not in the main file
import math
import numpy as np
from qutip import *

def arbstate(N,i,state): # return an arbitrary prepared state
	if i == 0:
		sigma0 = state
	else:
		sigma0 = qeye(2)
	if i != 0:
		for n in range(0,i-1):
			sigma0 = tensor(sigma0, qeye(2))
		sigma0 = tensor(sigma0, state)
		for n in range(i,N):
			sigma0 = tensor(sigma0, qeye(2))
		return sigma0
	else:
		for n in range(0,N):
			sigma0 = tensor(sigma0,qeye(2))
		return sigma0

def state0(N, i): #preparing the |0><0| state of the i-th element of a chain with length N
	if i == 0:
		sigma0 = Qobj(np.array([[1.,0.],[0.,0.]]))
	else:
		sigma0 = qeye(2)
	if i != 0:
		for n in range(0,i-1):
			sigma0 = tensor(sigma0, qeye(2))
		sigma0 = tensor(sigma0, Qobj(np.array([[1.,0.],[0.,0.]])))
		for n in range(i,N):
			sigma0 = tensor(sigma0, qeye(2))
		return sigma0
	else:
		for n in range(0,N):
			sigma0 = tensor(sigma0,qeye(2))
		return sigma0

def state1(N, i): #preparing the |1><1| state of the i-th element of a chain with length N
	if i == 0:
		sigma0 = Qobj(np.array([[0.,0.],[0.,1.]]))
	else:
		sigma0 = qeye(2)
	if i != 0:
		for n in range(0,i-1):
			sigma0 = tensor(sigma0, qeye(2))
		sigma0 = tensor(sigma0, Qobj(np.array([[0.,0.],[0.,1.]])))
		for n in range(i,N):
			sigma0 = tensor(sigma0, qeye(2))
		return sigma0
	else:
		for n in range(0,N):
			sigma0 = tensor(sigma0,qeye(2))
		return sigma0

def stateX(N, i): #preparing the X state of the i-th element of a chain with length N
	if i == 0:
		sigma0 = Qobj(np.array([[0.,1.],[1.,0.]]))
	else:
		sigma0 = qeye(2)
	if i != 0:
		for n in range(0,i-1):
			sigma0 = tensor(sigma0, qeye(2))
		sigma0 = tensor(sigma0, Qobj(np.array([[0.,1.],[1.,0.]])))
		for n in range(i,N):
			sigma0 = tensor(sigma0, qeye(2))
		return sigma0
	else:
		for n in range(0,N):
			sigma0 = tensor(sigma0,qeye(2))
		return sigma0


def stateY(N, i): #preparing the Y state of the i-th element of a chain with length N
	if i == 0:
		sigma0 = Qobj(np.array([[0.,-1j],[1j,0.]]))
	else:
		sigma0 = qeye(2)
	if i != 0:
		for n in range(0,i-1):
			sigma0 = tensor(sigma0, qeye(2))
		sigma0 = tensor(sigma0, Qobj(np.array([[0.,-1j],[1j,0.]])))
		for n in range(i,N):
			sigma0 = tensor(sigma0, qeye(2))
		return sigma0
	else:
		for n in range(0,N):
			sigma0 = tensor(sigma0,qeye(2))
		return sigma0


def stateZ(N, i): #preparing the Y state of the i-th element of a chain with length N
	if i == 0:
		sigma0 = Qobj(np.array([[1.,0.],[0.,-1.]]))
	else:
		sigma0 = qeye(2)
	if i != 0:
		for n in range(0,i-1):
			sigma0 = tensor(sigma0, qeye(2))
		sigma0 = tensor(sigma0, Qobj(np.array([[1.,0.],[0.,-1.]])))
		for n in range(i,N):
			sigma0 = tensor(sigma0, qeye(2))
		return sigma0
	else:
		for n in range(0,N):
			sigma0 = tensor(sigma0,qeye(2))
		return sigma0

def vector0(N, i): #preparing the |0> state of the i-th element of a chain with length N
	if i == 0:
		sigma0 = Qobj(np.array([[1.],[0.]]))
	else:
		sigma0 = qeye(2)
	if i != 0:
		for n in range(0,i-1):
			sigma0 = tensor(sigma0, qeye(2))
		sigma0 = tensor(sigma0, Qobj(np.array([[1.],[0.]])))
		for n in range(i,N):
			sigma0 = tensor(sigma0, qeye(2))
		return sigma0
	else:
		for n in range(0,N):
			sigma0 = tensor(sigma0,qeye(2))
		return sigma0

def vector1(N, i): #preparing the |1> state of the i-th element of a chain with length N
	if i == 0:
		sigma0 = Qobj(np.array([[0.],[1.]]))
	else:
		sigma0 = qeye(2)
	if i != 0:
		for n in range(0,i-1):
			sigma0 = tensor(sigma0, qeye(2))
		sigma0 = tensor(sigma0, Qobj(np.array([[0.],[1.]])))
		for n in range(i,N):
			sigma0 = tensor(sigma0, qeye(2))
		return sigma0
	else:
		for n in range(0,N):
			sigma0 = tensor(sigma0,qeye(2))
		return sigma0