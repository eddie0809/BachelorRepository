import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
import scipy

sigmaX = np.array([[0.,1.],[1.,0.]])
sigmaY = np.array([[0.,-1j],[1j,0.]])
sigmaZ = np.array([[1.,0.],[0.,-1.]])
sigmaPlus = .5 * (sigmaX + 1j*sigmaY)
sigmaMinus = .5 * (sigmaX - 1j*sigmaY)
kb = 8.617333262e-05

def prepState0(N, i): #preparing the X state of the i-th element of a chain with length N
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

def prepStateX(N, i): #preparing the X state of the i-th element of a chain with length N
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


def prepStateY(N, i): #preparing the Y state of the i-th element of a chain with length N
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


def prepStateZ(N, i): #preparing the Y state of the i-th element of a chain with length N
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

"""
#some examples to see if it works
ax = prepStateX(3, 3) #this should be the sigmaX matrix along the main diagonal, 16x16 matrix 
ay = prepStateY(3, 0) 
az = prepStateZ(3, 3)
x0 = prepStateX(3, 0) #
bx = prepStateX(5, 4) #a larger chain length to see if it scales

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


N = 3
states = {}
for i in range(N+1): #prepare every x state
	states["sigmaX", i] = prepStateX(N,i)
	states["sigmaY", i] = prepStateY(N,i)
	states["sigmaZ", i] = prepStateZ(N,i)
	
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

EigVa, EigVe = LA.eig(Hint)
print(np.real(EigVa))
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
first = prepState0(N,0)
rho = np.matmul(rho, first)
#plt.matshow(rho)
#plt.show()

# print(Z)


###
#
# Time evolution
#
###
from scipy import constants

deltat = 0
U = np.exp(-1.j * Hint)
print(U)
plt.matshow(np.real(U))
plt.colorbar()
plt.show()
rho1 = np.matmul(np.matrix.getH(U), np.matmul(rho, U))
#print(rho1)
#plt.matshow(np.real(rho1))
#plt.colorbar()
#plt.matshow(np.imag(rho1))
#plt.colorbar()
#plt.show()


