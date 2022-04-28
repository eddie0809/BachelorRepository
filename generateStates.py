import prepState
import numpy as np


N=4
myStates = {}
myBasis = {}
for i in range(N+1):
    myStates["sigmaX", i] = prepState.stateX(N,i)
    myStates["sigmaY", i] = prepState.stateY(N,i)
    myStates["sigmaZ", i] = prepState.stateZ(N,i)
    myBasis["zero",i] = prepState.state0(N,i)
    myBasis["one",i] = prepState.state1(N,i)


l = 1

Hint = 0
for i in range(N): # interaction hamiltonian, XY-Heisenberg
	Hint = Hint + l / 2 * np.sqrt((i+1)*(N-i)) * (myStates["sigmaX", i]* myStates["sigmaX", i+1] + myStates["sigmaY", i] * myStates["sigmaY", i+1])
J = 1/2
Hint = J * Hint # change interaction strength 

H0 = 0
for i in range(1,N+1):
		H0 = H0 - myStates["sigmaZ", i]