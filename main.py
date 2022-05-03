import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from numpy import linalg as LA
import prepState
from timeEvo import *
import myPlots
from generateStates import *
from qutip import *
import partialtrace

# for using tex formatting and font in plots

plt.rcParams.update({"text.usetex": True,}) 
mpl.rcParams['text.latex.preamble'] = [r'\usepackage[utf8]{inputenc}\usepackage[T1]{fontenc}\usepackage{lmodern}\inputencoding{utf8}\usepackage{amsmath}\usepackage{amssymb}\usepackage{dsfont}\usepackage{mathtools}\usepackage{physics}']
mpl.rcParams['font.family'] = ['serif']



beta = 1

#rho = np.exp(-beta * H0.diag())
#Z = np.sum(rho)
#rho = np.diag(rho / Z)
rho = (-beta * sigmaz()).expm()
everyRho = [(-beta * myStates["sigmaZ", n]).expm()/(-beta * myStates["sigmaZ", n]).expm().tr() for n in range(N+1)]
Z = rho.tr()
rho = rho/Z


firstState = Qobj([[1,0],[0,0]])
initSystem = firstState.unit()

for n in range(N):
	initSystem = tensor(initSystem, rho)


t = np.linspace(0,np.pi,200)

result = qutip.mesolve(Hint, initSystem, t, [], []).states

#myPlots.plotFidelity(initSystem, t, "Fidelity, $\\rho_0 = \dyad{1}$")#, 'Fidelityoneone')


#print(result[-1].ptrace(N))

# kullback leibler divergence
qKLdiv = [entropy_relative(Qobj([[0,0],[0,1]]), result[dt].ptrace(N)) for dt in range(len(t))]
myKLdiv = [entropy_relative(firstState.unit(), timeEvo(dt, initSystem, Hint).ptrace(N)) for dt in t]
print(myKLdiv[-1])
dumbHeat = [(myKLdiv[i+1] - myKLdiv[i])/beta for i in range(0,len(myKLdiv)-1)]
myHeat = [1.j * (commutator(Hint, timeEvo(dt, initSystem, Hint)) * Hges).tr() for dt in t]

fig, ax = plt.subplots()
ax.plot(t, myKLdiv, label="$\mathcal{D}_\\text{KL}(\\rho^1(0)||\\rho^N(t))$")
#plt.title("Kullback-Leibler Divergence from $\\rho^1(0)$ to $\\rho^N(t)$")
#plt.plot(t[:-1], deriv(t[2]-t[1], myKLdiv))
ax.set_xlabel("time")
#plt.ylabel("$\mathcal{D}_\\text{KL}$")
#plt.show()
ax.plot(t[:-1], deriv(t[2]-t[1], myKLdiv), label="$\dot{\mathcal{D}}_\\text{KL}$")
#plt.title("Information flow from $\\rho^1(0)$ to $\\rho^N(t)$")
#plt.xlabel("time")
#plt.ylabel("$\dot{\mathcal{D}}_\\text{KL}$")
#plt.show()
ax.plot(t, np.real(myHeat), label="$\dot{\mathcal{Q}}$")
ax.plot(t[:-1], dumbHeat)
ax.set_title("$\mathcal{D}_\\text{KL}, \dot{\mathcal{D}}_\\text{KL}$ and $\dot{\mathcal{Q}}$ from $\\rho^1(0)= \dyad{0}$ to $\\rho^N(t)$")
#plt.xlabel("time")
#plt.ylabel("$\dot{Q}_\\text{KL}$")
ax.legend()
ax.axhline(np.pi/3 * beta,color='grey')
ax.axhline(0,color='grey')
plt.show()