from unittest import result
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

#plt.rcParams.update({"text.usetex": True,}) 
#mpl.rcParams['text.latex.preamble'] = [r'\usepackage[utf8]{inputenc}\usepackage[T1]{fontenc}\usepackage{lmodern}\inputencoding{utf8}\usepackage{amsmath}\usepackage{amssymb}\usepackage{dsfont}\usepackage{mathtools}\usepackage{physics}']
#mpl.rcParams['font.family'] = ['serif']


def initSystem(alpha, firstState, secondState, N):
	astar = -alpha
	rhoSystem = tensor(firstState, secondState)+ alpha * tensor(Qobj([[0,1],[0,0]]), Qobj([[0,0],[1,0]])) + astar * tensor(Qobj([[0,0],[1,0]]), Qobj([[0,1],[0,0]]))
	for i in range(1,N):
		rhoSystem = tensor(rhoSystem, Qobj([[0,0],[0,1]]))
	return rhoSystem

beta = 1
rho = (-beta * sigmaz()).expm()
Z = rho.tr()
rho = rho/Z


alpha = 1.j * .25 * 1/(np.cosh(beta))**2

firstState = rho

rhoSystem = initSystem(alpha, firstState, rho, N)
rhoSystempoint5 = initSystem(.5*alpha, firstState, rho, N)
rhoSystempoint1 = initSystem(.1*alpha, firstState, rho, N)
rhoSystemzero = initSystem(.0*alpha, firstState, rho, N)

print("Spur der Dichtematrix des Systems: "+str(rhoSystem.tr()))

t_stop = np.pi
dt=1e-3

t = np.arange(0,t_stop,dt)

#result = qutip.mesolve(Hint, initSystem, t, [], []).states
bigZ = sigmaz()
"""
#myPlots.plotFidelity(initSystem, t, "State Transfer $F(t)=\mel{N}{U(t)}{1}$ of the first state to the last state,\\\\ \\ \\ \\ with $\\rho_0(0) = \\frac{1}{2}(\dyad{0}+\dyad{1})$ and $\\rho_{i\\neq 0}(0) = \dyad{0}$", 'FidelityONE')
egestest = [0 for dt in t]
colors=["#7aa0c4","#ca82e1" ,"#8bcd50","#e18882","#acb053"]
fig, axs = plt.subplots(5,1, sharex=True, sharey=True)#, figsize=(6,9))
axs = axs.ravel()
for i in [0,1,2,3,4]:
	expsigmaZ = []
	expsigmaZ = [(bigZ * result[dt].ptrace(i)).tr() for dt in range(len(t))]
	egestest = [egestest[n] + expsigmaZ[n] for n in range(len(expsigmaZ))]
	axs[i].plot(t, expsigmaZ, color = colors[i],label="$\expval{\sigma^z}_"+str(i)+"$")
	axs[i].set_ylabel("$\expval{\sigma^z}$")
	#axs[i].legend()
axs[0].set_title("$\expval{\sigma^z}$")
axs[4].set_xlabel("time")
fig.legend(loc=7)
fig.tight_layout()
fig.subplots_adjust(right=0.8)   
plt.show()
print("Gesamtenergie: "+str(sum(egestest)/len(egestest))+"\nStandardabweichung: " +str(np.std(egestest))+"\ntanh(beta) = "+str(np.tanh(beta)))
"""

maxEvo = qutip.mesolve(Hint, rhoSystem, t, [], []).states
point5Evo = qutip.mesolve(Hint, rhoSystempoint5, t, [], []).states
point1Evo = qutip.mesolve(Hint, rhoSystempoint1, t, [], []).states
zeroEvo = qutip.mesolve(Hint, rhoSystemzero, t, [], []).states


maxFid = [fidelity(firstState, maxEvo[dt].ptrace(N)) for dt in range(len(t))]
point5Fid = [fidelity(firstState, point5Evo[dt].ptrace(N)) for dt in range(len(t))]
point1Fid = [fidelity(firstState, point1Evo[dt].ptrace(N)) for dt in range(len(t))]
zeroFid = [fidelity(firstState, zeroEvo[dt].ptrace(N)) for dt in range(len(t))]

fig, (axs1, axs2, axs3, axs4) = plt.subplots(4,1,sharex=True, sharey=True, figsize=(12,6))


axs1.plot(t, maxFid, label="alpha = alpha max")
print('-')
axs2.plot(t, point5Fid, label="alpha = .5*alpha_max")
print('--')
axs3.plot(t, point1Fid, label="alpha = .1*alpha_max")
print('---')
axs4.plot(t, zeroFid, label="alpha = 0")
print('----')


fig.legend(loc=7)
fig.tight_layout()
fig.subplots_adjust(right=0.8)

plt.savefig("testingdingsplot.pdf")
plt.show()


# kullback leibler divergence

qEntmax = [entropy_vn(maxEvo[dt].ptrace(N), base=np.e) for dt in range(len(t))]
qEntpoint5 = [entropy_vn(point5Evo[dt].ptrace(N), base=np.e) for dt in range(len(t))]
qEntpoint1 = [entropy_vn(point1Evo[dt].ptrace(N), base=np.e) for dt in range(len(t))]
qEntzero = [entropy_vn(zeroEvo[dt].ptrace(N), base=np.e) for dt in range(len(t))]
print('#')
infoFlowmax = np.square(np.diff(qEntmax)/t[1])
infoFlowpoint5 = np.square(np.diff(qEntpoint5)/t[1])
infoFlowpoint1 = np.square(np.diff(qEntpoint1)/t[1])
infoFlowzero = np.square(np.diff(qEntzero)/t[1])
print('##')
maxHeatlim = [-1.j *(np.pi/3 *1 / (np.log(2)**2)) * (commutator(Hint, maxEvo[dt]).ptrace(N) * sigmaz()).tr() for dt in range(len(t))]
point5Heatlim = [-1.j *(np.pi/3 *1 / (np.log(2)**2)) * (commutator(Hint, point5Evo[dt]).ptrace(N) * sigmaz()).tr() for dt in range(len(t))]
point1Heatlim = [-1.j *(np.pi/3 *1 / (np.log(2)**2)) * (commutator(Hint, point1Evo[dt]).ptrace(N) * sigmaz()).tr() for dt in range(len(t))]
zeroHeatlim = [-1.j *(np.pi/3 *1 / (np.log(2)**2)) * (commutator(Hint, zeroEvo[dt]).ptrace(N) * sigmaz()).tr() for dt in range(len(t))]
print('###')
fig, (ax1, ax2, ax3, ax4) = plt.subplots(4,1,sharex=True, sharey=True, figsize=(12,6))
print('####')

ax1.plot(t, maxHeatlim, label="alpha = alpha max")
ax1.plot(t[:-1], infoFlowmax, label="alpha = alpha max")
print('-')
ax2.plot(t, point5Heatlim, label="alpha = .5*alpha_max")
ax2.plot(t[:-1], infoFlowpoint5, label="alpha = .5*alpha_max")
print('--')
ax3.plot(t, point1Heatlim, label="alpha = .1*alpha_max")
ax3.plot(t[:-1], infoFlowpoint1, label="alpha = .1*alpha_max")
print('---')
ax4.plot(t, zeroHeatlim, label="alpha = 0")
ax4.plot(t[:-1], infoFlowzero, label="alpha = 0")
print('----')

#ax1.axhline(0,color='grey')
#ax2.axhline(0,color='grey')
#ax3.axhline(0,color='grey')


fig.legend(loc=7)
fig.tight_layout()
fig.subplots_adjust(right=0.8)

plt.savefig("testplot.pdf")
plt.show()

"""
def myqEnt(N):
	qEnt = [entropy_vn(timeEvo(dt, initSystem, Hint).ptrace(N), base=np.e) for dt in t]
	myInfoflow = [(qEnt[i+1] - qEnt[i])/t[1] for i in range(0,len(qEnt)-1)]
	myHeatlim = [-1.j *(np.pi/3 *1 / (np.log(2)**2)) * (commutator(Hint, result[dt]).ptrace(N) * sigmaz()).tr() for dt in range(len(t))]
	#myHeatlim = [(np.pi/3 *1 / (np.log(2)**2)) * expect(timeEvo(dt, initSystem, Hint).ptrace(N), sigmaz()) for dt in t]
	#myHeatlim = deriv(t[1], myHeatlim)
	#myHeat = [-1.j * (commutator(timeEvo(dt, initSystem, Hint), Hint).ptrace(N) * -sigmaz()).tr() for dt in t]
	return qEnt, myInfoflow, myHeatlim#, myHeat
myKLdiv = [entropy_relative(firstState.unit(), timeEvo(dt, initSystem, Hint).ptrace(N)) for dt in t]
#print(myKLdiv[-1])

helper = myqEnt(N)[1] 
myInfoflow2 = [helper[n]**2 for n in range(0,199)]

myTest = [-1.j * (commutator(timeEvo(dt, initSystem, Hint), Hint).unit().ptrace(N) * sigmay()).tr() for dt in t]
#print(myTest)
myHeat = [-1.j * (commutator(timeEvo(dt, initSystem, Hint), Hint).unit().ptrace(N) * sigmaz()).tr() for dt in t]
#print(np.mean(myHeat))
#print([commutator(result[dt], Hint).ptrace(N) for dt in range(len(t))])
fig, ax = plt.sub(figsize=(32/3,6))
#ax.plot(t, myKLdiv, label="$\mathcal{D}_\\text{KL}(\\rho^1(0)||\\rho^N(t))$")
#ax.set_title("Kullback-Leibler Divergence (relative entropy) \\\\ from $\\rho^1(0)$ to $\\rho^N(t)$ with $\\rho^1(0) = e^{-\\beta\sigma^z}/\Tr[e^{-\\beta\sigma^z}]$")
# 
# #\\frac{e^{-\\beta\sigma^z}}{\Tr[e^{-\\beta\sigma^z}]}

#plt.plot(t[:-1], myInfoflow2, label="$\dot{\mathcal{I}}_N^2$")
#plt.plot(t, myqEnt(N)[2], label="$\dot{E}_N\pi/(3\ln^2(2))$")
#ax.plot(t, myHeat, label="$\dot{\mathcal{Q}}_0$")
ax.set_xlabel("time")
#ax.set_ylabel("$\dot{\mathcal{Q}}_0$")
#ax.set_title("Heat flow at Qubit 0, $\dot{\mathcal{Q}}_0$. The grey bar indicates its limit predicted by Pendry, 1983.")
#plt.show()
#ax.plot(t[:-1], deriv(t[2]-t[1], myKLdiv), label="$\dot{\mathcal{D}}_\\text{KL}$")
#plt.title("Information flow from $\\rho^1(0)$ to $\\rho^N(t)$")
#plt.xlabel("time")
#plt.ylabel("$\dot{\mathcal{D}}_\\text{KL}$")
#plt.show()
#ax.plot(t,qEnt, label="$\mathcal{S}_{vn}$")
#ax.plot(t, np.real(myHeat), label="$\dot{\mathcal{Q}}$")
#ax.plot(t[:-1], dumbHeat)
ax.set_title("$\dot{\mathcal{I}}^2$ and its upper bound, $\dot{E}\pi/(3\ln^2(2))$ [Pendry, 1983], with $\\rho^0(0) = e^{-\\beta\sigma^z}/\Tr[e^{-\\beta\sigma^z}]$")
#plt.xlabel("time")
#plt.ylabel("$\dot{Q}_\\text{KL}$")
ax.legend()
#ax.axhline(np.pi/3 /(beta**2),color='grey')
ax.axhline(0,color='grey')
#plt.savefig("lost.png")
#plt.show()
plt.close('all')
"""