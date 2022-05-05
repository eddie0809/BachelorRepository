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



beta = 1/100
rho = (-beta * sigmaz()).expm()
#everyRho = [(-beta * myStates["sigmaZ", n]).expm()/(-beta * myStates["sigmaZ", n]).expm().tr() for n in range(N+1)]
Z = rho.tr()
rho = rho/Z


alpha = 1.j * .25 * 1/(np.cosh(beta))**2
astar = -alpha
firstState = rho
initSystem = tensor(firstState, rho) #+ alpha * tensor(Qobj([[0,1],[0,0]]), Qobj([[0,0],[1,0]])) + astar * tensor(Qobj([[0,0],[1,0]]), Qobj([[0,1],[0,0]]))

for n in range(1,N):
	initSystem = tensor(initSystem, Qobj([[1,0],[0,0]]))

print("Spur der Dichtematrix des Systems: "+str(initSystem.tr()))

t = np.linspace(0,np.pi,200)

result = qutip.mesolve(Hint, initSystem, t, [], []).states
bigZ = sigmaz()
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




#print(result[-1].ptrace(N))

# kullback leibler divergence
"""
def myqEnt(N):
	qEnt = [entropy_vn(timeEvo(dt, initSystem, Hint).ptrace(N), base=np.e) for dt in t]
	myInfoflow = [(qEnt[i+1] - qEnt[i])/t[1] for i in range(0,len(qEnt)-1)]
	myHeatlim = [-1.j *(np.pi/3 *1 / (np.log(2)**2)) * (commutator(result[dt], Hint).unit().ptrace(N) * sigmaz()).tr() for dt in range(len(t))]
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
fig, ax = plt.subplots(figsize=(32/3,6))
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