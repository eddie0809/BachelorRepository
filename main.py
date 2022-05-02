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





beta = 1#/3000 # 1/(kb*300) ----> choose a value for beta then calculate the temperature from that, avoids overflow error (kb is small)
#rho = np.exp(-beta * H0.diag())
#Z = np.sum(rho)
#rho = np.diag(rho / Z)
rho = (-beta * sigmaz()).expm()
everyRho = [(-beta * myStates["sigmaZ", n]).expm()/(-beta * myStates["sigmaZ", n]).expm().tr() for n in range(N+1)]
Z = rho.tr()
rho = rho/Z


initSystem = qeye(2)

for n in range(N):
	initSystem = tensor(initSystem, (-beta * sigmaz()).expm())


t = np.linspace(0,np.pi,200)

result = qutip.mesolve(Hges, initSystem, t, [], []).states

#myPlots.plotFidelity(initSystem, t, "Fidelity, $\\rho_0 = \dyad{1}$")#, 'Fidelityoneone')




myEntropy=[entropy_vn(timeEvo(dt, initSystem, Hges)) for dt in t]
myIdot = deriv(t[1]-t[0], myEntropy)



plt.plot(t[:-1], myIdot)
plt.show()

plt.close('all')


