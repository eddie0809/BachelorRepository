import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from numpy import linalg as LA
import scipy
import prepState
from timeEvo import *
import myPlots
from generateStates import *
from qutip import *
import partialtrace

# for using tex formatting and font in plots
#
#plt.rcParams.update({"text.usetex": True,}) 
#mpl.rcParams['text.latex.preamble'] = [r'\usepackage[utf8]{inputenc}\usepackage[T1]{fontenc}\usepackage{lmodern}\inputencoding{utf8}\usepackage{amsmath}\usepackage{amssymb}\usepackage{mathtools}']
#mpl.rcParams['font.family'] = ['serif']


#sigmaX = np.array([[0.,1.],[1.,0.]])
#sigmaY = np.array([[0.,-1j],[1j,0.]])
#sigmaZ = np.array([[1.,0.],[0.,-1.]])
#sigmaPlus = .5 * (sigmaX + 1j*sigmaY)
#sigmaMinus = .5 * (sigmaX - 1j*sigmaY)
kb = 8.617333262e-05


beta = 1 # 1/(kb*300) ----> choose a value for beta then calculate the temperature from that, avoids overflow error (kb is small)
#rho = np.exp(-beta * H0.diag())
#Z = np.sum(rho)
#rho = np.diag(rho / Z)
rho = (-beta * H0).expm()
Z = rho.tr()
rho = rho/Z
print(rho)
print()
asdf = 0.9 * myBasis["one", 0] + 0.1 * myBasis["zero", 0] + 0 * myStates["sigmaX", 0] + 0 * myStates["sigmaY", 0] + 0 * myStates["sigmaZ", 0]

print(rho.tr())
print()

#first = what(0,0,0,0,1)
#first = first#/first.tr()

print(asdf)
thermal = asdf.unit()
thermal = (rho * thermal).unit()
print()
print(thermal)
#pt = partialtrace.PartialTrace(N)
#print(pt.get_first_state(rho))
helpme = thermal.ptrace(0)
print()
print(helpme)
print()
#print(asdf.unit().ptrace(0))

plotname = "Fidelity, Thermal state"

t = np.linspace(0,np.pi,100)

myPlots.plotFidelity(thermal, t, plotname)
S=[]
for dt in t:
	S.append(entropy_vn(timeEvo(dt,thermal,Hint)))
print(np.mean(S))
plt.plot(t, S, '+')
plt.show()

"""
pt = partialtrace.PartialTrace(N)
#rho_0 = pt.get_first_state(rho)
for dt in t:
	#rho_0 = pt.get_first_state(timeEvo(dt, rho, Hint))
	rho_F = timeEvo(dt, rho, Hint)
	#rho_N = pt.get_last_state(rho_F)
	fidelityList.append(fidelity(rho, rho_F))
fidelityList = np.real(fidelityList)
print(max(fidelityList))

for dt in t:
	rho_F = timeEvo(dt, rho, Hint)
	rho_N = pt.get_last_state(rho_F)
	fidelityList.append(fidelity(rho_0, rho_N))
fidelityList = np.real(fidelityList)

plt.plot(t, fidelityList)#, '+')
plt.ylim(0,1)
plt.show()#block=False)
#plt.pause(5)
"""
plt.close('all')


