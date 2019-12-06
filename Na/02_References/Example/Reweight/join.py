import numpy as np

temp=350
N=9
deltaFreeEnergy=np.zeros(N)
for j in range(N):
	name="analysis." + str(int(j)) + ".fes" 
	data=np.genfromtxt(name)
	beta=1./(0.0083144621*temp)
	y=data[:,1]
	y=np.exp(-beta*y)
	# Solid Basin
	division=500
	integral1=np.trapz(y[:division])
	freeEnergyLeft=-(1./beta)*np.log(integral1)
	# Liquid Basin
	integral2=np.trapz(y[division:])
	freeEnergyRight=-(1./beta)*np.log(integral2)
	deltaFreeEnergy[j]=freeEnergyLeft-freeEnergyRight
print(temp,np.mean(deltaFreeEnergy),np.std(deltaFreeEnergy),np.std(deltaFreeEnergy)/np.sqrt(N),np.amax(deltaFreeEnergy),np.amin(deltaFreeEnergy))
