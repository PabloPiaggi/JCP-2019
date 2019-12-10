import numpy as np

temp=np.linspace(350,450,11) # K
press=np.linspace(0,10000,11)  # bar

print("# vim:ft=plumed")
print("")
print("energy: READ FILE=COLVARtrim VALUES=energy  IGNORE_TIME")
print("b1: READ FILE=COLVARtrim VALUES=b1.bias  IGNORE_TIME")
print("uwall: READ FILE=COLVARtrim VALUES=uwall.bias IGNORE_TIME")
print("refcv: READ FILE=COLVARtrim VALUES=refcv.morethan  IGNORE_TIME")
print("voltmp: READ FILE=COLVARtrim VALUES=vol  IGNORE_TIME")
print("")
print("vol: COMBINE ARG=voltmp PARAMETERS=10.0 PERIODIC=NO")
print("renergy: COMBINE ARG=energy PARAMETERS=-25000 PERIODIC=NO")
print("weights: REWEIGHT_BIAS TEMP=400. ARG=b1.bias,uwall.bias")
print("")

for i in range(temp.shape[0]):
	for j in range(press.shape[0]):
		name=str(int(temp[i])) + "-" + str(int(press[j]))
		print("weights" + name  + ": REWEIGHT_TEMP_PRESS TEMP=400 PRESSURE=301.10704285 REWEIGHT_TEMP=" + str(temp[i]) + " REWEIGHT_PRESSURE=" + str(press[j]*0.06022140857) + " ENERGY=renergy VOLUME=vol")
		#print("weights" + name  + ": REWEIGHT_TEMP TEMP=375 REWEIGHT_TEMP=" + str(temp[i]) + " ARG=renergy")
		print("HISTOGRAM ...")
		print("  ARG=refcv.morethan")
		print("  GRID_MIN=-0.5")
		print("  GRID_MAX=250.5")
		print("  GRID_BIN=1000")
		print("  BANDWIDTH=1.")
		print("  LOGWEIGHTS=weights,weights" + name)
		print("  NORMALIZATION=true")
		print("  LABEL=hh"+ name)
		print("  CLEAR=10000")
		print("... HISTOGRAM")
		print("")
		print("DUMPGRID GRID=hh" + name + " FILE=histo" + name + " STRIDE=10000 FMT=%16.10f")
		print("ff" + name + ": CONVERT_TO_FES GRID=hh" + name + " TEMP=" + str(temp[i]) )
		print("DUMPGRID GRID=ff" + name + " FILE=fes" + name + " STRIDE=10000 FMT=%16.10f")
		print("")

print("")

print("PRINT ARG=* FILE=COLVAR STRIDE=1 FMT=%16.10f")
