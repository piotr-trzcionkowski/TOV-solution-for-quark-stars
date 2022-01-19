#!/usr/bin/python3.6

import numpy as np
from math import pi
from scipy.integrate import odeint
f=open("M(R).dat", "w")
f1=open("ro(R)1.dat", "w")
f2=open("ro(R)2.dat", "w")
f3=open("ro(R)3.dat", "w")

#stale fizyczne
c = 299792458. #m/s
G = 6.6743015 * 10** (-11) #m**3/kg/s**2

#warunek początkowy masy maksymalnej
Mmax=0

# model 
def model(z, P, n):
	r, m = z
	ro=9. / 14. * 10 ** 18 + 3 * P / (c ** 2)
	drdP = -1 * r * (c ** 2) * (r * (c ** 2) - 2 * G * m) / (G * (ro * (c **2) + P) * (m * (c ** 2) + 4 * pi * P * (r ** 3) ) )
	dmdP = drdP * 4 * pi * (r ** 2) * ro
	if (n == 28):
		f1.write(str(r / 1000) + '\t' + str(ro / 10 ** 18) + '\n')
	elif (n == 32):
		f2.write(str(r / 1000) + '\t' + str(ro / 10 ** 18) + '\n')
	elif (n == 36):
		f3.write(str(r / 1000) + '\t' + str(ro / 10 ** 18) + '\n')
	return [drdP, dmdP]

# solve ODE
for i in np.linspace(24, 40, 201):
	PC = 10. ** i
	z0 = [0.01, 16 / 3 * pi * 0.01 ** 3 * 9. / 14. * 10 ** 18 + 3 * PC / (c ** 2)]
	# span for next time step
	Pspan = np.linspace(PC,0, 10000)
	# solve for next step
	z = odeint(model, z0, Pspan, args=(i,))
	if (Mmax<(z[-1][1]/(2*10**30))):
		Mmax=z[-1][1]/(2*10**30)
	f.write(str(z[-1][0]/1000) + '\t' + str(z[-1][1]/(2*10**30)) + '\n')

print('Mmax = ', Mmax)
f.close()
f1.close()
f2.close()
f3.close()