from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.optimize import fsolve

import sys
import csv
import numpy as np
import math
import matplotlib.pyplot as plt

file_name = sys.argv[1]
x     = []
u     = []
u_ref = None
BOLTZMAN_CONST = 1.380649E-023

with open(file_name, 'r') as specie_file:
	reader = csv.DictReader(specie_file, delimiter=",")
	for row in reader:
		u.append(float(row['Viscosity (Pa*s)']))
		x.append(float(row['Temperature (K)']))

vhs_tref = float(input("Enter the desired reference temperature (K): "))
m        = float(input("Enter the molecular mass (Dalton): "))
t_lower  = float(input("Enter the fit lower temperature bound: "))
t_upper  = float(input("Enter the fit upper temperature bound: "))

m = m * 1.6605390666E-027
f2 = interp1d(x, u, kind='cubic')
u_ref = f2(vhs_tref)

def func(x, w):
   return u_ref * np.power(np.divide(x, vhs_tref), w)

def func2(temps, d, w):
	return (15.0 * np.sqrt(np.multiply( (1/math.pi) * m * BOLTZMAN_CONST, x))) / (d * d * 2 * (5-2 * w) * (7-2 * w))

bounds=np.where((np.array(x) <= t_upper) & (np.array(x) >= t_lower))[0]
w, pcov = curve_fit(func, x[bounds[0]:bounds[-1]], u[bounds[0]:bounds[-1]], bounds=(0.0,1.0))
d_ref = math.sqrt((15 * math.sqrt(m * BOLTZMAN_CONST * vhs_tref / math.pi)) / (2*(5-2*w)*(7-2*w)*u_ref))

print('======= Result: ========')
print("reference diameter (d_ref): " + str(d_ref))
print("temperature exponent (omega): " + str(w))

plt.figure()
plt.plot(x,u, label='reference viscosity data')
plt.plot(x,func(x, w), label='DSMC VHS model curve fit')
plt.legend()
plt.grid(True, which='both', axis='both')
plt.xlabel("Viscosity Pa*s")
plt.xlabel("Temperature K")
plt.savefig('vhs_' + file_name + '_curvefit.png')
plt.close()
