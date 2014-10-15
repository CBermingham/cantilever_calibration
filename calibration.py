import csv, math
import scipy.optimize as optimization
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import kv

def Sho_power_response(omega, omega_R, Q, A, B):
	return A + (B * omega_R**4) / ((omega**2 - omega_R**2)**2 + ((omega**2 * omega_R**2) / Q**2))

#def Sho_power_response(omega, omega_R, Q, B):
#	return (B * omega_R**4) / ((omega**2 - omega_R**2)**2 + ((omega**2 * omega_R**2) / Q**2))


width = 1E-6 # width / m
eta = 1.8205E-5 # dynamic viscosity of air / Pa.s
rho = 1.2047 # density of air / kg.m^-3
L = 50E-6 # length / m

f = open('test13.txt', 'rU')
lines=f.readlines()[1:]
f.close()

frequency = []
PSD = []
for l in lines:
 	b = l.split()
 	PSD.append(float(b[3]))
 	frequency.append(float(b[0]))

ang_freq = [i * 2 * math.pi for i in frequency]

#initial = [1000, 10, 1E-21]
initial = [1000, 0.04, 1E-21, 1E-21]
fit = optimization.curve_fit(Sho_power_response, ang_freq, PSD, initial)

fit_data = [Sho_power_response(f, fit[0][0], fit[0][1], fit[0][2], fit[0][3]) for f in ang_freq]
#fit_data = [Sho_power_response(f, fit[0][0], fit[0][1], fit[0][2]) for f in ang_freq]

res_ang_freq = fit[0][0]
quality_factor = fit[0][1]
white_noise_floor = fit[0][2]
Pdc = fit[0][3]
#Pdc = fit[0][2]
print 'Resonant frequency =', res_ang_freq / (2 * math.pi), 'Hz'
print 'Quality factor =', quality_factor
print 'White noise floor =', white_noise_floor, 'V^2 / Hz'
print 'DC power response of the cantilever measured from the photodetector =', Pdc, 'V^2 / Hz'

Re = (rho * res_ang_freq * width**2) / (4 * eta) # Reynold's number for the flow
tau = np.log10(Re) 
Tc = (1 + ((4 * 1j * kv(1, -1j * np.power(1j * Re, 0.5))) / (np.power(1j * Re, 0.5) * 
    kv(0, -1j * np.power(1j * Re, 0.5)))))
ohm = (0.91324 - 0.48274 * tau + 0.46842 * tau**2 - 0.12886 * tau**3 + 0.044055 * 
    tau**4 - 0.0035117 * tau**5 + 0.00069085 * tau**6) / (1 - 0.56964 * tau + 
    0.48690 * tau**2 - 0.13444 * tau**3 + 0.045155 * tau**4 - 0.0035862 * tau**5 + 
    0.00069085 * tau**6) + 1j * (-0.024134 - 0.029256 * tau + 0.016294 * 
    tau**2 - 0.00010961 * tau**3 + 0.000064577 * tau**4 - 0.000044510 * tau**5) / (1 
    - 0.59702 * tau + 0.55182 * tau**2 - 0.18357 * tau**3 + 0.079156 * tau**4 - 
    0.014369 * tau**5 + 0.0028361 * tau**6)
Hydrodynamic_function = ohm * Tc

k = 0.1906 * rho * width**2 * L * quality_factor * Hydrodynamic_function.imag * res_ang_freq**2
print 'Spring constant of cantilever =', k, 'N/m'

calibration = math.sqrt((2 * 1.3806488E-23 * 293.15)/(math.pi * k * res_ang_freq * 2 * math.pi * Pdc * quality_factor))
print 'Calibration =', calibration * 10E9, 'nm / V'

plt.scatter(frequency, PSD, color='g', s=1)
plt.plot(frequency, fit_data, color='r')
plt.xscale('log')
plt.yscale('log')
plt.ylim(10E-24, 10E-18)
plt.axvline(x=res_ang_freq / (2 * math.pi), ymin=0, ymax=1, linewidth=1, color='k')
plt.xlabel('Frequency / Hz')
plt.ylabel('PSD / V^2 / Hz')
plt.show()