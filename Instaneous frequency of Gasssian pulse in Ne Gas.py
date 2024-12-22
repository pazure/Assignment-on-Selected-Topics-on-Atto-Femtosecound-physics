# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 13:56:23 2024

@author: azure
"""

import numpy as np 
from matplotlib import pyplot as plt


tau = 50e-15
t0 = tau/np.sqrt(2*np.log(2))
wavelength = 800e-9
w_0 = 2*np.pi*3e8/wavelength
I0 = 10e12# W/cm^2
E0 = np.sqrt((2*I0)/((3e8)*(8.8541e-12))) #Electric field amplitude as a function in V/m
#print(E0) # For self check
Lenght_Ne = 1 
n2 =0.85e-19#cm^2/W
t = np.linspace(-100e-15, 100e-15, num=2000)
normalise_time = t/tau


n0 =1 #Non linear refractive index of Ne
n= n0 +n2*I0 #Non linear refractive index of Ne 


#FWHM funtion to find the pulse duration of each distance 
def FWHM (x,y):
     max_y = max(y)  # Find the maximum y value
     half_max_y = max_y/2
     dataPoints = np.stack((x,y), axis=-1)
     condition = dataPoints[:,1] > half_max_y
     
     xpoints = np.extract(condition, t)
     return (max(xpoints)-min(xpoints))


#Calculating intensity as a function time from Equation given in selected topics lecture notes
I_t = I0*np.exp(-(2*np.log(2)*(t**2)/tau**2))

# Field dependent refractive index 
n_t = n0 +n2*I_t

# Input Electric field as a function of time 
E_t = E0*np.exp(-2*np.log(2)*((t)**2/tau**2))*np.exp(1j*w_0*t)
Et_conjugate = np.conjugate(E_t)


#Output Electric field
E0_t =np.sqrt((2*I_t)/((3e8)*(8.8541e-12)))
E_lt = np.exp(-1j*2*np.pi*n/wavelength)*E_t # Efield at distance l and t 
Intensity_out = abs(E_lt**2)



plt.figure(1)
plt.plot(t, E_t, label = 'Intput Electric Field ')
plt.xlabel('Time',fontsize = 'large')
plt.ylabel('Electric field', fontsize = 'large')
plt.legend(loc = 'upper right',fontsize = 'large')
plt.show()


#Ploting the intensities  for input and output pulse

plt.figure(2)
plt.plot(t, abs(E_t**2), label = ("Input Temporal Intensity for input pulse" )+ (f" metres with FWHM = {FWHM(t, abs(E_t**2)):.2e}"))
plt.xlabel('Time',fontsize = 'large')
plt.ylabel('Electric field', fontsize = 'large')
plt.legend(loc = 'upper right',fontsize = 'large')
plt.show()

plt.figure(4)
#e
plt.plot(t, I_t, label = (f"Output Temporal Intensity for: { 1:.2e}" )+ (f" metres with FWHM = {FWHM(t, I_t):.2e}"))
plt.xlabel('Time(s)',fontsize = 'large')
plt.ylabel('Intensity (W/m)', fontsize = 'large')
plt.legend(loc = 'upper right',fontsize = 'large')
plt.show()


#%%








# Phase Modulation as a function of time
phi_t = w_0*t -(I_t*n*Lenght_Ne)/wavelength
# Slope phase modulation (Instantaneous frequency)
Instant_frequency = w_0+(4*np.pi*n2*Lenght_Ne)/(wavelength*tau**2)*t*I0

# Spectral change in frequency 
omega_t = w_0+(4*np.pi*n2*Lenght_Ne)/(wavelength*tau**2)*t*I_t

#Input electric field in Spectral domain 
E_fourier = E0*np.sqrt(np.pi/(2*np.log(2)))*t0* np.exp(-(((omega_t-w_0)**2)/8*np.log(2))*t0**2)

#out put electric field in the spectrall domain
E_out_fourier =np.exp(-1j*2*np.pi*n/wavelength)*E_fourier


#Ploting  spectrall intensities 
plt.figure(5)
plt.plot(omega_t,E_out_fourier, label ='Spectrall Amplitude of Output Pulse')
plt.plot(omega_t,E_fourier, label = 'Spectral amplitude for input pulse')
plt.xlabel('Time(s)',fontsize = 'large')
plt.ylabel('Amplitude', fontsize = 'large')
plt.legend(loc = 'upper right',fontsize = 'large')
plt.show()













#%%

# Create figure and axis
fig, ax1 = plt.subplots()

# Plot on the first y-axis (left)
ax1.plot(t, I_t, 'g-')
ax1.set_xlabel('Time (fs)',fontsize='large')
ax1.set_ylabel('Numalise Electric field', color='g',fontsize='large')

# Create another axis sharing the same x-axis
ax2 = ax1.twinx()
ax2.plot(t,omega_t, 'b-')
ax2.set_ylabel(' Frequency modulation', color='b',fontsize='large')
plt.legend(loc = 'upper right',fontsize = 'large')
plt.title('Angular frenquency modualtion and Electric field ', fontsize='large')
plt.show()


