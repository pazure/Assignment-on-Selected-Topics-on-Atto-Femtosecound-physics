# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 13:56:23 2024

@author: azure
"""

import numpy as np 
from matplotlib import pyplot as plt

# Constants
c = 3e8  # Speed of light (m/s)
epsilon_0 = 8.8541e-12  # Vacuum permittivity (F/m)

# Pulse parameters
tau = 50e-15  # Pulse duration (s)
t0 = tau / np.sqrt(2 * np.log(2))
wavelength = 800e-9  # Wavelength (m)
w_0 = 2 * np.pi * c / wavelength  # Central angular frequency (rad/s)
I0 = 1e13  # Peak intensity (W/m²) converted from W/cm²
E0 = np.sqrt((2 * I0) / (c * epsilon_0))  # Electric field amplitude (V/m)
Length_Ne = 1  # Length of the medium (m)
n2 = 0.85e-19  # Nonlinear refractive index (cm²/W)
n2 = n2 * 1e-4  # Convert n2 to m²/W




# Time array
t = np.linspace(-200e-15, 200e-15, num=4000)  # Time array (s)
dt = t[1] - t[0]  # Time step (s)

# Linear refractive index of Ne
n0 = 1

# Spectral change in frequency 
omega_t0 = w_0+(4*np.pi*n2*1)/(wavelength*tau**2)*t*I0

# Input electric field as a function of time
E_t = E0 * np.exp(-2 * np.log(2) * (t**2 / tau**2)) * np.exp(1j * w_0 * t)
I_t = np.abs(E_t)**2  # Input intensity as a function of time

# Nonlinear refractive index as a function of time
n_t = n0 + n2 * I_t

# Nonlinear phase shift
phi_NL = - (2 * np.pi / wavelength) * (n_t - n0) * Length_Ne  # Only nonlinear part

# Total phase
phi_total = w_0 * t + phi_NL

# Instantaneous angular frequency
omega_t = np.gradient(phi_total, dt)  # Numerical derivative of phase

# Output electric field after passing through the medium
E_out_t = E0 * np.exp(-2 * np.log(2) * (t**2 / tau**2)) * np.exp(1j * phi_total)
I_out_t = np.abs(E_out_t)**2  # Output intensity as a function of time

# Prepare frequency array for Fourier Transform
N = len(t)
# Frequency array centered around w_0 with an appropriate range
delta_omega = 5e14  # Frequency range (rad/s)
omega = np.linspace(w_0 - delta_omega, w_0 + delta_omega, N)  # Angular frequency array (rad/s)




plt.figure(1)
plt.plot(t, E_t, label = 'Intput Electric Field ')
plt.xlabel('Time',fontsize = 'large')
plt.ylabel('Electric field', fontsize = 'large')
plt.legend(loc = 'upper right',fontsize = 'large')
plt.show()

# Vectorized Numerical Fourier Transform Function
def numerical_fourier_transform(E_t, t, omega):
    exponent = -1j * np.outer(omega, t)
    integrand = E_t * np.exp(exponent)  # Shape: (N_omega, N)
    E_omega = np.trapz(integrand, t, axis=1)  # Integrate over time for each frequency
    return E_omega


#FWHM funtion to find the pulse duration of each distance 
def FWHM (x,y):
     max_y = max(y)  # Find the maximum y value
     half_max_y = max_y/2
     dataPoints = np.stack((x,y), axis=-1)
     condition = dataPoints[:,1] > half_max_y
     
     xpoints = np.extract(condition, t)
     return (max(xpoints)-min(xpoints))



# Compute Fourier Transform of the input electric field
E_in_omega = numerical_fourier_transform(E_t, t, omega)
# Compute Fourier Transform of the output electric field
E_out_omega = numerical_fourier_transform(E_out_t, t, omega)

# Normalize the Fourier Transforms
E_in_omega /= np.max(np.abs(E_in_omega))
E_out_omega /= np.max(np.abs(E_out_omega))

# Plotting the spectral amplitudes of input and output pulses on the same graph
plt.figure(figsize=(12, 6))
plt.plot(omega / 1e15, np.abs(E_in_omega), label='Input Pulse Spectrum'+(f" has  FWHM = {FWHM(omega / 1e15, np.abs(E_in_omega)):.5e}"))
plt.plot(omega / 1e15, np.abs(E_out_omega), label='Output Pulse Spectrum'+(f" has  FWHM = {FWHM(omega / 1e15, np.abs(E_in_omega)):.5e}"))
plt.xlabel('Angular Frequency $\\omega$ (rad/fs)', fontsize='large')
plt.ylabel('Normalized Spectral Amplitude', fontsize='large')
plt.title('Spectral Amplitudes of Input and Output Pulses', fontsize='large')
plt.legend(fontsize='large')
plt.grid(True)
plt.show()

# Plotting the instantaneous angular frequency omega_t
# Create figure and axis

plt.figure(figsize=(12, 6))
plt.plot(t * 1e15, omega_t / 1e15, label=' Angular Frequency Modulation $\\omega_t$')
#plt.plot(t * 1e15, omega_t0/1e15, label= 'Instaneoud frequency')
plt.xlabel('Time (fs)', fontsize='large')
plt.ylabel('Angular Frequency $\\omega_t$ (rad/fs)', fontsize='large')
plt.title('Instantaneous Angular Frequency', fontsize='large')
plt.legend(fontsize='large')
plt.grid(True)
plt.show()

# Plotting the input and output electric field intensities in the time domain
plt.figure(figsize=(12, 6))
plt.plot(t * 1e15, I_t / np.max(I_t), label='Input Pulse Intensity' + (f" has  FWHM = {FWHM(t * 1e15, I_t / np.max(I_t)):.2e}"))
plt.plot(t * 1e15, I_out_t / np.max(I_out_t), label='Output Pulse Intensity'+ (f" has  FWHM = {FWHM(t * 1e15, I_out_t / np.max(I_out_t)):.2e}"))
plt.xlabel('Time (fs)', fontsize='large')
plt.ylabel('Normalized Intensity', fontsize='large')
plt.title('Temporal Intensity Profiles of Input and Output Pulses', fontsize='large')
plt.legend(fontsize='large')
plt.grid(True)
plt.show()

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

#####
# Plotting the Electric field of input and output pulses on the same graph
plt.figure(figsize=(12, 6))
plt.plot(t, E_t, label='Input Electric field')
plt.plot(t, E_out_t, label='Output Electric Field')
plt.xlabel('T(fs)', fontsize='large')
plt.ylabel('Electric field', fontsize='large')
plt.title('Electric field of Input and Output Pulses', fontsize='large')
plt.legend(fontsize='large')
plt.grid(True)
plt.show()
