# -*- coding: utf-8 -*-
"""
Created on Sun Nov  3 13:32:35 2024

@author: azure
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar

plt.close('all')

# Constants and Parameters
wl = 800e-9           # Wavelength (m)
c = 3e8               # Speed of light (m/s)
intensity = 4e18      # Intensity (W/m^2)
e = 1.66e-19          # Electron charge (C)
eo = 8.854e-12        # Permittivity of free space (F/m)
m_e = 9.109e-31       # Electron mass (kg)
Ip_eV = 20            # Ionization potential in eV
Ip_joules = Ip_eV * e # Ionization potential in Joules

# Derived Parameters
T = wl / c            # Period (s)
w = 2 * np.pi * c / wl # Angular frequency (rad/s)
E0 = np.sqrt((2 * intensity) / (c * eo)) # Electric field amplitude (V/m)
Up = (e * E0) ** 2 / (4 * m_e * w ** 2) # Ponderomotive energy (J)
Up_ev = Up / e        # Ponderomotive energy in eV

# Time and Electric Field Calculation
t = np.arange(-0*T, T, 0.02 * T) # Time array
Et = np.cos(w * t)                  # Electric field over time
#Numerically solving the Return time that leads to recombination 

ti = 0
position0 = 0.5*(np.cos(w*t)-np.cos(w*ti) +np.sin(w*ti)*w*(t-ti))
plt.figure(2)
plt.plot(t/T,Et, 'r-', label = 'Electric field')
# plt.plot(t/T,position0,'b--', label = (' Electron trajectory for t\' = 0'))
plt. title('Electric Field vs time(1/period)')
plt.xlabel('Time(1/period)')
plt.ylabel('Electric Field')
plt.legend()
plt.grid()

plt.figure(3)
# plt.plot(t/T,Et, 'r-', label = 'Electric field')

tb  = np.arange(0.0*T,  T, 0.08 * T)
for ti in tb:
    position = 0.5*(np.cos(w*t)-np.cos(w*ti) +np.sin(w*ti)*w*(t-ti))
    plt.plot(t/T,position,'--', label=f'Born time = {ti/T:.2f} T')
plt. title('Electric Field vs time(1/period)')
plt.xlabel('Time(1/period)')
plt.ylabel('Electric Field')
plt.legend()
plt.grid()




#Return kinectic energy

def energy ( input_phase, outPutphase):
    
  
    return 2*(np.sin(outPutphase)-np.sin(input_phase))**2

# def energy2 (ti_valid, t_valid,w):
     
#     return 2 *( np.cos(3*np*arcsin(2/np.pi))


def outPutphase( input_phase):
     
    return np.pi/2 - 3*np.arcsin(2/(np.pi)*input_phase-1)


input_phase  = np.linspace(0, 0.25, 100)*np.pi*2

phase_R = outPutphase(input_phase)

energy = energy(input_phase, phase_R)



plt.figure(7)
plt.plot(input_phase, phase_R)
plt.xlabel("Emiision Phase ")
plt.ylabel('recombiantion Phase ')
plt.title('Recombination phase  vs Emission Phase ')
plt.legend(loc = 'best')
plt.grid()


plt.figure(8)
plt.plot(phase_R/(w*T), energy, 'r-' ,label = 'Recombination  Time ')
plt.plot(input_phase/(w*T), energy, 'b-',label = 'Emission Time ')
# plt.plot(t/T,Et, 'g-', label = 'Electric field')
plt.legend()
plt.grid()
plt.xlabel("Time(Optical cycle) ")
plt.ylabel('Kinectic Energy(Up) ')
plt.title('Energy(Up) vs Recombinaiton Time and Emission time')
plt.show()


plt.figure(9, figsize=(6,8) )
plt.subplot(211)
plt.plot(phase_R/(w*T), energy, 'r-' ,label = 'Recombination  Time ')
plt.plot(input_phase/(w*T), energy, 'b-',label = 'Emission Time ')
# plt.plot(t/T,Et, 'g-', label = 'Electric field')
plt.legend()
plt.grid()
plt.xlabel("Time(Optical cycle) ")
plt.ylabel('Kinectic Energy(Up) ')
plt.title('Energy(Up) vs Recombinaiton Time and Emission time')
plt.subplot(212)
# plt.plot(phase_R/(w*T), energy, 'r-' ,label = 'Recombination  Time ')
# plt.plot(input_phase/(w*T), energy, 'b-',label = 'Emission Time ')
plt.plot(t/T,Et, 'g-', label = 'Electric field')
plt.xlabel("Time(Optical cycle) ")
plt.ylabel('Kinectic Energy(Up) ')
plt.legend()
plt.grid()
plt.legend(loc='best')
plt.tight_layout()
plt.show()


plt.figure(10)
plt.plot(input_phase/(w*T), phase_R/(w*T))
plt.xlabel("Ionization time(Optical cycle)  ")
plt.ylabel('Recombination Time(optical cycle) ')
plt.title('Recombination time   vs Ionization time ')
plt.legend(loc = 'best')
plt.grid()




# plt.figure(8)
# plt.plot(phase_R, energy)
# plt.xlabel("Recombination Phase ")
# plt.ylabel('Kinectic Energy ')

# # Plot Electric Field and Electron Trajectories
# fig, ax1 = plt.subplots()
# ax1.plot(t , Et, 'r-', label='Electric Field')
# ax1.set_xlabel('Time (1/Period)')
# ax1.set_ylabel('Electric Field (a.u.)', color='red')
# ax1.tick_params(axis='y', labelcolor='red')

# # Electron trajectories
# born_time = np.arange(0.0*T, 0.3 * T, 0.01 * T)
# positions = calculate_trajectory(born_time,w, t)

# # Plotting each trajectory
# ax2 = ax1.twinx()
# for i, position in enumerate(positions):
#     ax2.plot(t , position, linestyle='--', label=f'Born time = {born_time[i]/T:.2f} T')
# ax2.set_ylabel('Normalized Electron Position (a.u.)')
# ax2.legend(loc='upper right')

# plt.title('Electric Field and Electron Trajectories for Different Born Times')
# plt.grid()

# plt.show()

