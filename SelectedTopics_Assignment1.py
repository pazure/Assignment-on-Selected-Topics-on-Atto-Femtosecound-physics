# -*- coding: utf-8 -*-
"""
Created on Sat Sep 21 09:07:54 2024

@author: azure
"""

import numpy as np 
from matplotlib import pyplot as plt


# Defining the  paramaters

t = np.linspace(-100e-15, 100e-15, num=10000)
w_0 = (2*3.142*3E8)/(800e-9)
t_p = 50e-15
tG_0 = t_p/ np.sqrt(2*np.log(2))
group_delay = 3.61e-26



a = 0




#FWHM funtion to find the pulse duration of each distance 
def FWHM (x,y):
     max_y = max(y)  # Find the maximum y value
     half_max_y = max_y/2
     dataPoints = np.stack((x,y), axis=-1)
     condition = dataPoints[:,1] > half_max_y
     
     xpoints = np.extract(condition, t)
     return (max(xpoints)-min(xpoints))

# input Pulse
z_in =0
x = ((tG_0)**2)/(4*(1+a**2))

y = (a*x)-((3.61e-26)*z_in/2)
w = ((x**2) + (y**2))
tG_z = np.sqrt(4*(w)/x)


E_0 = 1

#Calculating the phase
phase = -(y)/(4*(w))*(t**2)
    

E_in = E_0*np.exp(-(t/tG_0)**2)*(np.exp(1j*phase))
E_conjugate = np.conjugate(E_in)
Intensity = abs(np.multiply(E_in,E_conjugate))

plt.figure(1)
plt.plot(t,Intensity, label = (f"Input Temporal Intensity for: { z_in:.2e}" )+ (f" metres with FWHM = {FWHM(t, Intensity):.2e}"))
plt.legend(loc = 'upper right',fontsize='large')
plt.xlabel('Time(fs)',fontsize='large')
plt.ylabel('Intensity',fontsize='large')
plt.show()

print(f"Scientific notation: {group_delay:.2e}")




    

# Output pulse 
plt.figure(2)
Distance = np.array([0,20e-3,30e-3,40e-3])
for z in Distance:
    y = (a*x)-((3.61e-26)*z/2)
    w = ((x**2) + (y**2))
    tG_z = np.sqrt(4*(w)/x)

    E_out = np.exp(-(1+(1j*y/x))*(t/tG_z)**2)
    E_outConjugate= np.conjugate(E_out)
    Intensity_out = abs(np.multiply(E_out,E_outConjugate))
    plt.plot(t,Intensity_out, label = (f"Output Temporal Intensity for: { z:.2e}" )+ (f" metres with FWHM = {FWHM(t, Intensity_out):.2e}"))
    #plt.plot(t, Intensity_out, label = 'Output Temporal Intensity for '  + str(z) + ' metres with FWHM at '+ str(FWHM(t, Intensity_out)))

plt.legend(loc = 'upper right',fontsize='large')
plt.xlabel('Time(fs)',fontsize='large')
plt.ylabel('Intensity',fontsize='large')
plt.show()


# Numerically Estimating the length for which the pulse duration is double that 
# of thee input pulse 
plt.figure(3)
Distance = np.arange(0,80e-3,5e-3 )
for z in Distance:
    y = (a*x)-((3.61e-26)*z/2)
    w = ((x**2) + (y**2))
    tG_z = np.sqrt(4*(w)/x)

    E_out = np.exp(-(1+(1j*y/x))*(t/tG_z)**2)
    E_outConjugate= np.conjugate(E_out)
    Intensity_out = abs(np.multiply(E_out,E_outConjugate))
    plt.plot(t,Intensity_out, label = (f"Output Temporal Intensity for: { z:.2e}" )+ (f" metres with FWHM = {FWHM(t, Intensity_out):.2e}"))
    #plt.plot(t, Intensity_out, label = 'Output Temporal Intensity for '  + str(z) + ' metres with FWHM at '+ str(FWHM(t, Intensity_out)))

plt.legend(loc = 'upper right',  fontsize='large')
plt.xlabel('Time(fs)',fontsize='large')
plt.ylabel('Intensity',fontsize='large')
plt.show()






