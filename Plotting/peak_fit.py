# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 09:34:07 2024

@author: wbsst
Manual Peak fitting for low freq
"""
#imports
import numpy as np
from matplotlib import pyplot as plt

#opening experiamental data
infile = open(r"C:/Users/wbsst/Box/RugGroupShared/Raman Data/Will/MgTA_VT_10/78/data_normalized.txt", "r")
x_exp, y_exp = np.loadtxt(infile, unpack=True, delimiter=",")

#Peaks: frequency, intenisty, bandwidth, stick height
peaks = [ 
 [104,2.7,5,1],
 [119,2.5,2,1],
 [144,1.7,2,1],
 [197,3.8,4,1],
 [212.5,1.9,3.5,1]
 ]


#creating each list
y = []
stem_freqs = []
y = []
y_stem = []
x = np.linspace(0, 1000, 100000)


#create curve for each point
def function(F, A,B):
    for i in x:
        y = A * ( ( float(B) /np.pi ) / ((x - F)**2  + B**2  ) )
        return y


#Call function
for i in range(len(peaks)):
    y_temp = function(peaks[i][0],peaks[i][1],peaks[i][2])
    y_stem.append(peaks[i][3]/10)
    y.append(y_temp)
    stem_freqs.append(peaks[i][0])
y = sum(y)

#plot Theory
plt.plot(x, y, color='blue', label="Fit")
plt.stem(stem_freqs,y_stem, markerfmt=' ',linefmt = 'b',basefmt = ' ')
ax = plt.gca()
ax.set_xlim(50,250)
ax.set_ylim(0,0.4)
plt.title('MgTA')
plt.xlabel(r'Wavenumber (cm$^{-1}$)')
plt.ylabel('Intensity (Normalized Counts)')
#plt.grid(True)

#plot Experiamental
plt.plot(x_exp,y_exp-0.6, label="Raman 78", color='green')#, label="CuACAC 78K")

#plt.legend(["Theory","Experiamental"])
#show plot
plt.legend()
plt.savefig(r"C:\Users\wbsst\Desktop/temp.png", dpi=500)
plt.show()
