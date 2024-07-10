# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 13:10:20 2024

@author: wbsst
"""

import os
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import subprocess
import matplotlib.ticker as ticker
import numpy as np


#Import Unpolarized Data
infile_UnPOL = open(r"300.asc", "r")
x_UnPOL, y_UnPOL = np.loadtxt(infile_UnPOL, unpack=True, delimiter="\t")


directory = os.path.join(os.getcwd(), "data")

angles = []

# Gather the angles
for item in os.listdir(directory):
    item_path = os.path.join(directory, item)
    if os.path.isdir(item_path):
        try:
            angles.append(int(item))
        except ValueError:
            # If the directory name is not an integer, skip it
            continue
angles.sort()

# Grab X points
x = []
x_infile_path = os.path.join(directory, str(angles[0]), "data.asc")
with open(x_infile_path, "r") as x_infile:
    for line in x_infile:
        x.append(round(float(line.split()[0])))

# Extract the Y Data
y = [[] for _ in range(len(angles))]
y1 = [[] for _ in range(len(angles))]
for i, angle in enumerate(angles):
    y_infile_path = os.path.join(directory, str(angle), "data.asc")
    with open(y_infile_path, "r") as y_infile:
        for line in y_infile:
            y[i].append(float(line.split()[1]))
            y1[i].append(float(line.split()[2]))



# Create the DataFrame for y data
data_y = {'x': x}
for i, angle in enumerate(angles):
    data_y[f'{angle}'] = y[i]

df_y = pd.DataFrame(data_y)

# Filter and transpose y DataFrame
df_y = df_y[(df_y['x'] >= 10) & (df_y['x'] <= 300)]
df_y.set_index('x', inplace=True)
df_y = df_y.T

# Create the DataFrame for y1 data
data_y1 = {'x': x}
for i, angle in enumerate(angles):
    data_y1[f'{angle}'] = y1[i]

df_y1 = pd.DataFrame(data_y1)

# Filter and transpose y1 DataFrame
df_y1 = df_y1[(df_y1['x'] >= 10) & (df_y1['x'] <= 300)]
df_y1.set_index('x', inplace=True)
df_y1 = df_y1.T

# Determine common vmin and vmax for color scale
vmin = min(df_y.min().min(), df_y1.min().min())
vmax = max(df_y.max().max(), df_y1.max().max())

# Create a figure with two subplots
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 14))

#EXP Data
ax3.plot(x_UnPOL, y_UnPOL, color='b', label='300K')
ax3.set_title('Unpolarized Data')
ax3.set_xlabel('Raman Shift (cm$^{-1}$)')
ax3.set_ylabel('Intensity (Counts)')
ax3.set_xlim(10, 300)
ax3.set_ylim(14000, 31000)
ax3.legend()

# Plot heatmap for y data
sns.heatmap(df_y, annot=False, cmap='magma', fmt="", ax=ax1, vmin=vmin, vmax=vmax,cbar=False)
ax1.set_title('TMA-AZO Polarized Raman - A')
ax1.set_xlabel('')
ax1.set_ylabel(r'Angle ($\degree$)')
ax1.set_xticks([])

# Plot heatmap for y1 data
sns.heatmap(df_y1, annot=False, cmap='magma', fmt="", ax=ax2, vmin=vmin, vmax=vmax,cbar = False)
ax2.set_title('TMA-AZO Polarized Raman - B')
ax2.set_xlabel('')
ax2.set_ylabel(r'Angle ($\degree$)')
ax2.set_xticks([])



# Show the plots
plt.tight_layout()
plt.show()

subprocess.run("cat goodbye.txt", shell = True, executable="/bin/bash")
