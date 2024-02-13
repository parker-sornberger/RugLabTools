# GetCartesianCoords
A container class which parses the fractional coordinates and cartesian coordinates from a CRYSTAL frequency output file. 
## Capabilities 
+ Easily get arrays of fractional or cartesian coordinates
+ Create supercells
+ Get all lattice parameters of a crystal
+ Write out the cartesian coordinates to an XYZ file
# ParseFreq
Composed with the GetCartesianCoords class which does most of the heavy lifing for getting geometries. 
## Capabilities:
+ Get all vibrational modes found in an output file
  1. Eigenvalue of vibration
  2. Frequency of vibration in cm<sup>-1</sup> and THz
  3. Intensity of vibration
  4. Eigenvectors of vibrations
+ Shift molecules in a cell by a specific vibrational mode
+ Make XYZ animations of vibrations to visualize with VMD
+ Find distances between atoms, angles, and dihedrals
