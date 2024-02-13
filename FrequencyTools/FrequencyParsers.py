import numpy as np
from functools import cached_property
import pandas as pd
from typing import Union
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.weight'] = 'bold'
__all__ = ["GetCartesianCoords", "ParseFreq", ] #->will export these if used as module

class GetCartesianCoords:
    def __init__(self, file:str, upper_coord_block:bool=False, supercell:Union[tuple,None]=None, sig_figs:Union[int,None] = None):
        """
        Container for fractional and cartesian coordinates from a Frequency output file from CRYSTAL.
        

        Parameters
        ----------
        file : str
            Name of CRYSTAL Frequency output file
        upper_coord_block : bool, optional
            CRYSTAL has two sections of coordinates, one for the asymmetric unit, the other for the full cell. The upper is for the asymmetric unit. The default is False.
        supercell : tuple or None, optional
            If a supercell should be used, specify it in terms of INTEGER units of a,b,c. The default is None.
        sig_figs : int, optional
            If your CIF file used for the CRYSTAL input has a different number of sig-figs, specify these here. Otherwise, in Mercury, there will be different angles and distances than what is calculated here. The default is None.



        """
        self.alpha = 1
        self.beta = 1
        self.gamma = 1
        self.a = 1
        self.b = 1
        self.c = 1
        self.sig_figs = sig_figs
        self.unique_ids = set()
        self.neighbors = {} #dict of dicts or list of dicts?
        self.upper_coord_block = upper_coord_block
        with open(file, "r") as file:
            self.lines = file.readlines()
        self._frac_coords = []
        self.atoms = []
        self.atoms = []
        self.pulled_coords = []
        self._parse()
        self._is_supercell = False
        if supercell:
            self._is_supercell = True
            self._supercell(supercell)
        cot = lambda x : 1/ np.tan(np.radians(x))
        csc = lambda x : 1 / np.sin(np.radians(x))
        cos = lambda x : np.cos(np.radians(x))
        sin = lambda x : np.sin(np.radians(x))

        self.trans_matrix = np.array([[self.a *sin(self.beta)* np.sqrt(1 - (cot(self.alpha)*cot(self.beta) - csc(self.alpha)*csc(self.beta)*cos(self.gamma) )**2), 0, 0],
                                      [self.a * csc(self.alpha) * cos(self.gamma) - self.a * cot(self.alpha) * cos(self.beta), self.b * sin(self.alpha), 0],
                                      [self.a * cos(self.beta), self.b*cos(self.alpha), self.c]
                                      ])

    def _X(self,x,y,z, a,b,c, alpha, beta, gamma, vol):
        alpha, beta, gamma = list(map(np.radians, (alpha, beta, gamma)))    
        t1 = x/a
        t2 = y*np.cos(gamma)/(x*np.sin(gamma))
        t3 = b * c * z * (( np.cos(alpha)*np.cos(gamma) - np.cos(beta) )/vol *np.sin(gamma))

        return t1 - t2 + t3
    def _Y(self,x,y,z, a,b,c, alpha, beta, gamma, vol):
        alpha, beta, gamma = list(map(np.radians, (alpha, beta, gamma)))
        t1 = y/(b * np.sin(gamma))
        t2 = a * c* z * (( np.cos(beta) * np.cos(gamma) - np.cos(alpha)  )  /(np.sin(alpha) * vol))
        return t1 + t2
        
        
    def _Z(self,x,y,z, a,b,c, alpha, beta, gamma, vol):
        alpha, beta, gamma = list(map(np.radians, (alpha, beta, gamma))) 
        return a * b * z * (np.sin(gamma)/vol)
    def _xyz_to_abc(self,*args):
        return self._X(*args), self._Y(*args), self._Z(*args)
    def cart_to_new_frac(self, coords):
        new_frac = []
        for coord in coords:
            newcoord = self._xyz_to_abc(*coord, self.a, self.b, self.c, self.alpha, self.beta, self.gamma, self.volume)
            new_frac.append(newcoord)
        return new_frac

    def _transform(self, vec): 
        return self.trans_matrix @ vec

    @cached_property 
    def cartesian_coords(self):
        if self._is_supercell:
            return np.array(list(map(self._transform, self._frac_coords)))
        return self.pulled_coords

    @cached_property
    def fractional_coordinates(self):
        return np.array(self._frac_coords)

    def _get_keywords(self):
        abc_keyword = "A              B              C           ALPHA      BETA       GAMMA"
        frac_coord_begin_keyword = "ATOM                 X/A                 Y/B                 Z/C"
        end_frac_coord_keyword = " T = ATOM BELONGING TO THE ASYMMETRIC UNIT"
        if self.upper_coord_block:
            end_frac_coord_keyword = "TRANSFORMATION MATRIX PRIMITIVE-CRYSTALLOGRAPHIC CELL"
        neighbor_block_begin_keyword = "    ATOM  N     R/ANG      R/AU   NEIGHBORS (ATOM LABELS AND CELL INDICES)"
        neighbor_end_keyword = " SYMMETRY ALLOWED INTERNAL DEGREE(S) OF FREEDOM"
        cart_start_keyword = "*      ATOM          X(ANGSTROM)         Y(ANGSTROM)         Z(ANGSTROM)"
        cart_start = 0
        cart_end_keyword = "LOCAL ATOMIC FUNCTIONS BASIS SET"
        cart_end = 0
        neighbor_block_begin = 0
        neighbor_block_end = 0
        
        abc_index = 0
        begin_frac_coord_index = 0
        end_frac_coord_index = 0

        for i, line in enumerate(self.lines):
            if neighbor_block_begin_keyword in line:
                neighbor_block_begin = i
            if neighbor_end_keyword in line:
                if neighbor_block_begin and neighbor_block_end > neighbor_block_begin:
                    continue
                neighbor_block_end = i
            if self.upper_coord_block and end_frac_coord_index:
                continue
            if abc_keyword in line:
                abc_index = i
                
            if frac_coord_begin_keyword in line:
                begin_frac_coord_index = i
            if end_frac_coord_keyword in line:

                end_frac_coord_index = i
            if cart_start_keyword in line:
                cart_start = i
            if cart_end_keyword in line:
                cart_end = i

        return (abc_index,  begin_frac_coord_index, end_frac_coord_index, 
    neighbor_block_begin, neighbor_block_end, cart_start, cart_end)
    def _get_abc_and_angles_lazy(self, abc_index):
        abc_angles = self.lines[abc_index+1]

        a,b,c,alpha,beta,gamma = list(map(float, abc_angles.split()))
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        vol_line = self.lines[abc_index-1].split()
        try:
            vol_ind = vol_line.index("VOLUME=")
            self.volume = float(vol_line[vol_ind+1])
        except Exception:
            pass
        
        

        
    def _get_fractional_coords(self, begin, end):
        
        frac_coords = self.lines[begin+2:end-1]

        self.e = frac_coords
        for line in frac_coords:

            index, asym, number, sym, *coords = line.strip().split()
            
            self.atoms.append(sym)

            if self.sig_figs is not None:
                _round = lambda x : float(F"{float(x):.{self.sig_figs}g}")
                self._frac_coords.append(list(map(_round, coords)))
                
            else:
                self._frac_coords.append(list(map(float, coords)))
    def _get_neighbors_and_unique_ids(self, neighbor_block_begin, neighbor_block_end):
        chunk = self.lines[neighbor_block_begin: neighbor_block_end]
        _id = 0
        for line in chunk:
            if not len(line.split()):
                continue
            if not line.split()[0].isdigit():
                continue
            contents = line.strip().split()
            #print(contents)
            _id = int(contents[0])-1
            if _id not in self.unique_ids:
                info = {}

                self.unique_ids.add(_id)
                symbol = contents[1]
                info = {"ID":_id, "symbol":symbol, 
                        "Neighbors R away":[] , "Dist (Ang)":[], "Dist (AU)":[],
                        "Neighbor label":[], "Neighbor ID":set(), "Supercell":[]}
            info["Neighbors R away"].append(int(contents[2]))
            info["Dist (Ang)"].append(float(contents[3]))
            info["Dist (AU)"].append(float(contents[4]))
            info["Neighbor ID"].add(int(contents[5])-1)
            info["Neighbor label"].append(F"{contents[5]}{contents[6]}")

            self.neighbors[_id] = info
    def _pull_cart_coords(self, cart_start, cart_end):
        chunk = self.lines[cart_start:cart_end]
        
        for line in chunk:
            if "*" in line or not line.strip():
                continue
            _id, a_num, symbol, *coords = line.strip().split()
            coords = list(map(float, coords))
            self.pulled_coords.append(coords)
        self.pulled_coords = np.array(self.pulled_coords)
        
            
        
    def _parse(self):
        abc_index,  begin_frac_coord_index,end_frac_coord_index,neighbor_block_begin, neighbor_block_end, cart_start, cart_end = self._get_keywords()

        self._get_abc_and_angles_lazy(abc_index)
        self._get_fractional_coords(begin_frac_coord_index, end_frac_coord_index)
        self._get_neighbors_and_unique_ids(neighbor_block_begin, neighbor_block_end)
        self._orig_atoms = self.atoms.copy()
        self._pull_cart_coords(cart_start, cart_end)


    def _supercell(self, supercell):

        self._frac_coords = np.array(self._frac_coords) 
        atoms = []
        assert len(supercell) == 3, "supercell must be 3D"
        a, b, c = supercell
        A = self._frac_coords[:,0].T
        B = self._frac_coords[:,1].T
        C = self._frac_coords[:,2].T
        i = j = k = 0
        coords = [self._frac_coords.copy()]

        for i in range(-a//2, a//2):


            for j in range(-b//2, b//2):


                for k in range( -c//2, c//2):

                    if i==j==k==0:
                        continue
                    coords.append(np.array([A+i, B+j, C+k]).T) 

        [atoms.extend(F"{atom}{i}" for atom in self.atoms) for i, _ in enumerate(range(a*b*c))]
        
        self.atoms = atoms

        self._frac_coords = np.vstack(coords)

    def write_xyz(self, filename = 'freqxyz',atoms = None, coords = None, IDs= None,labels = False, pulled_coords = False):
        num_atoms = len(self.atoms) if atoms is None else len(atoms)
        coords = self.cartesian_coords if coords is None else coords
        coords = self.pulled_coords if pulled_coords else coords 
        atoms = self.atoms if atoms is None else atoms
        if IDs:
            num_atoms = len(IDs)
            coords = self.cartesian_coords[list(IDs)]
            #print(coords)
            atoms = [atoms[_id] for _id in IDs]
        with open(F"{filename}.xyz", "w") as file:
            file.write(F"{num_atoms}\n\n")
            for i, (atom, coord) in enumerate(zip(atoms, coords )):
                if not labels:
                    file.write(F"{atom} {' '.join(str(c) for c in coord)}\n")
                else:
                    file.write(F"{atom} {' '.join(str(c) for c in coord)} {i+1}\n")
    

    

            
        
    
class ParseFreq:
    
    
    def __init__(self, file:str, upper_coord_block:bool=False, supercell:Union[tuple,None]=None, sig_figs:Union[int,None] = None):
        """
        Same __init__ params as GetCartesianCoords since ParseFreq is composed with this class. 

        Parameters
        ----------
        file : TYPE
            DESCRIPTION.
        supercell : TYPE, optional
            DESCRIPTION. The default is None.
        upper_coord_block : TYPE, optional
            DESCRIPTION. The default is False.
        sig_figs : TYPE, optional
            DESCRIPTION. The default is None.

        Key attribute of this class is "modes"
        The modes dictionary contains 'EigenValue', 'WaveNumber', 'THz', 'ActiveIR', 'Intensity', 'ActiveRaman', and 'EigenVector' as keys.
        Eigenvalue is the eigenvalue of the vibrational eigenvector --> float
        Wavenumber is the energy in cm^-1 --> float
        THz is the energy in THz --> float
        ActiveIR is a flag if a vibration is IR active --> bool
        Intensity if the intensity of this vibration --> float
        ActiveRaman is a flag if a vibration is Raman active --> bool
        EigenVector is an array of the eigenvector of a vibration for each atom in a unit cell --> np.ndarray
            For supercells, eigenvectors are simply vertically stacked. 
        

        """
        self.GetCoords = GetCartesianCoords(file, upper_coord_block=upper_coord_block, supercell=supercell, sig_figs = sig_figs)
        self.supercell = supercell
        self.lines = self.GetCoords.lines
        self.cart_coords = self.GetCoords.cartesian_coords
        self.unique_atoms = self.GetCoords.unique_ids
        self.num_atoms = len(self.GetCoords.atoms) if not supercell else len(self.GetCoords._orig_atoms)
        self.atoms = self.GetCoords.atoms if not supercell else self.GetCoords._orig_atoms
        self.modes = {}
        self.vectors = {}
        self._parse()
    def _get_intens_indices(self):
        begin_key = "(HARTREE**2)   (CM**-1)     (THZ)             (KM/MOL)"
        end_key = "NORMAL MODES NORMALIZED TO CLASSICAL AMPLITUDES (IN BOHR)"
        eigen_key = "FREQ(CM**-1)"
        cap = ' *******************************************************************************\n'
        raman_freq_end = "<RAMAN><RAMAN><RAMAN><RAMAN><RAMAN><RAMAN><RAMAN><RAMAN><RAMAN><RAMAN><RAMAN>"
        bindex = 0
        eindex = 0
        eigenindex = []
        is_raman_freq = False
        for i, line in enumerate(self.lines):
            if begin_key in line:
                bindex = i
            if (end_key in line or raman_freq_end in line) and not is_raman_freq:
                if raman_freq_end in line:
                    is_raman_freq = True
                eindex = i
            if eigen_key in line:
                eigenindex.append(i)

        return bindex, eindex, eigenindex
    def _parse_intense(self, bindex, eindex):
        block = self.lines[bindex+1: eindex-1]

        for line in block:
            _, mode,  eigv, wavnum, thz, irrep, ir, intens, raman = line.strip().replace("(","").replace(")","").split()
            self.modes[int(mode)] = dict(EigenValue=float(eigv), WaveNumber=float(wavnum), THz=float(thz), 
                                    ActiveIR=False if ir =="I" else True, Intensity=float(intens), 
                                    ActiveRaman=False if raman =="I" else True,EigenVector=None)
    def _parse_vectors(self, eigenindex):

        pairs = []
        for i in range(1,len(eigenindex)):
            pair = (eigenindex[i-1]+2, eigenindex[i]-1)         
            pairs.append(pair)
        mode_index = 1
        for upper_chunk, lower_chunk in pairs:
            vector_block = self.lines[upper_chunk:lower_chunk]
            vectors = []
            conv = lambda x : list(map(float, x))
            for line in vector_block:
                line = conv(line.strip().split()[4:]) if "AT" in line else conv(line.strip().split()[1:])
                vectors.append(line)
            vectors = np.array(vectors)
            for mode in range(vectors.shape[1]):
                self.modes[mode_index]["EigenVector"] = vectors[:, mode].copy().reshape(self.num_atoms, 3)
                mode_index+=1

    def _parse(self):
        bindex, eindex, eigenindex = self._get_intens_indices()
        self.e = eigenindex
        self._parse_intense(bindex, eindex)
        self._parse_vectors(eigenindex)
        
    def shift_molecule(self, mode, d=1, write = "", other_coords = None,eigen_slices = None, energy_cutoff=2000, mode_2 = None, d2=1):
        vector = self.modes[mode]["EigenVector"]

        atoms = None
        coords = None if other_coords is None else other_coords
        if vector is None or self.modes[mode]['WaveNumber'] > energy_cutoff :
            return 
        if self.supercell:
            copies = int(self.cart_coords.shape[0] // vector.shape[0])

            vector_copies = [vector for _ in range(copies-1)]
            vector = np.vstack([vector, *vector_copies])
            atoms = [self.GetCoords.atoms[i] for i in eigen_slices] if eigen_slices else []
        if eigen_slices is not None:
            vector = vector[np.ix_(tuple(eigen_slices))]
        new_mol = self.cart_coords + d * vector if other_coords is None else other_coords + d * vector
        if mode_2:
            vect2 = self.modes[mode_2]["EigenVector"] # we can't always say this exists...
            new_mol = new_mol + d2 * vect2
        if write:
            self.GetCoords.write_xyz(filename=write, atoms = atoms, coords = coords)
        return new_mol
    def make_animation(self, mode, max_d, steps=10, file_stem="step", other_coords = None, eigen_slices = None):
        dsteps = np.linspace(0, max_d, steps)
        dsteps = list(dsteps) + list(reversed(dsteps))
        if other_coords is not None:
            num_atoms = len(other_coords)
            atoms = [self.GetCoords.atoms[ind] for ind in eigen_slices]
        else:
            num_atoms = self.num_atoms if not self.supercell else len(self.GetCoords.atoms)
            atoms = self.atoms if not self.supercell else self.GetCoords.atoms
        mols = []
        for d in dsteps:
            mols.append(self.shift_molecule(mode, d, other_coords = other_coords,eigen_slices = eigen_slices ))
        with open(F"{file_stem}.xyz", "w") as file:
            file.write(F"{num_atoms }\n")
            for i, (step,mol) in enumerate(zip(dsteps, mols)):
                file.write(F"{step} TIMESTEP: {i}\n")
                for atom, coord in zip(atoms, mol):
                    file.write(F"{atom} {' '.join(str(c) for c in coord)}\n")
                file.write(F"{num_atoms}\n") if i+1 < len(mols) else None
    def get_dist_between_atoms(self, atomID1, atomID2, 
                               displace_along_mode:int=None, displacement_scale:float = 1, relative_change=False):
        atom1_coords = self.cart_coords[atomID1]
        atom2_coords = self.cart_coords[atomID2]
        dist = np.linalg.norm(atom1_coords-atom2_coords)
        if displace_along_mode is None:
            return dist
        shifted_mol = self.shift_molecule(displace_along_mode, d=displacement_scale)
        shift_atom1 = shifted_mol[atomID1]
        shift_atom2 = shifted_mol[atomID2]
        shift_dist = np.linalg.norm(shift_atom1-shift_atom2)
        change_percent = (abs(shift_dist-dist)/dist)*100
        if relative_change:
            return change_percent
        return abs(shift_dist-dist)
    
    def pi_pi_dist(self, atomID1, atomID2,
                               displace_along_mode:int=None, displacement_scale:float = 1, relative_change=False):

        one_mol_over = atomID2 + len(self.atoms)
        return self.min_dist(atomID1, one_mol_over, displace_along_mode=displace_along_mode, 
                                           displacement_scale=displacement_scale, relative_change=relative_change, do_not_search=True)
        
        
    def min_dist(self, atomID1, atomID2,
                 displace_along_mode:int=None, displacement_scale:float = 1, relative_change=False, do_not_search=False):
        """
        Use when one atom may be part of the asymmetric unit. Otherwise, use get_dist_between_atoms

        """
        eigen = None
        if self.supercell:
            atoms = self.GetCoords.atoms #make this change in other functions too 
        else:
            atoms = self.atoms
        if atomID1 not in self.unique_atoms and atomID2 not in self.unique_atoms or do_not_search:
            return self.get_dist_between_atoms(atomID1, atomID2, displace_along_mode=displace_along_mode, 
                                               displacement_scale=displacement_scale, relative_change=relative_change)
        unique_atom,other = (atomID1, atomID2) if atomID1 not in self.unique_atoms else (atomID2, atomID1)
        other_type = atoms[other]
        like_other = list(filter(lambda x : x[1]==other_type and x[0] != unique_atom, enumerate(atoms)))
        #print(like_other)
        atom1_coords = self.cart_coords[unique_atom]
        dist = lambda x_: np.linalg.norm(atom1_coords-self.cart_coords[x_[0]])
        neighbor = min(like_other, key = dist)[0]
        atom2_coords = self.cart_coords[neighbor]
        dist_ = np.linalg.norm(atom1_coords-atom2_coords)

        if displace_along_mode is None:
            return dist_
        shifted_mol = self.shift_molecule(displace_along_mode, d=displacement_scale, eigen_slices=eigen)
        shift_atom1 = shifted_mol[unique_atom]
        shift_atom2 = shifted_mol[neighbor]
        shift_dist = np.linalg.norm(shift_atom1-shift_atom2)
        change_percent = (abs(shift_dist-dist_)/dist_)*100
        if relative_change:
            return change_percent
        return abs(shift_dist-dist_)
        
        
    def get_angle_between_atoms(self, atomID1, atomID2, atomID3,
                                displace_along_mode:int=None, displacement_scale:float = 1,relative_change=False):
        atom1_coords = self.cart_coords[atomID1]
        atom2_coords = self.cart_coords[atomID2]
        atom3_coords = self.cart_coords[atomID3]
        
        ab = atom1_coords - atom2_coords 
        cb = atom3_coords - atom2_coords
        
        
        
        dot = ab @ cb
        
        mag1 = np.linalg.norm(ab)
        mag2 = np.linalg.norm(cb)
        return np.degrees(np.arccos(dot/(mag1*mag2)))
        
        
        
    def get_dihedral_angle(self, atomID1, atomID2, atomID3, atomID4,
                           displace_along_mode:int=None, displacement_scale:float = 1,relative_change=False):
        
        def get_coords(coords, *aids):
            return tuple(coords[aid] for aid in aids)


        def dihed(u1, u2, u3, u4):
            a1 = u2 - u1
            a2 = u3 - u2
            a3 = u4 - u3
        
            v1 = np.cross(a1, a2)
            v1 = v1 / (v1 * v1).sum(-1)**0.5
            v2 = np.cross(a2, a3)
            v2 = v2 / (v2 * v2).sum(-1)**0.5
            porm = np.sign((v1 * a3).sum(-1))
            rad = np.arccos((v1*v2).sum(-1) / ((v1**2).sum(-1) * (v2**2).sum(-1))**0.5)
            if not porm == 0:
                #print(4)
                rad = rad * porm
        
            return np.degrees(rad)


        u1, u2, u3, u4 = get_coords(self.cart_coords, atomID1, atomID2, atomID3, atomID4)
        ori_dihed = dihed(u1, u2, u3, u4)
        if displace_along_mode is None:
            return ori_dihed
        shifted_mol = self.shift_molecule(displace_along_mode, d=displacement_scale, )
        u1, u2, u3, u4 = get_coords(shifted_mol, atomID1, atomID2, atomID3, atomID4)
        shift_dihed = dihed(u1, u2, u3, u4)
        diff = abs(ori_dihed-shift_dihed)
        percent = (diff/ori_dihed) *100
        if relative_change:
            return percent
        return diff





def plot_mag_contributions(df:pd.DataFrame, rows:int = 10, name:str = ""):
    """
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame of x, y, and z contributions to the magnitude of an eigenvector. The columns should be x, y z and look something like:
                        x         y         z
            0    0.930328  0.000000  0.069672

    rows : int, optional
        How many eigenvectors should be plotted. The default is 10.
    name : TYPE, optional
        To save the plot, give it a name. Without one, no plot is saved. Saved as a PNG with 300 DPI. The default is "".


    """
    rcParams['xtick.bottom'] = rcParams['xtick.labelbottom'] = False
    rcParams['xtick.top'] = rcParams['xtick.labeltop'] = True

    colors = ["blue", "orange", "grey"]
    plt.figure(figsize = (rows//2 + 1 if rows %2 else rows //2 , 6))
    width = 0.4
    for index, row in df.iterrows():
        if index >= rows:
            break
        plt.bar(index+1, row[0], width = width, color = colors[0])
        plt.bar(index+1, row[1],width = width, bottom = row[0],  color= colors[1])
        plt.bar(index+1, row[2], width = width,bottom = row[0]+row[1],  color = colors[2])
    plt.xticks(range(1, rows+1), )
    plt.legend(["x", "y", "z"], bbox_to_anchor = (0.6,-0.1), ncol = 3, frameon=False)
    plt.xlim([0, rows+1])
    plt.tight_layout()
    if name:
        plt.savefig(F"{name}.png", dpi=300)
    
    plt.show()  

def mode_versus_relative_change(d:dict, rows:int = 10, name:str = ""):
    """
    

    Parameters
    ----------
    d : dict
        Dictionary with form {'Name of bond/dihedral/angle':{mode:dict(Quantity=quantity, DisplacedQuantity=disp_quantitiy, RelativeChange=rel_disp)  }} for all available modes.
    rows : int, optional
        How many modes to plot the relative change of. The default is 10.
    name : str, optional
        To save the plot, give it a name. Without one, no plot is saved. Saved as a PNG with 300 DPI. The default is "".



    """
    

    fs = 30
    rel_change = {key:[v["RelativeChange"] for k, v in val.items()] for key, val in d.items()}
    print(rel_change.keys())
    for _name, rcs in rel_change.items():
        if any(x < 0 for x in rcs[:rows]):
            plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
            plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True
        else:
            plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = True
            plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = False
        plt.figure(figsize =(rows//2 + 1 if rows %2 else rows //2 , 6) )
        plt.margins(x=0)
        modes = list(range(1,len(rcs)+1))
        plt.xticks(range(1, rows+1))

        plt.bar(modes[:rows], rcs[:rows])
        plt.title(_name,fontsize=fs+5)
        plt.ylabel("% Change", fontsize=fs)
        plt.xlabel("Mode",fontsize=fs)
        if name:
            plt.savefig(F"{name}_{_name}.png", dpi=300, bbox_inches = "tight")
        plt.show() 


   
      
if __name__ == "__main__":
    
    import os
    
    path = r"Directory of output file"
    
    file = "Your output file"
    file = os.path.join(path, file)

    C = GetCartesianCoords(F"{file}.out", upper_coord_block=0, )
    p = ParseFreq(F"{file}.out", supercell = (2,2,2), upper_coord_block=False,sig_figs =None )
   
    

    info = {"H129-C125-C105-C117":{}}
    for mode in p.modes:
        if any(item is None for item in p.modes[mode]) or p.modes[mode]["WaveNumber"]>2000:
            continue
        for key, dihed in zip(info.keys(),[(297,292,104,116)]):
            
            angle = p.get_dihedral_angle(*dihed)
            disp_angle = p.get_dihedral_angle(*dihed, displace_along_mode=mode)
            rel_disp = p.get_dihedral_angle(*dihed, displace_along_mode=mode, relative_change = 1)
            info[key][mode] = dict(Dihedral=angle, DisplacedDihedral=disp_angle, RelativeChange=rel_disp)
    dict_of_df = {k: pd.DataFrame(v).T for k,v in info.items()}
    df = pd.concat(dict_of_df, axis=1)
     

    modes = []
    
    for index, mode in p.modes.items():
        mode = mode["EigenVector"]
        if mode is None:
            continue
        x, y, z = np.sqrt(np.power(mode, 2).sum(axis=0))
        
        modes.append(dict(x=x, y=y, z=z))
    df = pd.DataFrame(modes)
    df = (df.T/df.T.sum()).T
    df.to_csv(F"{file}_all_modes.csv")
    
