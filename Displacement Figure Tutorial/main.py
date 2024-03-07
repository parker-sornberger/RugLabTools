# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 17:23:34 2024

@author: wbsst + MTR
"""
import os
import re
import numpy as np

loc = os.getcwd()

print(loc)

#READ IN CIF CONTAINING DIMER
infile0 = open("CIF_FROM_MERCURY.cif", "r")
outfile = open("dimer.txt",'w')
infile = infile0.readlines()
for i in range(37,len(infile)-2):
    infile[i]=infile[i].split()
    outfile.write(infile[i][0] + " " + infile[i][2] + " " + infile[i][3] + " " + infile[i][4] + "\n")
infile0.close()
outfile.close()


CRYSTALOUT = 'output.txt'
Supercell = [3, 3, 3]
dispamp = 10
NumberOfSteps = 10
# atoms_dimer = 82
# atoms_fragment = int(atoms_dimer/2)


# dimer_label_xyz = 'labeled_dimer.xyz'
# f1_label_xyz = 'labeled_f1.xyz'
# f2_label_xyz = 'labeled_f2.xyz'



MovieFull = True



ImagFreq = 0
FreqScale = 1



for i in range(4,5,1):
# for i in range(8,9,1):
    ModeNumber = 15
    
    atomicpos = open('atompos.txt', 'w')
    
    with open(CRYSTALOUT, 'r') as output:
        CellTag=False
        PositionTag = False
        for line in output:
            if CellTag==True:
                #cellparam = '26.09087278     7.07930800    13.55956828    90.000000  90.000000  90.000000'
                cellparam = line
                #print(line)
                CellTag = False
            ###This grabs atomic positions
            if PositionTag==True:
                if not line.strip():
                    PositionTag=False
    
                else:
                    #print(atomicpos, line)
                    atomicpos.write(line)
                    splity = line.split()
                    atomspercell = splity[0]
    
            #if 'CRYSTALLOGRAPHIC CELL (VOLUME= ' in line:
            if ' PRIMITIVE CELL - CENTRING CODE' in line:
                output.readline()
                CellTag = True
    
            if '      1 T' in line:
                if '1 THZ    ' in line:
                    pass
                else:
                    PositionTag=True
                    #print(atomicpos, line)
                    atomicpos.write(line)
    atomicpos.close()
    
    cellsplit = cellparam.split()
    #print(cellsplit)
    ###UNIT CELL VALUES####
    a=float(cellsplit[0])
    b=float(cellsplit[1])
    c=float(cellsplit[2])
    alp=np.radians(float(cellsplit[3]))
    beta=np.radians(float(cellsplit[4]))
    gam=np.radians(float(cellsplit[5]))
    print(a)
    print(b)
    print(c)
    print(alp)
    print(beta)
    print(gam)
    ##########################
    
    v = a*b*c*np.sqrt(1-np.cos(alp)*np.cos(alp)-np.cos(beta)*np.cos(beta)-np.cos(gam)*np.cos(gam)+2*np.cos(alp)*np.cos(beta)*np.cos(gam))
    xv = np.array([a, b*np.cos(gam), c*np.cos(beta)])
    yv = np.array([0, b*np.sin(gam), c*((np.cos(alp)-(np.cos(beta)*np.cos(gam)))/np.sin(gam))])
    zv = np.array([[0,0,v/(a*b*np.sin(gam))]])
    
    
    ###MAKE XYZ and SUPERCELL
    
    supercount = 0
    superpos = open('supercell.txt', 'w')
    superxyz = open('supercell_xyz.xyz', 'w')
    #print(atomspercell)
    superxyz.write(str(int(atomspercell)*(int(Supercell[0]))*(int(Supercell[1]))*(int(Supercell[2]))))
    superxyz.write('\n')
    
    #print >> superxyz, atomspercell
    #print >> superxyz, 'here is an xyz'
    superxyz.write('sample text')
    superxyz.write('\n')
    
    
    for x in range(0, Supercell[0]):
        for y in range(0, Supercell[1]):
            for z in range(0, Supercell[2]):
                with open('atompos.txt', 'r') as pos:
    
                    for line in pos:
                        stripped = line.split()
                        adx = stripped[0]+stripped[3]
    
                        #Coord in Fractional
                        coord = np.array([x+float(stripped[4]), y+float(stripped[5]), z+float(stripped[6])])
                        #Convert to Cartesian
                        a1 = np.array([[a, b*np.cos(gam), c*np.cos(beta)] , [0, b*np.sin(gam), c*((np.cos(alp)-(np.cos(beta)*np.cos(gam)))/np.sin(gam))] , [0,0,v/(a*b*np.sin(gam))]])
                        cart = np.dot(a1, coord)
                        # print cart
                        posstring = str(cart).lstrip('[').rstrip(']')
                        #Make XYZ and a temp file maintaining indicies
                        #print >> superpos, adx + ' ' + posstring
                        #superpos.write(atomspercell)
                        #superpos.write('Sample Text')
                        superpos.write(adx + ' ' + posstring + '\n')
                    
                        superxyz.write(stripped[3] + ' ' + posstring)
                        superxyz.write('\n')
                        #print >> superxyz, stripped[3] + ' ' + posstring
                        supercount+=1
    #print(coord)
    superpos.close()
    superxyz.close()
    #print(supercount)
    
    
    
    
    EigList = []
    FreqTag1 = False
    FreqTag2 = False
    LabelList = []
    EigFile = open('Eigs.txt', 'w')
    
    ModeCount = 0
    
    with open(CRYSTALOUT, 'r') as freqout:
        for line in freqout:
            if FreqTag2 == True:
                #CLEAR ALL VARIABLES HERE, write to eigenvector file, start over.
    
                trimmedEigs = []
                for item in EigList:
                    splititem = item.split()
                    #TrimmedEigs are the 6x1 eigenvectors (columns correspond to mode)
                    trimmedEigs.append(splititem[-6:])
                    if 'AT.' in splititem:
                        alabel = splititem[1]+splititem[2]
                        LabelList.append(alabel)
    
                for idx, freq in enumerate(freqlist):
                    #print 'Index: ' + str(idx) + ' Freq: ' + freq
                    if 'FreqScale' in locals():
                        #print(freq)
                        freqs = FreqScale*float(freq)
                        #EigFile.write('\n')
                        EigFile.write('----MODE: ' + str(ModeCount*6+idx+1-ImagFreq) + ' ---FREQ: ' + str(freqs) + ' ---')
                        EigFile.write('\n')
                        #print >> EigFile, '----MODE: ' + str(ModeCount*6+idx+1-ImagFreq) + ' ---FREQ: ' + str(freqs) + ' ---'
                    else:
                        #print >> EigFile, '----MODE: ' + str(ModeCount*6+idx+1-ImagFreq) + ' ---FREQ: ' + freq + ' ---'
                        EigFile.write('----MODE: ' + str(ModeCount*6+idx+1-ImagFreq) + ' ---FREQ: ' + freq + ' ---' + '\n')
                    for idx2, atom in enumerate(LabelList):
                        eigrowx = trimmedEigs[idx2*3]
                        eigrowy = trimmedEigs[idx2*3+1]
                        eigrowz = trimmedEigs[idx2*3+2]
                        #eigrowx = str(float(trimmedEigs[idx2*3])*1.5)
                        #eigrowy = str(float(trimmedEigs[idx2*3+1])*1.5)
                        #eigrowz = str(float(trimmedEigs[idx2*3+2])*1.5)
                        EigFile.write(atom + ' ' + str(1.5*float(eigrowx[idx])) + ' ' + str(1.5*float(eigrowy[idx])) + ' ' + str(1.5*float(eigrowz[idx])))
                        EigFile.write('\n')
                        #print >> EigFile, atom + ' ' + str(1.5*float(eigrowx[idx])) + ' ' + str(1.5*float(eigrowy[idx])) + ' ' + str(1.5*float(eigrowz[idx]))
                        #print >> EigFile, atom + ' ' + eigrowx[idx] + ' ' + eigrowy[idx] + ' ' + eigrowz[idx]
    
                ModeCount += 1
    
    
                LabelList = []
                EigList = []
                FreqTag1 = True
                FreqTag2 = False
    
                freqline = line.split()
                freqlist = freqline[1:7]
                if 'ESTIMATED NUMBER OF IMAGINARY FREQS' in line:
                    break
                if '*********' in line:
                    #print(line)
                    break
            if FreqTag1 == True:
                if 'FREQ(CM**-1) ' in line:
                    #print line
                    freqline = line.split()
                    freqlist = freqline[1:7]
    
                    freqout.readline()
                    continue
                if not line.strip():
                    FreqTag1 = False
                    FreqTag2 = True
                else:
                    #print line
                    EigList.append(line)
                    #print line
            if 'NORMAL MODES NORMALIZED TO CLASSICAL AMPLITUDES' in line:
                FreqTag1 = True
                freqout.readline()
    
    EigFile.close()
    #print(ModeNumber)
    ModeFlag = False
    modedict = {}
    with open('Eigs.txt', 'r') as Eigs:
         for line in Eigs:
             if ModeFlag==True:
                 if '----MODE:' in line:
                     break
                 splity = line.split()
                 label = splity[0]
                 #print([label])
                 xyzlist = [float(splity[1]), float(splity[2]), float(splity[3])]
                 #\print(xyzlist)
                 
                 modedict[label]=xyzlist
             if '----MODE:' in line:
                 splity = line.split()
                 if splity[1]==str(ModeNumber):
                     ModeFlag=True
    
    
    
    
    ##############################################################################################
    ##############################################################################################
    ##############################################################################################
    ############################ Dimer Projection File Generation ################################
    ##############################################################################################
    ##############################################################################################
    ##############################################################################################
    ##############################################################################################
    
    # def displace_frag(label_xyz,num_atoms_xyz,sign,label):
    #     if sign == 1:
    #         sign_label = 'p'
    #     if sign == -1:
    #         sign_label = 'n'
    #     test_out = open(str(label)+'_'+str(ModeNumber)+'_'+str(sign_label)+'disp%'+str(dispamp)+'%'+'.xyz','w')
    #     test_out.write(str(num_atoms_xyz)+'\n')
    #     test_out.write('Displaced mode '+str(ModeNumber)+'\n')
    #     counta = 0  
    #     with open(str(label_xyz), 'r') as pos:
    #              for line in pos:
    #                  split = line.split()
    #                  label1 = split[0]
    #                  atom = split[1]
    #                  Xi = float(split[2])
    #                  Yi = float(split[3])
    #                  Zi = float(split[4])
    #                  displist = modedict[label1]
    #                  Xf = Xi + int(sign)*dispamp*displist[0]
    #                  Yf = Yi + int(sign)*dispamp*displist[1]
    #                  Zf = Zi + int(sign)*dispamp*displist[2]
    #                  counta += 1
    #                  test_out.write(str(atom)+' '+ str(Xf) + ' ' + str(Yf) + ' ' + str(Zf) +'\n')
    #     test_out.close()
        
        
    # # dimer positive, negative
    # displace_frag(dimer_label_xyz,atoms_dimer,1,'ab')
    # displace_frag(dimer_label_xyz,atoms_dimer,-1,'ab')
    # # f1 positive, negative
    # displace_frag(f1_label_xyz,atoms_fragment,1,'a')
    # displace_frag(f1_label_xyz,atoms_fragment,-1,'a')
    # # f2 positive, negative
    # displace_frag(f2_label_xyz,atoms_fragment,1,'b')
    # displace_frag(f2_label_xyz,atoms_fragment,-1,'b')
    
    ##############################################################################################
    ##############################################################################################
    ##############################################################################################
    ##############################################################################################
    
    
    
    
    # def moviemaker(modefreq):
    #     FreqNumber = 109
    #     ModeFlag = False
    #     modedict = {}
    #     with open('Eigs.txt', 'r') as Eigs:
    #         for line in Eigs:
    #             if ModeFlag==True:
    #                 if '----MODE:' in line:
    #                     break
    #                 splity = line.split()
    #                 label = splity[0]
    #                 xyzlist = [float(splity[1]), float(splity[2]), float(splity[3])]
    #                 modedict[label]=xyzlist
    #             if '----MODE:' in line:
    #                 #print('I am working')
    #                 splity = line.split()
    #                 if splity[1]==str(modefreq):
    #                     FreqNumber = int(float(splity[3]))
    #                     ModeFlag=True
    
    #     totalatoms = str(int(atomspercell)*int(Supercell[0])*int(Supercell[1])*int(Supercell[2]))
    #     #print(totalatoms)
    #     totalatoms = supercount
    #     stepsize = dispamp/NumberOfSteps
    #     reversesteps = NumberOfSteps*-1
    #     totalnewatoms = atoms_dimer
    #     fname = 'NEW_dimer_mode_' + str(modefreq) + '_Freq_' + str(FreqNumber) + '.xyz'
    
    #     MOVIExyz = open(fname, 'w')
    
    #     for x in range(0, NumberOfSteps,1):
    #         dispamnt = x*stepsize
    #         MOVIExyz.write(str(totalnewatoms))
    #         MOVIExyz.write('\n')
    #         MOVIExyz.write('Displacement Amount: ' + str(dispamnt))
    #         MOVIExyz.write('\n')
    #         #print >> MOVIExyz, totalatoms
    #         #print >> MOVIExyz, 'Displacement Amount: ' + str(dispamnt)
    #         #with open('supercell.txt', 'r') as pos:
    #         with open(dimer_label_xyz,'r') as pos:
    #             for line in pos:
    #                 split = line.split()
    #                 label1 = split[0]
    #                 #print(label1)
    #                 #Xi = float(split[1])
    #                 #Yi = float(split[2])
    #                 #Zi = float(split[3])
    #                 Xi = float(split[2])
    #                 Yi = float(split[3])
    #                 Zi = float(split[4])
    #                 #print(label1)
    #                 displist = modedict[label1]
    #                 #print(label1)
    #                 Xf = Xi + dispamnt*displist[0]
    #                 Yf = Yi + dispamnt*displist[1]
    #                 Zf = Zi + dispamnt*displist[2]
    #                 atomsymbol = re.match(r"([0-9]+)([a-z]+)", label1, re.I)
    #                 if atomsymbol:
    #                     asymb = atomsymbol.groups()
    #                     asymb1 = asymb[1]
    #                     #print asymb1
    #                 #print >> MOVIExyz, asymb1 + ' ' + str(Xf) + ' ' + str(Yf) + ' ' + str(Zf)
    #                 #MOVIExyz.write('\n')
    #                 MOVIExyz.write(asymb1 + ' ' + str(Xf) + ' ' + str(Yf) + ' ' + str(Zf))
    #                 MOVIExyz.write('\n')
    #     for x in range(NumberOfSteps, reversesteps, -1):
    #         #print x
    #         dispamnt = x*stepsize
    #         #MOVIExyz.write('\n')
    #         MOVIExyz.write(str(totalnewatoms))
    #         MOVIExyz.write('\n')
    #         #print >> MOVIExyz, totalatoms
    #         MOVIExyz.write('Displacement Amount: ' + str(dispamnt))
    #         MOVIExyz.write('\n')
    #         #print >> MOVIExyz, 'Displacement Amount: ' + str(dispamnt)
    #         #with open('supercell.txt', 'r') as pos:
    #         with open(dimer_label_xyz,'r') as pos:
    #             for line in pos:
    #                 split = line.split()
    #                 label1 = split[0]
    # ##                Xi = float(split[1])
    # ##                Yi = float(split[2])
    # ##                Zi = float(split[3])
    #                 Xi = float(split[2])
    #                 Yi = float(split[3])
    #                 Zi = float(split[4])
    
    #                 displist = modedict[label1]
    #                 Xf = Xi + dispamnt*displist[0]
    #                 Yf = Yi + dispamnt*displist[1]
    #                 Zf = Zi + dispamnt*displist[2]
    #                 atomsymbol = re.match(r"([0-9]+)([a-z]+)", label1, re.I)
    #                 if atomsymbol:
    #                     asymb = atomsymbol.groups()
    #                     asymb1 = asymb[1]
    #                     #print asymb1
    #                # print >> MOVIExyz, asymb1 + ' ' + str(Xf) + ' ' + str(Yf) + ' ' + str(Zf)
    #                 #MOVIExyz.write('\n')
    #                 MOVIExyz.write(asymb1 + ' ' + str(Xf) + ' ' + str(Yf) + ' ' + str(Zf))
    #                 MOVIExyz.write('\n')
    #         #with open('DNTT_disp_new.xyz', 'w') as disp:
    #     for x in range(reversesteps, 0):
    #         dispamnt = x*stepsize
    #         #MOVIExyz.write('\n')
    #         MOVIExyz.write(str(totalnewatoms))
    #         MOVIExyz.write('\n')
    #         MOVIExyz.write('Displacement Amount: ' + str(dispamnt))
    #         MOVIExyz.write('\n')
    #         #print >> MOVIExyz, totalatoms
    #         #print >> MOVIExyz, 'Displacement Amount: ' + str(dispamnt)
    #         #with open('supercell.txt', 'r') as pos:
    #         with open(dimer_label_xyz,'r') as pos:
    #             for line in pos:
    #                 split = line.split()
    #                 label1 = split[0]
    # ##                Xi = float(split[1])
    # ##                Yi = float(split[2])
    # ##                Zi = float(split[3])
    #                 Xi = float(split[2])
    #                 Yi = float(split[3])
    #                 Zi = float(split[4])
    #                 displist = modedict[label1]
    #                 Xf = Xi + dispamnt*displist[0]
    #                 Yf = Yi + dispamnt*displist[1]
    #                 Zf = Zi + dispamnt*displist[2]
    #                 atomsymbol = re.match(r"([0-9]+)([a-z]+)", label1, re.I)
    #                 if atomsymbol:
    #                     asymb = atomsymbol.groups()
    #                     asymb1 = asymb[1]
    #                     #print asymb1
    #                 #print >> MOVIExyz, asymb1 + ' ' + str(Xf) + ' ' + str(Yf) + ' ' + str(Zf)
    #                 #MOVIExyz.write('\n')
    #                 MOVIExyz.write(asymb1 + ' ' + str(Xf) + ' ' + str(Yf) + ' ' + str(Zf))
    #                 MOVIExyz.write('\n')
    #     MOVIExyz.write(str(totalnewatoms))
    #     MOVIExyz.write('\n')
    #     MOVIExyz.write('Displacement Amount: 0')
    #     MOVIExyz.write('\n')
    #     #print >> MOVIExyz, totalatoms
    #     #print >> MOVIExyz, 'Displacement Amount: 0'
    #     #with open('supercell.txt', 'r') as pos:
    #     with open(dimer_label_xyz,'r') as pos:
    #         for line in pos:
    #             split = line.split()
    #             label1 = split[0]
    # ##            Xi = float(split[1])
    # ##            Yi = float(split[2])
    # ##            Zi = float(split[3])
    #             Xi = float(split[2])
    #             Yi = float(split[3])
    #             Zi = float(split[4])
    #             atomsymbol = re.match(r"([0-9]+)([a-z]+)", label1, re.I)
    #             if atomsymbol:
    #                 asymb = atomsymbol.groups()
    #                 asymb1 = asymb[1]
    #                 #print asymb1
    #             #print >> MOVIExyz, asymb1 + ' ' + str(Xi) + ' ' + str(Yi) + ' ' + str(Zi)
    #             #MOVIExyz.write('\n')    
    #             MOVIExyz.write(asymb1 + ' ' + str(Xi) + ' ' + str(Yi) + ' ' + str(Zi))
    #             MOVIExyz.write('\n')
    #     MOVIExyz.close()
    
    
    ##################
    
    # if MovieFull == True:
    #     moviemaker(ModeNumber)
    
    ###Mode analysis
    def moviemaker_periodic(modefreq):
        FreqNumber = 69
        ModeFlag = False
        modedict = {}
        with open('Eigs.txt', 'r') as Eigs:
            for line in Eigs:
                if ModeFlag==True:
                    if '----MODE:' in line:
                        break
                    splity = line.split()
                    label = splity[0]
                    xyzlist = [float(splity[1]), float(splity[2]), float(splity[3])]
                    modedict[label]=xyzlist
                if '----MODE:' in line:
                    #print('I am working')
                    splity = line.split()
                    if splity[1]==str(modefreq):
                        FreqNumber = int(float(splity[3]))
                        ModeFlag=True
    
        totalatoms = str(int(atomspercell)*int(Supercell[0])*int(Supercell[1])*int(Supercell[2]))
        #print(totalatoms)
        totalatoms = supercount
        stepsize = dispamp/NumberOfSteps
        reversesteps = NumberOfSteps*-1
    
        fname = 'Mode_' + str(modefreq) + '_Freq_' + str(FreqNumber) + '.xyz'
    
        MOVIExyz = open(fname, 'w')
    
        for x in range(0, NumberOfSteps,1):
            dispamnt = x*stepsize
            MOVIExyz.write(str(totalatoms))
            MOVIExyz.write('\n')
            MOVIExyz.write('Displacement Amount: ' + str(dispamnt))
            MOVIExyz.write('\n')
            #print >> MOVIExyz, totalatoms
            #print >> MOVIExyz, 'Displacement Amount: ' + str(dispamnt)
            with open('supercell.txt', 'r') as pos:
                for line in pos:
                    split = line.split()
                    label1 = split[0]
                    #print(label1)
                    Xi = float(split[1])
                    Yi = float(split[2])
                    Zi = float(split[3])
                    #print(label1)
                    displist = modedict[label1]
                    #print(label1)
                    Xf = Xi + dispamnt*displist[0]
                    Yf = Yi + dispamnt*displist[1]
                    Zf = Zi + dispamnt*displist[2]
                    atomsymbol = re.match(r"([0-9]+)([a-z]+)", label1, re.I)
                    if atomsymbol:
                        asymb = atomsymbol.groups()
                        asymb1 = asymb[1]
                        #print asymb1
                    #print >> MOVIExyz, asymb1 + ' ' + str(Xf) + ' ' + str(Yf) + ' ' + str(Zf)
                    #MOVIExyz.write('\n')
                    MOVIExyz.write(asymb1 + ' ' + str(Xf) + ' ' + str(Yf) + ' ' + str(Zf))
                    MOVIExyz.write('\n')
        for x in range(NumberOfSteps, reversesteps, -1):
            #print x
            dispamnt = x*stepsize
            #MOVIExyz.write('\n')
            MOVIExyz.write(str(totalatoms))
            MOVIExyz.write('\n')
            #print >> MOVIExyz, totalatoms
            MOVIExyz.write('Displacement Amount: ' + str(dispamnt))
            MOVIExyz.write('\n')
            #print >> MOVIExyz, 'Displacement Amount: ' + str(dispamnt)
            with open('supercell.txt', 'r') as pos:
                for line in pos:
                    split = line.split()
                    label1 = split[0]
                    Xi = float(split[1])
                    Yi = float(split[2])
                    Zi = float(split[3])
    
                    displist = modedict[label1]
                    Xf = Xi + dispamnt*displist[0]
                    Yf = Yi + dispamnt*displist[1]
                    Zf = Zi + dispamnt*displist[2]
                    atomsymbol = re.match(r"([0-9]+)([a-z]+)", label1, re.I)
                    if atomsymbol:
                        asymb = atomsymbol.groups()
                        asymb1 = asymb[1]
                        #print asymb1
                   # print >> MOVIExyz, asymb1 + ' ' + str(Xf) + ' ' + str(Yf) + ' ' + str(Zf)
                    #MOVIExyz.write('\n')
                    MOVIExyz.write(asymb1 + ' ' + str(Xf) + ' ' + str(Yf) + ' ' + str(Zf))
                    MOVIExyz.write('\n')
            #with open('DNTT_disp_new.xyz', 'w') as disp:
        for x in range(reversesteps, 0):
            dispamnt = x*stepsize
            #MOVIExyz.write('\n')
            MOVIExyz.write(str(totalatoms))
            MOVIExyz.write('\n')
            MOVIExyz.write('Displacement Amount: ' + str(dispamnt))
            MOVIExyz.write('\n')
            #print >> MOVIExyz, totalatoms
            #print >> MOVIExyz, 'Displacement Amount: ' + str(dispamnt)
            with open('supercell.txt', 'r') as pos:
                for line in pos:
                    split = line.split()
                    label1 = split[0]
                    Xi = float(split[1])
                    Yi = float(split[2])
                    Zi = float(split[3])
    
                    displist = modedict[label1]
                    Xf = Xi + dispamnt*displist[0]
                    Yf = Yi + dispamnt*displist[1]
                    Zf = Zi + dispamnt*displist[2]
                    atomsymbol = re.match(r"([0-9]+)([a-z]+)", label1, re.I)
                    if atomsymbol:
                        asymb = atomsymbol.groups()
                        asymb1 = asymb[1]
                        #print asymb1
                    #print >> MOVIExyz, asymb1 + ' ' + str(Xf) + ' ' + str(Yf) + ' ' + str(Zf)
                    #MOVIExyz.write('\n')
                    MOVIExyz.write(asymb1 + ' ' + str(Xf) + ' ' + str(Yf) + ' ' + str(Zf))
                    MOVIExyz.write('\n')
        MOVIExyz.write(str(totalatoms))
        MOVIExyz.write('\n')
        MOVIExyz.write('Displacement Amount: 0')
        MOVIExyz.write('\n')
        #print >> MOVIExyz, totalatoms
        #print >> MOVIExyz, 'Displacement Amount: 0'
        with open('supercell.txt', 'r') as pos:
            for line in pos:
                split = line.split()
                label1 = split[0]
                Xi = float(split[1])
                Yi = float(split[2])
                Zi = float(split[3])
                atomsymbol = re.match(r"([0-9]+)([a-z]+)", label1, re.I)
                if atomsymbol:
                    asymb = atomsymbol.groups()
                    asymb1 = asymb[1]
                    #print asymb1
                #print >> MOVIExyz, asymb1 + ' ' + str(Xi) + ' ' + str(Yi) + ' ' + str(Zi)
                #MOVIExyz.write('\n')    
                MOVIExyz.write(asymb1 + ' ' + str(Xi) + ' ' + str(Yi) + ' ' + str(Zi))
                MOVIExyz.write('\n')
        MOVIExyz.close()
    
    MODES = open('modesvmd.nmm', 'w')
    counta = 0
    with open('dimer.txt', 'r') as pos:
         for line in pos:
             split = line.split()
             label1 = split[0]
             #print(label1)
             Xi = float(split[1])
             Yi = float(split[2])
             Zi = float(split[3])

         #print(modedict)
         #print(displist)
                 
             displist = modedict[label1]    
             Xf = Xi + dispamp*displist[0]
             Yf = Yi + dispamp*displist[1]
             Zf = Zi + dispamp*displist[2]
        

             counta += 1

             MODES.write('ATOM ' + str(counta) + ' ' + str(displist[0]) + ' ' + str(displist[1]) + ' ' + str(displist[2]) + '\n')


    MODES.close()
    
    
    
    ##################
    
    if MovieFull == True:
        # moviemaker(ModeNumber)
        moviemaker_periodic(ModeNumber)
    
