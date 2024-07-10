atomicpos = open('atompos.txt', 'w')
CRYSTALOUT = 'output.out'
cellparam = []
with open(CRYSTALOUT, 'r') as output:
    CellTag=False
    PositionTag = False
    cellparam = []
    for line in output:
        if CellTag==True:
            cellparam.append(line)
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


#cellparam = cellparam.split()

infile0 = open("atompos.txt", "r")
outfile = open("CIF_FOR_MERCURY.cif",'w')
infile = infile0.readlines()
for i in range(0,len(infile)):
    infile[i]=infile[i].split()
for i in range(0,len(cellparam)):
    cellparam[i] = cellparam[i].split()
    outfile.write(f'''
    data__optstep{i+1}
    _symmetry_cell_setting           triclinic
    _symmetry_space_group_name_H-M   'P1'
    _symmetry_Int_Tables_number      1
    _cell_length_a                   {cellparam[i][0]}
    _cell_length_b                   {cellparam[i][1]}
    _cell_length_c                   {cellparam[i][2]}
    _cell_angle_alpha                {cellparam[i][3]}
    _cell_angle_beta                 {cellparam[i][4]}
    _cell_angle_gamma                {cellparam[i][5]}
    loop_
    _atom_site_label
    _atom_site_type_symbol
    _atom_site_fract_x
    _atom_site_fract_y
    _atom_site_fract_z
    ''')
    for j in range((i*int(splity[0])),int(splity[0])*(i+1)):
        outfile.write(infile[j][0] + infile[j][3] + " " +infile[j][3] + " " + infile[j][4] + " " + infile[j][5] + " " + infile[j][6] + "\n")
    
outfile.write("#END")
outfile.close()
infile0.close()
