atomicpos = open('atompos.txt', 'w')
CRYSTALOUT = 'output.txt'
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

cellparam = cellparam.split()

infile0 = open("atompos.txt", "r")
outfile = open("CIF_FOR_MERCURY.cif",'w')
infile = infile0.readlines()
outfile.write(f'''
data__ref
_symmetry_cell_setting           monoclinic
_symmetry_space_group_name_H-M   'P1'
_symmetry_Int_Tables_number      1
_cell_length_a                   {cellparam[0]}
_cell_length_b                   {cellparam[1]}
_cell_length_c                   {cellparam[2]}
_cell_angle_alpha                {cellparam[3]}
_cell_angle_beta                 {cellparam[4]}
_cell_angle_gamma                {cellparam[5]}
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
''')
for i in range(0,len(infile)):
    infile[i]=infile[i].split()
    outfile.write(infile[i][0] + infile[i][3] + " " +infile[i][3] + " " + infile[i][4] + " " + infile[i][5] + " " + infile[i][6] + "\n")
    
outfile.write("#END")
outfile.close()
infile0.close()