from operator import itemgetter

# This script is made with lot of assumption to the file format. Please check them before using this script.
# Cell length information should present with 2nd row, and that line must contains lo in xlo, ylo, zlo coordinates and hi in xhi, yhi, zhi coordinates.

target_atom = 3 # Check based on atom representation in lammps
low_angstrom = 2.0 # Distance criteria between first atom and target atom. more than which we will select for performing PMF                              #
high_angstrom = 3.5 # Distance criteria between first atom and target atom. more than which we will select for performing PMF                            #
bond_pair_distance = {'H_O':1.2, 'O_H':1.2, 'Si_O':2.0, 'O_Si':2.0, 'O_O':3, 'H_H':0.0, 'Si_H':0.0, 'H_Si':0.0, 'Si_Si':3.5}
lammps_atom_code = {1:'O', 2:'H', 3:'Si'}


import re

def dist_calculator_pbc(set_1, set_2, cell_length):
    zlx = float(cell_length[0])
    zly = float(cell_length[1])
    zlz = float(cell_length[2])
    zlx2 = float(cell_length[0])/2
    zly2 = float(cell_length[1])/2
    zlz2 = float(cell_length[2])/2
    x1 = float(set_1.split()[3])
    y1 = float(set_1.split()[4])
    z1 = float(set_1.split()[5])
    x2 = float(set_2.split()[3])
    y2 = float(set_2.split()[4])
    z2 = float(set_2.split()[5])
    dx = x1 - x2
    if dx > zlx2:
        dx = dx - zlx
    if dx < -zlx2:
        dx = dx + zlx
    dy = y1 - y2
    if dy > zly2:
        dy = dy - zly
    if dy < -zly2:
        dy = dy + zly
    dz = z1 - z2
    if dz > zlz2:
        dz = dz - zlz
    if dz < -zlz2:
        dz = dz + zlz
    distance = (dx*dx+dy*dy+dz*dz) ** 0.5

    return distance

def bond_checker(first_atom_info_list, second_atom_info_list, cell_length): # Both these info are to be given in dlpoly config processed format
    atom_type_1 = lammps_atom_code.get(int(first_atom_info_list.strip().split()[1]))
    atom_type_2 = lammps_atom_code.get(int(second_atom_info_list.strip().split()[1]))
    bond = str(atom_type_1) + '_' + str(atom_type_2)
    bond_length_criteria = float(bond_pair_distance.get(bond))
    dist_bet_two_atoms = dist_calculator_pbc(first_atom_info_list, second_atom_info_list, cell_length)
    if dist_bet_two_atoms <= bond_length_criteria:
        return 'yes'
    else:
        return 'No'


# x1_1 = open('input_structure_2.data').readlines()
x1_1 = open('input_structure_2.data').readlines()

cell_length = []

position_of_b4_atom_coords = []

for lines1_1 in x1_1:
    r1 = re.search(r'hi', lines1_1)
    r2 = re.search(r'lo', lines1_1)
    if r1 and r2:
        cell_length.append(float(lines1_1.strip().split()[1]) - float(lines1_1.strip().split()[0]))
#         cell_length.append(lines1_1.strip().split()[1])
    r3 = re.search(r'^Atoms', lines1_1)
    if r3:
        position_of_b4_atom_coords.append(x1_1.index(lines1_1))

target_atom_coordinates = []
oxygen_belonging_full_water_coordinates = []
all_atom_coordinates = []
for lines2_1 in x1_1[position_of_b4_atom_coords[0]+1:]:
    if lines2_1 != '\n':
        atom_type = (lines2_1.strip().split()[1])
        if int(atom_type) == target_atom:
            target_atom_coordinates.append(lines2_1)

# # Splitting up of oxygen belonging to water and oxygen belonging to glass
pos_of_water_until = int(target_atom_coordinates[0].split()[0])
all_water_molecules = []
position_of_all_o = []
for lines2_2 in x1_1[position_of_b4_atom_coords[0]+1:position_of_b4_atom_coords[0]+1+pos_of_water_until]:
    if lines2_2 != '\n':
        all_water_molecules.append(lines2_2)
        if int(lines2_2.split()[1]) == 1:
            position_of_all_o.append(all_water_molecules.index(lines2_2))

# Verification if the water molecule is full or not
for lines2_3 in position_of_all_o:
    oxygen = all_water_molecules[lines2_3]
    hydrogen1 = all_water_molecules[lines2_3+1]
    hydrogen2 = all_water_molecules[lines2_3+2]
    path1 = bond_checker(oxygen, hydrogen1, cell_length)
    path2 = bond_checker(oxygen, hydrogen2, cell_length)
    if path1 and path2 == 'yes':
        oxygen_belonging_full_water_coordinates.append(oxygen)

# # # # Same water close to multiple Si. We stack all Si based on water molecules. We chose the water vs Si pair with shortest distance among them.
out_1 = open('atom_pairs.txt', 'w')
selected_atom_pairs = []
for lines3_1 in target_atom_coordinates:
    for lines3_2 in oxygen_belonging_full_water_coordinates:
        x2_3 = dist_calculator_pbc(lines3_1, lines3_2, cell_length)
        # if x2_3 <= 3.5:
        #     print(lines3_1, lines3_2, x2_3)
        if x2_3 <= high_angstrom and x2_3 >= low_angstrom:
            sel = ('O\t%s\tSi\t%s\t%.2f' %(lines3_2.split()[0], lines3_1.split()[0], x2_3))
            selected_atom_pairs.append('O\t%s\tSi\t%s\t%.2f' %(lines3_2.split()[0], lines3_1.split()[0], x2_3))



first_layer_filter = {}
for lines4_1 in selected_atom_pairs:
    splitted_sap = lines4_1.split('\t')
    water_base = splitted_sap[0] + '_' + splitted_sap[1]
    first_layer_filter.setdefault(water_base, []).append(splitted_sap)


for lines4_2, lines4_3 in first_layer_filter.items():
    if len(lines4_3) > 1:
        sorted_1 = sorted(lines4_3, key=itemgetter(4))
        for lines4_4 in sorted_1[:1]:
            out_1.write('\t'.join(lines4_4) + '\n')
    if len(lines4_3) == 1:
        for lines4_5 in lines4_3:
            out_1.write('\t'.join(lines4_5) + '\n')


