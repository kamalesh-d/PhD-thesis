from operator import itemgetter

# This script is made with lot of assumption to the file format. Please check them before using this script.
# Cell length information should present with 2nd row, and that line must contains lo in xlo, ylo, zlo coordinates and hi in xhi, yhi, zhi coordinates.

target_atom = 3 # Check based on atom representation in lammps
low_silanol_cutoff = 3.0 # Distance criteria between first atom and target atom. more than which we will select for performing PMF                              #
high_silanol_cutoff = 3.75 # Distance criteria between first atom and target atom. more than which we will select for performing PMF                            #
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
all_oxygen_in_system = []
all_hydrogen_in_system = []
for lines2_1 in x1_1[position_of_b4_atom_coords[0]+1:]:
    if lines2_1 != '\n':
        all_atom_coordinates.append(lines2_1.strip())
        atom_type = (lines2_1.strip().split()[1])
        if int(atom_type) == target_atom:
            target_atom_coordinates.append(lines2_1.strip())
        if lammps_atom_code.get(int(atom_type)) == 'O':
            all_oxygen_in_system.append(lines2_1.strip())
        if lammps_atom_code.get(int(atom_type)) == 'H':
            all_hydrogen_in_system.append(lines2_1.strip())


# # # # Oxygen atom connected to the target Si atom:
all_o_connected_to_si = []
all_si_connected_to_o = []
for lines3_1 in target_atom_coordinates:
    for lines3_2 in all_oxygen_in_system:
        si_o_bond = (bond_checker(lines3_1, lines3_2, cell_length))
        if si_o_bond == 'yes':
            if lines3_2 not in all_o_connected_to_si:
                all_o_connected_to_si.append(lines3_2)
            if lines3_1 not in all_si_connected_to_o:
                all_si_connected_to_o.append(lines3_1)

all_h_of_silanol = []
all_o_of_silanol = []
for lines3_3 in all_o_connected_to_si:
    for lines3_4 in all_hydrogen_in_system:
        si_o_h_bond = bond_checker(lines3_3, lines3_4, cell_length)
        if si_o_h_bond == 'yes':
            all_h_of_silanol.append(lines3_4)
            all_o_of_silanol.append(lines3_3)

all_si_of_silanol = []
for lines3_5 in all_si_connected_to_o:
    for lines3_6 in all_o_of_silanol:
        if bond_checker(lines3_5, lines3_6, cell_length) == 'yes':
            if lines3_5 not in all_si_of_silanol:
                all_si_of_silanol.append(lines3_5)

selected_pair_1 = []
for lines3_7, lines3_8 in zip(all_o_of_silanol, all_si_of_silanol):
    si_o_dist = (dist_calculator_pbc(lines3_7, lines3_8, cell_length))
    if si_o_dist > 3 and si_o_dist < 3.75:
        oxy_splitted = lines3_7.strip().split()
        selected_pair_1.append(lammps_atom_code.get(int(lines3_7.strip().split()[1])) + '\t' + lines3_7.strip().split()[0] + '\t' + lammps_atom_code.get(int(lines3_8.split()[1])) + '\t' + lines3_8.strip().split()[0] + '\t' + str('%.2f' %(si_o_dist)) + '\n' )


out_1 = open('atom_pairs.txt', 'w')
for lines4_1 in selected_pair_1:
    out_1.write(lines4_1)




