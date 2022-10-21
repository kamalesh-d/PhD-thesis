import os
from operator import itemgetter

os.system('cp CONFIG SOURCE_CONFIG')

###########################################################################################################################################################
# This python script creates CONFIG only if the water molecules are full and real. If target oxygen is already dissociated, such borken water will not be created
num_line_con=4 #Each atom is defined by number of lines in config file has to be given here                                                               #
input_file_name = 'SOURCE_CONFIG' # Type the file name to be processed                                                                                    #
first_atom = 'O'  # Note: This must be specified with water, because the script is written to grab the neighbouring h2 atoms along with oxygen atom       #
target_atom = 'Si' # User can change depending on the target towards which we will move the water molecule                                                #
atom_next_to_water = 'Si' # User has to report the first atom placed very after the Water atom                                                            #
low_angstrom = 2.0 # Distance criteria between first atom and target atom. more than which we will select for performing PMF                              #
high_angstrom = 3.5 # Distance criteria between first atom and target atom. more than which we will select for performing PMF                            #
bond_pair_distance = {'H_O':1.2, 'O_H':1.2, 'Si_O':2.0, 'O_Si':2.0, 'O_O':3, 'H_H':1.5, 'Si_H':2.0, 'H_Si':2.0, 'Si_Si':3.5}
# You will get the output as CONFIG_*    where * represent numbers from 1 to n depending on number of atoms falling under the criteria we described       #
# Note: This python script is created assuming first 5 lines in input file has not started the first atom element co-ordinate description                 #
# Note: This python script is updated with config to impose PMF only for O and Si                                                                         #
###########################################################################################################################################################

def config_input_processor(input_file):
    x1_1 = open('%s' %(input_file)).readlines()
    compiled = []
    list_1 = []
    final_list = []
    for lines1_1 in x1_1[5:]:
        list_1.append(lines1_1.strip())
    x1_2 = range(int(len(list_1)/num_line_con))
    for lines1_2 in x1_2:
        sub_list = []
        i = num_line_con*lines1_2
        while i < (lines1_2+1)*num_line_con:
            sub_list.append(list_1[i])
            i = i + 1
        final_list.append(sub_list)
    for lines1_3 in final_list:
        atom_type = (lines1_3[0].split()[0])
        atom_number = (lines1_3[0].split()[1])
        x_coordinate = lines1_3[1].split()[0]
        y_coordinate = lines1_3[1].split()[1]
        z_coordinate = lines1_3[1].split()[2]
        compile = atom_type + '\t' + atom_number + '\t' + x_coordinate + '\t' + y_coordinate + '\t' + z_coordinate
        compiled.append(compile.strip().split('\t'))
    return compiled

def normalize_num(list_of_list):
    l = len(list_of_list)
    aligned_1 = []
    for lines6_1, lines6_2 in zip(list_of_list, range(0,l)):
        x6_1 = list_of_list[lines6_2][1]
        mod_atom_num = lines6_2 + 1
        lines6_1[1] = str(mod_atom_num)
        aligned_1.append(lines6_1)
    return aligned_1


def dist_calculator_pbc(set_1, set_2, cell_length):
    zlx = float(cell_length[0])
    zly = float(cell_length[1])
    zlz = float(cell_length[2])
    zlx2 = float(cell_length[0])/2
    zly2 = float(cell_length[1])/2
    zlz2 = float(cell_length[2])/2
    x1 = float(set_1[2])
    y1 = float(set_1[3])
    z1 = float(set_1[4])
    x2 = float(set_2[2])
    y2 = float(set_2[3])
    z2 = float(set_2[4])
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



def dist_calculator(set_1, set_2):
    x1 = float(set_1[2])
    y1 = float(set_1[3])
    z1 = float(set_1[4])
    x2 = float(set_2[2])
    y2 = float(set_2[3])
    z2 = float(set_2[4])
    a1 = (x2 - x1) ** 2
    b1 = (y2 - y1) ** 2
    c1 = (z2 - z1) ** 2
    c2 = a1 + b1 + c1
    distance = float(c2 ** 0.5)
    return distance

def bond_checker(first_atom_info_list, second_atom_info_list, cell_length): # Both these info are to be given in dlpoly config processed format
    bc_first_atom_type = str(first_atom_info_list[0])
    bc_second_atom_type = str(second_atom_info_list[0])
    bond = bc_first_atom_type + '_' + bc_second_atom_type
    bond_length_criteria = float(bond_pair_distance.get(bond))
    dist_bet_two_atoms = dist_calculator(first_atom_info_list, second_atom_info_list)
    if dist_bet_two_atoms <= bond_length_criteria:
        return 'yes'
    else:
        return 'No'

# Distance between H2O and Si atom less than 3 angstrom are computed

x1_1 = input_file_name
x0_1 = open(x1_1).readlines()[2:5]
x_cell_length = x0_1[0].split()[0]
y_cell_length = x0_1[1].split()[1]
z_cell_length = x0_1[2].split()[2]
cell_length = [x_cell_length, y_cell_length, z_cell_length]

x1_2 = config_input_processor(x1_1)
atom_type_num_info = []
water_atom_coords = []
target_atom_coords = []
water_to_grab = []
silicon_to_grab = []
all_atoms = []
for lines1_1 in x1_2:
    atom_type = lines1_1[0]
    all_atoms.append(atom_type)
    if atom_type == target_atom:
        target_atom_coords.append(lines1_1)


water_pos = all_atoms.index('%s' %(atom_next_to_water))
for lines1_2 in x1_2[0:water_pos]:
    atom_type_1 = lines1_2[0]
    if atom_type_1 == first_atom:
        water_atom_coords.append(lines1_2)
distance_information = []
pmf_to_perform = []
for lines2_1 in target_atom_coords:
    for lines2_2 in water_atom_coords:
        x2_3 = dist_calculator_pbc(lines2_1,lines2_2, cell_length)
        if x2_3 <= high_angstrom and x2_3 >= low_angstrom:
            silicon_to_grab.append(lines2_1)
            water_to_grab.append(lines2_2)
            atom_type_num_info.append(lines2_2[0] + '_' + lines2_2[1] + '_' + lines2_1[0] + '_' + lines2_1[1])
            distance_information.append(str(x2_3))

out_2 = open('pre_dist_info.txt', 'w')
# The complete water molecule and silicon atoms are pushed to the end of the list
for lines4_1, lines4_2, i_2, lines4_3 in zip(water_to_grab, silicon_to_grab, atom_type_num_info, distance_information):
    x4_1 = input_file_name
    x4_2 = config_input_processor(x4_1)
    for lines4_3_2 in x4_2:
        if lines4_3_2 == lines4_2:
            pos_silicon = x4_2.index(lines4_3_2)
            x4_2.append(x4_2.pop(pos_silicon))
    full_water_verified = []
    for lines4_3_1 in x4_2:
        if lines4_3_1 == lines4_1:
            pos_water = x4_2.index(lines4_3_1)
            atom_type_1 = (x4_2[pos_water])
            atom_type_2 = (x4_2[pos_water+1])
            atom_type_3 = (x4_2[pos_water+2])
            if bond_checker(atom_type_1, atom_type_2, cell_length) is 'yes':
                if bond_checker(atom_type_1, atom_type_3, cell_length) is 'yes':
                    x4_2.append(x4_2.pop(pos_water))
                    x4_2.append(x4_2.pop(pos_water))
                    x4_2.append(x4_2.pop(pos_water))
                    full_water_verified.append('yes')
    if len(full_water_verified) != 0:
        if full_water_verified[0] is 'yes':
            out_1 = open('CONFIG_%s' % (i_2), 'w')
            i_3 = 0
            x6_1 = open('%s' %(x4_1)).readlines()
            for lines_1 in x6_1[:1]:
                out_1.write(lines_1)
            for lines_3 in x6_1[1:2]:
                x_3 = (lines_3.split())
                x_3[0] = '         0'
                out_1.write('         '.join(x_3) + '\n')
            for lines_4 in x6_1[2:5]:
                out_1.write((lines_4))

            x5_1 = normalize_num(x4_2)
            for lines5_1 in x5_1:
                out = '    '.join(lines5_1[0:2]) + '\n' + '    '.join(lines5_1[2:]) + '\n'
                out_1.write(out)
            out_2.write('\t'.join(i_2.split('_')) + '\t' + lines4_3 + '\n')

