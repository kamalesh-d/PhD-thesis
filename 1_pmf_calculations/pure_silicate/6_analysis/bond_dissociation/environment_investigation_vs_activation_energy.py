import os
import re
import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import seaborn as sns



bond_pair_distance = {'H_O':1.2, 'O_H':1.2, 'Si_O':2.0, 'O_Si':2.0, 'O_O':3, 'H_H':0.0, 'Si_H':0.0, 'H_Si':0.0, 'Si_Si':3.5, 'Al_O':2.4, 'O_Al':2.4, 'Ca_O':3.24, 'O_Ca':3.24, 'Al_Al':0, 'Ca_Al':0, 'H_Al':0, 'Si_Al':0, 'Al_Si':0, 'Ca_Si':0}
major_cations_in_glass = ['Si', 'Al']

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

def arrange_based_column(list_of_list, column_to_sort):
    sorted_list = []
    l = len(list_of_list)
    for linesd3_1 in range(0,l):
        for linesd3_2 in range(0, l-linesd3_1-1):
            if int(list_of_list[linesd3_2][column_to_sort-1]) > int(list_of_list[linesd3_2+1][column_to_sort-1]):
                tempo = list_of_list[linesd3_2]
                list_of_list[linesd3_2] = list_of_list[linesd3_2+1]
                list_of_list[linesd3_2+1] = tempo
    return list_of_list


def history_to_config_converter(input_file, timestep_to_convert):
    x1_1 = open('%s' %(input_file)).readlines()
    compiled_lines = []
    to_filter_lines = []
    for lines1_2 in x1_1[2:]:
        to_filter_lines.append(lines1_2)
    position_of_his_1000 = []
    total_num_of_atoms = []
    for lines1_3 in to_filter_lines:
        r1 = re.search(r'^timestep', lines1_3)
        if r1:
            a1_1 = lines1_3.split()[1]
            if int(a1_1) == timestep_to_convert:
                a1_2 = to_filter_lines.index(lines1_3)
                position_of_his_1000.append(str(a1_2))
                total_num_of_atoms.append((lines1_3.split()[2]))
    timestep_line_num = int(''.join(position_of_his_1000))
    filter_line_list = []
    for lines1_1 in x1_1[:2]:
        f1_1 = (lines1_1.split()[0:2])
        f1_1.append(''.join(total_num_of_atoms))
    x1_2 = range(int(len(to_filter_lines[timestep_line_num+4:timestep_line_num+int(''.join(total_num_of_atoms))*4])/4+1))
    for lines1_2 in x1_2:
        sub_list = []
        i = 4*lines1_2
        while i < (lines1_2+1)*4:
            sub_list.append(to_filter_lines[timestep_line_num+4:][i])
            i = i + 1
        filter_line_list.append(sub_list)
    for lines1_6 in filter_line_list:
        new_list = ((['%s   %s' %(lines1_6[0].split()[0], lines1_6[0].split()[1]), lines1_6[1].strip(), lines1_6[2].strip(), lines1_6[3].strip()]))
        compiled_lines.append(new_list)
    assembled_list = []
    for linesd1_1 in compiled_lines:
        atom_type = linesd1_1[0].split()[0]
        atom_num = linesd1_1[0].split()[1]
        atom_coord = linesd1_1[1].split()
        pre_assembled = [atom_type, atom_num, atom_coord[0], atom_coord[1], atom_coord[2]]
        assembled_list.append(pre_assembled)
    return assembled_list

def bond_atoms_finder(first_atom_info_list, second_atom_info_list, cell_length):
    bc_first_atom_type = str(first_atom_info_list[0])
    bc_second_atom_type = str(second_atom_info_list[0])
    bond = bc_first_atom_type + '_' + bc_second_atom_type
    bond_length_criteria = float(bond_pair_distance.get(bond))
    dist_bet_two_atoms = dist_calculator_pbc(first_atom_info_list, second_atom_info_list, cell_length)
    if dist_bet_two_atoms <= bond_length_criteria:
        return 'yes'
    else:
        return 'No'



def broken_bond_atoms_finder(initial_bonded_config, final_broken_config, cell_length):
    initial_target_al = initial_bonded_config[-4]
    oxy_connected_to_initial_al = []
    for lines1_1 in initial_bonded_config:
        for lines1_2 in initial_bonded_config[-4:-3]:
            x2_1 = bond_atoms_finder(lines1_1, lines1_2, cell_length)
            if x2_1 == 'yes':
                if lines1_1[0] == 'O':
                    oxy_connected_to_initial_al.append(lines1_1)
    # Now we will find all the atoms connected to the oxygen connected with Si from initial structure
    initial_three_atom_pairs = []
    for lines1_3 in initial_bonded_config:
        for lines1_4 in oxy_connected_to_initial_al:
            x3_1 = bond_atoms_finder(lines1_3, lines1_4, cell_length)
            if x3_1 == 'yes':
                if lines1_3[0] == 'Si' or lines1_3[0] == 'Al':
                    if lines1_3 != initial_target_al:
                        initial_three_atom_pairs.append(initial_target_al[0] + '_' + initial_target_al[1] + '_' + lines1_4[0] + '_' + lines1_4[1] + '_' + lines1_3[0] + '_' + lines1_3[1])

    final_target_al = final_broken_config[-4]
    oxy_connected_to_final_al = []
    for lines2_1 in final_broken_config:
        for lines2_2 in final_broken_config[-4:-3]:
            x4_1 = bond_atoms_finder(lines2_1, lines2_2, cell_length)
            if x4_1 == 'yes':
                if lines2_1[0] == 'O':
                    oxy_connected_to_final_al.append(lines2_1)

    # Now we will find all the atoms connected to the oxygen connected with Si from final structure
    final_three_atom_pairs = []
    for lines2_3 in final_broken_config:
        for lines2_4 in oxy_connected_to_final_al:
            x5_1 = bond_atoms_finder(lines2_3, lines2_4, cell_length)
            if x5_1 == 'yes':
                if lines2_3[0] == 'Si' or lines2_3[0] == 'Al':
                    if lines2_3 != final_target_al:
                        final_three_atom_pairs.append(final_target_al[0] + '_' + final_target_al[1] + '_' + lines2_4[0] + '_' + lines2_4[1] + '_' + lines2_3[0] + '_' + lines2_3[1])
    # print(len(initial_three_atom_pairs), len(final_three_atom_pairs))
    broken_bond_atoms = []
    for lines3_1 in initial_three_atom_pairs:
        if lines3_1 not in final_three_atom_pairs:
            broken_bond_atoms.append(lines3_1)
    return broken_bond_atoms

def target_atoms_bonded(final_config, cell_length):
    final_target_al = final_config[-4]
    oxy_connected_to_final_al = []
    for lines2_1 in final_config:
        for lines2_2 in final_config[-4:-3]:
            x4_1 = bond_atoms_finder(lines2_1, lines2_2, cell_length)
            if x4_1 == 'yes':
                if lines2_1[0] == 'O':
                    oxy_connected_to_final_al.append(lines2_1)

    # Now we will find all the atoms connected to the oxygen connected with Si from final structure
    final_three_atom_pairs = []
    for lines2_3 in final_config:
        for lines2_4 in oxy_connected_to_final_al:
            x5_1 = bond_atoms_finder(lines2_3, lines2_4, cell_length)
            if x5_1 == 'yes':
                if lines2_3[0] == 'Si' or lines2_3[0] == 'Al':
                    if lines2_3 != final_target_al:
                        final_three_atom_pairs.append(final_target_al[0] + '_' + final_target_al[1] + '_' + lines2_4[0] + '_' + lines2_4[1] + '_' + lines2_3[0] + '_' + lines2_3[1])
    return final_three_atom_pairs

# def broken_atom_neighbours_finder(history_file, broken_atom_type, broken_atom_number, cell_length):
#     collect_neighbours_cations_5_ang = []
#     for linesd1_2 in history_file:
#         for linesd1_1 in history_file:
#             if linesd1_1[0] == broken_atom_type and linesd1_1[1] == broken_atom_number:
#                 if linesd1_2 != linesd1_1:
#                     if linesd1_2[0] in major_cations_in_glass:
#                         if dist_calculator_pbc(linesd1_1, linesd1_2, cell_length) <= 5.0:
#                             collect_neighbours_cations_5_ang.append(linesd1_2)
#     return len(collect_neighbours_cations_5_ang)

def broken_atom_depth_finder(history_file, broken_atom_type, broken_atom_number, cell_length):
    # collect_broken_atom_depth = []
    for linesd1_2 in history_file:
        if linesd1_2[0] == broken_atom_type and linesd1_2[1] == broken_atom_number:
            z_coordinate_broken_atom = float(linesd1_2[4])
            z_coordinate_silicate_starting = float(cell_length[0])/2
            depth_of_broke_atom = z_coordinate_silicate_starting - abs(z_coordinate_broken_atom)
            print(depth_of_broke_atom)
            return depth_of_broke_atom

def normalize_the_list(list):
    largest_num = max(list)
    normalized = []
    for linesd15_1 in list:
        normalized.append(linesd15_1/largest_num)
    return normalized

def pbc_diff_finder(a_atom_coordinates, b_atom_coordinates, cell_length_list):
    zlx = float(cell_length_list[0])
    zly = float(cell_length_list[1])
    zlz = float(cell_length_list[2])
    zlx2 = float(cell_length_list[0])/2
    zly2 = float(cell_length_list[1])/2
    zlz2 = float(cell_length_list[2])/2
    x1 = float(a_atom_coordinates[0])
    y1 = float(a_atom_coordinates[1])
    z1 = float(a_atom_coordinates[2])
    x2 = float(b_atom_coordinates[0])
    y2 = float(b_atom_coordinates[1])
    z2 = float(b_atom_coordinates[2])
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
    return [dx, dy, dz]


def bond_angle_finder(history_file_config, list_of_broken_bond_atoms, cell_length):
    all_angles = []
    for linesd1_1 in list_of_broken_bond_atoms:
        three_points = []
        for linesd18_3 in history_file_config:
            if linesd18_3[1] == linesd1_1.split('_')[1]:
                three_points.append(linesd18_3)

        for linesd18_4 in history_file_config:
            if linesd18_4[1] == linesd1_1.split('_')[3]:
                three_points.append(linesd18_4)

        for linesd18_5 in history_file_config:
            if linesd18_5[1] == linesd1_1.split('_')[5]:
                three_points.append(linesd18_5)

        a = np.array([float(three_points[0][2]), float(three_points[0][3]), float(three_points[0][4])])
        b = np.array([float(three_points[1][2]), float(three_points[1][3]), float(three_points[1][4])])
        c = np.array([float(three_points[2][2]), float(three_points[2][3]), float(three_points[2][4])])
        ba = pbc_diff_finder(a, b, cell_length)
        bc = pbc_diff_finder(c, b, cell_length)
        cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
        angle = np.arccos(cosine_angle)
        all_angles.append(np.degrees(angle))
    return (all_angles)


# # # # Difference in stress or shear between broken vs bonded atoms
coordination_table = open('compiled_coordination_activation_table_updated_1.txt').readlines()
new_coordination_table = []
current_directory = os.getcwd()
all_x = []
all_y = []
out_1 = open('compiled_coordination_activation_table_updated_2.txt', 'w')
out_1.write(coordination_table[0].strip() + '\t' + 'Broken_bond\'s_depth' + '\n')
all_diff_angle = []
for lines12_1 in coordination_table[1:]:
    target_directory = ('/'.join(current_directory.split('/')[:-2]))
    glass_directory = ('%s/%s' %(target_directory, lines12_1.split('\t')[0]))
    pmf_directory = ('%s/5_pmf/%s' %(glass_directory, lines12_1.split('\t')[1]))
    history_file = '%s/HISTORY_%s_1.50' %(pmf_directory, lines12_1.split('\t')[1])
    initial_config = history_to_config_converter(history_file, 0)
    final_config = history_to_config_converter(history_file, 10000)
