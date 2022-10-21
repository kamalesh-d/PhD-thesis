import glob
import re
import os
import numpy as np
import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy.stats import linregress
import statistics



# # # # Distance and history file converter is completely modified in this

# Feature: Any particular timestep of the history can be chosen to fetch the exact number of atoms within that particular timestep

distance_bet_si_oh = 1.6 # This is the distance between Si and O of H2O, whose history file you wish to convert
name_water_against_target_atom = 'Si'
timestep_to_convert_history = 10000
diff_cut_off = 0.05
pmf_atoms_num = 4 # This is to define number of atoms are pmf in this script. Only the last number of atoms (That you define) will be considered
neighbour_grab_cutoff = 5.0
bond_pair_distance = {'H_O':1.2, 'O_H':1.2, 'Si_O':2.0, 'O_Si':2.0, 'O_O':3, 'H_H':0.0, 'Si_H':0.0, 'H_Si':0.0, 'Si_Si':3.5, 'Al_O':2.4, 'O_Al':2.4, 'Ca_O':3.24, 'O_Ca':3.24, 'Al_Al':0, 'Ca_Al':0, 'H_Al':0, 'Si_Al':0, 'Al_Si':0, 'Ca_Si':0}
num_line_con = 4

def remove_tail_dot_zeros(a):
    x10_1 = str(a)
    if (x10_1[-1]) == '0':
        return ('%.1f' %(float(x10_1)))
    if (x10_1[-1]) != '0':
        return ('%.2f' %(float(x10_1)))

def dist_calculator_pbc(set_1, set_2, cell_length):
    zlx = float(cell_length[0])
    zly = float(cell_length[1])
    zlz = float(cell_length[2])
    zlx2 = float(cell_length[0])/2
    zly2 = float(cell_length[1])/2
    zlz2 = float(cell_length[2])/2
    x1 = float(set_1[5])
    y1 = float(set_1[6])
    z1 = float(set_1[7])
    x2 = float(set_2[5])
    y2 = float(set_2[6])
    z2 = float(set_2[7])
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
        new_list = ((['%s   %s   %s   %s   %s' %(lines1_6[0].split()[0], lines1_6[0].split()[1], lines1_6[0].split()[2], lines1_6[0].split()[3], lines1_6[0].split()[4]), lines1_6[1].strip(), lines1_6[2].strip(), lines1_6[3].strip()]))
        compiled_lines.append(new_list)
    assembled_list = []
    for linesd1_1 in compiled_lines:
        atom_type = linesd1_1[0].split()[0]
        atom_num = linesd1_1[0].split()[1]
        atom_add_info = linesd1_1[0].split()
        atom_coord = linesd1_1[1].split()
        third_line = linesd1_1[2].split()
        force = linesd1_1[3].split()
        pre_assembled = [atom_type, atom_num, atom_add_info[2], atom_add_info[3], atom_add_info[4], atom_coord[0], atom_coord[1], atom_coord[2], third_line[0], third_line[1], third_line[2], force[0], force[1], force[2]]
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
        return 'no'



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

def target_atom_environment(initial_bonded_config, final_broken_config, cell_length):
    initial_target_atom = initial_bonded_config[-4]
    oxy_connected_to_bonded_target_atom = []
    for lines1_1 in initial_bonded_config:
        for lines1_2 in initial_bonded_config[-4:-3]:
            x2_1 = bond_atoms_finder(lines1_1, lines1_2, cell_length)
            if x2_1 == 'yes':
                if lines1_1[0] == 'O':
                    oxy_connected_to_bonded_target_atom.append(lines1_1)

    # Now we will find all the atoms connected to the oxygen connected with Si from initial structure
    initial_bonded_h = []
    initial_bonded_si = []
    initial_bonded_al = []

    for lines1_3 in initial_bonded_config:
        for lines1_4 in oxy_connected_to_bonded_target_atom:
            x3_1 = bond_atoms_finder(lines1_3, lines1_4, cell_length)
            if x3_1 == 'yes':
                if lines1_3 != initial_target_atom:
                    if lines1_3[0] == 'Si':
                        initial_bonded_si.append(lines1_3)
                    if lines1_3[0] == 'H':
                        initial_bonded_h.append(lines1_3)
                    if lines1_3[0] == 'Al':
                        initial_bonded_al.append(lines1_3)

    final_target_atom = final_broken_config[-4]
    oxy_connected_to_final_target_atom = []
    for lines2_1 in final_broken_config:
        for lines2_2 in final_broken_config[-4:-3]:
            x4_1 = bond_atoms_finder(lines2_1, lines2_2, cell_length)
            if x4_1 == 'yes':
                if lines2_1[0] == 'O':
                    oxy_connected_to_final_target_atom.append(lines2_1)

    # Now we will find all the atoms connected to the oxygen connected with Si from final structure
    final_bonded_h = []
    final_bonded_si = []
    final_bonded_al = []
    for lines2_3 in final_broken_config:
        for lines2_4 in oxy_connected_to_final_target_atom:
            x5_1 = bond_atoms_finder(lines2_3, lines2_4, cell_length)
            if x5_1 == 'yes':
                if lines2_3 != final_target_atom:
                    if lines2_3[0] == 'Si':
                        final_bonded_si.append(lines2_3)
                    if lines2_3[0] == 'H':
                        final_bonded_h.append(lines2_3)
                    if lines2_3[0] == 'Al':
                        final_bonded_al.append(lines2_3)


    initial_coordination_pattern = 'O%s ' %(len(oxy_connected_to_bonded_target_atom)) + 'Al%s ' %(len(initial_bonded_al)) + 'Si%s ' %(len(initial_bonded_si)) + 'H%s ' %(len(initial_bonded_h))
    final_coordination_pattern = 'O%s ' %(len(oxy_connected_to_final_target_atom)) + 'Al%s ' %(len(final_bonded_al)) + 'Si%s ' %(len(final_bonded_si)) + 'H%s ' %(len(final_bonded_h))

    return (initial_coordination_pattern + '\t' + final_coordination_pattern)


    # print(len(initial_three_atom_pairs), len(final_three_atom_pairs))
    # broken_bond_atoms = []
    # for lines3_1 in initial_three_atom_pairs:
    #     if lines3_1 not in final_three_atom_pairs:
    #         broken_bond_atoms.append(lines3_1)
    # # print(broken_bond_atoms)
    # # # return broken_bond_atoms

def integrate_for_pmf(diff_distance, force):
    a = np.array(force)
    m = np.r_[False, a > 0, False]
    idx = np.flatnonzero(m[:-1] != m[1:])
    I = (idx[1::2] - idx[::2]).argmax()
    elems = a[idx[2 * I]:idx[2 * I + 1]]
    d_x1_1 = sum(elems)
    d_x1_2 = d_x1_1*diff_distance
    return d_x1_2

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

def find_pattern_and_replace(input_file, pattern_in_target_line, replace_line, file_directory_location, output_file_name):
    xd2_1 = open('%s/%s' %(file_directory_location, input_file)).readlines()
    outd2_1 = open('%s/%s' %(file_directory_location, output_file_name), 'w')
    replace_line_pos_1 = []
    for linesd2_1 in xd2_1:
        rd2 = re.search(r'%s' %(pattern_in_target_line), linesd2_1)
        if rd2:
            replace_line_pos_1.append(str(xd2_1.index(linesd2_1)))
    if replace_line_pos_1 != []:
        xd2_1[int(''.join(replace_line_pos_1))] = replace_line + '      %s' %(pattern_in_target_line) + '\n'
        for linesd2_2 in xd2_1:
            outd2_1.write(linesd2_2)



def bond_angle_finder(config_file_name, list_of_broken_bond_atoms):
    config_processed = config_input_processor(config_file_name)
    three_points = []
    baf1_1 = open('%s' %(config_file_name)).readlines()
    x_length = (baf1_1[2].split()[0])
    y_length = (baf1_1[3].split()[1])
    z_length = (baf1_1[4].split()[2])
    cell_length = [x_length, y_length, z_length]
    for lines1_1 in config_processed:
        if (lines1_1[1]) == (list_of_broken_bond_atoms[3:6]):
            three_points.append(lines1_1)
        if (lines1_1[1]) == (list_of_broken_bond_atoms[9:12]):
            three_points.append(lines1_1)
        if (lines1_1[1]) == (list_of_broken_bond_atoms[16:19]):
            three_points.append(lines1_1)


    # https://stackoverflow.com/questions/35176451/python-code-to-calculate-angle-between-three-point-using-their-3d-coordinates
    # https://www.youtube.com/watch?time_continue=200&v=AN4a53rlkhM&feature=emb_logo
    a = np.array([float(three_points[0][2]), float(three_points[0][3]), float(three_points[0][4])])
    b = np.array([float(three_points[1][2]), float(three_points[1][3]), float(three_points[1][4])])
    c = np.array([float(three_points[2][2]), float(three_points[2][3]), float(three_points[2][4])])
    # ba = a - b
    # bc = c - b
    ba = pbc_diff_finder(a, b, cell_length)
    bc = pbc_diff_finder(c, b, cell_length)
    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)
    return (np.degrees(angle))

def list_to_sort(list_of_list, column_to_sort, type_of_sort):
    if type_of_sort == 'ascending':
        list_of_list.sort(key=lambda x: x[column_to_sort-1])
    if type_of_sort == 'descending':
        list_of_list.sort(key=lambda x: x[column_to_sort-1], reverse=True)
    return list_of_list

def config_realigned(config_file):
    sorted_config = []
    only_si = []
    only_o = []
    only_h = []
    out_1 = open('%s/analysis/stress_analysis/number_mapping' %(glass_directory), 'w')
    for lines2_1 in config_file:
        if lines2_1[0] == 'Si':
            sorted_config.append(lines2_1)
            only_si.append(lines2_1)
    for lines2_2 in config_file:
        if lines2_2[0] == 'O':
            sorted_config.append(lines2_2)
            only_o.append(lines2_2)
    for lines2_3 in config_file:
        if lines2_3[0] == 'H':
            sorted_config.append(lines2_3)
            only_h.append(lines2_3)
    new_num_range = range(1, len(sorted_config) + 1)
    all_old_number = []
    all_new_number = []
    for lines3_1, lines3_2 in zip(sorted_config, new_num_range):
        old_number = (lines3_1[1])
        all_old_number.append(str(old_number))
        lines3_1[1] = str(lines3_2)
        all_new_number.append(lines3_2)
    out_1.write('Old_number' + '\t' + 'New_number' + '\n')
    for lines4_1, lines4_2 in zip(all_old_number, all_new_number):
        out_1.write(str(lines4_1) + '\t' + str(lines4_2) + '\n')
    find_pattern_and_replace('FIELD', '#number_of_silicon', 'NUMMOLS %s' %(len(only_si)), '%s/analysis/stress_analysis' %(glass_directory), 'FIELD')
    find_pattern_and_replace('FIELD', '#number_of_oxygen', 'NUMMOLS %s' %(len(only_o)), '%s/analysis/stress_analysis' %(glass_directory), 'FIELD')
    find_pattern_and_replace('FIELD', '#number_of_hydrogen', 'NUMMOLS %s' %(len(only_h)), '%s/analysis/stress_analysis' %(glass_directory), 'FIELD')
    return sorted_config

"""

directory = os.getcwd()
#for lines0_1 in range(19,20):
for lines0_1 in [0, 2, 3, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]:
    target_directory = ('/'.join(directory.split('/')[:-2]))
    glass_directory = '%s/glass_num_%s/5_pmf' %(target_directory, lines0_1)
    os.system('cd %s/analysis && rm -r stress_analysis' %(glass_directory))
    os.system('cd %s/analysis && mkdir stress_analysis' %(glass_directory))
    os.system('cp FIELD TABLE LOCAL_STRESS_KD_0620_V1.f %s/analysis/stress_analysis' %(glass_directory))
    os.system('cd %s/glass_num_%s/4_equilibrate_sol_glass_2 && cp HISTORY_V4 CONTROL ../5_pmf/analysis/stress_analysis' %(target_directory, lines0_1))
    out_2 = open('%s/analysis/stress_analysis/HISTORY' %(glass_directory), 'w')
    timesteps = [0, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000]
    #
    realigned_history = []
    history_file_read = open('%s/analysis/stress_analysis/HISTORY_V4' %(glass_directory)).readlines()
    for lines2_1 in history_file_read[0:2]:
        realigned_history.append(lines2_1)

    all_timestep_contents = []
    all_sorted_config = []

    for lines1_1 in timesteps:
        history_file = history_to_config_converter('%s/analysis/stress_analysis/HISTORY_V4' %(glass_directory), lines1_1)
        sorted_config = config_realigned(history_file)
        all_sorted_config.append(sorted_config)
        timestep_contents = []
        for lines2_2 in history_file_read:
            r2 = re.search(r'^timestep +%s ' %(lines1_1), lines2_2)
            if r2:
                timestep_position_found = history_file_read.index(lines2_2)
                timestep_contents.append(history_file_read[timestep_position_found])
                timestep_contents.append(history_file_read[timestep_position_found + 1])
                timestep_contents.append(history_file_read[timestep_position_found + 2])
                timestep_contents.append(history_file_read[timestep_position_found + 3])
        all_timestep_contents.append(timestep_contents)

    for lines3_1, lines3_2 in zip(all_timestep_contents, all_sorted_config):
        realigned_history.append(lines3_1[0])
        realigned_history.append(lines3_1[1])
        realigned_history.append(lines3_1[2])
        realigned_history.append(lines3_1[3])
        for lines3_3 in lines3_2:
            realigned_history.append('%s                %s   %s   %s   %s\n' %(lines3_3[0], lines3_3[1], lines3_3[2], lines3_3[3], lines3_3[4]))
            realigned_history.append('     %s        %s         %s\n' %(lines3_3[5], lines3_3[6], lines3_3[7]))
            realigned_history.append('     %s        %s         %s\n' %(lines3_3[8], lines3_3[9], lines3_3[10]))
            realigned_history.append('     %s        %s         %s\n' %(lines3_3[11], lines3_3[12], lines3_3[13]))
    for lines4_1 in realigned_history:
        out_2.write(lines4_1)
    os.system('cd %s/analysis/stress_analysis && gfortran LOCAL_STRESS_KD_0620_V1.f && ./a.out -filepos HISTORY -filefield FIELD -filecontrol CONTROL' %(glass_directory))

    number_mapping = open('%s/analysis/stress_analysis/number_mapping' %(glass_directory)).readlines()
    pressure = open('%s/analysis/stress_analysis/Pressure.dat' %(glass_directory)).readlines()
    out_3 = open('%s/analysis/stress_analysis/number_mapped_pressure' %(glass_directory), 'w')
    out_3.write('glass_number   old_number   new_number   stress   shear\n')
    for lines5_1, lines5_2 in zip(number_mapping[1:], pressure):
        num_map_splitted = lines5_1.split()
        press_splitted = lines5_2.split()
        compiled_num_press = ['glass_num_%s' %(lines0_1), num_map_splitted[0], num_map_splitted[1], press_splitted[1], press_splitted[2]]
        out_3.write('   '.join(compiled_num_press) + '\n')


current_directory = os.getcwd()
number_of_water_in_glass = []
coordination_table_contents = []
out_1 = open('compiled_coordination_activation_table', 'w')
for lines0_1 in [0, 2, 3, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]:
    target_directory = ('/'.join(current_directory.split('/')[:-2]))
    coordination_table = open('%s/glass_num_%s/5_pmf/analysis/coordination_activation_table' %(target_directory, lines0_1)).readlines()
    for lines0_2 in coordination_table:
        if lines0_2.strip().split('\t') not in coordination_table_contents:
            coordination_table_contents.append(lines0_2.strip().split('\t'))

for lines1_2 in coordination_table_contents[0:1]:
    out_1.write('\t'.join(lines1_2) + '\n')

for lines1_1 in (list_to_sort(coordination_table_contents[1:], 7, 'descending')):
    out_1.write('\t'.join(lines1_1) + '\n')


"""


# # # Appending stress and shear information in coordination_table
coordination_table = open('compiled_coordination_activation_table').readlines()
new_coordination_table = []
current_directory = os.getcwd()
for lines0_1 in [0, 2, 3, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]:

    target_directory = ('/'.join(current_directory.split('/')[:-2]))
    glass_directory = open('%s/glass_num_%s/5_pmf/analysis/stress_analysis/number_mapped_pressure' %(target_directory, lines0_1)).readlines()


    for lines1_1 in glass_directory:
        for lines1_2 in coordination_table:
            glass_number_table = lines1_2.split('\t')[0]
            target_atom_number_table = lines1_2.split('\t')[1].split('_')[-1]
            glass_number_pressure = lines1_1.split()[0]
            old_atom_number_pressure = (lines1_1.split()[1])
            if glass_number_table == glass_number_pressure and target_atom_number_table == old_atom_number_pressure:
                temp_list = lines1_2.strip().split('\t')
                temp_list.append(lines1_1.split()[3])
                temp_list.append(lines1_1.split()[4])
                new_coordination_table.append(temp_list)
out_1 = open('compiled_coordination_activation_table_updated', 'w')
out_1.write('Glass_number' + '\t' + 'Filename' + '\t' + 'Before_breaking' + '\t' + 'After_breaking' + '\t' + 'x_coordinate' + '\t' + 'y_coordinate' + '\t' + 'z_coordinate' + '\t'  + 'Activation_barrier' + '\t' + '#_broken_bonds' + '\t' + 'List_of_broken_bond' + '\t' + 'Bond_angle' + '\t' +  'distance_water_broken_atom' + '\t' + 'Bond_angle_water_target_si_broken_si' + '\t'  + 'Dissociation_mechanism'  + '\t'  + 'Stress' + '\t' + 'Shear' + '\n')
for lines3_1 in new_coordination_table:
    out_1.write('\t'.join(lines3_1) + '\n')



