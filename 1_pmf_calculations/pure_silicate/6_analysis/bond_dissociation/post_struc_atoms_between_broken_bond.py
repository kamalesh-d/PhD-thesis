import glob
import re
import os
import numpy as np

# Feature: Any particular timestep of the history can be chosen to fetch the exact number of atoms within that particular timestep

distance_bet_si_oh = 1.6 # This is the distance between Si and O of H2O, whose history file you wish to convert
name_water_against_target_atom = 'Si'
timestep_to_convert_history = 10000
diff_cut_off = 0.05
pmf_atoms_num = 4 # This is to define number of atoms are pmf in this script. Only the last number of atoms (That you define) will be considered
neighbour_grab_cutoff = 5.0
bond_pair_distance = {'H_O':1.2, 'O_H':1.2, 'Si_O':2.0, 'O_Si':2.0, 'O_O':3, 'H_H':0.0, 'Si_H':0.0, 'H_Si':0.0, 'Si_Si':3.5, 'Al_O':2.4, 'O_Al':2.4, 'Ca_O':3.24, 'O_Ca':3.24, 'Al_Al':0, 'Ca_Al':0, 'H_Al':0, 'Si_Al':0, 'Al_Si':0, 'Ca_Si':0}
num_line_con = 2

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


def target_atom_environment(initial_bonded_config, final_broken_config, cell_length, target_atom_number):
    initial_target_atom = []
    for lines4_1 in initial_bonded_config:
        if lines4_1[1] == '%s' %(target_atom_number):
            for lines4_1_1 in lines4_1:
                initial_target_atom.append(lines4_1_1)

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
                    if lines1_3 != initial_bonded_config[-4]:
                        if lines1_3[0] == 'Si':
                            initial_bonded_si.append(lines1_3)
                        if lines1_3[0] == 'H':
                            initial_bonded_h.append(lines1_3)
                        if lines1_3[0] == 'Al':
                            initial_bonded_al.append(lines1_3)

    final_target_atom = []
    for lines4_2 in final_broken_config:
        if lines4_2[1] == '%s' %(target_atom_number):
            for lines4_2_1 in lines4_2:
                final_target_atom.append(lines4_2_1)

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
                    if lines2_3 != final_broken_config[-4]:
                        if lines2_3[0] == 'Si':
                            final_bonded_si.append(lines2_3)
                        if lines2_3[0] == 'H':
                            final_bonded_h.append(lines2_3)
                        if lines2_3[0] == 'Al':
                            final_bonded_al.append(lines2_3)


    initial_coordination_pattern = 'O%s ' %(len(oxy_connected_to_bonded_target_atom)) + 'Al%s ' %(len(initial_bonded_al)) + 'Si%s ' %(len(initial_bonded_si)) + 'H%s ' %(len(initial_bonded_h))
    final_coordination_pattern = 'O%s ' %(len(oxy_connected_to_final_target_atom)) + 'Al%s ' %(len(final_bonded_al)) + 'Si%s ' %(len(final_bonded_si)) + 'H%s ' %(len(final_bonded_h))

    # print(initial_coordination_pattern)
    return (initial_coordination_pattern + '\t' + final_coordination_pattern)



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


def bond_angle_finder(config_file_name, list_of_three_atoms):
    config_processed = config_input_processor(config_file_name)
    three_points = []
    baf1_1 = open('%s' %(config_file_name)).readlines()
    x_length = (baf1_1[2].split()[0])
    y_length = (baf1_1[3].split()[1])
    z_length = (baf1_1[4].split()[2])
    cell_length = [x_length, y_length, z_length]
    for lines1_1 in config_processed:
        if (lines1_1[1]) == (list_of_three_atoms[0]):
            three_points.append(lines1_1)
        if (lines1_1[1]) == (list_of_three_atoms[1]):
            three_points.append(lines1_1)
        if (lines1_1[1]) == (list_of_three_atoms[2]):
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

def dissociation_mechanism_finder(bonded_conf, broken_conf, list_of_broken_bonds, cell_length):
    broke_bond_pattern = []
    for linesd1_21 in list_of_broken_bonds:
        broken_atom_environment = target_atom_environment(bonded_conf, broken_conf, cell_length, linesd1_21.split('_')[5])
        position_of_h = []
        for linesd1_24 in range(len(broken_atom_environment.split('\t')[0].strip().split())):
            if (broken_atom_environment.split('\t')[0].strip().split()[linesd1_24][0]) == 'H':
                position_of_h.append(linesd1_24)
        initial_broken_atom_env = broken_atom_environment.split('\t')[0].strip().split()[position_of_h[0]]
        final_broken_atom_env = broken_atom_environment.split('\t')[1].strip().split()[position_of_h[0]]
        broke_bond_pattern.append(int(final_broken_atom_env[-1]) - int(initial_broken_atom_env[-1]))
    if len(broke_bond_pattern) == 1:
        if broke_bond_pattern[0] == 0:
            return 'only_indirect_mech'
        if broke_bond_pattern[0] >= 1:
            return 'only_direct_mech'
        else:
            return 'strange_mech'
    if len(broke_bond_pattern) > 1:
        number_of_broken_bonds = len(broke_bond_pattern)
        dir = []
        indir = []
        for linesd1_22 in broke_bond_pattern:
            if int(linesd1_22) >= 1:
                dir.append(linesd1_22)
        if len(dir) > 0 and len(indir) > 0:
            return 'direct_and_indirect_mech'
        if len(dir) > 0 and len(indir) == 0:
            return 'only_direct_mech'
        if len(indir) > 0 and len(dir) == 0:
            return 'only_indirect_mech'



x1_1 = open('all_distance_information.txt').readlines()
out_1 = open('coordination_activation_table', 'w')
list_broken_bond_info = []
out_1.write('Glass_number' + '\t' + 'Filename' + '\t' + 'Before_breaking' + '\t' + 'After_breaking' + '\t' + 'x_coordinate' + '\t' + 'y_coordinate' + '\t' + 'z_coordinate' + '\t'  + 'Activation_barrier' + '\t' + '#_broken_bonds' + '\t' + 'List_of_broken_bond' + '\t' + 'Bond_angle' + '\t' + 'distance_water_broken_atom' + '\t' + 'Bond_angle_water_target_si_broken_si' + '\t' + 'Dissociation_mechanism' + '\n')
for lines1_1 in x1_1:
    x1_2 = ('_'.join(lines1_1.split()[:-1]))
    os.system('cd %s && cp HISTORY*1.50 ..' %(x1_2))
    initial_filename = 'HISTORY_%s_1.50' %(x1_2)
    final_filename = 'HISTORY_%s_1.50' %(x1_2)

    x2_1 = open(initial_filename).readlines()
    cell_length = []
    cell_length.append(x2_1[3].split()[0])
    cell_length.append(x2_1[4].split()[1])
    cell_length.append(x2_1[5].split()[2])

    broken_conf = history_to_config_converter(final_filename, 10000)
    bonded_conf = history_to_config_converter(initial_filename, 0)
    broken_conf_list = []
    x6_1 = broken_bond_atoms_finder(bonded_conf, broken_conf, cell_length)
    x7_1 = target_atom_environment(bonded_conf, broken_conf, cell_length, bonded_conf[-4][1])

    if len(x6_1) >= 1:
        dissociation_mechanism = dissociation_mechanism_finder(bonded_conf, broken_conf, x6_1, cell_length)
    if len(x6_1) == 0:
        dissociation_mechanism = 'N/A'

    activation_barrier_file = open('all_results/result_%s.txt' %(x1_2)).readlines()
    force = []
    for lines41_1 in activation_barrier_file:
        force.append(float(lines41_1.split('\t')[1].strip()))
    activation_barrier = integrate_for_pmf(diff_cut_off, force)



    broken_atom_number = []
    initial_oxy_of_water = bonded_conf[-3]
    target_si = (bonded_conf[-4])
    initial_broken_atom = []
    water_target_si_broken_si = []
    if len(x6_1) >= 1:
        for lines3_1 in x6_1:
            broken_atom_number.append(lines3_1.split('_')[5])
        for lines3_2 in bonded_conf:
            if lines3_2[1] == broken_atom_number[0]:
                initial_broken_atom.append(lines3_2)
        dist_water_broken_si = dist_calculator_pbc(initial_oxy_of_water, initial_broken_atom[0], cell_length)
        water_target_si_broken_si.append(initial_oxy_of_water[1])
        water_target_si_broken_si.append(target_si[1])
        water_target_si_broken_si.append(initial_broken_atom[0][1])
    else:
        dist_water_broken_si = 'N/A'



    x8_1 = (str(initial_filename) + '\t' + str('\t'.join(x6_1)) + '\n')
    file_info = x8_1.strip().split('\t')[0][8:-5]
    os.system('cd %s && cp CONFIG_%s ..' %(file_info, file_info))
    input_filename_1 = 'CONFIG_%s' %(file_info)
    broken_bond_1 = (x8_1.strip().split('\t'))

    if len(broken_bond_1) >= 2:
        broken_bond_atoms = broken_bond_1[1]
        bond_angle_water_broken_si = bond_angle_finder(input_filename_1, water_target_si_broken_si)
        bond_angle = bond_angle_finder(input_filename_1, [broken_bond_atoms.split('_')[1], broken_bond_atoms.split('_')[3], broken_bond_atoms.split('_')[5]])
        current_dir = os.getcwd().split('/')
        glass_num = (current_dir[-2])
        out_1.write('%s' %(glass_num) + '\t' + str(initial_filename[8:-5]) + '\t' + x7_1 + '\t' + str(bonded_conf[-4][2]) + '\t' + str(bonded_conf[-4][3]) + '\t' + str(bonded_conf[-4][4]) + '\t' + str(activation_barrier) + '\t' + str(len(x6_1)) + '\t' + str(x6_1) + '\t' + str(bond_angle) + '\t' + str(dist_water_broken_si) + '\t' + str(bond_angle_water_broken_si) + '\t' + str(dissociation_mechanism) + '\n')


