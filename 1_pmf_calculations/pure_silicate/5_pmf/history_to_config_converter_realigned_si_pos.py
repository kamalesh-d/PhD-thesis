import numpy as np
import re


distance_bet_si_oh = 1.6 # This is the distance between Si and O of H2O, whose history file you wish to convert
timestep_to_convert_history = 10000
num_line_con = 4


def remove_tail_dot_zeros(a):
    x10_1 = str(a)
    if (x10_1[-1]) == '0':
        return ('%.1f' %(float(x10_1)))
    if (x10_1[-1]) != '0':
        return ('%.2f' %(float(x10_1)))


def normalize_num(list_of_list):
    l = len(list_of_list)
    aligned_1 = []
    for lines6_1, lines6_2 in zip(list_of_list, range(0,l)):
        x6_1 = list_of_list[lines6_2][1]
        mod_atom_num = lines6_2 + 1
        lines6_1[1] = str(mod_atom_num)
        aligned_1.append(lines6_1)
    return aligned_1

def config_input_processor(input_file):
    x1_1 = input_file
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


def history_to_config_converter(input_file):
    x1_1 = open('%s' %(input_file)).readlines()
    out_1 = input_file.split('_')
    out_2 = open('CONFIG_' + str(out_1[1]) + '_' + str(out_1[2]) + '_' + str(out_1[3]) + '_' + str(out_1[4] + '_bond_reformation'), 'w')
    compiled_lines = []
    to_filter_lines = []
    for lines1_1 in x1_1[:2]:
        f1_1 = (lines1_1.split()[0:2])
        compiled_lines.append('    '.join(f1_1) + '   \n')
    for lines1_2 in x1_1[2:]:
        to_filter_lines.append(lines1_2)
    position_of_his_1000 = []
    for lines1_3 in to_filter_lines:
        r1 = re.search(r'^timestep', lines1_3)
        if r1:
            a1_1 = lines1_3.split()[1]
            if int(a1_1) == timestep_to_convert_history:
                a1_2 = to_filter_lines.index(lines1_3)
                position_of_his_1000.append(str(a1_2))
    timestep_line_num = int(''.join(position_of_his_1000))
    for lines1_5 in to_filter_lines[timestep_line_num+1:timestep_line_num+4]:
        compiled_lines.append(lines1_5)
    for lines1_6 in to_filter_lines[timestep_line_num+4:-16]:
        compiled_lines.append(lines1_6)
    for lines1_6_1 in to_filter_lines[-12:]:
        compiled_lines.append(lines1_6_1)
    for lines1_6 in to_filter_lines[-16:-12]:
        compiled_lines.append(lines1_6)
    compiled_lines_processed = config_input_processor(compiled_lines)
    compiled_lines_processed_aligned = normalize_num(compiled_lines_processed)
    new_compiled = compiled_lines[:1]
    (new_compiled.append('0  3  \n'))
    for lines3_1 in compiled_lines[2:5]:
        new_compiled.append(lines3_1)
    for lines3_2 in compiled_lines_processed_aligned:
        config_format = lines3_2[0] + '       ' + lines3_2[1] + '\n' + '   ' + lines3_2[2] + '   ' + lines3_2[3] + '   ' + lines3_2[4] + '\n'
        new_compiled.append(config_format)
    for lines1_7 in new_compiled:
        out_2.write(lines1_7)
    return compiled_lines



x1_1 = open('bond_open_dist_info').readlines()
all_vir_pmf = []
for lines1_1 in x1_1:
    x1_2 = lines1_1.split('\t')
    atom_type_num_info = x1_2[0] + '_' + x1_2[1] + '_' + x1_2[2] + '_' + x1_2[3]
    # x5_1 = remove_tail_dot_zeros('%.2f' %(lines2_1))
    filename = ('HISTORY_' + str(atom_type_num_info) + '_%.2f' %(float(x1_2[4])))
    x6_1 = history_to_config_converter(filename)

    
