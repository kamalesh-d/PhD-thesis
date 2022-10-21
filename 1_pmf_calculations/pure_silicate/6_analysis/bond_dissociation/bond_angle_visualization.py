import numpy as np
import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import glob
from statistics import mean
import os
from scipy.stats import linregress
import seaborn as sns


num_line_con=2 #Each atom is defined by number of lines in config file has to be given here                                                               #
diff_cut_off = 0.05
lower_cut_off= 1.9
higher_cut_off = 3.1


def get_mean(list_values):
    mean_value = sum(list_values)/len(list_values)
    return mean_value


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


def integrate_for_pmf(diff_distance, force):
    a = np.array(force)
    m = np.r_[False, a > 0, False]
    idx = np.flatnonzero(m[:-1] != m[1:])
    I = (idx[1::2] - idx[::2]).argmax()
    elems = a[idx[2 * I]:idx[2 * I + 1]]
    d_x1_1 = sum(elems)
    d_x1_2 = d_x1_1*diff_distance
    return d_x1_2

# def Sort_based_on_column(sub_li):
def list_to_sort(list_of_list, column_to_sort, type_of_sort):
    if type_of_sort == 'ascending':
        list_of_list.sort(key=lambda x: x[column_to_sort-1])
    if type_of_sort == 'descending':
        list_of_list.sort(key=lambda x: x[column_to_sort-1], reverse=True)
    return list_of_list

def best_fit_slope(xs,ys):
    m = (((mean(xs)*mean(ys)) - mean(xs*ys)) /
         ((mean(xs)**2) - mean(xs**2)))
    return m






# compiled_coordination_table
current_directory = os.getcwd()
number_of_water_in_glass = []
coordination_table_contents = []
out_1 = open('compiled_coordination_activation_table', 'w')
#for lines0_1 in range(0,19):
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


