import re
import numpy as np
# import statistics
# import matplotlib
# matplotlib.use('TkAgg')
# import matplotlib.pyplot as plt


# Vir con being calculated here
diff_cut_off = 0.05
lower_cut_off= 1.5
higher_cut_off = 3.8



def remove_tail_dot_zeros(a):
    x10_1 = str(a)
    if (x10_1[-1]) == '0':
        return ('%.1f' %(float(x10_1)))
    if (x10_1[-1]) != '0':
        return ('%.2f' %(float(x10_1)))



def num_after_point(x):
    s = str(x)
    if not '.' in s:
        return 0
    return len(s) - s.index('.') - 1

def except_zero(x):
    x10_1 = str(x)
    end_num = []
    for lines10_1 in x10_1:
        if lines10_1 != '0' and lines10_1 != '.':
            end_num.append(lines10_1)
    if int(''.join(end_num)) == 5:
        return 2
    if int(''.join(end_num)) == 25:
        return 4

def output_file_processor(filename):
    x1_3_d = open('%s' %(filename)).readlines()
    i_1_d = 0
    lines_after_termination_d = []
    for lines2_1_d in x1_3_d:
        r1_d = re.search(r'run terminated after\b', lines2_1_d)
        if r1_d:
            i_1_d += 1
        if 0 < i_1_d:
            lines_after_termination_d.append(lines2_1_d)
    vir_pmf_d = lines_after_termination_d[12].split()[8]
    return float(vir_pmf_d)

def integrate_for_pmf(diff_distance, force):
    a = np.array(force)
    m = np.r_[False, a > 0, False]
    idx = np.flatnonzero(m[:-1] != m[1:])
    I = (idx[1::2] - idx[::2]).argmax()
    elems = a[idx[2 * I]:idx[2 * I + 1]]
    d_x1_1 = sum(elems)
    d_x1_2 = d_x1_1*diff_distance
    return d_x1_2



x1_1 = open('dist_info.txt').readlines()
out_1 = open('result.txt', 'w')
all_vir_pmf = []
for lines1_1 in x1_1:
    x1_2 = lines1_1.split('\t')
    atom_type_num_info = x1_2[0] + '_' + x1_2[1] + '_' + x1_2[2] + '_' + x1_2[3]
    x2_1 = np.arange(lower_cut_off, higher_cut_off+diff_cut_off, diff_cut_off)
    all_vir_pmf_bet_two_atoms = {}
    file_name = ('Figure_'+str(atom_type_num_info))
    x_axis = np.arange(lower_cut_off, higher_cut_off+diff_cut_off, diff_cut_off)
    y_axis = []

# Individual figures for every pmf method
    for lines2_1 in (x2_1):
        x5_1 = ('%.2f' %(lines2_1))
        x2_3 = output_file_processor('OUTPUT_'+str(atom_type_num_info) + '_' + x5_1)
        x2_4 = (-1)*(x2_3)/float(x5_1)
        y_axis.append(x2_4)
    for lines4_1, lines4_2 in zip(x_axis, y_axis):
        out_1.write(str('%.2f' %(lines4_1)) + '\t' + str(lines4_2) + '\n')
