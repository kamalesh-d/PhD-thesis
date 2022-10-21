import os
import re


#Note: Files required in the directory to do PMF: CONFIG_O_*_Si_* , FIELD_O_*_Si_* , ./DLPOLY.Y, dist_info.txt, TABLE, CONTROL, PMF_1_master_glass.py, pmf_4_output_analyzer, run_python_script (A shell script for submitting jobs in ceres)

lower_cut_off = 1.50 # This defines the maximum closest distance to be taken between water and target atom
higher_cut_off = 3.80 # This defines the maximum longest distance to be taken between water and target atom
diff_cut_off = 0.05 # If you change diff_cut_off verify def update field %.2f


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
    if int(''.join(end_num)) == 1:
        return 2
    if int(''.join(end_num)) == 2:
        return 2



def update_field(field_file, distance_to_update):
    x10_1 = open('%s' %(field_file)).readlines()
    x10_2 = open('FIELD', 'w')
    search_line = []
    rep_line = []
    for lines10_1 in x10_1:
        r1 = re.search(r'^pmf', lines10_1)
        if r1:
            r2 = (lines10_1.strip().split())
            if len(r2) == 2:
                search_line.append(lines10_1)
                updated_value = float(distance_to_update)
                replace_line = 'pmf      %.2f\n' %(updated_value)
                rep_line.append(replace_line)
        modified_lines = (lines10_1.replace('%s' % (''.join(search_line)), '%s' % (''.join(rep_line))))
        x10_2.write(modified_lines)


x1_1 = open('dist_info.txt').readlines()
for lines1_1 in x1_1:
    x1_2 = lines1_1.split('\t')
    atom_type_num_info = x1_2[0] + '_' + x1_2[1] + '_' + x1_2[2] + '_' + x1_2[3]
    os.system('cp CONFIG_%s SOURCE_CONFIG' %(atom_type_num_info))
    os.system('cp FIELD_%s SOURCE_FIELD' %(atom_type_num_info))
    os.system('cp CONFIG_%s CONFIG' %(atom_type_num_info))
    os.system('cp FIELD_%s FIELD' %(atom_type_num_info))
    x2_1 = '%.4f' % (float(x1_2[4]))
    x2_2 = (num_after_point(diff_cut_off))
    x2_3 = except_zero(diff_cut_off)
    dist = round(float(x2_1) * x2_3, x2_2-1) / x2_3
    start_pos = round(float(dist+diff_cut_off) * x2_3, x2_2-1) / x2_3
    x2_4_pos = [x / 100.0 for x in range(int(start_pos*100), int(higher_cut_off*100+diff_cut_off*100), int(diff_cut_off*100))]
    x2_4_neg_1 = [x / 100.0 for x in range(int(lower_cut_off*100), int(dist*100+diff_cut_off*100), int(diff_cut_off*100))]
    x2_4_neg = reversed(x2_4_neg_1)


    for lines2_1 in x2_4_neg:
        field_info_to_change_n = '%.2f' %(lines2_1)
        x3_1 = atom_type_num_info + '_' + str(field_info_to_change_n)
        os.system('module load dl-poly/4.09_MahaGaro-intel')
        os.system('DLPOLY.Z')
        # os.system('cp REVCON CONFIG')
        os.system('cp FIELD FIELD_%s' %(x3_1))
        os.system('mv REVCON REVCON_%s' %(x3_1))
        os.system('mv HISTORY HISTORY_%s' %(x3_1))
        os.system('mv OUTPUT OUTPUT_%s' %(x3_1))
        os.system('mv RDFDAT RDFDAT_%s' %(x3_1))
        os.system('mv REVIVE REVIVE_%s' %(x3_1))
        os.system('mv STATIS STATIS_%s' %(x3_1))
        update_field('FIELD', float(field_info_to_change_n)-diff_cut_off)


    print('Completed_negative')
    os.system('cp SOURCE_CONFIG CONFIG')
    os.system('cp SOURCE_FIELD FIELD')
    os.system('cp CONFIG_%s CONFIG' %(atom_type_num_info))
    os.system('cp FIELD_%s FIELD' %(atom_type_num_info))
    print('Files copied')


    for lines2_2 in x2_4_pos:
        field_info_to_change_p = '%.2f' %(lines2_2)
        x3_1_p = atom_type_num_info + '_' + str(field_info_to_change_p)
        os.system('module load dl-poly/4.09_MahaGaro-intel')
        os.system('DLPOLY.Z')
        # os.system('cp REVCON CONFIG')
        os.system('cp FIELD FIELD_%s' %(x3_1_p))
        os.system('mv REVCON REVCON_%s' %(x3_1_p))
        os.system('mv HISTORY HISTORY_%s' %(x3_1_p))
        os.system('mv OUTPUT OUTPUT_%s' %(x3_1_p))
        os.system('mv RDFDAT RDFDAT_%s' %(x3_1_p))
        os.system('mv REVIVE REVIVE_%s' %(x3_1_p))
        os.system('mv STATIS STATIS_%s' %(x3_1_p))
        update_field('FIELD', float(field_info_to_change_p)+diff_cut_off)

