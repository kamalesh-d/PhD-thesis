import os
import re

os.system('cp FIELD SOURCE_FIELD')
diff_cut_off = 0.05

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



x1_1 = open('all_distance_information.txt').readlines()
dist_info_field = []
for lines1_1 in x1_1:
    x1_2 = (lines1_1.split('\t'))
    x1_3 = x1_2[4]
    x1_4 = x1_2[0] + '_' + x1_2[1] + '_' + x1_2[2] + '_' + x1_2[3]
    dist_info_field.append(x1_3)
    x2_1 = open('FIELD').readlines()
    search_line = []
    rep_line = []
    x3_2 = open('FIELD_%s' %(x1_4), 'w')
    for lines2_1 in (x2_1):
        r1 = re.search(r'^pmf', lines2_1)
        if r1:
            x2_2 = (lines2_1.strip().split())
            if len(x2_2) == 2:
                search_line.append(lines2_1)
                x2_3 = '%.4f' %(float(x1_3))
                x2_4 = (num_after_point(diff_cut_off))
                x2_5 = except_zero(diff_cut_off)
                x2_6 = round(float(x2_3)*x2_5, x2_4-1) / x2_5
                replace_line = 'pmf        %s\n' %(x2_6)
                rep_line.append(replace_line)
    for lines2_2 in x2_1:
        modified_lines = (lines2_2.replace('%s' % (''.join(search_line)), '%s' % (''.join(rep_line))))
        x3_2.write(modified_lines)




