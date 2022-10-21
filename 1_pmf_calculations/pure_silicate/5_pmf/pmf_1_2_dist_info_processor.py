from operator import itemgetter
import os

# lowest_cut_off = 2.0 #Only above this cut-off, we will start considering the water molecule is full

x1_1 = open('pre_dist_info.txt').readlines()
dict1_1 = {}
to_be_removed = []
for lines1_1 in x1_1:
    splitted = lines1_1.strip().split('\t')
    water = splitted[0] + '_' + splitted[1]
    dict1_1.setdefault(water, []).append(splitted)

for lines2_1, lines2_2 in dict1_1.items():
    if len(lines2_2) >= 2:
        sorted_1 = sorted(lines2_2, key=itemgetter(4))
        for lines3_1 in sorted_1[1:]:
            to_be_removed.append('\t'.join(lines3_1))

out_1 = open('all_distance_information.txt', 'w')
for lines4_1, lines4_2 in dict1_1.items():
    if len(lines4_2) == 1:
        for lines4_3 in lines4_2:
            out_1.write('\t'.join(lines4_3)+'\n')
    if len(lines4_2) >= 2:
        sorted_2 = sorted(lines4_2, key=itemgetter(4))
        for lines5_1 in sorted_2[:1]:
            out_1.write('\t'.join(lines5_1)+'\n')

