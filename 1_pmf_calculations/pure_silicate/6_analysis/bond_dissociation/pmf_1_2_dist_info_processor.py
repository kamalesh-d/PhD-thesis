from operator import itemgetter
import os

# lowest_cut_off = 2.0 #Only above this cut-off, we will start considering the water molecule is full

x1_1 = open('pre_dist_info.txt').readlines()
dict1_1 = {}
for lines1_1 in x1_1:
    splitted = lines1_1.strip().split('\t')
    water = splitted[0] + '_' + splitted[1]
    dict1_1.setdefault(water, []).append(splitted)


first_curation = []
for lines4_1, lines4_2 in dict1_1.items():
    if len(lines4_2) == 1:
        for lines4_3 in lines4_2:
            first_curation.append('\t'.join(lines4_3)+'\n')

    if len(lines4_2) >= 2:
        sorted_2 = sorted(lines4_2, key=itemgetter(4))
        for lines5_1 in sorted_2[:1]:
            first_curation.append('\t'.join(lines5_1)+'\n')

dict5_1 = {}
for lines5_1 in first_curation:
    splitted5 = lines5_1.strip().split('\t')
    water5 = splitted5[2] + '_' + splitted5[3]
    dict5_1.setdefault(water5, []).append(splitted5)


out_1 = open('all_distance_information.txt', 'w')
for lines6_1, lines6_2 in dict5_1.items():
    if len(lines6_2) == 1:
        for lines6_3 in lines6_2:
            out_1.write('\t'.join(lines6_3)+'\n')
    if len(lines6_2) >= 2:
        sorted6_2 = sorted(lines6_2, key=itemgetter(4))
        for lines6_1 in sorted6_2[:1]:
            out_1.write('\t'.join(lines6_1)+'\n')

