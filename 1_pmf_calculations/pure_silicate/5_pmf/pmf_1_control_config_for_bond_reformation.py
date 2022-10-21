import os
import glob

# x1_1 = open('all_distance_information.txt').readlines()
# for lines1_1 in x1_1:
#     x1_2 = lines1_1.strip().split('\t')
#     assembled = x1_2[0] + '_' + x1_2[1] + '_' + x1_2[2] + '_' + x1_2[3]
#     os.system('cp history_to_config_converter.py %s' %(assembled))
#     os.system('cd %s && python history_to_config_converter.py' %(assembled))
#     os.system('mkdir config_bond_reformation')
#     os.system('cd %s && mv CONFIG*_bond_reformation ../config_bond_reformation' %(assembled))
#     os.system('tar -czf config_bond_reformation.tar.gz config_bond_reformation')



os.system('mkdir config_bond_reformation')
x1_1 = open('bond_opening_distance_information').readlines()
for lines1_1 in x1_1:
#    out_1 = open('bond_open_dist_info', 'w')
    x1_2 = lines1_1.strip().split('\t')
    assembled = x1_2[0] + '_' + x1_2[1] + '_' + x1_2[2] + '_' + x1_2[3]
#    out_1.write(lines1_1)
#    os.system('mv bond_open_dist_info %s' %(assembled))
    os.system('cp history_to_config_converter.py %s' %(assembled))
    os.system('cd %s && python history_to_config_converter.py' %(assembled))
    os.system('cd %s && mv CONFIG*_bond_reformation ../config_bond_reformation' %(assembled))


x2_1 = glob.glob('./config_bond_reformation/*_bond_reformation')
for lines2_1 in x2_1:
    splitted2_0 = lines2_1.split('/')
    filename = splitted2_0[2]
    splitted2_1 = filename.split('_bond')
    new_filename = splitted2_1[0]
    old_filename = splitted2_1[0] + '_bond_reformation'
    os.system('cd config_bond_reformation && mv %s %s' %(old_filename, new_filename))

os.system('tar -czf config_bond_reformation.tar.gz config_bond_reformation')




