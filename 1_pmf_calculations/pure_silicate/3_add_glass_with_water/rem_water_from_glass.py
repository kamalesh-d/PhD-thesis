import os
#import statistics

# os.system('cp CONFIG SOURCE_CONFIG')

###########################################################################################################################################################
num_line_con=2 #Each atom is defined by number of lines in config file has to be given here                                                               #
input_file_glass_water_name = 'CONFIG_GLASS_WATER' # Type the glass config file name with empty space for water addition                                  #
small_cut_off = 0.00                                                                                                                                      #
# first_atom = 'OW'  # Note: This must be specified with water, because the script is written to grab the neighbouring h2 atoms along with oxygen atom    #
# target_atom = 'Si' # User can change depending on the target towards which we will move the water molecule                                              #
# angstrom = 2.54 # Distance criteria between first atom and target atom. less than which we will select for performing PMF                               #
# You will get the output as CONFIG_*    where * represent numbers from 1 to n depending on number of atoms falling under the criteria we described       #
# Note: This python script is created assuming first 5 lines in input file has not started the first atom element co-ordinate description                 #
###########################################################################################################################################################

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

def config_output_processor(list_of_list_of_all):
    # List of list of all should have atom type, atom number, x_coord, y_coord, z_coord
    out_put_ready = []
    for lines30_1 in list_of_list_of_all:
        atom_type_o = lines30_1[0]
        atom_number_o = lines30_1[1]
        x_coord_o = lines30_1[2]
        y_coord_o = lines30_1[3]
        z_coord_o = lines30_1[4]
        compile_o = atom_type_o + '   ' + atom_number_o + '    ' + '\n  ' + x_coord_o + '     ' + y_coord_o + '    ' + z_coord_o + '    \n'
        out_put_ready.append(compile_o)
    return out_put_ready

def normalize_num(list_of_list):
    l = len(list_of_list)
    aligned_1 = []
    for lines6_1, lines6_2 in zip(list_of_list, range(0,l)):
        x6_1 = list_of_list[lines6_2][1]
        mod_atom_num = lines6_2 + 1
        lines6_1[1] = str(mod_atom_num)
        aligned_1.append(lines6_1)
    return aligned_1

def dist_calculator(set_1, set_2):
    x1 = float(set_1[2])
    y1 = float(set_1[3])
    z1 = float(set_1[4])
    x2 = float(set_2[2])
    y2 = float(set_2[3])
    z2 = float(set_2[4])
    a1 = (x2 - x1) ** 2
    b1 = (y2 - y1) ** 2
    c1 = (z2 - z1) ** 2
    c2 = a1 + b1 + c1
    distance = float(c2 ** 0.5)
    return distance

def lattice_finder(input_file):
    x20_1 = open(input_file).readlines()
    i_20 = 0
    lattice_info = []
    for lines20_1 in x20_1:
        i_20 = i_20 + 1
        if i_20 > 2 and i_20 < 6:
            lattice_info.append(lines20_1)
    lattice_info_x = float(lattice_info[0].split()[0])
    lattice_info_y = float(lattice_info[1].split()[1])
    lattice_info_z = float(lattice_info[2].split()[2])
    return [lattice_info_x, lattice_info_y, lattice_info_z]


x2_1 = input_file_glass_water_name
lat_water = lattice_finder(x2_1)
lat_cutoff_water_to_remove = lat_water[0]/2

x2_2 = config_input_processor(x2_1)
ow_atom_num_to_remove = []
for lines1_1 in x2_2:
    atom_type_4 = str(lines1_1[0])
    if atom_type_4 == 'OW':
        z_coord = abs(float(lines1_1[4]))
        if z_coord < lat_cutoff_water_to_remove + small_cut_off:
            ow_atom_num_to_remove.append(lines1_1)

h2o_atom_num_to_remove = []
for lines5_1 in ow_atom_num_to_remove:
    x5_1 = int(lines5_1[1])
    h2o_atom_num_to_remove.append(x5_1)
    h2o_atom_num_to_remove.append(x5_1+1)
    h2o_atom_num_to_remove.append(x5_1+2)


selected_h2o_glass = []
for lines4_1 in x2_2:
    atom_num_of_water_1 = int(lines4_1[1])
    if atom_num_of_water_1 not in h2o_atom_num_to_remove:
        selected_h2o_glass.append(lines4_1)

for lines6_1 in selected_h2o_glass:
    if lines6_1[0] == 'OW':
        lines6_1[0] = 'O'
    if lines6_1[0] == 'HW':
        lines6_1[0] = 'H'

x7_1 = open('CONFIG_GLASS_WATER_REFINED', 'w')
x7_2 = (normalize_num(selected_h2o_glass))
x7_3 = config_output_processor(x7_2)

i_1 = 0
for lines8_1 in open(x2_1).readlines():
    i_1 = i_1 + 1
    if i_1 < 6:
        x7_1.write(lines8_1)
for lines8_2 in x7_3:
    x7_1.write(lines8_2)

num_h2o = []
for lines10_1 in selected_h2o_glass:
    atom_type_10 = lines10_1[0]
    if atom_type_10 == 'H':
        num_h2o.append(atom_type_10)

print(len(num_h2o)/2)

