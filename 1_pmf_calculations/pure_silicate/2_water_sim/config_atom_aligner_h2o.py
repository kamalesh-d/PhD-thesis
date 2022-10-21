
############
num_line_con=4 #Each atom is defined by number of lines in config file has to be given here
# list_of_atoms = ['Si', 'B', 'O', 'Na'] #List of atoms to be specified in the format you need to obtain in output
input_file_name = 'CONFIG_H2O' # Type the file name to be processed
# You will get the output as pro_config
# Note: This python script is created assuming first 5 lines in input file has not started the first atom element co-ordinate description
############




def config_processor(input_file):
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
        # atom_mass = lines1_3[0].split()[2]
        x_coordinate = lines1_3[1].split()[0]
        y_coordinate = lines1_3[1].split()[1]
        z_coordinate = lines1_3[1].split()[2]
        compile = atom_type + '   ' + atom_number + '   ' + '\n'  + x_coordinate + '   ' + y_coordinate + '   ' + z_coordinate
        compiled.append(compile.strip().split('   '))
    return compiled
def self_sort(list_of_list):
    l = len(list_of_list)
    aligned = []
    aligned_sorted = []
    # for lines3_1 in list_of_list:
    first_set = []
    h_atom_set = []
    second_set = []
    third_set = []
    i_3 = 0
    for lines3_1 in list_of_list:
        i_3 = i_3 +1
        if i_3 <= (l/3):
            first_set.append(lines3_1)
        if i_3 > (l/3):
            h_atom_set.append(lines3_1)
    h_o = range(0,len(h_atom_set),2)
    h_e = range(1,len(h_atom_set),2)
    for lines3_2 in h_o:
        second_set.append(h_atom_set[lines3_2])
    for lines3_3 in h_e:
        third_set.append(h_atom_set[lines3_3])
    for lines4_1, lines4_2, lines4_3 in zip(first_set, second_set, third_set):
        aligned.append(lines4_1)
        aligned.append(lines4_2)
        aligned.append(lines4_3)
    return aligned


x1_0 = open(input_file_name).readlines()
x1_1 = config_processor(input_file_name)
x3_1 = self_sort(x1_1)
l = len(x1_1)
x2_1 = open('CONFIG', 'w')
for lines3_1 in x1_0[:5]:
    x2_1.write(lines3_1)
for lines5_2, lines5_3 in zip(x3_1, range(1, l + 1)):
    (lines5_2[1]) = str(lines5_3)
    x2_1.write('   '.join(lines5_2) + '\n')
