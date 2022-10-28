
num_line_con=4 #Each atom is defined by number of lines in config file has to be given here

atom_code = {'H': 2, 'O': 1, 'Si':3}
atom_mass = {'H': 1, 'O': 15.999, 'Si':28.085}


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

def curve_fitting(temperature, phi_value, viscousity):
    a = 3
    b = 2
    t0 = 2
    ced = (temperature * phi_value)/viscousity
    ced_1 = (ced + a) * t0
    return ced_1

def cell_length_finder(file):
    xclf1_1 = open(file).readlines()
    # for linesclf1_1 in (xclf1_1[2:5]):
        # print(linesclf1_1.strip().split())
    x_cell_length = xclf1_1[2].strip().split()[0]
    y_cell_length = xclf1_1[3].strip().split()[1]
    z_cell_length = xclf1_1[4].strip().split()[2]
    return [x_cell_length, y_cell_length, z_cell_length]


def config_to_lammps_format_input(coordinates_list, cell_length, output_filename):
    out_1 = open(output_filename, 'w')
    out_1.write('LAMMPS input data prepared from DLPOLY\n\n')
    out_1.write('%s atoms' %(len(coordinates_list)) + '\n\n')
    x1_1 = coordinates_list
    all_atom_types = []
    for linesd1_1 in x1_1:
        if linesd1_1[0] not in all_atom_types:
            all_atom_types.append(linesd1_1[0])

    out_1.write('%s atom types\n\n' %(len(all_atom_types)))
    out_1.write('0.0 %s xlo xhi\n' %(cell_length[0]))
    out_1.write('0.0 %s ylo yhi\n' %(cell_length[1]))
    out_1.write('0.0 %s zlo zhi\n' %(cell_length[2]))
    out_1.write('\n')
    out_1.write('Masses\n\n')

    for linesctl1_1 in (all_atom_types):
        atomic_mass = (atom_mass.get(linesctl1_1))
        atomic_code = atom_code.get(linesctl1_1)
        out_1.write('   %s  %s\n' %(atomic_code, atomic_mass))
    out_1.write('\n')
    out_1.write('Atoms # Charge\n\n')
    for lines1_1, lines1_2 in zip(x1_1, range(1, len(x1_1)+1)):
        out_1.write(str(lines1_2) + '       ' + str(atom_code.get(lines1_1[0])) + '       ' + str(0.0)  + '       ' + str(lines1_1[2]) + '       ' + str(lines1_1[3])  + '       ' + str(lines1_1[4]) + '\n')



x1_1 = config_input_processor('REVCON_V4')
cell_length = cell_length_finder('REVCON_V4')
config_to_lammps_format_input(x1_1, cell_length, 'input_structure.data')














