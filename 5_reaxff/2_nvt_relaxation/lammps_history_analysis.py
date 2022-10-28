from operator import itemgetter


def lammps_history_split_basedon_timestep(filename):
    xd1_1 = open(filename).readlines()
    item_timestep_positions = []
    for linesd1_1 in range(len(xd1_1)):
        if xd1_1[linesd1_1] == 'ITEM: TIMESTEP\n':
            item_timestep_positions.append(linesd1_1)
    list_of_list_file = []
    for linesa1_1 in range(len(item_timestep_positions)):
        try:
            x2_3 = xd1_1[item_timestep_positions[linesa1_1]:item_timestep_positions[linesa1_1+1]]
            list_of_list_file.append(x2_3)
        except:
            x2_4 = xd1_1[item_timestep_positions[linesa1_1]:]
            list_of_list_file.append(x2_4)
    return list_of_list_file


def extract_info_from_lammps_history_list_of_list(list_of_list_history_lammps):
    dict_history_informations = {}

    for lines1_1 in list_of_list_history_lammps:
        dict_history_informations['Timestep'] = int(lines1_1[1].strip())
        dict_history_informations['Number_of_atoms'] = int(lines1_1[3].strip())
        dict_history_informations['X_cell_length'] = lines1_1[5].strip()
        dict_history_informations['Y_cell_length'] = lines1_1[6].strip()
        dict_history_informations['Z_cell_length'] = lines1_1[7].strip()
        dict_history_informations['All_coordinates'] = list_of_list_history_lammps[0][9:]

    return dict_history_informations

def rearrange_list(list_to_rearrange, lines_of_list_seperated_by):
    splitted_list = []
    for linesdr1_1 in list_to_rearrange:
        if lines_of_list_seperated_by == '':
            splitted_list.append([int(linesdr1_1.split()[0]), linesdr1_1.split()[1], linesdr1_1.split()[2], linesdr1_1.split()[3], linesdr1_1.split()[4], linesdr1_1.split()[5], linesdr1_1.split()[6]])
        else:
            splitted_list.append([int(linesdr1_1.split(lines_of_list_seperated_by)[0]), linesdr1_1.split(lines_of_list_seperated_by)[1], linesdr1_1.split(lines_of_list_seperated_by)[2], linesdr1_1.split(lines_of_list_seperated_by)[3], linesdr1_1.split(lines_of_list_seperated_by)[4], linesdr1_1.split(lines_of_list_seperated_by)[5], linesdr1_1.split(lines_of_list_seperated_by)[6]])

    first_item = itemgetter(0)
    return (sorted(splitted_list, key = first_item))


def last_history_to_lammps_input(filename):
    x2_1 = lammps_history_split_basedon_timestep(filename)
    list_of_history_to_be_converted = x2_1[-1]
    dict_lammps_history = extract_info_from_lammps_history_list_of_list([list_of_history_to_be_converted])
    input_for_lammps = []
    input_for_lammps.append('LAMMPS input structure data file written by python\n\n')
    input_for_lammps.append(str(dict_lammps_history.get('Number_of_atoms')) + ' atoms\n\n')


    # # Computing number of atom types
    all_atom_types_with_mass = {}
    for linesd2_1 in dict_lammps_history.get('All_coordinates'):
        all_atom_types_with_mass[int(linesd2_1.split()[1])] = linesd2_1.split()[6]

    all_atom_types_with_mass_sorted = {k: all_atom_types_with_mass[k] for k in sorted(all_atom_types_with_mass)}
    input_for_lammps.append('%s atom types\n\n' %(len(all_atom_types_with_mass_sorted)))

    input_for_lammps.append('%s xlo xhi\n' %(dict_lammps_history.get('X_cell_length')))
    input_for_lammps.append('%s ylo yhi\n' %(dict_lammps_history.get('Y_cell_length')))
    input_for_lammps.append('%s zlo zhi\n\n' %(dict_lammps_history.get('Z_cell_length')))

    input_for_lammps.append('Masses\n\n')

    for linesd2_2, linesd2_3 in all_atom_types_with_mass_sorted.items():
        input_for_lammps.append('  %s  %s\n' %(linesd2_2, linesd2_3))

    input_for_lammps.append('\n')

    input_for_lammps.append('Atoms # charge\n\n')
    for splitted2_4 in (rearrange_list(dict_lammps_history.get('All_coordinates'), '')):
        input_for_lammps.append('%s       %s       %s       %s       %s       %s\n' %(splitted2_4[0], splitted2_4[1], splitted2_4[2], splitted2_4[3], splitted2_4[4], splitted2_4[5]))
    return input_for_lammps

x1_1 = last_history_to_lammps_input('analysis_nvt.dump')
out_1 = open('input_structure_2.data', 'w')
for lines1_1 in x1_1:
    out_1.write(lines1_1)

