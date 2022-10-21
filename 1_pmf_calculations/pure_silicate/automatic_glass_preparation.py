import os
import re
import subprocess
import random

"""
This script is a master script (pipeline), which will be used to control multiple scripts step by step to perform the PMF calculations.
For instance, when we choose 'glass_preparation' in 'type_of_run', it will compile all the scripts to prepare the glass. 
All options for selection in type_of_run are: glass_preparation, double_the_box, water_simulation, add_glass_with_water, relax_glass_water, run_pmf, result_analysis, bond_reformation, aluminum_run_pmf, neural_network_input_preparation, neural_network_input_compiled
"""

density = {'water': 0.9970, 'silicate': 2.2}
number_of_glasses = 2

def find_and_replace(input_file, starting_word_in_line, replace_line, file_directory_location, output_file_name):
    xd1_1 = open('%s/%s' %(file_directory_location, input_file)).readlines()
    out_1 = open('%s/%s' %(file_directory_location, output_file_name), 'w')
    rep_line_pos = []
    for lines1_1 in xd1_1:
        r1 = re.search(r'\b%s' %(starting_word_in_line), lines1_1)
        if r1:
            rep_line_pos.append(str(xd1_1.index(lines1_1)))
    xd1_1[int(''.join(rep_line_pos))] = replace_line + '\n'
    for lines2_1 in xd1_1:
        out_1.write(lines2_1)

def find_pattern_and_replace(input_file, pattern_in_target_line, replace_line, file_directory_location, output_file_name):
    xd2_1 = open('%s/%s' %(file_directory_location, input_file)).readlines()
    outd2_1 = open('%s/%s' %(file_directory_location, output_file_name), 'w')
    replace_line_pos_1 = []
    for linesd2_1 in xd2_1:
        rd2 = re.search(r'%s' %(pattern_in_target_line), linesd2_1)
        if rd2:
            replace_line_pos_1.append(str(xd2_1.index(linesd2_1)))
    if replace_line_pos_1 != []:
        xd2_1[int(''.join(replace_line_pos_1))] = replace_line + '      %s' %(pattern_in_target_line) + '\n'
        for linesd2_2 in xd2_1:
            outd2_1.write(linesd2_2)

def get_cell_length(input_file):
    xd3_1 = open('%s' %(input_file)).readlines()
    cell_length = []
    cell_length.append(xd3_1[2].split()[0])
    cell_length.append(xd3_1[3].split()[1])
    cell_length.append(xd3_1[4].split()[2])
    return cell_length

def number_of_molecules_finder_from_volume(cell_length_list, molecule_type):
    volume_of_system = float(cell_length_list[0]) * float(cell_length_list[1]) * float(cell_length_list[2])
    mass = 2.99118047
    density_of_molecule = density.get(molecule_type)
    num_of_molecule = (volume_of_system * density_of_molecule * 0.1)/mass
    return int(num_of_molecule)

def fortran_script_glass_water_mixer(input_file_name, output_file_name, directory_location, cut_off_updated):
    xgw1_1 = open('%s/%s' %(directory_location, input_file_name)).readlines()
    outgw1 = open('%s/%s' %(directory_location, output_file_name), 'w')
    gw_line_to_find = []
    gw_line_to_replace = []
    for linesgw1 in xgw1_1:
        gwr1 = re.search('cut_off_changing_location', linesgw1)
        if gwr1:
            gw_line_to_find.append(linesgw1)
            gwr1_splitted = linesgw1.split()
            gwr1_splitted[0] = '      rcut=%sd0' %(cut_off_updated)
            gw_line_to_replace.append('    '.join(gwr1_splitted) + '\n')

    for linesgw_2, linesgw_3 in zip(gw_line_to_find, gw_line_to_replace):
        for gwn, gwi in enumerate(xgw1_1):
            if gwi == linesgw_2:
                xgw1_1[gwn] = linesgw_3
    for linesgw_4 in xgw1_1:
        outgw1.write(linesgw_4)


def fortran_script_water_simulation_modifier(input_file_name, output_file_name, directory_location):
    x4_1 = open('%s/%s' %(directory_location, input_file_name)).readlines()
    out4_1 = open('%s/%s' %(directory_location, output_file_name), 'w')
    line_to_find = []
    line_to_replace = []
    x5_1 = number_of_molecules_finder_from_volume(box_refined_cell_length, 'water')
    for lines4_1 in x4_1:
        r2 = re.search(r'x_coordinate_to_be_modified', lines4_1)
        if r2:
            r2_splitted = (lines4_1.split())
            r2_splitted[0] = '        ZL_x=%s' %(box_refined_cell_length[0])
            line_to_find.append(lines4_1)
            line_to_replace.append('    '.join(r2_splitted) + '\n')
        r3 = re.search(r'y_coordinate_to_be_modified', lines4_1)
        if r3:
            r3_splitted = (lines4_1.split())
            r3_splitted[0] = '        ZL_y=%s' % (box_refined_cell_length[1])
            line_to_find.append(lines4_1)
            line_to_replace.append('    '.join(r3_splitted) + '\n')
        r4 = re.search(r'z_coordinate_to_be_modified', lines4_1)
        if r4:
            r4_splitted = (lines4_1.split())
            r4_splitted[0] = '        ZL_z=%s' % (box_refined_cell_length[2])
            line_to_find.append(lines4_1)
            line_to_replace.append('    '.join(r4_splitted) + '\n')
        r5 = re.search(r'number_of_oxygen_atoms', lines4_1)
        if r5:
            r5_splitted = lines4_1.split()
            r5_splitted[0] = '        NA(3)=%s' %(x5_1)
            line_to_find.append(lines4_1)
            line_to_replace.append('    '.join(r5_splitted) + '\n')
        r6 = re.search(r'number_of_hydrogen_atoms', lines4_1)
        if r6:
            r6_splitted = lines4_1.split()
            r6_splitted[0] = '        NA(4)=%s' %(x5_1*2)
            line_to_find.append(lines4_1)
            line_to_replace.append('    '.join(r6_splitted) + '\n')
        r7 = re.search(r'total_number_of_water_molecules', lines4_1)
        if r7:
            r7_splitted = lines4_1.split()
            r7_splitted[0] = '      PARAMETER(IM=%s)' %(x5_1*3)
            line_to_find.append(lines4_1)
            line_to_replace.append('    '.join(r7_splitted) + '\n')

    for lines4_2, lines4_3 in zip(line_to_find, line_to_replace):
        for n, i in enumerate(x4_1):
            if i == lines4_2:
                x4_1[n] = lines4_3

    for lines4_4 in x4_1:
        out4_1.write(lines4_4)

def get_num_water_in_glass(glass_number, file_name_location):
    xd6_1 = open(file_name_location).readlines() 
    for linesd6_1 in xd6_1:
        if (linesd6_1.strip().split('\t')[0] == glass_number):
            return (linesd6_1.strip().split('\t')[1])

type_of_run = 'neural_network_input_compiled' #keywords - glass_preparation, double_the_box, water_simulation, add_glass_with_water, relax_glass_water, run_pmf, result_analysis, aluminum_run_pmf, neural_network_input_preparation, neural_network_input_compiled
number_of_glasses_to_make = 20

list_nnic = []
out21_1 = open('glass_water_mixing_problem', 'w')
for lines0_1 in range(1, number_of_glasses_to_make+1):
    os.system('module load python/3.6.10')
    os.system('module load dl-poly/4.09_MahaGaro-intel')
    glass_name_number = 'glass_num_%s' % (lines0_1)
    directory = os.getcwd()
    number_of_water_in_glass = []
    current_glass_directory = '%s/%s' % (directory, glass_name_number)
    if type_of_run == 'glass_preparation':
        #Glass preparation
        os.system('mkdir %s' %(glass_name_number))
        x32_1 = random.random()*50000
        os.system('cp config_aleatoire.f %s' %(current_glass_directory))
        find_pattern_and_replace('config_aleatoire.f', '!number_to_be_changed_in_config_aleatoire' , '       DO K=1,%s' %(int(x32_1)), current_glass_directory, 'config_aleatoire.f')
        os.system('cp -r 1_box_relax 2_water_sim 3_add_glass_with_water 4_equilibrate_sol_glass_2 FIELD glasscript.sh TABL* TEMP* %s' %(glass_name_number))
        os.system('cd %s && gfortran config_aleatoire.f && ./a.out' %(current_glass_directory))
        os.system('cd %s && cp CONFIG_CJ1 CONFIG' %(current_glass_directory))
        os.system('cd %s && chmod +x glasscript.sh' %(current_glass_directory))
        os.system('cd %s && sbatch -n 1 ./glasscript.sh' %(current_glass_directory))


    elif type_of_run == 'double_the_box':
        #Double the size of the box along z-axis
        x2_1 = open('%s/REVCON_V4' %(current_glass_directory)).readlines()
        out2_1 = open('%s/1_box_relax/CONFIG' %(current_glass_directory), 'w')
        z_dimension = []
        for lines2_1 in x2_1[4:5]:
            z_dimension.append('        0.0000000000')
            z_dimension.append('0.0000000000')
            z_dimension.append(str(float(lines2_1.split()[2])*2))
        box_shape_info_updated = []
        for lines2_2 in x2_1[1:2]:
            x2_2 = (lines2_2.split())
            x2_2[1] = '3'
            box_shape_info_updated.append('       '.join(x2_2))
        x2_1[1] = ('    '.join(box_shape_info_updated) + '\n')
        x2_1[4] = ('       '.join(z_dimension) + '\n')
        for lines3_1 in x2_1:
            out2_1.write(lines3_1)

        os.system('cd %s && cp FIELD TABLE 1_box_relax' %(current_glass_directory))
        os.system('cd %s/1_box_relax && chmod +x glasscript_P1.sh' %(current_glass_directory))
        os.system('cd %s/1_box_relax && sbatch glasscript_P1.sh' %(current_glass_directory))

    elif type_of_run == 'water_simulation':
        #2_water_simulation
        os.system('cd %s/1_box_relax && cp REVCON_V4 CONFIG_GLASS' %(current_glass_directory))
        os.system('cd %s/1_box_relax && cp CONFIG_GLASS ../3_add_glass_with_water' %(current_glass_directory))
        box_refined_cell_length = get_cell_length('%s/1_box_relax/REVCON_V4' %(current_glass_directory))
        number_of_pure_water_molecules = number_of_molecules_finder_from_volume(box_refined_cell_length, 'water')
        find_pattern_and_replace('FIELD_pre', '#number_of_water_molecules', 'NUMMOLS %s' %(number_of_pure_water_molecules), '%s/2_water_sim' %(current_glass_directory), 'FIELD')
        x25_1 = fortran_script_water_simulation_modifier('config_aleatoire_H2O_220319.f', 'config_aleatoire_H2O_updated.f', '%s/2_water_sim/' %(current_glass_directory))
        os.system('cd %s/2_water_sim && gfortran config_aleatoire_H2O_updated.f' %(current_glass_directory))
        os.system('cd %s/2_water_sim && ./a.out && python config_atom_aligner_h2o.py' %(current_glass_directory))
        os.system('cd %s/2_water_sim && chmod +x glasscript_P1.sh' %(current_glass_directory))
        os.system('cd %s/2_water_sim && sbatch ./glasscript_P1.sh' %(current_glass_directory))


    elif type_of_run == 'add_glass_with_water':
        #3_add_glass_with_water
        os.system('cp -r 3_add_glass_with_water 4_equilibrate_sol_glass_2 %s' %(current_glass_directory))
        os.system('cd %s/2_water_sim && cp REVCON_V4 ../3_add_glass_with_water/CONFIG_WATER' %(current_glass_directory))
        os.system('cd %s/1_box_relax && cp REVCON_V4 ../3_add_glass_with_water/CONFIG_GLASS' %(current_glass_directory))
        fortran_script_glass_water_mixer('placement_H2O.f', 'placement_H2O_updated.f', '%s/3_add_glass_with_water' %(current_glass_directory), 0.1)
        os.system('cd %s/3_add_glass_with_water && gfortran placement_H2O_updated.f' %(current_glass_directory))
        os.system('cd %s/3_add_glass_with_water && ./a.out -filepos1 CONFIG_GLASS -filepos2 CONFIG_WATER' %(current_glass_directory))
        os.system('cd %s/3_add_glass_with_water && python rem_water_from_glass.py' %(current_glass_directory))
        output = subprocess.Popen('cd %s/3_add_glass_with_water && python rem_water_from_glass.py' %(current_glass_directory), stdout=subprocess.PIPE, shell=True)
        (num_water_in_glass, err) = output.communicate()
        box_refined_cell_length = get_cell_length('%s/1_box_relax/REVCON_V4' %(current_glass_directory))
        number_of_pure_water_molecules = number_of_molecules_finder_from_volume(box_refined_cell_length, 'water')
        if (int(float(num_water_in_glass))) <= (number_of_pure_water_molecules/2) - 2:
            out21_1.write(glass_name_number + '\t%s' %(int(float(num_water_in_glass))) + '\t%s' %(number_of_pure_water_molecules) + '= yes problem\n')
        else:
            out21_1.write(glass_name_number + '\t%s' %(int(float(num_water_in_glass))) + '\t%s' %(number_of_pure_water_molecules) + '= no problem\n')


    #4 Relax glass added with water  ! Addition of field file has to be taken care here
#    elif type_of_run == 'relax_glass_water':
        find_pattern_and_replace('FIELD', '#number_of_water_in_glass', 'NUMMOL %s' %(int(float(num_water_in_glass))), '%s/4_equilibrate_sol_glass_2' %(current_glass_directory), 'FIELD')
        os.system('cd %s/3_add_glass_with_water && cp CONFIG_GLASS_WATER_REFINED ../4_equilibrate_sol_glass_2' %(current_glass_directory))
        os.system('cd %s/4_equilibrate_sol_glass_2 && cp CONFIG_GLASS_WATER_REFINED CONFIG' %(current_glass_directory))
        os.system('cd %s/4_equilibrate_sol_glass_2 && chmod +x glasscript_P1.sh &&  sbatch glasscript_P1.sh' %(current_glass_directory))

    #5_pmf_calculations
    elif type_of_run == 'run_pmf':
        number_of_water_in_glass_pmf = get_num_water_in_glass(glass_name_number, '%s/glass_water_mixing_problems/glass_water_mixing_problem' %(directory))
        os.system('cp -r 5_pmf %s' %(current_glass_directory))
        find_pattern_and_replace('FIELD', '#number_of_water_in_glass', 'NUMMOL %s' %(int(float(number_of_water_in_glass_pmf))-1), '%s/5_pmf' %(current_glass_directory), 'FIELD')
        os.system('cd %s/4_equilibrate_sol_glass_2 && cp REVCON_V4 ../5_pmf/CONFIG' %(current_glass_directory))
        os.system('cd %s/5_pmf && python 1_1_config_replicator_only_full_water_acc_to_pmf_pbc.py' %(current_glass_directory))
        os.system('cd %s/5_pmf && python pmf_1_2_dist_info_processor.py' %(current_glass_directory))
        os.system('cd %s/5_pmf && python field_replicator_according_to_pmf.py' %(current_glass_directory))
        os.system('cd %s/5_pmf && python pmf_1_control_file.py' %(current_glass_directory))

    elif type_of_run == 'aluminum_run_pmf':
        number_of_water_in_glass_pmf = get_num_water_in_glass(glass_name_number, '%s/glass_water_mixing_problems/glass_water_mixing_problem' %(directory))
        os.system('cp -r 5_pmf %s' % (current_glass_directory))
        os.system('cd %s/4_equilibrate_sol_glass_2 && cp REVCON_V4 ../5_pmf/al/CONFIG' % (current_glass_directory))
        os.system('cd %s/4_equilibrate_sol_glass_2 && cp REVCON_V4 ../5_pmf/si/CONFIG' % (current_glass_directory))

        #       find_pattern_and_replace('FIELD', '#number_of_water_in_glass', 'NUMMOL %s' %(int(float(number_of_water_in_glass_pmf))-1), '%s/5_pmf/si' %(current_glass_directory), 'FIELD')
        #        find_pattern_and_replace('FIELD', '#number_of_water_in_glass', 'NUMMOL %s' %(int(float(number_of_water_in_glass_pmf))-1), '%s/5_pmf/al' %(current_glass_directory), 'FIELD')

        os.system('cd %s/5_pmf/al && python 1_1_config_replicator_only_full_water_acc_to_pmf_pbc.py' % (current_glass_directory))
        os.system('cd %s/5_pmf/si && python 1_1_config_replicator_only_full_water_acc_to_pmf_pbc.py' % (current_glass_directory))

        os.system('cd %s/5_pmf/al && python pmf_1_2_dist_info_processor.py' % (current_glass_directory))
        os.system('cd %s/5_pmf/si && python pmf_1_2_dist_info_processor.py' % (current_glass_directory))

        os.system('cd %s/5_pmf/al && python field_replicator_according_to_pmf.py' % (current_glass_directory))
        os.system('cd %s/5_pmf/si && python field_replicator_according_to_pmf.py' % (current_glass_directory))

        os.system('cd %s/5_pmf/al && python pmf_1_control_file.py' % (current_glass_directory))
        os.system('cd %s/5_pmf/si && python pmf_1_control_file.py' % (current_glass_directory))

    elif type_of_run == 'neural_network_input_preparation':
        print(current_glass_directory)
        # # Si in aluminosilicate glass
        # os.system('cd 6_analysis/si && cp neural_network_data_gen_si_in_alumino.py FIELD TABLE LOCAL_STRESS_KD_0620_V1.f %s/5_pmf/si' %(current_glass_directory))
        # os.system('cd %s/5_pmf/si && python3 pmf_5_results_to_tar.py' %(current_glass_directory))
        # os.system('cd %s/5_pmf/si && python3 neural_network_data_gen_si_in_alumino.py' %(current_glass_directory))

        # # Al in aluminosilicate glass
        # os.system('cd 6_analysis/al && cp neural_network_data_gen_al_in_alumino.py FIELD TABLE LOCAL_STRESS_KD_0620_V1.f %s/5_pmf/al' %(current_glass_directory))
        # os.system('cd %s/5_pmf/al && python3 pmf_5_results_to_tar.py' %(current_glass_directory))
        # os.system('cd %s/5_pmf/al && python3 neural_network_data_gen_al_in_alumino.py' %(current_glass_directory))

        # # Si in pure silicate glass
        os.system('cd 6_analysis/bond_dissociation && cp pmf_4_output_analyzer_id.py pmf_5_results_to_tar.py neural_network_data_gen.py FIELD TABLE LOCAL_STRESS_KD_0620_V1.f %s/5_pmf' % (current_glass_directory))
        os.system('cd %s/5_pmf && python3 pmf_5_results_to_tar.py' % (current_glass_directory))
        os.system('cd %s/5_pmf && python3 neural_network_data_gen.py' % (current_glass_directory))


    elif type_of_run == 'neural_network_input_compiled':
        os.system('cp %s/4_equilibrate_sol_glass_2/REVCON_V4 6_analysis/glass_structures/REVCON_V4_%s' %(current_glass_directory, lines0_1))

        # # Si in aluminosilicate glass
        # out_3_1 = open('6_analysis/si/neural_network_data_compiled', 'w')
        # xnnic_1 = open('%s/5_pmf/si/neural_input_preparation_table' %(current_glass_directory)).readlines()
        # for linesnnic1 in xnnic_1:
        #     if linesnnic1 not in list_nnic:
        #         list_nnic.append(linesnnic1)
        # for linesnnic2 in list_nnic:
        #     out_3_1.write(linesnnic2)

    # # Si in aluminosilicate glass
    #     out_3_1 = open('6_analysis/al/neural_network_data_compiled', 'w')
    #     xnnic_1 = open('%s/5_pmf/al/neural_input_preparation_table' %(current_glass_directory)).readlines()
    #     for linesnnic1 in xnnic_1:
    #         if linesnnic1 not in list_nnic:
    #             list_nnic.append(linesnnic1)
    #     for linesnnic2 in list_nnic:
    #         out_3_1.write(linesnnic2)

    # # Si in puresilicate glass
        out_3_1 = open('6_analysis/bond_dissociation/neural_network_data_compiled', 'w')
        xnnic_1 = open('%s/5_pmf/neural_input_preparation_table' % (current_glass_directory)).readlines()
        for linesnnic1 in xnnic_1:
            if linesnnic1 not in list_nnic:
                list_nnic.append(linesnnic1)
        for linesnnic2 in list_nnic:
            out_3_1.write(linesnnic2)

