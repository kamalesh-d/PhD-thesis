import os
import re




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


type_of_run = '5_bond_dissociation_analysis' #keywords - 1_npt_relaxation, 2_nvt_relaxation, 3_pmf_calculation
current_directory = os.getcwd()
for lines0_1 in [0, 2, 3, 4, 6, 8, 9]:
    glass_directory = '%s/glass_num_%s' %(current_directory, lines0_1)
    if type_of_run == '1_npt_relaxation':
        os.system('cd 1_npt_relaxation && rm REVCON_V4 && cp REVCON_V4_%s REVCON_V4' %(lines0_1))
        os.system('mkdir glass_num_%s' %(lines0_1))
        os.system('cp -r 1_npt_relaxation %s' %(glass_directory))
        os.system('cd %s/1_npt_relaxation && python converter.py' %(glass_directory))
        os.system('cd %s/1_npt_relaxation && chmod +x run_lammps_npt.sh && sbatch ./run_lammps_npt.sh' %(glass_directory))

    if type_of_run == '2_nvt_relaxation':
        os.system('cp -r 2_nvt_relaxation %s' %(glass_directory))
        os.system('cd %s/1_npt_relaxation && python lammps_history_analysis.py && cp input_structure_2.data ../2_nvt_relaxation' %(glass_directory))
        os.system('cd %s/2_nvt_relaxation && chmod +x run_lammps_nvt.sh && sbatch ./run_lammps_nvt.sh' %(glass_directory))

    if type_of_run == '3_pmf_calculation':
        os.system('cp -r 3_pmf_calculation %s' %(glass_directory))
        os.system('cd %s/2_nvt_relaxation && python lammps_history_analysis.py && cp input_structure_2.data ../3_pmf_calculation' %(glass_directory))
        os.system('cd %s/3_pmf_calculation && python water_close_to_si_search.py' %(glass_directory))
        x1_1 = open('%s/3_pmf_calculation/atom_pairs.txt' %(glass_directory)).readlines()
        for lines1_1 in x1_1:
            splitted_1 = lines1_1.strip().split('\t')
            file_name = splitted_1[0] + '_' + splitted_1[1] + '_' + splitted_1[2] + '_' + splitted_1[3]
            os.system('cd %s/3_pmf_calculation && mkdir %s' %(glass_directory, file_name))
            out_1 = open('%s/3_pmf_calculation/%s/dist_info.txt' %(glass_directory, file_name), 'w')
            out_1.write(lines1_1)
            os.system('cd %s/3_pmf_calculation && cp in.interface input_structure_2.data run_lammps_pmf.sh ffield_new %s' %(glass_directory, file_name))
            x1_3 = find_pattern_and_replace('in.interface', '#oxygen_belonging_to_water', 'group A1 id %s ' %(splitted_1[1]), '%s/3_pmf_calculation/%s' %(glass_directory, file_name), 'in.interface')
            x1_4 = find_pattern_and_replace('in.interface', '#silica_belonging_to_glass', 'group A2 id %s ' %(splitted_1[3]), '%s/3_pmf_calculation/%s' %(glass_directory, file_name), 'in.interface')
            os.system('cd %s/3_pmf_calculation/%s && chmod +x run_lammps_pmf.sh && sbatch ./run_lammps_pmf.sh' %(glass_directory, file_name))



    if type_of_run == '4_bond_reformation':

#        os.system('cp -r 4_bond_reformation %s' %(glass_directory))
 #       os.system('cd %s/2_nvt_relaxation && python lammps_history_analysis.py && cp input_structure_2.data ../4_bond_reformation' %(glass_directory))
  #      os.system('cd %s/4_bond_reformation && python search_for_close_silanol_groups.py' %(glass_directory))
        x1_1 = open('%s/4_bond_reformation/atom_pairs.txt' %(glass_directory)).readlines()
        for lines1_1 in x1_1:
            splitted_1 = lines1_1.strip().split('\t')
            file_name = splitted_1[0] + '_' + splitted_1[1] + '_' + splitted_1[2] + '_' + splitted_1[3]
            os.system('cd %s/4_bond_reformation && mkdir %s' %(glass_directory, file_name))
            out_1 = open('%s/4_bond_reformation/%s/dist_info.txt' %(glass_directory, file_name), 'w')
            out_1.write(lines1_1)
            os.system('cd %s/4_bond_reformation && cp in.interface input_structure_2.data run_lammps_pmf.sh ffield_new %s' %(glass_directory, file_name))
            x1_3 = find_pattern_and_replace('in.interface', '#oxygen_belonging_to_water', 'group A1 id %s ' %(splitted_1[1]), '%s/4_bond_reformation/%s' %(glass_directory, file_name), 'in.interface')
            x1_4 = find_pattern_and_replace('in.interface', '#silica_belonging_to_glass', 'group A2 id %s ' %(splitted_1[3]), '%s/4_bond_reformation/%s' %(glass_directory, file_name), 'in.interface')
            os.system('cd %s/4_bond_reformation/%s && chmod +x run_lammps_pmf.sh && sbatch ./run_lammps_pmf.sh' %(glass_directory, file_name))



    if type_of_run == '5_bond_dissociation_analysis':
        x3_1 = open('%s/3_pmf_calculation/atom_pairs.txt' %(glass_directory)).readlines()
        for lines3_1 in x3_1:
            splitted_3 = lines3_1.strip().split('\t')
            file_name_3 = splitted_3[0] + '_' + splitted_3[1] + '_' + splitted_3[2] + '_' + splitted_3[3]
            os.system('cp %s/3_pmf_calculation/%s/hn.data 5_results/bond_dissociation/hn_data/hn_%s.data' %(glass_directory, file_name_3, file_name_3))
            os.system('cp %s/3_pmf_calculation/%s/pmf_calculation.dump 5_results/bond_dissociation/pmf_calculation/pmf_calculation_%s.dump' %(glass_directory, file_name_3, file_name_3))





