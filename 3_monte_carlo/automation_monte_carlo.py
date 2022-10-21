import re
import os
import glob
# import matplotlib.pyplot as plt
import statistics

os.system('module load python/3.6.10')

# Notes:
# If you are running this script, the initial dissolution rate of Si is calculated based on num of si at 500 and si saturation level is calculated based on steps between 20,000 to 30,000 time steps.
position_for_initial_dissolution = 500
print('Note: initial dissolution is calculated based on num of si at 500 and si saturated is calculated by averaging 20,000 to 30,000 timesteps')

#
glass_composition = {'SiO2':65, 'B2O3':17.3, 'Na2O':13.6, 'Al2O3':4.1}
atomic_mass = {'Si': 28.085, 'O':15.999, 'B':10.81, 'Na':22.989, 'Al':26.981}
avogadro_number = 6.023 * 10 ** 23
# Place the size of the monte carlo glass in x and y directions.
num_monte_glass_nodes_x = 50
num_monte_glass_nodes_y = 50
monte_nodes_distance = 3.5
monte_carlo_file_name = 'monte_carlo_3D_V14b_SBNA2_P25.f'
print('Note: Check if the monte carlo glass has only 50 nodes on x-axis and 50 nodes on y-axis and each nodes are seperated by 3.5 ang distance')
surf_over_vol = 50
print('Note: Check the surface area by volume ratio used for the different glass during the experiments. Look to keep this dynamic for each 6 glasses')

consecutive_number_to_verify_slope = 7
threshold_of_slope = 0.05


# wbreak1 = [100, 200, 300, 400, 500, 600, 700, 800, 900]



def find_pattern_and_replace(input_file, pattern_in_target_line, replace_line, file_directory_location, output_file_name):
    xd2_1 = open('%s/%s' % (file_directory_location, input_file)).readlines()
    outd2_1 = open('%s/%s' % (file_directory_location, output_file_name), 'w')
    replace_line_pos_1 = []
    for linesd2_1 in xd2_1:
        rd2 = re.search(r'%s' % (pattern_in_target_line), linesd2_1)
        if rd2:
            replace_line_pos_1.append(str(xd2_1.index(linesd2_1)))
    if replace_line_pos_1 != []:
        xd2_1[int(''.join(replace_line_pos_1))] = replace_line + '\n'
        for linesd2_2 in xd2_1:
            outd2_1.write(linesd2_2)

def grab_from_output(list_of_contents, parameter_to_fetch):
    for linesd5_1 in list_of_contents:
        r3 = re.search(r'^%s' %(parameter_to_fetch), linesd5_1)
        if r3:
            return linesd5_1.strip().split(':')[1]



p_1 = glob.glob('trial_*')
directories_num = []
for linesd1_1 in p_1:
    directories_num.append(int(linesd1_1.split('_')[1]))

create_directory_from = sorted(directories_num)[-1] + 1


#create_directory_from = 0

type_of_run = 'restart_calculation'
if type_of_run == 'wbreak_pilot_study':
    y1_1 = open('wbreak.txt').readlines()
    for lines1_0, lines1_1 in zip(range(create_directory_from, create_directory_from+len(y1_1)), y1_1):
        os.system('mkdir trial_%s' %(lines1_0))
        os.system('cp script_V11_p07 %s trial_%s' %(monte_carlo_file_name, lines1_0))
        current_trial_directory = 'trial_%s' % (lines1_0)

        wbreak1 = '%s.d0' %(lines1_1.strip())
        wbreak2 = '%s.d0' %(lines1_1.strip())
        wbreak3 = '%s.d0' %(lines1_1.strip())

        find_pattern_and_replace('%s' %(monte_carlo_file_name), 'wbreak1 =', '      wbreak1 =  %s ! Fréquence d ouverture des liaisons Si-Si sur 1000 (333 = 1 chance sur 3)' % (wbreak1), current_trial_directory, '%s' %(monte_carlo_file_name))
        find_pattern_and_replace('%s' %(monte_carlo_file_name), 'wbreak2 =', '      wbreak2 =  %s ! Fréquence d ouverture des liaisons Si-Si sur 1000 (333 = 1 chance sur 3)' % (wbreak2), current_trial_directory, '%s' %(monte_carlo_file_name))
        find_pattern_and_replace('%s' %(monte_carlo_file_name), 'wbreak3 =', '      wbreak3 =  %s ! Fréquence d ouverture des liaisons Si-Si sur 1000 (333 = 1 chance sur 3)' % (wbreak3), current_trial_directory, '%s' %(monte_carlo_file_name))

if type_of_run == 'wred_pilot_study':
    y1_2 = open('wred.txt').readlines()
    for lines1_2, lines1_3 in zip(range(create_directory_from, create_directory_from+len(y1_2)), y1_2):
        os.system('mkdir trial_%s' % (lines1_2))
        os.system('cp script_V11_p07 %s trial_%s' % (monte_carlo_file_name, lines1_2))
        current_trial_directory = 'trial_%s' % (lines1_2)

        # # Wred changes will be made here
        wred = '%s.d0' %(lines1_3.strip())
        wredal = '%s.d0' %(lines1_3.strip())

        find_pattern_and_replace('%s' %(monte_carlo_file_name), 'wred =', '      wred =  %s/1.d0! Fréquence d ouverture des liaisons Si-Si sur 1000 (333 = 1 chance sur 3)' %(wred), current_trial_directory, '%s' %(monte_carlo_file_name))
        find_pattern_and_replace('%s' %(monte_carlo_file_name), 'wredal =', '      wredal =  %s/1.d0! Fréquence d ouverture des liaisons Si-Si sur 1000 (333 = 1 chance sur 3)' %(wredal), current_trial_directory, '%s' %(monte_carlo_file_name))
        os.system('cd trial_%s && sbatch -n 1 ./script_V11_p07' %(lines1_2))

if type_of_run == 'ncvoisw_pilot_study':
    y1_3 = open('nc_voisw.txt').readlines()
    for lines1_4, lines1_5 in zip(range(create_directory_from, create_directory_from+len(y1_3)), y1_3):
        os.system('mkdir trial_%s' % (lines1_4))
        os.system('cp script_V11_p07 %s trial_%s' % (monte_carlo_file_name, lines1_4))
        current_trial_directory = 'trial_%s' % (lines1_4)

        # # # nc_voisw
        nc_voisw = str(lines1_4)
        find_pattern_and_replace('%s' %(monte_carlo_file_name), '      nc_voisw', '      nc_voisw = %s ! Fréquence d ouverture des liaisons Si-Si sur 1000 (333 = 1 chance sur 3)' %(nc_voisw), current_trial_directory, '%s' %(monte_carlo_file_name))
        os.system('cd trial_%s && sbatch -n 1 ./script_V11_p07' %(lines1_5))


if type_of_run == 'wbreak_wred_ncvoisw_pilot_study':
    wbreak_1 = []
    wred_1 = []
    nc_voisw_1 = []
    y1_1 = open('wbreak1.txt').readlines()
    y1_2 = open('wred.txt').readlines()
    y1_3 = open('nc_voisw.txt').readlines()

    for lines3_1 in y1_1:
        for lines3_2 in y1_2:
            for lines3_3 in y1_3:

                os.system('mkdir trial_%s' %(create_directory_from))
                os.system('cp script_V11_p07 %s trial_%s' %(monte_carlo_file_name, create_directory_from))
                current_trial_directory = 'trial_%s' % (create_directory_from)
                create_directory_from = create_directory_from + 1

                wbreak1 = '%s.d0' % (lines3_1.strip())
                wbreak2 = '%s.d0' % (lines3_1.strip())
                wbreak3 = '%s.d0' % (lines3_1.strip())

                # # Wred changes will be made here
                wred = '%sd0' % (lines3_2.strip())
                wredal = '%sd0' % (lines3_2.strip())

                # # # nc_voisw
                nc_voisw = str(lines3_3)
                # wbreak1
                find_pattern_and_replace('%s' %(monte_carlo_file_name), 'wbreak1 =', '      wbreak1 =  %s ! Fréquence d ouverture des liaisons Si-Si sur 1000 (333 = 1 chance sur 3)' % (wbreak1), current_trial_directory, '%s' %(monte_carlo_file_name))
                find_pattern_and_replace('%s' %(monte_carlo_file_name), 'wbreak2 =', '      wbreak2 =  %s ! Fréquence d ouverture des liaisons Si-Si sur 1000 (333 = 1 chance sur 3)' % (wbreak2), current_trial_directory, '%s' %(monte_carlo_file_name))
                find_pattern_and_replace('%s' %(monte_carlo_file_name), 'wbreak3 =', '      wbreak3 =  %s ! Fréquence d ouverture des liaisons Si-Si sur 1000 (333 = 1 chance sur 3)' % (wbreak3), current_trial_directory, '%s' %(monte_carlo_file_name))
                # wred
                find_pattern_and_replace('%s' %(monte_carlo_file_name), 'wred =', '      wred =  %s/1.d0! Fréquence d ouverture des liaisons Si-Si sur 1000 (333 = 1 chance sur 3)' %(wred), current_trial_directory, '%s' %(monte_carlo_file_name))
                find_pattern_and_replace('%s' %(monte_carlo_file_name), 'wredal =', '      wredal =  %s/1.d0! Fréquence d ouverture des liaisons Si-Si sur 1000 (333 = 1 chance sur 3)' %(wredal), current_trial_directory, '%s' %(monte_carlo_file_name))
                # nc_voisw
                find_pattern_and_replace('%s' %(monte_carlo_file_name), '      nc_voisw', '      nc_voisw = %s ! Fréquence d ouverture des liaisons Si-Si sur 1000 (333 = 1 chance sur 3)' %(nc_voisw), current_trial_directory, '%s' %(monte_carlo_file_name))
                os.system('cd %s && sbatch -n 1 ./script_V11_p07' %(current_trial_directory))


if type_of_run == 'nsit_pilot_study':
    # You have to change only nsit, nsits will be automatically added with 50.
    nsit = []
    y1_5 = open('nsit.txt').readlines()

    for lines3_5 in y1_5:
        os.system('mkdir trial_%s' %(create_directory_from))
        os.system('cp script_V11_p07 %s trial_%s' %(monte_carlo_file_name, create_directory_from))
        current_trial_directory = 'trial_%s' % (create_directory_from)
        create_directory_from = create_directory_from + 1

        # nsit changed
        find_pattern_and_replace('%s' %(monte_carlo_file_name), '      parameter \(nsit=', '      parameter (nsit=%s) !Total size' % (lines3_5.strip()), current_trial_directory, '%s' %(monte_carlo_file_name))
        find_pattern_and_replace('%s' %(monte_carlo_file_name), '      parameter \(nsits=', '      parameter (nsits=%s) !Size of the glass' % (int(lines3_5.strip())-50), current_trial_directory, '%s' %(monte_carlo_file_name))
        os.system('cd %s && sbatch -n 1 ./script_V11_p07' %(current_trial_directory))


if type_of_run == 'wbreak_wred_wsaut_pilot_study':
    wbreak_1 = []
    wred_1 = []
    nc_voisw_1 = []
    y1_1 = open('wbreak1.txt').readlines()
    y1_2 = open('wred.txt').readlines()
    y1_3 = open('wsaut.txt').readlines()

    for lines3_1 in y1_1:
        for lines3_2 in y1_2:
            for lines3_3 in y1_3:

                os.system('mkdir trial_%s' %(create_directory_from))
                os.system('cp script_V11_p07 %s trial_%s' %(monte_carlo_file_name, create_directory_from))
                current_trial_directory = 'trial_%s' % (create_directory_from)
                create_directory_from = create_directory_from + 1

                wbreak1 = '%s.d0' % (lines3_1.strip())
                wbreak2 = '%s.d0' % (lines3_1.strip())
                wbreak3 = '%s.d0' % (lines3_1.strip())

                # # Wred changes will be made here
                wred = '%sd0' % (lines3_2.strip())
                wredal = '%sd0' % (lines3_2.strip())

                # # # wsaut
                wsaut = str(lines3_3)
                # wbreak1
                find_pattern_and_replace('%s' %(monte_carlo_file_name), 'wbreak1 =', '      wbreak1 =  %s ! Fréquence d ouverture des liaisons Si-Si sur 1000 (333 = 1 chance sur 3)' % (wbreak1), current_trial_directory, '%s' %(monte_carlo_file_name))
                find_pattern_and_replace('%s' %(monte_carlo_file_name), 'wbreak2 =', '      wbreak2 =  %s ! Fréquence d ouverture des liaisons Si-Si sur 1000 (333 = 1 chance sur 3)' % (wbreak2), current_trial_directory, '%s' %(monte_carlo_file_name))
                find_pattern_and_replace('%s' %(monte_carlo_file_name), 'wbreak3 =', '      wbreak3 =  %s ! Fréquence d ouverture des liaisons Si-Si sur 1000 (333 = 1 chance sur 3)' % (wbreak3), current_trial_directory, '%s' %(monte_carlo_file_name))
                # wred
                find_pattern_and_replace('%s' %(monte_carlo_file_name), 'wred =', '      wred =  %s/1.d0! Fréquence d ouverture des liaisons Si-Si sur 1000 (333 = 1 chance sur 3)' %(wred), current_trial_directory, '%s' %(monte_carlo_file_name))
                find_pattern_and_replace('%s' %(monte_carlo_file_name), 'wredal =', '      wredal =  %s/1.d0! Fréquence d ouverture des liaisons Si-Si sur 1000 (333 = 1 chance sur 3)' %(wredal), current_trial_directory, '%s' %(monte_carlo_file_name))
                # wsaut
                find_pattern_and_replace('%s' %(monte_carlo_file_name), '      wsaut =', '      wsaut = %s ! FrÃ©quence de sauts des W' %(wsaut.strip()), current_trial_directory, '%s' %(monte_carlo_file_name))
                find_pattern_and_replace('%s' %(monte_carlo_file_name), '      wsaut34 =', '      wsaut34 = %s ! FrÃ©quence de sauts des W' %(wsaut.strip()), current_trial_directory, '%s' %(monte_carlo_file_name))
                find_pattern_and_replace('%s' %(monte_carlo_file_name), '      wsaut56 =', '      wsaut56 = %s ! FrÃ©quence de sauts des W' %(wsaut.strip()), current_trial_directory, '%s' %(monte_carlo_file_name))


                os.system('cd %s && sbatch -n 1 ./script_V11_p07' %(current_trial_directory))



calculations_directory_completed = create_directory_from - 1
output_analyze_from_to = [210, 217]
if type_of_run == 'output_analysis':

    os.system('mkdir results')
    # out_1 = open('results/analysis.txt', 'w')
    # out_1.write('Number of Si predicted in saturated place' + '\t\t\t' + 'Number of B predicted in saturated place' + '\t\t\t' + 'B/Si ratio' + '\n')

    for lines6_21 in range(output_analyze_from_to[0], output_analyze_from_to[1]+1):
        os.system('cp trial_%s/output_V11_p07 results/output_%s' %(lines6_21, lines6_21))
        os.system('cp trial_%s/structure_den25 results/structure_den25_%s' %(lines6_21, lines6_21))
        os.system('cd results && grep nredept output_%s > temp_%s' %(lines6_21, lines6_21))
    os.system('tar -czf results.tar.gz results')


restart_calculation_from_to = [219, 219]
if type_of_run == 'restart_calculation':
    for lines7_21 in range(restart_calculation_from_to[0], restart_calculation_from_to[1]+1):
        current_trial_directory_to_restart = 'trial_%s' %(lines7_21)
        os.system("cd %s && mkdir old_calculations && cp * old_calculations" %(current_trial_directory_to_restart))
        find_pattern_and_replace('%s' % (monte_carlo_file_name), '      irepri=', '      irepri=1 !irepri=1 pour reprendre un calcul', current_trial_directory_to_restart, '%s' % (monte_carlo_file_name))
        os.system('cd %s && sbatch -n 1 ./script_V11_p07' % (current_trial_directory_to_restart))
        print("Restart completed for %s" %(current_trial_directory_to_restart))

