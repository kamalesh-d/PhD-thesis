
import matplotlib.pyplot as plt
import glob
import re
from scipy.stats import linregress
import matplotlib as mpl
import numpy as np



mpl.rc('figure', max_open_warning = 0)

consecutive_number_to_verify_slope = 7
threshold_of_slope = 0.05


def calculate_b_over_si_ratio(filename, timescale_observed):
    x2_1 = open(filename).readlines()
    num_si_in_sol = {}
    col7 = {}
    col12 = {}
    for lines2_1 in x2_1:
        r2 = re.search(r'^ nredept', lines2_1)
        if r2:
            splitted = lines2_1.strip().split()
            num_si_in_sol[int(splitted[1])] = int(splitted[2])
            col7[int(splitted[1])] = int(splitted[6])
            col12[int(splitted[1])] = int(splitted[11])

    fig, axs = plt.subplots(2)
    fig.suptitle('MC si and B')
    time_step = list(col7.keys())
    time_wise_num_of_si = list(col7.values())
    time_wise_num_of_b = list(col12.values())
    total_number_of_steps = len(time_step)

    starting_position_predicted = time_step.index(timescale_observed)

    axs[0].plot(col7.keys(), col7.values(), label="Si")
    axs[0].legend()
    axs[0].set_ylabel('Number of Si in the solution')
    axs[1].plot(col12.keys(), col12.values(), label="B")
    axs[1].legend()
    axs[1].set_xlabel('Monte Carlo Steps')
    axs[1].set_ylabel('Number of B in the solution')
    axs[0].axvline(list(col7.keys())[starting_position_predicted])
    axs[1].axvline(list(col7.keys())[starting_position_predicted])
    num_of_si_at_start = list(col7.values())[starting_position_predicted]
    num_of_b_at_start = list(col12.values())[starting_position_predicted]

    with open("analysis.txt", "r") as f:
        x3_1 = f.readlines()
    with open("analysis.txt", "w") as f:
        for lines3_1 in x3_1:
            r1 = re.search('^%s' %(filename), lines3_1)
            if r1:
                f.write(filename + '\t' + str(num_of_si_at_start) + '\t' + str(num_of_b_at_start) + '\t' + str(num_of_b_at_start/num_of_si_at_start) + '\n')
            else:
                f.write(lines3_1)

    plt.savefig('image_%s.png' %(filename))


def find_si_saturation_point(timestep, num_of_si_released):
    max_num_si = num_of_si_released[0]
    for i in range(0, len(num_of_si_released)):
        if num_of_si_released[i] > max_num_si:
            max_num_si = num_of_si_released[i]
    max_timestep = timestep[num_of_si_released.index(max_num_si)]
    return {'tsat': max_timestep, 'max_num_of_si': max_num_si, 'tsat_index': num_of_si_released.index(max_num_si)}


# calculate_b_over_si_ratio('output_78', 200000)

# Options = 'mc_data_analysis' and 'mc_exp_comparison', 'mc_and_mc_exp_comparison'
# if type_of_run == 'mc_and_mc_exp_comparison':

type_of_run = 'mc_exp_comparison'

if type_of_run == 'mc_data_analysis':
    x1_1 = glob.glob('output_*')
    out_1 = open('analysis.txt', 'w')
    out_1.write('filename' + '\t' + 'Number of Si predicted in saturated place' + '\t' + 'Number of B predicted in saturated place' + '\t' + 'B/Si ratio' + '\t' + 'si_saturation_timestep' + '\n')
    error = open('redo.txt', 'w')

    for lines1_1 in x1_1:
        x2_1 = open(lines1_1).readlines()
        num_si_in_sol = {}
        col3 = {}
        col7 = {}
        col8 = {}
        col9 = {}
        col10 = {}
        col11 = {}
        col12 = {}
        for lines2_1 in x2_1:
            r2 = re.search(r'^ nredept', lines2_1)
            if r2:
                splitted = lines2_1.strip().split()
                num_si_in_sol[int(splitted[1])] = int(splitted[2])
                col3[int(splitted[1])] = int(splitted[2])
                col7[int(splitted[1])] = int(splitted[6])
                col8[int(splitted[1])] = int(splitted[7])
                col9[int(splitted[1])] = int(splitted[8])
                col10[int(splitted[1])] = int(splitted[9])
                col11[int(splitted[1])] = int(splitted[10])
                col12[int(splitted[1])] = int(splitted[11])

        fig, axs = plt.subplots(2)
        fig.suptitle('MC si and B')
        time_step = list(col7.keys())
        time_wise_num_of_si = list(col7.values())
        time_wise_num_of_b = list(col12.values())
        total_number_of_steps = len(time_step)
        all_slope = []

        for lines6_1 in range(0, total_number_of_steps-7):
            x1 = time_step[lines6_1]
            x2 = time_step[lines6_1+5]
            y1 = time_wise_num_of_si[lines6_1]
            y2 = time_wise_num_of_si[lines6_1+5]
        saturation_data = find_si_saturation_point(time_step, time_wise_num_of_si)

        starting_position_predicted = saturation_data.get('tsat_index')
        timestep_max_si_released = saturation_data.get('tsat')

        axs[0].plot(col7.keys(), col7.values(), label="Si")
        axs[0].legend()
        axs[0].set_ylabel('Number of Si in the solution')
        axs[1].plot(col12.keys(), col12.values(), label="B")
        axs[1].legend()
        axs[1].set_xlabel('Monte Carlo Steps')
        axs[1].set_ylabel('Number of B in the solution')
        axs[0].axvline(timestep_max_si_released)
        axs[1].axvline(timestep_max_si_released)
        num_of_si_at_start = list(col7.values())[starting_position_predicted]
        num_of_b_at_start = list(col12.values())[starting_position_predicted]

        out_1.write(lines1_1 + '\t' + str(num_of_si_at_start) + '\t' + str(num_of_b_at_start) + '\t' + str(num_of_b_at_start/num_of_si_at_start) + '\t' + str(timestep_max_si_released) + '\n')
        plt.savefig('image_%s.png' %(lines1_1))


def list_filename_tsat(filename):
    x3_0 = open(filename).readlines()
    dict_file_tsat = {}
    for lines3_0 in x3_0[1:]:
        splitted3_0 = lines3_0.strip().split('\t')
        filename = (splitted3_0[0])
        si_sat_timestep = splitted3_0[-1]
        dict_file_tsat[filename] = si_sat_timestep
    return dict_file_tsat

dict_si_saturation_timestep = list_filename_tsat('analysis.txt')

"""
Input parameters to be given starts here
"""

si_saturation_experimental = 190.91*3600 # In seconds
name_of_the_glass = 'SBNA2'
glass_size = 500 * 3.5 # Each atoms placed minimum of distance 3.5 angstrom
#Glass composition
number_of_si = 755108
number_of_b = 239194 + 120048
number_of_al = 135650
total_atoms = number_of_si + number_of_b + number_of_al

"""
Input parameters to be given ends here
"""

if type_of_run == 'mc_exp_comparison':
    x10_1 = glob.glob('output_83')
    for lines10_1 in x10_1:
        si_saturation_timestep = int(dict_si_saturation_timestep.get(lines10_1))
        montecarlo_time_in_seconds = si_saturation_experimental / si_saturation_timestep

        filename_splitted = lines10_1.strip().split('_')
        x10_2 = open(lines10_1).readlines()
        all_si_in_sol = []
        all_timesteps = []
        all_b_in_sol = []

        fig, axs = plt.subplots(2)
        fig.suptitle('Fitting of monte-carlo with experimental data \n for %s' %(name_of_the_glass))

        for lines10_2 in x10_2:
            r2 = re.search(r'nredept', lines10_2)
            if r2:
                splitted10_1 = lines10_2.strip().split()
                all_timesteps.append(int(splitted10_1[1]))
                all_si_in_sol.append(int(splitted10_1[6]))
                all_b_in_sol.append(int(splitted10_1[11]))

        # Comparison of experimental vs MC simulation results
        temps_s = (np.array(all_timesteps) * montecarlo_time_in_seconds)
        temps_h = temps_s/3600
        mc_racine_temps = np.sqrt(temps_h)
        mc_equiv_thick_b = (np.array(all_b_in_sol)/(number_of_b)) * glass_size * 0.1
        mc_equiv_thick_si = (np.array(all_si_in_sol)/(number_of_si)) * glass_size * 0.1
        #
        #
        # # Experimental data visualization
        exp_racine_temps = []
        exp_equiv_thick_b = []
        exp_equiv_thick_si = []
        x12_1 = open('experimental_data.txt').readlines()
        for lines12_1 in x12_1[1:]:
            splitted12_1 = lines12_1.strip().split()
            exp_racine_temps.append(float(splitted12_1[0]))
            exp_equiv_thick_b.append(float(splitted12_1[1]))
            exp_equiv_thick_si.append(float(splitted12_1[2]))

        axs[0].plot(exp_racine_temps, exp_equiv_thick_si, label='Experimental Si in solution')
        axs[0].plot(mc_racine_temps, mc_equiv_thick_si, label='Monte-Carlo Si in solution')
        axs[1].plot(exp_racine_temps, exp_equiv_thick_b, label='Experimental B in solution')
        axs[1].plot(mc_racine_temps, mc_equiv_thick_b, label='Monte-Carlo B in solution')
        axs[0].set_ylabel('Equivalent thickness \n of Si (nm)')
        axs[1].set_ylabel('Equivalent thickness \n of B (nm)')
        axs[1].set_xlabel('Square root of time (Hour)')
        axs[0].legend()
        axs[1].legend()
        plt.savefig('exp_vs_monte_carlo_%s.png' %(filename_splitted[1]), dpi=300)
        # plt.show()

