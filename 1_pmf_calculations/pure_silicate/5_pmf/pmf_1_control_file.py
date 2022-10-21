import os



x1_1 = open('all_distance_information.txt').readlines()
for lines1_1 in x1_1:
    x1_2 = lines1_1.strip().split('\t')
    assembled = x1_2[0] + '_' + x1_2[1] + '_' + x1_2[2] + '_' + x1_2[3]

    os.system('mkdir %s' %(assembled))

    out_1 = open('%s/dist_info.txt' %(assembled), 'w')
    out_1.write(lines1_1)
    os.system('mv CONFIG_%s %s' %(assembled, assembled))
    os.system('mv FIELD_%s %s' %(assembled, assembled))
    os.system('cp CONTROL %s' %(assembled))
    os.system('cp TABLE %s' %(assembled))
    os.system('cp pmf_1_master_glass.py %s' %(assembled))
    os.system('cp pmf_4_output_analyzer_id.py %s' %(assembled))
    os.system('cp run_python_script %s' %(assembled))
    ## os.system('cd %s && sbatch -n 1 python3 pmf_1_master_glass.py &' %(assembled))   # How do we submit the python script in sbatch format
    ## os.system('source /softs/_environnement/dl_poly-4.0.6.env')   # How do we submit the python script in sbatch format
    os.system('cd %s && sbatch run_python_script' %(assembled))   # How do we submit the python script in sbatch format

