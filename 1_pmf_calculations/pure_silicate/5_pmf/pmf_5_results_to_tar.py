import os



x1_1 = open('all_distance_information.txt').readlines()
os.system('mkdir all_results')
for lines1_1 in x1_1:
    x1_2 = lines1_1.strip().split('\t')
    assembled = x1_2[0] + '_' + x1_2[1] + '_' + x1_2[2] + '_' + x1_2[3]
    os.system('cp pmf_4_output_analyzer_id.py %s' %(assembled))
    os.system('cd %s && python pmf_4_output_analyzer_id.py' %(assembled))
    os.system('cd %s && cp result.txt result_%s.txt' %(assembled, assembled))
    os.system('cd %s && cp result_%s.txt ..' %(assembled, assembled))
    os.system('mv result_%s.txt all_results' %(assembled))
    os.system('tar -czf all_results.tar.gz all_results')
