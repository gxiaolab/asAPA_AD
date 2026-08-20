#!/usr/bin/python
import sys
import re
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats
##########

## python chisquare_test.py rbp_file_with_counts
## Performs chisquare test on rbp and outputs p-value

##########

infile = sys.argv[1]
out_path = sys.argv[2]

with open(infile) as rbp_output:
    all_lines = rbp_output.readlines()
    for i in range(len(all_lines)):
        line = all_lines[i]
        peak = line.strip()
        info = peak.split()
        rbp = info[0]
        file_name = f"{rbp}.jpg"
        file_path = os.path.join(out_path, file_name)
        test_count = int(info[1])
        ctrl_count = int(info[2])
        total_test = int(info[3])
        total_control = int(info[4])
        true_prop_test = test_count / total_test
        true_prop_ctrl = ctrl_count / total_control
        prop_diff = true_prop_test - true_prop_ctrl
        prop_test = (test_count + 1) / (total_test + 1)
        prop_ctrl = (ctrl_count + 1) / (total_control + 1)
        ratio = [prop_test, prop_ctrl]
        fold = np.log2(ratio[0]/ratio[1])
        ctrl_model = ["H"]*ctrl_count + ["N"]*(total_control - ctrl_count) ## H stands for hit. N stands for Not
        Sim = np.array(ctrl_model)
        results=[]
        for j in range(10000):
            p_ctrlmodel = np.random.choice(Sim,total_test)
            result = np.sum(p_ctrlmodel=="H")
            results.append(result/total_test)
        results = np.array(results)
        stdev = np.std(results)
        mean = np.mean(results)
        zscore = (true_prop_test - mean) / stdev
        pval = scipy.stats.norm.sf(abs(zscore))
        ##if true_prop_test >= true_prop_ctrl:
        ##    pval = np.sum(results >= true_prop_test) / 1000
        ##else:
        ##    pval = np.sum(results <= true_prop_test) / 1000
        ##if pval == 0:
        ##    pval = 0.001
        print(f"{rbp}\t{test_count}\t{total_test}\t{true_prop_test}\t{ctrl_count}\t{total_control}\t{true_prop_ctrl}\t{prop_diff}\t{fold}\t{pval}")
        #if i <= 3:
        #    p = sns.displot(data=results)
        #    p.set(xlabel="Proportion of RBP Hits", ylabel="Number of Simulations")
        #    plt.savefig(file_path, dpi=300)
