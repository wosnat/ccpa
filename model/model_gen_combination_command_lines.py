#!/usr/bin/env python
# coding: utf-8


import model.basic_model_lib as md
from itertools import combinations
from pprint import pprint

pro_param_list = [
    'gamma_n_p',
    #'gamma_c_p',
    #'gamma_refractory_n_p',
    #'gamma_refractory_c_p',
    'excretion_n_p',
    #'excretion_c_p',
    'v_n_max_p',
    #'v_c_max_p',
    'mu_inf_p',
    'k_n_p',
    #'k_c_p',
    'q_n_min_p',
    #'q_c_min_p',
    'q_n_max_p',
    #'q_c_max_p',
    'mortality_p',
    #'r0_p',
    #'b_p'
]

alt_param_list = [
    'gamma_n_a',
    'gamma_c_a',
    #'gamma_refractory_n_a',
    #'gamma_refractory_c_a',
    'excretion_n_a',
    'excretion_c_a',
    'v_in_max_a',
    'v_n_max_a',
    'v_c_max_a',
    'mu_inf_a',
    'k_n_a',
    'k_in_a',
    'k_c_a',
    'q_n_min_a',
    'q_c_min_a',
    'q_n_max_a',
    'q_c_max_a',
    'mortality_a',
    #'r0_a',
    #'b_a'
]

print(len(alt_param_list))
# m = md.ModelProALT()
# params = m.parameters.keys()
# #
# # pro_params = [p for p in params if p.endswith('_p')]
# #
# alt_params = [p for p in params if p.endswith('_a')]
# # pprint(pro_params)
# # print(len(pro_params))
# pprint(alt_params)


# pro_sensitivity
def print_pro_slist():
    param_combs = list()
    for i in range (1,len(pro_param_list)+1):
    #for i in range(1, 2):
        print(i, len(param_combs))
        param_combs.extend(combinations(pro_param_list, i))
        #pprint(param_combs)
    sbatch_str = 'sbatch --cpus-per-task 10 -p hive7d,preempt7d --wrap '
    cmd_str = '/data/home/snir/oweissber/model_work/ccpa/model/basic_model_lib.py  --optimize --params_opt {params_list} --ref_csv /data/home/snir/oweissber/model_work/ccpa/model/med4_1a3_ref.csv --disable_alt --ref_pro_col PRO --outfile /data/home/snir/oweissber/model_work/run_optimize/pro_sensitivity/{model_name}.csv'

    with open('sensitivity_pro_list', 'w') as fp:
        for i in param_combs:
            model_name = '_'.join(i)
            params_list = ' '.join(i)
            print (sbatch_str + " '" +  cmd_str.format(params_list=params_list, model_name=model_name) + "'", file=fp)


def print_alt_slist():
    param_combs = list()
    alt_param_list = [i.replace("_p", "_a") for i in pro_param_list]
    print(alt_param_list)
    for i in range(1, len(alt_param_list) + 1):
        # for i in range(1, 2):
        print(i, len(param_combs))
        param_combs.extend(combinations(alt_param_list, i))
        # pprint(param_combs)
    sbatch_str = 'sbatch --cpus-per-task 10 -p hive7d,preempt7d --wrap '
    cmd_str = '/data/home/snir/oweissber/model_work/ccpa/model/basic_model_lib.py  --optimize --params_opt {params_list} --ref_csv /data/home/snir/oweissber/model_work/ccpa/model/med4_1a3_ref.csv --disable_pro --ref_alt_col ALT --outfile /data/home/snir/oweissber/model_work/run_optimize/alt_sensitivity/{model_name}.csv'

    with open('sensitivity_alt_list', 'w') as fp:
        for i in param_combs:
            model_name = '_'.join(i)
            params_list = ' '.join(i)
            print(sbatch_str + " '" + cmd_str.format(params_list=params_list, model_name=model_name) + "'", file=fp)


#print_alt_slist()
