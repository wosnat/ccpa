#!/usr/bin/env python
# coding: utf-8



import sys, os
import math
from sympy import *

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import metrics


# param - name, symbol, value
# variable: name, symbol, value, formula
# values over time
# comparison to reference

# init default values for all params & variables
# override params & variables
# generate simplified formulas (lambdify)
# variables: name, symbol, value, formula, function
# simulate
# - capture values
# -


class ModelParam:
    def __init__(self, name, symbolname, initval):
        self.name = name
        self.symbol = symbols(symbolname)
        self.value = initval

    def override_value(self, value):
        self.value = value


class ModelVariable:
    def __init__(self, name, symbolname, initval, organism = None):
        self.name = name
        self.initval = initval
        self.symbol = symbols(symbolname)
        self.lambda_symbol = symbols(symbolname.replace('^','_'))
        self.value = initval
        self.nextvalue = None
        self.formula = None
        self.lambdafun = None
        self.organism = organism

    def run_lambdify(self, param_values, variables_list, intermediate_variables, intermediate_evaluation_order):
        s = self.formula
        for n in reversed(intermediate_evaluation_order):
            v = intermediate_variables[n]
            #print(s, v.symbol, v.formula)
            s = s.subs({v.symbol: v.formula})
            #print('after:', s)

        subs_expr = s.subs(param_values)
        subs_expr2 = subs_expr.subs(variables_list)
        #print('before:', self.formula)
        #print('after', subs_expr2)
        var_list_for_lambda = [v[1] for v in variables_list]

        self.lambdafun = lambdify(var_list_for_lambda, subs_expr2)

    def evaluate(self, variables_values):
        self.nextvalue = self.lambdafun(*variables_values)

    def copy_from_nextvalue(self):
        self.value = self.nextvalue

    def reset_initial_val(self):
        self.value = self.initval

    def override_initial_value(self, value):
        self.value = value
        self.initval = value

class ModelProALT:

    def __init__(self):
        self._disable_organism = None
        self.reference_days = list()
        self.init_params()
        self.init_intermediate_variables()
        self.init_variables()
        self.init_symbols()
        self.init_formulas()
        self.init_var_list()


    def init_params(self):
        """ parameters and initial values  - from hard coded list """
        # constants
        # Redfield ratio
        redfield_C_to_N = 6.625
        # ???
        pro_radius = 0.3628;  # "MED4" = 9312
        alt_radius = 1;  # "MED4" = 9312
        pro_vol = (4 / 3) * math.pi * pro_radius ** 3;
        alt_vol = (4 / 3) * math.pi * alt_radius ** 3;
        pro_alt_vol_ratio = alt_vol / pro_vol
        


        parameters = [
            # (name, symbolname, initval)
            # gamma=fraction of heterotroph mort/resp to inorganic form
            ('gamma_n_p', 'gamma_N^P', 0.04),
            ('gamma_c_p', 'gamma_C^P', 0.04),
            ('gamma_n_a', 'gamma_N^A', 0.04),
            ('gamma_c_a', 'gamma_C^A', 0.04),

            # gamma refractory =fraction of heterotroph mort/resp to refractory inorganic form
            ('gamma_refractory_n_p', 'gamma_refractory_N^P', 0.1),
            ('gamma_refractory_c_p', 'gamma_refractory_C^P', 0.1),
            ('gamma_refractory_n_a', 'gamma_refractory_N^A', 0.1),
            ('gamma_refractory_c_a', 'gamma_refractory_C^A', 0.1),

            # excretion
            # excretion fraction - % day-1
            ('excretion_n_p', 'Excretion_N^P', 0.01),
            ('excretion_c_p', 'Excretion_C^P', 0.01),
            ('excretion_n_a', 'Excretion_N^A', 0.01),
            ('excretion_c_a', 'Excretion_C^A', 0.01),

            # uptake (umol N cell-1 day-1), /86400 for per sec
            ('v_n_max_p', 'V_N_max^P', 1.9e-9),
            ('v_c_max_p', 'V_C_max^P', 8e-9),
            ('v_in_max_a', 'V_IN_max^A', 1.9e-9 * pro_alt_vol_ratio**0.67 / 10 ),
            ('v_n_max_a', 'V_N_max^A', 1.9e-9 * pro_alt_vol_ratio**0.67 ),
            ('v_c_max_a', 'V_C_max^A', 8e-9* pro_alt_vol_ratio**0.67),
            # ('v_n_max_p', 'V_N_max^P', 9.1e-9 * vol**0.67),
            # ('v_c_max_p', 'V_C_max^P', 9.1e-9 * vol**0.67),
            # ('v_n_max_a', 'V_N_max^A', 9.1e-9 * vol**0.67),
            # ('v_c_max_a', 'V_C_max^A', 9.1e-9 * vol**0.67),

            # growth rates. day-1 (was 1e-5 sec-1),
            ('mu_inf_p', 'mu_inf^P', 0.86),
            ('mu_inf_a', 'mu_inf^A', 0.86*6),

            # k - half saturation umol Liter-1

            ('k_n_p', 'k_N^P', 0.17 * pro_vol**0.27),
            ('k_c_p', 'k_C^P', 0.17 * pro_vol**0.27),
            ('k_n_a', 'k_ON^A', 0.17 * alt_vol**0.27),
            ('k_in_a', 'k_N^A', 0.17 * alt_vol**0.27 / 10),
            ('k_c_a', 'k_OC^A', 0.17 * alt_vol**0.27),


            # minimum N quota (umol N cell-1 )
            # quota
            ('q_n_min_p', 'Q_N_min^P', 7e-10),
            ('q_c_min_p', 'Q_C_min^P', 3e-11),
            ('q_n_min_a', 'Q_N_min^A', 7e-10 * pro_alt_vol_ratio** 0.77),
            ('q_c_min_a', 'Q_C_min^A', 3e-11 * pro_alt_vol_ratio** 0.77),
            ('q_n_max_p', 'Q_N_max^P', 1.5e-9),
            ('q_c_max_p', 'Q_C_max^P', 6e-9),
            ('q_n_max_a', 'Q_N_max^A', 1.5e-9 * pro_alt_vol_ratio** 0.77),
            ('q_c_max_a', 'Q_C_max^A', 6e-9 * pro_alt_vol_ratio** 0.77),

            # ('q_n_min_p', 'Q_N_min^P', (1.36e-9) * vol ** 0.77),
            # ('q_c_min_p', 'Q_C_min^P', (1.36e-9) * redfield_C_to_N * vol ** 0.77),
            # ('q_n_min_a', 'Q_N_min^A', (1.36e-9) * vol ** 0.77),
            # ('q_c_min_a', 'Q_C_min^A', (1.36e-9) * redfield_C_to_N * vol ** 0.77),
            # ('q_n_max_p', 'Q_N_max^P', (1.36e-9) * 5 * vol ** 0.77),
            # ('q_c_max_p', 'Q_C_max^P', (1.36e-9) * 3 * redfield_C_to_N * vol ** 0.77),
            # ('q_n_max_a', 'Q_N_max^A', (1.36e-9) * 5 * vol ** 0.77),
            # ('q_c_max_a', 'Q_C_max^A', (1.36e-9) * 3 * redfield_C_to_N * vol ** 0.77),

            # mortality - day-1
            ('mortality_p', 'mortality^P', 0.1),
            ('mortality_a', 'mortality^A', 0.1),

            # CO2 parameters
            ('r0_p', 'r0^P', 0.18),   # dark respiration, sec-1 = 0.18 d-1, Geider & Osborne 1989
            ('r0_a', 'r0^A', 0.18),   # dark respiration, sec-1 = 0.18 d-1, Geider & Osborne 1989
            ('b_p', 'b^P', 0.01),   # respiration coefficient, no units, Geider & Osborne 1989
            ('b_a', 'b^A', 0.01),   # respiration coefficient, no units, Geider & Osborne 1989
            # seconds per day
            ('delta_t', 'Delta_t', 1/(24*3600)),
        ]

        self.parameters = {i[0] : ModelParam(*i) for i in parameters}

    def get_param_val(self, param_name):
        return self.parameters[param_name].value

    def get_variable_val(self, var_name):
        return self.variables[var_name].value

    def init_variables(self):
        redfield_C_to_N = 6.625

        variables = [

            # biomass (umol N l-1 ),
            ('b_n_p', 'B_N^P', self.get_param_val('q_n_min_p')*1e9, 'PRO'),
            #('b_c_p', 'B_C^P', self.get_param_val('q_c_min_p')*1e9, 'PRO'),
            ('b_n_a', 'B_N^A', self.get_param_val('q_n_min_a')*1e10, 'ALT'),
            #('b_c_a', 'B_C^A', self.get_param_val('q_c_min_a')*1e10, 'ALT'),

            # ('b_n_p', 'B_N^P', 0.5148),
            # ('b_c_p', 'B_C^P', 0.5148 * redfield_C_to_N),
            # ('b_n_a', 'B_N^A', 0.5148),
            # ('b_c_a', 'B_C^A', 0.5148 * redfield_C_to_N),

            # start cell density (cell L-1),
            ('x_p', 'X^P', 1e9, 'PRO'),
            ('x_a', 'X^A', 1e10, 'ALT'),

            # Nutrients (umol N l-1 ),
            ('n', 'N', 100),
            #('c', 'C', 2000),
            ('on', 'ON', 20),
            #('oc', 'OC', 60),
            ('on_refractory', 'ON_refractory', 0),
            #('oc_refractory', 'OC_refractory', 0),

        ]

        self.variables = {i[0] : ModelVariable(*i) for i in variables}

    def init_intermediate_variables(self):
        variables = [
            ('q_n_p', 'Q_N^P', None),
            #('q_c_p', 'Q_C^P', None),
            ('q_n_a', 'Q_N^A', None),
            #('q_c_a', 'Q_C^A', None),
            ('mu_p', 'mu^P', None),
            ('mu_a', 'mu^A', None),
            ('delta_b_n_p_uptake', 'B_N_uptake^P', None ),
            #('delta_b_on_p_uptake', 'B_ON_uptake^P', None ),
            ('delta_b_n_p_mortality', 'B_N_mortality^P', None  ),
            ('delta_b_n_p_excretion', 'B_N_excretion^P', None  ),
            ('delta_b_n_a_uptake', 'B_ON_uptake^A', None),
            ('delta_b_in_a_uptake', 'B_N_uptake^A', None),
            ('delta_b_n_a_mortality', 'B_N_mortality^A', None),
            ('delta_b_n_a_excretion', 'B_N_excretion^A', None),

            #('delta_b_c_p_uptake', 'B_C_uptake^P', None),
            #('delta_b_c_p_mortality', 'B_C_mortality^P', None),
            #('delta_b_c_p_excretion', 'B_C_excretion^P', None),
            #('delta_b_c_p_respiration', 'B_C_respiration^P', None),
            #('delta_b_c_a_uptake', 'B_OC_uptake^A', None),
            #('delta_b_c_a_mortality', 'B_C_mortality^A', None),
            #('delta_b_c_a_excretion', 'B_C_excretion^A', None),
            #('delta_b_c_a_respiration', 'B_C_respiration^A', None),
        ]
        self.intermediate_evaluation_order = [i[0] for i in variables]
        self.intermediate_variables = {i[0] : ModelVariable(*i) for i in variables}


    def init_symbols(self):
        self._symbols = {i: v.symbol for i,v in self.parameters.items()}
        self._symbols.update({i: v.symbol for i,v in self.variables.items()})
        self._symbols.update({i: v.symbol for i,v in self.intermediate_variables.items()})

    def symbol(self, name):
        return self._symbols[name]

    def add_formula(self, name, formula):
        if name in self.variables:
            self.variables[name].formula = formula
        else:
            self.intermediate_variables[name].formula = formula

    def init_formulas(self):

        s = lambda x : self.symbol(x)
        # formulas
        for n in ('n'):
            for i in ('p','a'):
                delta_quota = (s(f'q_{n}_max_{i}') - s(f'q_{n}_{i}')) / (s(f'q_{n}_max_{i}') - s(f'q_{n}_min_{i}'))
                delta_nutrient = s(n) / (s(n) + s(f'k_{n}_{i}'))
                # alt uptakes Organic carbon and nitrogen
                if (i == 'a'):
                    delta_nutrient = s(f'o{n}') / (s(f'o{n}') + s(f'k_{n}_{i}'))
                        
                self.add_formula(f'delta_b_{n}_{i}_uptake', s(f'v_{n}_max_{i}') * delta_nutrient * delta_quota * s(f'x_{i}'))
                self.add_formula(f'delta_b_{n}_{i}_mortality', s(f'mortality_{i}') * s(f'b_{n}_{i}'))
                self.add_formula(f'delta_b_{n}_{i}_excretion', s(f'excretion_{n}_{i}') *s (f'b_{n}_{i}'))
                self.add_formula(f'q_{n}_{i}', s(f'b_{n}_{i}') / s(f'x_{i}'))
        # inorganic N uptake by alt
        n = 'in'
        i = 'a'
        delta_quota = (s(f'q_n_max_{i}') - s(f'q_n_{i}')) / (s(f'q_n_max_{i}') - s(f'q_n_min_{i}'))
        delta_nutrient = s('n') / (s('n') + s(f'k_{n}_{i}'))                
        self.add_formula(f'delta_b_{n}_{i}_uptake', s(f'v_{n}_max_{i}') * delta_nutrient * delta_quota * s(f'x_{i}'))
        
        # cells
        for i in ('p', 'a'):
            #limit = max(q_n_min / q_n, q_c_min / q_c);
            #mu = mu_inf * (1 - limit);
            self.add_formula(
                f'mu_{i}',
                (1 - (s(f'q_n_min_{i}') / s(f'q_n_{i}'))) * s(f'mu_inf_{i}')
            )
            self.add_formula(
                f'x_{i}',
                s(f'x_{i}') + (s(f'x_{i}') * s(f'mu_{i}') - s(f'x_{i}') * s(f'mortality_{i}')) * s(f'delta_t')
            )


        # Nitrogen
        n = 'n'
        i = 'p'
        self.add_formula(
            f'b_{n}_{i}',
            s(f'b_{n}_{i}') + (s(f'delta_b_{n}_{i}_uptake') - s(f'delta_b_{n}_{i}_mortality') - s(f'delta_b_{n}_{i}_excretion'))  * s(f'delta_t'))
        i = 'a'
        self.add_formula(
            f'b_{n}_{i}',
            s(f'b_{n}_{i}') + (s(f'delta_b_{n}_{i}_uptake') + s(f'delta_b_i{n}_{i}_uptake') 
                                - s(f'delta_b_{n}_{i}_mortality') - s(f'delta_b_{n}_{i}_excretion'))  * s(f'delta_t'))

        self.add_formula(
            'n',
            s('n') + (
                    - s(f'delta_b_n_p_uptake')
                    - s(f'delta_b_in_a_uptake')
                    + s(f'gamma_n_a') * s(f'delta_b_n_a_excretion')
                    + s(f'gamma_n_p') * s(f'delta_b_n_p_excretion')
                    + s(f'gamma_n_p') * s(f'delta_b_n_p_mortality')
                    + s(f'gamma_n_a') * s(f'delta_b_n_a_mortality')
            ) * s(f'delta_t')
        )
        #dON_dt = (1-gamma_n) * m_star * b_n + epsilon_n*b_n;
        self.add_formula(
            'on',
            s('on') + (
                    - s(f'delta_b_n_a_uptake')
                    + (1 - s(f'gamma_n_p')) * s(f'delta_b_n_p_excretion') * (1 - s(f'gamma_refractory_n_p')) 
                    + (1 - s(f'gamma_n_a')) * s(f'delta_b_n_a_excretion') * (1 - s(f'gamma_refractory_n_a')) 
                    + (1 - s(f'gamma_n_p')) * s(f'delta_b_n_p_mortality') * (1 - s(f'gamma_refractory_n_p'))
                    + (1 - s(f'gamma_n_a')) * s(f'delta_b_n_a_mortality') * (1 - s(f'gamma_refractory_n_a'))
            ) * s(f'delta_t')
        )

        self.add_formula(
            'on_refractory',
            s('on_refractory') + (
                    + s(f'delta_b_n_p_excretion') * (s(f'gamma_refractory_n_p')) 
                    + (1 - s(f'gamma_n_p')) * s(f'delta_b_n_p_excretion') * (s(f'gamma_refractory_n_p')) 
                    + (1 - s(f'gamma_n_a')) * s(f'delta_b_n_a_excretion') * (s(f'gamma_refractory_n_a')) 
                    + (1 - s(f'gamma_n_p')) * s(f'delta_b_n_p_mortality') * (s(f'gamma_refractory_n_p'))
                    + (1 - s(f'gamma_n_a')) * s(f'delta_b_n_a_mortality') * (s(f'gamma_refractory_n_a'))
            ) * s(f'delta_t')
        )


    def init_var_list(self):
        if self._disable_organism is not None:
            self.variables_list = [v for v in self.variables.values() if v.organism != self._disable_organism]
        else:
            self.variables_list = list(self.variables.values())
            #list(self.intermediate_variables.values()) + list(self.variables.values())

    def lambadify(self):
        param_values = {i.symbol: i.value for i in self.parameters.values()}
        if self._disable_organism is not None:
            param_values.update({i.symbol : i.value for i in self.variables.values() if i.organism == self._disable_organism})
        var_symbols = [(v.symbol, v.lambda_symbol) for v in self.variables_list]
        #for i in self.intermediate_variables.values():
        #    i.run_lambdify(param_values, var_symbols)

        for i in self.variables_list:
            i.run_lambdify(param_values, var_symbols, self.intermediate_variables, self.intermediate_evaluation_order)

    # dxa_dt= (mu- m_star)*x;
    #i = delta_b_c_a_excretion.subs({b_c_a'): 5, b_c_p'): 3})


    def evaluate(self):
        #for n in self.intermediate_evaluation_order:
        #    i = self.intermediate_variables[n]
        #    var_values = [v.value for v in self.variables_list]
        #    i.evaluate(var_values)
        #    i.copy_from_nextvalue()

        # update the values computed for intermediate vars
        var_values = [v.value for v in self.variables_list]
        for i in self.variables_list:
            i.evaluate(var_values)

        # separate loop - so next values will not impact computations
        for i in self.variables_list:
            i.copy_from_nextvalue()

    def current_values(self, current_iteration):
        res = {v.name : v.value for v in self.variables.values()}
        #res.update({v.name : v.value for v in self.intermediate_variables.values()})
        res['day'] = current_iteration * self.get_param_val('delta_t')
        return res

    def disable_organism(self, organism):
        assert organism in ['PRO', 'ALT']
        self._disable_organism  =  organism

    def simulate(self, num_iterations, collect_every=3600*2):
        print('.', end='')
        self._set_initial_biomass()
        self._update_ref_iters()
        self.init_var_list()

        self.lambadify()
        results = list()
        ref_results = list()
        for i in range(num_iterations):
            if (collect_every is not None) and (i % collect_every  == 0):
                results.append(self.current_values(i))
                print('.', end='')

            if i in self.reference_iterations:
                ref_results.append(self.current_values(i))

            self.evaluate()
        return results, ref_results

    def set_referece_times(self, reference_times):
        """ set reference times (times where there is a reference value
        :param reference_values: list of numbers (reference time in days)
        """
        self.reference_days = reference_times

    def _update_ref_iters(self):
        delta_t = self.get_param_val('delta_t')
        self.reference_iterations = {int(i / delta_t) for i in self.reference_days}
        #print(self.reference_iterations)



    def override_param_values(self, values):
        for n,val in values.items():
            self.parameters[n].value = val

    def override_initial_values(self, values):
        for n,val in values.items():
            self.variables[n].override_initial_value(val)

    def _set_initial_biomass(self):
        # initial biomass is a function of quota and number of cells
        self.variables['b_n_p'].override_initial_value(self.get_param_val('q_n_min_p') * self.get_variable_val('x_p'))
        self.variables['b_c_p'].override_initial_value(self.get_param_val('q_c_min_p') * self.get_variable_val('x_p'))
        self.variables['b_n_a'].override_initial_value(self.get_param_val('q_n_min_a') * self.get_variable_val('x_a'))
        self.variables['b_c_a'].override_initial_value(self.get_param_val('q_c_min_a') * self.get_variable_val('x_a'))


    def reset_initial_values(self):
        for v in self.variables.values():
            v.reset_initial_val()


    def print_formulas(self):
        init_printing()
        for i in self.intermediate_variables:
            print(i, self.intermediate_variables[i].formula)
            print()

        for i in self.variables:
            print(i, self.variables[i].formula)
            print()

    def print_initial_values(self):
        for i in self.variables:
            print(i, self.variables[i].initval)

    def print_param_values(self):
        for i in self.parameters:
            print(i, self.parameters[i].value)


class DisplayModel:
    PRO_COLOR = 'MediumSeaGreen'
    PRO_FL_COLOR = 'DarkGreen'
    PRO_FCM_COLOR = 'MediumSpringGreen'

    ALT_COLOR = 'Gold'
    ALT_FCM_COLOR = 'PaleGoldenrod'
    N_COLOR = 'royalblue'
    ON_COLOR = 'PowderBlue'
    C_COLOR = 'FireBrick'
    OC_COLOR = 'LightCoral'
    QUOTA_COLOR='red'

    mcolors = {
        'b_n_p' : PRO_COLOR, 
        'b_c_p': PRO_COLOR, 
        'b_n_a': ALT_COLOR , 
        'b_c_a' : ALT_COLOR, 
        'x_p': PRO_COLOR, 
        'x_a': ALT_COLOR, 
        'n' : N_COLOR, 
        'c': C_COLOR, 
        'on': ON_COLOR, 
        'oc':OC_COLOR,
        'on_refractory': ON_COLOR, 
        'oc_refractory':OC_COLOR,
        'q_n_p' :PRO_COLOR, 
        'q_c_p': PRO_COLOR, 
        'q_n_a' :ALT_COLOR, 
        'q_c_a': ALT_COLOR, 
    }

    def __init__(self, res_df, model, model_name, reference_FL_df, reference_FCM_df, pro_only, max_day):
        self.res_df = self.apply_max_day(res_df, max_day)
        self.model=model
        self.model_name = model_name
        self.reference_FL_df = self.apply_max_day(reference_FL_df, max_day)
        self.reference_FCM_df = self.apply_max_day(reference_FCM_df, max_day)
        self.pro_only = pro_only
        
    def apply_max_day(self, df, max_day):
        if max_day is None:
            return df
        return df.loc[df['day'] < max_day]
        
    def display_quota(self, type, nutrient):
        i = type
        n = nutrient
        q = f'q_{n}_{i}'
        self.res_df[q] = self.res_df[f'b_{n}_{i}'] / self.res_df[f'x_{i}']
        sns.lineplot(data=self.res_df, x='day', y=q, label=q, color=self.mcolors[q], lw=5)

        for q in [
             f'q_{n}_min_{i}', f'q_{n}_max_{i}', 
        ]:
            sns.lineplot(data=self.res_df, x='day', y=self.model.get_param_val(q), label=q, color=self.QUOTA_COLOR)
        plt.legend(bbox_to_anchor=(1,1))
        plt.ylabel(f'umol {n.upper()}/cell')
        plt.title(f'{self.model_name} - Quota q_{n}_{i}')
        
        plt.show()

    def display_biomass(self, nutrient):
        n = nutrient
        for i in [
            f'b_{n}_p', f'b_{n}_a',  f'{n}', f'o{n}', f'o{n}_refractory'
           ]:
            if self.pro_only and i == f'b_{n}_a':
                # skip
                continue
            
            sns.lineplot(data=self.res_df, x='day', y=i, label=i, lw=5, color=self.mcolors[i])
        plt.legend(bbox_to_anchor=(1,1))
        plt.ylabel(f'umol {n.upper()}/L')
        plt.title(f'{self.model_name} - {n.upper()} Concentration')
        plt.show()
        
    def display_cells(self, logscale=False):
        xlst = ['x_p', 'x_a',]
        if self.pro_only:
            xlst = ['x_p']
        
        for i in xlst:
            sns.lineplot(data=self.res_df, x='day', y=i, label=i, lw=5, color=self.mcolors[i])

        sns.scatterplot(x=self.reference_FL_df.day, y=self.reference_FL_df.cells*1000, 
                        label='PRO, ref FL', color=self.PRO_FL_COLOR, s=50)
        sns.scatterplot(x=self.reference_FCM_df.day, y=self.reference_FCM_df['PRO.1']*1000, 
                        label='PRO, ref FCM', color=self.PRO_FCM_COLOR, s=100)
        if not self.pro_only:
            sns.scatterplot(x=self.reference_FCM_df.day, y=self.reference_FCM_df['ALT.1']*1000, 
                            label='ALT, ref FCM', color=self.ALT_FCM_COLOR, s=100)
            
        plt.legend(bbox_to_anchor=(1,1))
        plt.ylabel(f'cells/L')

        if logscale:
            plt.yscale('log')
        plt.title(f'{self.model_name} - X (cells/L)')
        plt.show()

    
    

def display_simulation_results(res_df, model, model_name, reference_FL_df=None, reference_FCM_df=None, pro_only=False, max_day=None):
    d = DisplayModel(res_df, model, model_name, reference_FL_df, reference_FCM_df, pro_only, max_day)
    # cells
    d.display_cells()
    d.display_cells(logscale=True)

    #biomass
    d.display_biomass('n')
    #d.display_biomass('c')

    # quotas
    d.display_quota('p', 'n')
    #d.display_quota('p', 'c')
    if not pro_only:
        d.display_quota('a', 'n')
        #d.display_quota('a', 'c')
    

    




def model_optimize_PRO(param_values, ref_df, param_names):
    m = ModelProALT()
    m.disable_organism('ALT')
    m.override_initial_values({
        'b_n_a': 0,
        'b_c_a': 0,
        'x_a': 0,
    })
    m.override_param_values({i : k for i, k in zip(param_names, param_values)})
    reference_days = ref_df['day'].unique().tolist()
    max_day = int(ref_df['day'].max()) + 1
    m.set_referece_times(reference_days)
    _, ref_res = m.simulate(num_iterations=max_day*3600*24, collect_every=None)
    ref_res_df = pd.DataFrame(ref_res)
    t = pd.merge(ref_res_df[['day', 'x_p']], ref_df[['day', 'VALUE']], on='day')
    #print(t)
    if (~np.isfinite(t)).any().any() or t.isna().any().any() or (t.shape[0] == 0):
        # bad solution
        print('bad')
        return 100000
    res = metrics.mean_squared_log_error(t['VALUE'], t['x_p'])
    print(res)
    return res


if __name__ == '__main__':
    m = ModelProALT()
    m.disable_organism('ALT')
    m.override_initial_values({
        'b_n_a': 0,
        'b_c_a': 0,
        'x_a': 0,
    })
    m.set_referece_times([0,0.5,0.6])
    #m.print_formulas()
    #m.lambadify()
    #m.evaluate()
    res, ref_res = m.simulate(num_iterations=1*3600*24)
    import pprint
    pprint.pprint(res)
    pprint.pprint(ref_res)