#!/usr/bin/env python
# coding: utf-8



import sys, os
import math
from sympy import *


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
    def __init__(self, name, symbolname, initval):
        self.name = name
        self.initval = initval
        self.symbol = symbols(symbolname)
        self.lambda_symbol = symbols(symbolname.replace('^','_'))
        self.value = initval
        self.nextvalue = None
        self.formula = None
        self.lambdafun = None

    def run_lambdify(self, param_values, variables_list, intermediate_variables, intermediate_evaluation_order):
        s = self.formula
        for n in reversed(intermediate_evaluation_order):
            v = intermediate_variables[n]
            #print(s, v.symbol, v.formula)
            s = s.subs({v.symbol: v.formula})
            #print('after:', s)

        subs_expr = s.subs(param_values)
        subs_expr2 = subs_expr.subs(variables_list)
        print('before:', self.formula)
        print('after', subs_expr2)
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
        radius = 0.3628;  # "MED4" = 9312
        vol = (4 / 3) * math.pi * radius ** 3;

        parameters = [
            # (name, symbolname, initval)
            # gamma=fraction of heterotroph mort/resp to inorganic form
            ('gamma_n_p', 'gamma_N^P', 0.04),
            ('gamma_c_p', 'gamma_C^P', 0.04),
            ('gamma_n_a', 'gamma_N^A', 0.04),
            ('gamma_c_a', 'gamma_C^A', 0.04),

            # excretion
            # excretion fraction - % day-1
            ('excretion_n_p', 'Excretion_N^P', 0.01),
            ('excretion_c_p', 'Excretion_C^P', 0.01),
            ('excretion_n_a', 'Excretion_N^A', 0.01),
            ('excretion_c_a', 'Excretion_C^A', 0.01),

            # uptake (umol N cell-1 day-1), /86400 for per sec
            ('v_n_max_p', 'V_N_max^P', 1.9e-9),
            ('v_c_max_p', 'V_C_max^P', 8e-9),
            ('v_n_max_a', 'V_N_max^A', 1.9e-9),
            ('v_c_max_a', 'V_C_max^A', 8e-9),
            # ('v_n_max_p', 'V_N_max^P', 9.1e-9 * vol**0.67),
            # ('v_c_max_p', 'V_C_max^P', 9.1e-9 * vol**0.67),
            # ('v_n_max_a', 'V_N_max^A', 9.1e-9 * vol**0.67),
            # ('v_c_max_a', 'V_C_max^A', 9.1e-9 * vol**0.67),

            # growth rates. day-1 (was 1e-5 sec-1),
            ('mu_inf_p', 'mu_inf^P', 0.86),
            ('mu_inf_a', 'mu_inf^A', 0.86),

            # k - half saturation umol Liter-1

            ('k_n_p', 'k_N^P', 0.17 * vol**0.27),
            ('k_c_p', 'k_C^P', 0.17 * vol**0.27),
            ('k_n_a', 'k_N^A', 0.17 * vol**0.27),
            ('k_c_a', 'k_C^A', 0.17 * vol**0.27),


            # minimum N quota (umol N cell-1 )
            # quota
            ('q_n_min_p', 'Q_N_min^P', 7e-10),
            ('q_c_min_p', 'Q_C_min^P', 3e-11),
            ('q_n_min_a', 'Q_N_min^A', 7e-10),
            ('q_c_min_a', 'Q_C_min^A', 3e-11),
            ('q_n_max_p', 'Q_N_max^P', 1.5e-9),
            ('q_c_max_p', 'Q_C_max^P', 6e-9),
            ('q_n_max_a', 'Q_N_max^A', 1.5e-9),
            ('q_c_max_a', 'Q_C_max^A', 6e-9),

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
            ('b_n_p', 'B_N^P', self.get_param_val('q_n_min_p')*1e6),
            ('b_c_p', 'B_C^P', self.get_param_val('q_c_min_p')*1e6),
            ('b_n_a', 'B_N^A', self.get_param_val('q_n_min_a')*1e7),
            ('b_c_a', 'B_C^A', self.get_param_val('q_c_min_a')*1e7),

            # ('b_n_p', 'B_N^P', 0.5148),
            # ('b_c_p', 'B_C^P', 0.5148 * redfield_C_to_N),
            # ('b_n_a', 'B_N^A', 0.5148),
            # ('b_c_a', 'B_C^A', 0.5148 * redfield_C_to_N),

            # start cell density (cell L-1),
            ('x_p', 'X^P', 1e6),
            ('x_a', 'X^A', 1e7),

            # Nutrients (umol N l-1 ),
            ('n', 'N', 800/8),
            ('c', 'C', 3000),
            ('on', 'ON', 20),
            ('oc', 'OC', 20 * redfield_C_to_N ),

        ]

        self.variables = {i[0] : ModelVariable(*i) for i in variables}

    def init_intermediate_variables(self):
        variables = [
            ('q_n_p', 'Q_N^P', None),
            ('q_c_p', 'Q_C^P', None),
            ('q_n_a', 'Q_N^A', None),
            ('q_c_a', 'Q_C^A', None),
            ('mu_p', 'mu^P', None),
            ('mu_a', 'mu^A', None),
            ('delta_b_n_p_uptake', 'B_N_uptake^P', None ),
            ('delta_b_n_p_mortality', 'B_N_mortality^P', None  ),
            ('delta_b_n_p_excretion', 'B_N_excretion^P', None  ),
            ('delta_b_n_a_uptake', 'B_N_uptake^A', None),
            ('delta_b_n_a_mortality', 'B_N_mortality^A', None),
            ('delta_b_n_a_excretion', 'B_N_excretion^A', None),

            ('delta_b_c_p_uptake', 'B_C_uptake^P', None),
            ('delta_b_c_p_mortality', 'B_C_mortality^P', None),
            ('delta_b_c_p_excretion', 'B_C_excretion^P', None),
            ('delta_b_c_p_respiration', 'B_C_respiration^P', None),
            ('delta_b_c_a_uptake', 'B_C_uptake^A', None),
            ('delta_b_c_a_mortality', 'B_C_mortality^A', None),
            ('delta_b_c_a_excretion', 'B_C_excretion^A', None),
            ('delta_b_c_a_respiration', 'B_C_respiration^A', None),
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
        for n in ('n','c'):
            for i in ('p','a'):
                delta_quota = (s(f'q_{n}_max_{i}') - s(f'q_{n}_{i}')) / (s(f'q_{n}_max_{i}') - s(f'q_{n}_min_{i}'))
                delta_nutrient = s(n) / (s(n) + s(f'k_{n}_{i}'))
                # alt uptakes Organic carbon
                if (n == 'c') and (i == 'a'):
                    delta_nutrient = s('oc') / (s('oc') + s(f'k_{n}_{i}'))
                self.add_formula(f'delta_b_{n}_{i}_uptake', s(f'v_{n}_max_{i}') * delta_nutrient * delta_quota * s(f'x_{i}'))
                self.add_formula(f'delta_b_{n}_{i}_mortality', s(f'mortality_{i}') * s(f'b_{n}_{i}'))
                self.add_formula(f'delta_b_{n}_{i}_excretion', s(f'excretion_{n}_{i}') *s (f'b_{n}_{i}'))
                self.add_formula(f'q_{n}_{i}', s(f'b_{n}_{i}') / s(f'x_{i}'))
        # cells
        for i in ('p', 'a'):
            #limit = max(q_n_min / q_n, q_c_min / q_c);
            #mu = mu_inf * (1 - limit);
            self.add_formula(
                f'mu_{i}',
                (1 - Max((s(f'q_n_min_{i}') / s(f'q_n_{i}')), (s(f'q_c_min_{i}') / s(f'q_c_{i}')))) * s(f'mu_inf_{i}')
            )
            self.add_formula(
                f'x_{i}',
                s(f'x_{i}') + (s(f'x_{i}') * (s(f'mu_{i}') - s(f'mortality_{i}')) * s(f'delta_t'))
            )


        # Nitrogen
        n = 'n'
        for i in ('p', 'a'):
            self.add_formula(
                f'b_{n}_{i}',
                s(f'b_{n}_{i}') + (s(f'delta_b_{n}_{i}_uptake') - s(f'delta_b_{n}_{i}_mortality') - s(f'delta_b_{n}_{i}_excretion'))  * s(f'delta_t'))

        self.add_formula(
            'n',
            s('n') + (
                    - s(f'delta_b_n_a_uptake')
                    - s(f'delta_b_n_p_uptake')
                    + s(f'gamma_n_p') * s(f'delta_b_n_p_mortality')
                    + s(f'gamma_n_a') * s(f'delta_b_n_a_mortality')
            ) * s(f'delta_t')
        )
        #dON_dt = (1-gamma_n) * m_star * b_n + epsilon_n*b_n;
        self.add_formula(
            'on',
            s('on') + (
                      s(f'delta_b_n_p_excretion')
                    + s(f'delta_b_n_p_excretion')
                    + (1 - s(f'gamma_n_p')) * s(f'delta_b_n_p_mortality')
                    + (1 - s(f'gamma_n_a')) * s(f'delta_b_n_a_mortality')
            ) * s(f'delta_t')
        )

        # Carbon
        n ='c'
        for i in ('p', 'a'):
            self.add_formula(f'delta_b_c_{i}_respiration', (s(f'r0_{i}') + s(f'b_{i}') * s(f'mu_{i}')) * s(f'b_c_{i}'))
            self.add_formula(
                f'b_{n}_{i}',
                s(f'b_{n}_{i}') + (
                        s(f'delta_b_{n}_{i}_uptake')
                        - s(f'delta_b_{n}_{i}_mortality')
                        - s(f'delta_b_{n}_{i}_excretion')
                        - s(f'delta_b_{n}_{i}_respiration')
                )  * s(f'delta_t'))

        # dC_dt = -(V_c_max * C * (q_c_max - q_c) * x) / ((C + k_c) * (q_c_max - q_c_min)) + ...
        #    gamma_c * m * b_c + (r0 + b * mu) * b_c - (C - Csat) / ((h * C) / (Kg * B * CO2_dis));
        # TODO: (C - Csat) / ((h * C) / (Kg * B * CO2_dis));
        self.add_formula(
            'c',
            s('c') + (
                    - s(f'delta_b_c_p_uptake')
                    + s(f'gamma_c_p') * s(f'delta_b_c_p_mortality')
                    + s(f'gamma_c_a') * s(f'delta_b_c_a_mortality')
                    + s(f'delta_b_c_a_excretion')
                    + s(f'delta_b_c_p_respiration')
                    + s(f'delta_b_c_a_respiration')
            ) * s(f'delta_t')
        )
        self.add_formula(
            'oc',
            s('oc') + (
                    - s(f'delta_b_c_a_uptake')
                    + s(f'delta_b_c_p_excretion')
                    + (1 - s(f'gamma_c_p')) * s(f'delta_b_c_p_mortality')
                    + (1 - s(f'gamma_c_a')) * s(f'delta_b_c_a_mortality')
            ) * s(f'delta_t')
        )

    def init_var_list(self):
        self.variables_list = list(self.variables.values()) #list(self.intermediate_variables.values()) + list(self.variables.values())

    def lambadify(self):
        param_values = {i.symbol: i.value for i in self.parameters.values()}
        var_symbols = [(v.symbol, v.lambda_symbol) for v in self.variables_list]
        #for i in self.intermediate_variables.values():
        #    i.run_lambdify(param_values, var_symbols)

        for i in self.variables.values():
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
        for i in self.variables.values():
            i.evaluate(var_values)

        # separate loop - so next values will not impact computations
        for i in self.variables.values():
            i.copy_from_nextvalue()

    def current_values(self, current_iteration):
        res = {v.name : v.value for v in self.variables.values()}
        #res.update({v.name : v.value for v in self.intermediate_variables.values()})
        res['day'] = current_iteration / self.get_param_val('delta_t')
        return res

    def simulate(self, num_iterations, collect_every=3600*2):
        self._set_initial_biomass()
        self._update_ref_iters()
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
        self.reference_iterations = {int(i * delta_t) for i in self.reference_days}



    def override_param_values(self, values):
        for n,val in values.items():
            self.parameters[n].value = val

    def override_initial_values(self, values):
        for n,val in values.items():
            self.variables[n].override_initial_val(val)

    def _set_initial_biomass(self):
        # initial biomass is a function of quota and number of cells
        self.variables['b_n_p'].override_initial_val(self.get_param_val('q_n_min_p') * self.get_variable_val('x_p'))
        self.variables['b_c_p'].override_initial_val(self.get_param_val('q_c_min_p') * self.get_variable_val('x_p'))
        self.variables['b_n_a'].override_initial_val(self.get_param_val('q_n_min_a') * self.get_variable_val('x_a'))
        self.variables['b_c_a'].override_initial_val(self.get_param_val('q_c_min_a') * self.get_variable_val('x_a'))


    def reset_initial_values(self):
        for v in self.variables_list:
            v.reset_initial_val()


    def print_formulas(self):
        init_printing()
        for i in self.intermediate_variables:
            print(i, self.intermediate_variables[i].formula)
            print()

        for i in self.variables:
            print(i, self.variables[i].formula)
            print()

if __name__ == '__main__':
    m = ModelProALT()
    #m.print_formulas()
    #m.lambadify()
    #m.evaluate()
    res = m.simulate(num_iterations=1*3600*24)
    import pprint
    pprint.pprint(res)