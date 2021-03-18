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
from matplotlib.backends.backend_pdf import PdfPages


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
    def __init__(self, name, symbolname, initval, organism = None, nutrient=None):
        self.name = name
        self.initval = initval
        self.symbol = symbols(symbolname)
        self.lambda_symbol = symbols(symbolname.replace('^','_'))
        self.value = initval
        self.nextvalue = None
        self.formula = None
        self.lambdafun = None
        self.organism = organism
        self.nutrient = nutrient

    def run_lambdify(self, param_values, variables_list, intermediate_variables, intermediate_evaluation_order):
        s = self.formula
        for n in reversed(intermediate_evaluation_order):
            v = intermediate_variables[n]
            #print(s, v.symbol, v.formula)
            s = s.subs({v.symbol: v.formula})
            #print('after:', s)

        #print(s, param_values)
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
        self._disable_nutrient = None
        self.reference_days = list()
        self.reference_iterations = set()
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
            # uptake (umol N cell-1 day-1), /86400 for per sec
            # k - half saturation umol Liter-1
            ('k_N_A', 'k_N^A', 0.17 * alt_vol**0.27),
            ('k_N_P', 'k_N^P',0.17 * pro_vol**0.27),
            ('v_N_A', 'v_N^A', 1.9e-9 * pro_alt_vol_ratio**0.67),
            ('v_N_P', 'v_N^P',1.9e-9),
            
            ('k_T_A', 'k_T^A', 0.05),
            ('k_T_P', 'k_T^P', 0.05),
            ('v_T_A', 'v_T^A', 1.9e-9 * pro_alt_vol_ratio**0.67),
            ('v_T_P', 'v_T^P',1.9e-9),

            ('m_A', 'm_A', 0.1),
            ('m_P', 'm_P', 0.1),
        
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
            ('P', 'P', 50),
            ('A', 'A', 50),
            ('N', 'N', 100),
            ('T', 'T', 0),
        ]

        self.variables = {i[0] : ModelVariable(*i) for i in variables}

    def init_intermediate_variables(self):
        variables = [
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
        # P uptake
        delta_uptake = dict()
        delta_toxin = dict()
        delta_mortility = dict()
        
        for i in ('P','A'):
            delta_nutrient = s('N') / (s('N') + s(f'k_N_{i}'))
            delta_toxin[i] = s('T') / (s('T') + s(f'k_T_{i}'))
            delta_uptake[i] = s(f'v_N_{i}') * delta_nutrient * s(i)
            
            delta_mortility[i] = s(f'm_{i}') * s(i)
            if i == 'P':
                delta_mortility[i] = s(f'm_{i}') * s(i) + s(f'm_{i}') * s(i) *delta_toxin[i]
            
            self.add_formula(i, s(i) + (delta_uptake[i] - delta_mortility[i])  * s(f'delta_t'))

            
        self.add_formula('N', s('N') + (delta_mortility['P'] + delta_mortility['A'] -  delta_uptake['P'] -  delta_uptake['A'])  * s(f'delta_t'))
        self.add_formula('T', s('T') + (s('v_T_P')* s('P') - delta_toxin['A'] * s('v_T_A') *s('A'))  * s(f'delta_t') )

    def init_var_list(self):
        self.variables_list = list(self.variables.values())

    def lambadify(self):
        param_values = {i.symbol: i.value for i in self.parameters.values()}
        var_symbols = [(v.symbol, v.lambda_symbol) for v in self.variables_list]
        for i in self.variables_list:
            i.run_lambdify(param_values, var_symbols, self.intermediate_variables, self.intermediate_evaluation_order)


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
        assert organism in ['PRO', 'ALT', None]
        self._disable_organism  =  organism

    def disable_nutrient(self, nutrient):
        assert nutrient in ['n', 'c', None]
        self._disable_nutrient  =  nutrient
        # override mu formula to take into account only one nutrient
        s = lambda x : self.symbol(x)
        n = 'n' if self._disable_nutrient == 'c' else 'c'
        for i in ['a', 'p']:
            self.add_formula(
                f'mu_{i}',
                (1 - (s(f'q_{n}_min_{i}') / s(f'q_{n}_{i}'))) * s(f'mu_inf_{i}')
            )


    def simulate(self, num_iterations, collect_every=3600*2):
        print('.', end='')
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
    REF_ON_COLOR = 'DeepSkyBlue'
    C_COLOR = 'FireBrick'
    OC_COLOR = 'LightCoral'
    REF_OC_COLOR = 'LightSalmon'
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
        'on_refractory': REF_ON_COLOR,
        'oc_refractory': REF_OC_COLOR,
        'q_n_p' :PRO_COLOR, 
        'q_c_p': PRO_COLOR, 
        'q_n_a' :ALT_COLOR, 
        'q_c_a': ALT_COLOR, 
    }

    def __init__(self, res_df, model, model_name, reference_FL_df, reference_FCM_df, pro_only, max_day,
                 ref_pro_col, ref_alt_col):
        self.res_df = self.apply_max_day(res_df, max_day)
        self.model=model
        self.model_name = model_name
        self.reference_FL_df = self.apply_max_day(reference_FL_df, max_day)
        self.reference_FCM_df = self.apply_max_day(reference_FCM_df, max_day)
        self.pro_only = pro_only
        self.ref_pro_col = ref_pro_col
        self.ref_alt_col = ref_alt_col
        
    def apply_max_day(self, df, max_day):
        if max_day is None:
            return df
        return df.loc[df['day'] < max_day]
        
    def display_quota(self, type, nutrient):
        if self.model is not None:
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

    def display_cells(self, logscale=False):
        xlst = ['x_p', 'x_a',]
        if self.pro_only:
            xlst = ['x_p']
        ax = None
        for i in xlst:
            ax = sns.lineplot(data=self.res_df, x='day', y=i, label=i, lw=5, color=self.mcolors[i], ax=ax)
        if self.reference_FL_df is not None:
            sns.scatterplot(x=self.reference_FL_df.day, y=self.reference_FL_df.cells*1000,
                            label='PRO, ref FL', color=self.PRO_FL_COLOR, s=50, ax=ax)
        if self.reference_FCM_df is not None:
            sns.scatterplot(x=self.reference_FCM_df.day, y=self.reference_FCM_df[self.ref_pro_col]*1000,
                            label='PRO, ref FCM', color=self.PRO_FCM_COLOR, s=100, ax=ax)
            if not self.pro_only:
                sns.scatterplot(x=self.reference_FCM_df.day, y=self.reference_FCM_df[self.ref_alt_col]*1000,
                                label='ALT, ref FCM', color=self.ALT_FCM_COLOR, s=100, ax=ax)
            
        plt.legend(bbox_to_anchor=(1,1))
        plt.ylabel(f'cells/L')

        if logscale:
            plt.yscale('log')
        plt.title(f'{self.model_name} - X (cells/L)')
        return ax

    
    

def display_simulation_results(res_df, model, model_name, reference_FL_df=None, reference_FCM_df=None,
                               pro_only=False, n_only=False, max_day=None,
                               ref_pro_col = 'PRO', ref_alt_col = 'ALT', ):
    d = DisplayModel(res_df, model, model_name, reference_FL_df, reference_FCM_df, pro_only, max_day,
                     ref_pro_col, ref_alt_col)

    sns.set(context='paper', style='white')
    # cells
    d.display_cells()
    plt.show()
    d.display_cells(logscale=True)
    plt.show()

    #biomass
    d.display_biomass('n')
    plt.show()
    if not n_only:
        d.display_biomass('c')
        plt.show()

    # quotas
    d.display_quota('p', 'n')
    plt.show()
    if not n_only:
        d.display_quota('p', 'c')
        plt.show()
    if not pro_only:
        d.display_quota('a', 'n')
        plt.show()
        if not n_only:
            d.display_quota('a', 'c')
            plt.show()


def _savefig(pdf):
    lgd = plt.legend()
    pdf.savefig(bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close()


def display_simulation_results_to_pdf(res_df, model, model_name, pdf_fpath, reference_FL_df=None, reference_FCM_df=None,
                               pro_only=False, n_only=False, max_day=None,
                                      ref_pro_col='PRO', ref_alt_col='ALT', ):
    d = DisplayModel(res_df, model, model_name, reference_FL_df, reference_FCM_df, pro_only, max_day,
                     ref_pro_col, ref_alt_col)
    with PdfPages(pdf_fpath) as pdf:

        # cells
        d.display_cells()
        _savefig(pdf)
        d.display_cells(logscale=True)
        _savefig(pdf)

        #biomass
        d.display_biomass('n')
        _savefig(pdf)
        if not n_only:
            d.display_biomass('c')
            _savefig(pdf)

        # quotas
        d.display_quota('p', 'n')
        _savefig(pdf)
        if not n_only:
            d.display_quota('p', 'c')
            _savefig(pdf)
        if not pro_only:
            d.display_quota('a', 'n')
            _savefig(pdf)
            if not n_only:
                d.display_quota('a', 'c')
                _savefig(pdf)


from scipy.optimize import differential_evolution


def genetic_optimization(param_names, ref_df, disable_organism, disable_nutrient, ref_pro_col, ref_alt_col,
                         workers, add_init, optimize_all):
    m = ModelProALT()
    reference_days = ref_df['day'].unique().tolist()
    max_day = int(ref_df['day'].max()) + 1
    num_iterations = max_day * 3600 * 24
    if add_init:
        init = [m.get_param_val(i) for i in param_names]
    else:
        init = 'latinhypercube'
    # optimize_params = dict(ref_df=ref_df, param_names=param_names, reference_days=reference_days, num_iterations=num_iterations,
    #                                          disable_organism=disable_organism, disable_nutrient=disable_nutrient,
    #                                          ref_pro_col=ref_pro_col, ref_alt_col=ref_alt_col)
    optimize_params = (ref_df, param_names, reference_days, num_iterations,
                                             disable_organism, disable_nutrient,
                                             ref_pro_col, ref_alt_col, optimize_all)

    # opt_func = lambda x : model_optimize(x, ref_df, param_names, reference_days, num_iterations,
    #                                          disable_organism=disable_organism, disable_nutrient=disable_nutrient,
    #                                          ref_pro_col=ref_pro_col, ref_alt_col=ref_alt_col, )

    def compute_bounds(param_name):
        i = m.get_param_val(param_name)
        if i < 0.1:
            return (i/10, i*10)
            #return (0, 0.2)
        else:
            return (i/10, i*2)
        
    
    bounds = [compute_bounds(i) for i in param_names]
    result = differential_evolution(model_optimize_all, bounds, disp=True, workers=workers, #updating='deferred',
                                    args=optimize_params, init=init)
    return result

def model_optimize_all(param_values, ref_df, param_names, reference_days, num_iterations, disable_organism,
                       disable_nutrient, ref_pro_col, ref_alt_col, optimize_all):
    if optimize_all:
        result = 0
        result += model_optimize(
            param_values, ref_df, param_names, reference_days, num_iterations, disable_nutrient=disable_nutrient,
            disable_organism=None, ref_pro_col='PRO_CO', ref_alt_col='ALT_CO'
        )
        result += model_optimize(
            param_values, ref_df, param_names, reference_days, num_iterations, disable_nutrient=disable_nutrient,
            disable_organism='PRO', ref_pro_col='PRO', ref_alt_col='ALT'
        )
        result += model_optimize(
            param_values, ref_df, param_names, reference_days, num_iterations, disable_nutrient=disable_nutrient,
            disable_organism='ALT', ref_pro_col='PRO', ref_alt_col='ALT'
        )
        return result
    else:
        return model_optimize(param_values, ref_df, param_names, reference_days, num_iterations, disable_organism,
                           disable_nutrient, ref_pro_col, ref_alt_col)


def model_optimize(param_values, ref_df, param_names, reference_days, num_iterations, disable_organism,
                       disable_nutrient, ref_pro_col, ref_alt_col):
    """

    :rtype: int
    """
    collect_every = None
    init_x_p, init_x_a = compute_x_init(ref_df, ref_pro_col, ref_alt_col, disable_organism)
            
    m, res, ref_res = run_model(param_values, param_names, reference_days, num_iterations, collect_every,
                                disable_organism, disable_nutrient, init_x_p, init_x_a)
    ref_res_df = pd.DataFrame(ref_res)
    if ((~np.isfinite(ref_res_df)).any().any() or
            ref_res_df.isna().any().any() or
            (ref_res_df.shape[0] == 0) or
            (ref_res_df.lt(0).any().any())
    ):
        # bad solution
        print('bad')
        return 100000

    ref_col_list = ['day']
    if disable_organism != 'ALT':
        ref_col_list.append(ref_alt_col)
    if disable_organism != 'PRO':
        ref_col_list.append(ref_pro_col)
        
    t = pd.merge(ref_res_df[['day', 'x_p', 'x_a']], ref_df[ref_col_list], on='day')
    #print(t)
    res = 0
    if disable_organism != 'PRO':
        res += metrics.mean_squared_log_error(t[ref_pro_col] * 1000, t['x_p'] )
    if disable_organism != 'ALT':
        res += metrics.mean_squared_log_error(t[ref_alt_col] * 1000, t['x_a'])
    #res = metrics.mean_squared_log_error(t['VALUE'], t['x_p'])
    print(res)
    return res


def run_model(param_values, param_names, reference_days, num_iterations, collect_every,
              disable_organism, disable_nutrient, init_x_p, init_x_a):

    m = ModelProALT()
    m.disable_organism(disable_organism)
    m.disable_nutrient(disable_nutrient)
    m.override_initial_values({
        'x_a': init_x_a,
        'x_p': init_x_p,
    })

    m.override_param_values({i: k for i, k in zip(param_names, param_values)})
    m.set_referece_times(reference_days)
    res, ref_res = m.simulate(num_iterations=num_iterations, collect_every=collect_every)
    return m, res, ref_res

def compute_x_init(ref_df, ref_pro_col, ref_alt_col, disable_organism):
    init_x_a = 0
    init_x_p = 0
    init_row = ref_df.loc[ref_df["day"] == 0]
    if disable_organism != 'ALT':
        init_x_a = 1e10
        if init_row.shape[0] > 0:
            init_x_a = init_row[ref_alt_col].min() * 1000
    if disable_organism != 'PRO':
        init_x_p = 1e9
        if init_row.shape[0] > 0:
            init_x_p = init_row[ref_pro_col].min() * 1000
    return init_x_p, init_x_a


if __name__ == '__main__':
    import argparse
    import json
    import pprint

    parser = argparse.ArgumentParser(description='Run models.')
    parser.add_argument("--optimize", help="run optimization (default: run model",
                        action="store_true")
    parser.add_argument("--opt_add_init", help="specify init values in optimization",
                        action="store_true")
    parser.add_argument('--params_opt', nargs='*', help='optimization params', default=[])
    parser.add_argument('--param_values', nargs='*', help='override param values', default=[])
    parser.add_argument('--param_names', nargs='*', help='override param names', default=[])
    parser.add_argument('--ref_csv', help='reference data (FCM cells/ml), csv file', required=True)

    parser.add_argument("--optimize_all", help="optimize all options", action="store_true")
    parser.add_argument("--disable_c", help="disable C", action="store_true")
    parser.add_argument("--disable_alt", help="disable ALT", action="store_true")
    parser.add_argument("--disable_pro", help="disable PRO", action="store_true")
    parser.add_argument('--ref_alt_col', help='name of ALT col', default='ALT')
    parser.add_argument('--ref_pro_col', help='name of PRO col', default='PRO')
    parser.add_argument('--workers', help='number of workers', type=int, default=-1)
    parser.add_argument('--maxday', help='max day of simulation', type=int, default=-1)


    #parser.add_argument("--outdir", help="output dir", default='.')
    parser.add_argument("--outfile", help="output filename", required=True)

    args = parser.parse_args()
    dpath = os.path.dirname(args.outfile)
    if dpath != '':
        os.makedirs(dpath, exist_ok=True)
    ref_df = pd.read_csv(args.ref_csv)
    disable_organism = 'ALT' if args.disable_alt else None
    if args.disable_pro:
        disable_organism = 'PRO'
    disable_nutrient = 'c' if args.disable_c else None
    model_name = os.path.splitext(os.path.basename(args.outfile))[0]

    if args.optimize:
        add_init=False
        assert(len(args.params_opt))
        result = genetic_optimization(args.params_opt, ref_df,
                                      disable_organism=disable_organism, disable_nutrient=disable_nutrient,
                                      ref_pro_col=args.ref_pro_col, ref_alt_col=args.ref_alt_col,
                                      workers= args.workers,
                                      add_init=args.opt_add_init, optimize_all=args.optimize_all)
        pprint.pprint(result)
        success = result.success
        res_dict = {
            'fun': result.fun,
             'message': result.message,
             'nfev': result.nfev,
             'nit': result.nit,
             'success': result.success,
            'model_name' : model_name,
        }
        res_dict.update({i:v for i,v in zip (args.params_opt, result.x)})
        pd.DataFrame([res_dict]).to_csv(args.outfile)

    else:
        # simulate (no optimization)
        reference_days = ref_df['day'].unique().tolist()
        max_day = args.maxday
        if max_day == -1:
            max_day = int(ref_df['day'].max()) + 1
        num_iterations = max_day * 3600 * 24
        init_x_p, init_x_a = compute_x_init(ref_df, args.ref_pro_col, args.ref_alt_col, disable_organism)
        collect_every = 3600*4
        m, res, ref_res = run_model(args.param_values, args.param_names, reference_days, num_iterations, collect_every,
                                    disable_organism, disable_nutrient, init_x_p, init_x_a)
        res_df = pd.DataFrame(res)
        pdf_fpath = f'{os.path.splitext(args.outfile)[0]}.pdf'
        display_simulation_results_to_pdf(res_df, m, model_name, pdf_fpath,
                                          #reference_FL_df=ref_fl_df,
                                          reference_FCM_df=ref_df,
                                          pro_only=args.disable_alt,
                                          n_only=args.disable_c,
                                          ref_pro_col=args.ref_pro_col, ref_alt_col=args.ref_alt_col, )
        res_df.to_csv(args.outfile)

# m = ModelProALT()
    # m.disable_organism('ALT')
    # m.override_initial_values({
    #     'b_n_a': 0,
    #     'b_c_a': 0,
    #     'x_a': 0,
    # })
    # m.set_referece_times([0,0.5,0.6])
    # #m.print_formulas()
    # #m.lambadify()
    # #m.evaluate()
    # res, ref_res = m.simulate(num_iterations=1*3600*24)
    # import pprint
    # pprint.pprint(res)
    # pprint.pprint(ref_res)