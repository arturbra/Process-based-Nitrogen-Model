import random

import numpy as np
import pandas as pd
import configparser

def water_flow_comparison_test(original_water_flow_results_file, GP, PZ, USZ, SZ):
    water_flow_results = pd.read_csv(original_water_flow_results_file)
    water_flow_results = water_flow_results.round(decimals=4)
    df_dict = {"t": range(0, len(GP.rain_inflow)),
               "Qin": GP.inflow,
               "Qet": PZ.evapotranspiration_overall,
               "hpEND": PZ.height_after,
               "Qpf": PZ.infiltration_to_filter_material,
               "Qover": PZ.overflow,
               "Qfs": USZ.infiltration_to_sz,
               "Qet1": PZ.evapotranspiration,
               "Qhc": USZ.capillary_rise,
               "Qpipe": SZ.pipe_outflow,
               "Qet2": USZ.evapotranspiration,
               "teta_usz": USZ.theta,
               "teta_sz": SZ.theta,
               "Qrain": GP.rain_inflow,
               "Qinfp": PZ.infiltration_to_surround,
               "Qinfsz": SZ.infiltration_to_surround,
               "hp": PZ.height,
               "s": USZ.wilting_point_moisture,
               "husz": USZ.height,
               "hsz": SZ.height,
               "nsz": SZ.porosity,
               "nusz": USZ.porosity,
               "hszEST": SZ.height_estimated
               }

    df = pd.DataFrame(df_dict)
    df = df.round(decimals=4)

    errors = []
    for c in water_flow_results.columns:
        comparison = np.where(df[c] == water_flow_results[c], True, False)
        if False in comparison:
            errors.append(c)
    return errors

def water_quality_comparison_test(file, nutrient):
    results = pd.read_csv(file)
    results = results.round(decimals=4)
    nutrient = nutrient.round(decimals=4)
    errors = []
    for c in results.columns:
        comparison = np.where(nutrient[c] == results[c], True, False)
        if False in comparison:
            errors.append(c)
    return errors

def read_variables_ini_file(ini_file):
    setup = configparser.ConfigParser()
    setup.read(ini_file)
    ini_keys = list(setup.keys())[1:]
    D = [dict(setup.items(k)) for k in ini_keys]
    return ini_keys, D

def change_variables_at_ini(ini_file, output_file_name):
    ini_groups, ini_variables = read_variables_ini_file(ini_file)
    banned_list = ['flagp', 'flagsz', 'dt', 'psz', 'show_summary', 'obs_file_ecoli']
    for k in range(len(ini_variables)):
        keys_variables = list(ini_variables[k].keys())
        for l in keys_variables:
            if l not in banned_list:
                ini_variables[k][l] = float(ini_variables[k][l]) * random.uniform(0.8, 1.2)

    setup = configparser.ConfigParser()
    for i in range(len(ini_groups)):
        setup[ini_groups[i]] = ini_variables[i]

    with open(output_file_name, 'w') as setupfile:
        setup.write(setupfile)

def generate_ini_to_test(ini_quanti, ini_quali):
    for i in range(5):
        change_variables_at_ini(ini_quanti, f'parameters_quanti_{i}.ini')
        change_variables_at_ini(ini_quali, f'parameters_quali_{i}.ini')


