import numpy as np
import pandas as pd


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



