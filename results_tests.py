import numpy as np
import pandas as pd

def water_flow_comparison_test(original_water_flow_results_file, WFR):
    water_flow_results = pd.read_csv(original_water_flow_results_file)
    water_flow_results = water_flow_results.round(decimals=4)
    water_flow_att = [att for att in dir(WFR) if not att.startswith("__")]
    a = [getattr(WFR, att) for att in water_flow_att]
    b = pd.DataFrame(a)
    b = b.transpose()
    b = b.round(decimals=4)
    b.columns = water_flow_att
    errors = []
    for c in water_flow_results.columns:
        name_t = "t" + c
        comparison = np.where(b[name_t] == water_flow_results[c], True, False)
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


