#!/usr/bin/env python

import parameters
import results_tests
import datetime


def water_flow_module(GP, USZ, PZ, SZ):
    for time in range(len(GP.rain_inflow)):
        # PZ
        PZ.overflow.append(PZ.f_weir_overflow(time, GP))
        PZ.height.append(PZ.f_height(time, GP))
        PZ.infiltration_to_surround.append(PZ.f_infiltration_to_surround(time, USZ))
        PZ.infiltration_to_filter_material.append(PZ.f_infiltration_to_filter_material(time, USZ, GP))
        PZ.height_after_now = PZ.f_height_after(time, GP)

        # USZ
        USZ.wilting_point_estimated.append(USZ.f_wilting_point_moisture_estimation(time, PZ, GP))
        USZ.capillary_rise.append(USZ.f_capillary_rise(time, GP))
        USZ.infiltration_to_sz.append(USZ.f_infiltration_to_sz(time, PZ, GP))
        USZ.wilting_point_estimated_after.append(USZ.f_wilting_point_moisture_next_interaction_estimation(time, SZ))

        # EVT
        PZ.evapotranspiration_overall.append(PZ.f_evapotranspiration_overall(time, USZ, GP))
        PZ.evapotranspiration.append(PZ.f_evapotranspiration(time, USZ, SZ))
        USZ.evapotranspiration.append(USZ.f_evapotranspiration(time, PZ))

        # SZ
        SZ.height_estimated.append(SZ.f_height_estimation(time, GP, USZ))
        SZ.infiltration_to_surround.append(SZ.f_infiltration_to_surround(time, PZ, USZ))
        SZ.pipe_outflow.append(SZ.f_underdrain_flow(time, GP, USZ))
        SZ.height_now = SZ.f_height(time, GP, USZ)
        USZ.height_now = USZ.f_height(GP, SZ)

        # Porosity#
        SZ.porosity_now = SZ.f_porosity(GP)
        USZ.porosity_now = USZ.f_porosity(GP, SZ)

        if time != 0:
            USZ.height_before = USZ.height[time - 1]
            USZ.porosity_before = USZ.porosity[time - 1]
            USZ.wilting_point_moisture_before = USZ.wilting_point_moisture[time - 1]

        USZ.wilting_point_moisture_now = USZ.f_wilting_point_moisture(time, GP, PZ)

        USZ.theta.append(USZ.wilting_point_moisture_now * USZ.porosity_now)
        SZ.theta.append(SZ.porosity_now)
        PZ.height_after.append(PZ.height_after_now)
        SZ.height.append(SZ.height_now)
        USZ.height.append(USZ.height_now)
        USZ.wilting_point_moisture.append(USZ.wilting_point_moisture_now)
        USZ.porosity.append(USZ.porosity_now)
        SZ.porosity.append(SZ.porosity_now)


# **4. Model routine**
def water_quality_module(GP, USZ, PZ, SZ, SOIL_PLANT, NH4, NO3, O2, DOC):
    for time in range(len(GP.rain_inflow) - 1):
        # PZ
        PZ.height_now = PZ.height_after[time]
        if PZ.height_now < 0.001:
            PZ.height_now = 0

        if time != 0:
            PZ.height_before = PZ.height_after[time - 1]
            SZ.pipe_outflow_now = SZ.pipe_outflow[time]

        if time < (len(GP.rain_inflow) - 1):
            USZ.theta_after = USZ.theta[time + 1]
            SZ.theta_after = SZ.theta[time + 1]
        else:
            USZ.theta_after = USZ.theta[time]
            SZ.theta_after = SZ.theta[time]

        O2.reaction_rate_pz_now = O2.f_reaction_pz()
        NH4.reaction_rate_pz_now = NH4.f_reaction_pz()
        NO3.reaction_rate_pz_now = NO3.f_reaction_pz()
        DOC.reaction_rate_pz_now = DOC.f_reaction_pz()

        if PZ.height_now == 0:
            O2.concentration_pz_now = 0
            NH4.concentration_pz_now = 0
            NO3.concentration_pz_now = 0
            DOC.concentration_pz_now = 0

        else:
            O2.concentration_pz_now = PZ.f_concentration(time, GP, O2)
            NH4.concentration_pz_now = PZ.f_concentration(time, GP, NH4)
            NO3.concentration_pz_now = PZ.f_concentration(time, GP, NO3)
            DOC.concentration_pz_now = PZ.f_concentration(time, GP, DOC)

        O2.concentration_pz.append(O2.concentration_pz_now)
        NH4.concentration_pz.append(NH4.concentration_pz_now)
        NO3.concentration_pz.append(NO3.concentration_pz_now)
        DOC.concentration_pz.append(DOC.concentration_pz_now)

        O2.concentration_pz_before = O2.concentration_pz[-1]
        NH4.concentration_pz_before = NH4.concentration_pz[-1]
        NO3.concentration_pz_before = NO3.concentration_pz[-1]
        DOC.concentration_pz_before = DOC.concentration_pz[-1]

        # USZ
        O2.concentration_usz_layers = O2.concentration_usz[time].copy()
        NH4.concentration_usz_layers = NH4.concentration_usz[time].copy()
        NO3.concentration_usz_layers = NO3.concentration_usz[time].copy()
        DOC.concentration_usz_layers = DOC.concentration_usz[time].copy()

        # SZ
        if GP.hpipe > 0:
            O2.concentration_sz_layers = O2.concentration_sz[time].copy()
            NH4.concentration_sz_layers = NH4.concentration_sz[time].copy()
            NO3.concentration_sz_layers = NO3.concentration_sz[time].copy()
            DOC.concentration_sz_layers = DOC.concentration_sz[time].copy()

        # USZ
        if USZ.m_usz != 0:

            O2.concentration_usz_layers_now = []
            NH4.concentration_usz_layers_now = []
            NO3.concentration_usz_layers_now = []
            DOC.concentration_usz_layers_now = []

            O2.reaction_rate_usz_layers = []
            NH4.reaction_rate_usz_layers = []
            NO3.reaction_rate_usz_layers = []
            DOC.reaction_rate_usz_layers = []

            NH4.concentration_soil_usz_now = []
            DOC.concentration_soil_usz_now = []

            # Predictive step
            for usz_layer in range(USZ.m_usz):
                if usz_layer < (USZ.m_usz - 1):
                    O2.concentration_usz_next_layer = O2.concentration_usz_layers[usz_layer + 1]
                    NH4.concentration_usz_next_layer = NH4.concentration_usz_layers[usz_layer + 1]
                    NO3.concentration_usz_next_layer = NO3.concentration_usz_layers[usz_layer + 1]
                    DOC.concentration_usz_next_layer = DOC.concentration_usz_layers[usz_layer + 1]
                else:
                    O2.concentration_usz_next_layer = 0
                    NH4.concentration_usz_next_layer = 0
                    NO3.concentration_usz_next_layer = 0
                    DOC.concentration_usz_next_layer = 0

                O2.concentration_soil_usz_layer = O2.f_concentration_soil()
                NO3.concentration_soil_usz_layer = NO3.f_concentration_soil()

                NH4.concentration_soil_usz_layer = NH4.concentration_soil_usz[time][usz_layer]
                NH4.concentration_soil_usz_now.append(
                    NH4.f_concentration_soil(NH4.concentration_soil_usz_layer, USZ.theta[time], NH4.kads,
                                             NH4.concentration_usz_layers[usz_layer], NH4.kdes, NH4.k_mb, SOIL_PLANT.ro,
                                             GP.dt))

                DOC.concentration_soil_usz_before = DOC.concentration_soil_usz[time][usz_layer]
                DOC.concentration_soil_usz_layer = DOC.f_concentration_soil(DOC.concentration_soil_usz_before,
                                                                            USZ.theta[time], DOC.kads,
                                                                            DOC.concentration_usz_layers[usz_layer],
                                                                            DOC.kdes, DOC.k_mb, SOIL_PLANT.ro, GP.dt)
                DOC.concentration_soil_usz_now.append(DOC.concentration_soil_usz_layer)

                USZ.unit_flux_now = USZ.f_unit_flux(time, usz_layer, GP, PZ, SZ)
                USZ.unit_flux.append(USZ.unit_flux_now)

                O2.reaction_rate_usz_now = O2.f_reaction_usz(usz_layer, GP, NH4) + O2.f_plant_uptake_usz(time,
                                                                                                         usz_layer, USZ,
                                                                                                         SOIL_PLANT)
                NH4.reaction_rate_usz_now = NH4.f_reaction_usz(usz_layer, GP) + NH4.f_plant_uptake_usz(time, usz_layer,
                                                                                                       USZ, SOIL_PLANT)
                NO3.reaction_rate_usz_now = NO3.f_reaction_usz(usz_layer, GP, NH4) + NO3.f_plant_uptake_usz(time,
                                                                                                            usz_layer,
                                                                                                            USZ,
                                                                                                            SOIL_PLANT)
                DOC.reaction_rate_usz_now = DOC.f_reaction_usz(usz_layer)

                O2.reaction_rate_usz_layers.append(O2.f_reaction_conversion_usz(GP, USZ))
                NH4.reaction_rate_usz_layers.append(NH4.f_reaction_conversion_usz(GP, USZ))
                NO3.reaction_rate_usz_layers.append(NO3.f_reaction_conversion_usz(GP, USZ))
                DOC.reaction_rate_usz_layers.append(DOC.f_reaction_conversion_usz(GP, USZ))

                O2.peclet = USZ.f_peclet(GP, O2)
                O2.delta_concentration = O2.f_delta_concentration_usz(time, usz_layer, GP, USZ, SOIL_PLANT)
                O2.concentration_usz_layers_now.append(O2.delta_concentration)

                NH4.peclet = USZ.f_peclet(GP, NH4)
                NH4.delta_concentration = NH4.f_delta_concentration_usz(time, usz_layer, GP, USZ, SOIL_PLANT)
                NH4.concentration_usz_layers_now.append(NH4.delta_concentration)

                NO3.peclet = USZ.f_peclet(GP, NO3)
                NO3.delta_concentration = NO3.f_delta_concentration_usz(time, usz_layer, GP, USZ, SOIL_PLANT)
                NO3.concentration_usz_layers_now.append(NO3.delta_concentration)

                DOC.peclet = USZ.f_peclet(GP, DOC)
                DOC.delta_concentration = DOC.f_delta_concentration_usz(time, usz_layer, GP, USZ, SOIL_PLANT)
                DOC.concentration_usz_layers_now.append(DOC.delta_concentration)

        #####   SZ   #####

        O2.concentration_sz_layers_now = []
        NH4.concentration_sz_layers_now = []
        NO3.concentration_sz_layers_now = []
        DOC.concentration_sz_layers_now = []

        O2.reaction_rate_sz_layers = []
        NH4.reaction_rate_sz_layers = []
        NO3.reaction_rate_sz_layers = []
        DOC.reaction_rate_sz_layers = []

        NH4.concentration_soil_sz_now = []
        DOC.concentration_soil_sz_now = []

        for sz_layer in range(SZ.m_sz):
            if sz_layer < (SZ.m_sz - 1):
                O2.concentration_sz_next_layer = O2.concentration_sz_layers[sz_layer + 1]
                NH4.concentration_sz_next_layer = NH4.concentration_sz_layers[sz_layer + 1]
                NO3.concentration_sz_next_layer = NO3.concentration_sz_layers[sz_layer + 1]
                DOC.concentration_sz_next_layer = DOC.concentration_sz_layers[sz_layer + 1]

            else:
                O2.concentration_sz_next_layer = 0
                NH4.concentration_sz_next_layer = 0
                NO3.concentration_sz_next_layer = 0
                DOC.concentration_sz_next_layer = 0

            O2.concentration_soil_sz_layer = O2.f_concentration_soil()
            NO3.concentration_soil_sz_layer = NO3.f_concentration_soil()

            # See the comments at USZ f_concentration_soil code
            NH4.concentration_soil_sz_before = NH4.concentration_soil_sz[time][sz_layer]
            NH4.concentration_soil_sz_layer = NH4.f_concentration_soil(NH4.concentration_soil_sz_before,
                                                                       SZ.theta[time], NH4.kads2,
                                                                       NH4.concentration_sz_layers[sz_layer], NH4.kdes2,
                                                                       NH4.k_mb, SOIL_PLANT.ro, GP.dt)
            NH4.concentration_soil_sz_now.append(NH4.concentration_soil_sz_layer)

            DOC.concentration_soil_sz_before = DOC.concentration_soil_sz[time][sz_layer]
            DOC.concentration_soil_sz_layer = DOC.f_concentration_soil(DOC.concentration_soil_sz_before,
                                                                       SZ.theta[time], DOC.kads2,
                                                                       DOC.concentration_sz_layers[sz_layer], DOC.kdes2,
                                                                       DOC.k_mb, SOIL_PLANT.ro, GP.dt)
            DOC.concentration_soil_sz_now.append(DOC.concentration_soil_sz_layer)

            SZ.unit_flux_now = SZ.f_unit_flux(time, sz_layer, PZ, USZ)
            SZ.unit_flux.append(SZ.unit_flux_now)

            O2.reaction_rate_sz_now = O2.f_reaction_sz(sz_layer, GP, NH4) + O2.f_plant_uptake_sz(time, sz_layer, SZ,
                                                                                                 SOIL_PLANT)
            NH4.reaction_rate_sz_now = NH4.f_reaction_sz() + NH4.f_plant_uptake_sz(time, sz_layer, SZ, SOIL_PLANT)
            NO3.reaction_rate_sz_now = NO3.f_reaction_sz(sz_layer, GP, O2, DOC) + NO3.f_plant_uptake_sz(time, sz_layer,
                                                                                                        SZ, SOIL_PLANT)
            DOC.reaction_rate_sz_now = DOC.f_reaction_sz(sz_layer)

            O2.reaction_rate_sz_layers.append(O2.f_reaction_conversion_sz(GP, SZ))
            NH4.reaction_rate_sz_layers.append(NH4.f_reaction_conversion_sz(GP, SZ))
            NO3.reaction_rate_sz_layers.append(NO3.f_reaction_conversion_sz(GP, SZ))
            DOC.reaction_rate_sz_layers.append(DOC.f_reaction_conversion_sz(GP, SZ))

            # Oxygen
            O2.peclet = SZ.f_peclet(GP, O2)
            O2.delta_concentration = O2.f_delta_concentration_sz(time, sz_layer, GP, USZ, SZ, SOIL_PLANT)
            O2.concentration_sz_layers_now.append(O2.delta_concentration)

            # Ammonia
            NH4.peclet = SZ.f_peclet(GP, NH4)
            NH4.delta_concentration = NH4.f_delta_concentration_sz(time, sz_layer, GP, USZ, SZ, SOIL_PLANT)
            NH4.concentration_sz_layers_now.append(NH4.delta_concentration)

            # Nitrate
            NO3.peclet = SZ.f_peclet(GP, NO3)
            NO3.delta_concentration = NO3.f_delta_concentration_sz(time, sz_layer, GP, USZ, SZ, SOIL_PLANT)
            NO3.concentration_sz_layers_now.append(NO3.delta_concentration)

            # DOC
            DOC.peclet = SZ.f_peclet(GP, DOC)
            DOC.delta_concentration = DOC.f_delta_concentration_sz(time, sz_layer, GP, USZ, SZ, SOIL_PLANT)
            DOC.concentration_sz_layers_now.append(DOC.delta_concentration)

        # Corrective step
        O2.mass_balance.append(O2.f_mass_balance_stormwater(time, USZ, SZ, PZ, GP))
        O2.mass_stormwater_now = O2.f_mass_stormwater(time, GP, PZ, USZ, SZ)
        O2.f_mass_balance_concentration(time, PZ, USZ)

        NH4.mass_balance.append(NH4.f_mass_balance_stormwater(time, USZ, SZ, PZ, GP, SOIL_PLANT))
        NH4.mass_stormwater_now = NH4.f_mass_stormwater(time, GP, PZ, USZ, SZ)
        NH4.f_mass_balance_concentration(time, PZ, USZ)

        NO3.mass_balance.append(NO3.f_mass_balance_stormwater(time, USZ, SZ, PZ, GP))
        NO3.mass_stormwater_now = NO3.f_mass_stormwater(time, GP, PZ, USZ, SZ)
        NO3.f_mass_balance_concentration(time, PZ, USZ)

        DOC.mass_balance.append(DOC.f_mass_balance_stormwater(time, USZ, SZ, PZ, GP, SOIL_PLANT))
        DOC.mass_stormwater_now = DOC.f_mass_stormwater(time, GP, PZ, USZ, SZ)
        DOC.f_mass_balance_concentration(time, PZ, USZ)

        # adding layers of USZ in time
        O2.concentration_usz.append(O2.concentration_usz_layers_now)
        O2.reaction_rate_usz.append(O2.reaction_rate_usz_layers)

        NH4.concentration_usz.append(NH4.concentration_usz_layers_now)
        NH4.concentration_soil_usz.append(NH4.concentration_soil_usz_now)
        NH4.reaction_rate_usz.append(NH4.reaction_rate_usz_layers)

        NO3.concentration_usz.append(NO3.concentration_usz_layers_now)
        NO3.reaction_rate_usz.append(NO3.reaction_rate_usz_layers)

        DOC.concentration_usz.append(DOC.concentration_usz_layers_now)
        DOC.concentration_soil_usz.append(DOC.concentration_soil_usz_now)
        DOC.reaction_rate_usz.append(DOC.reaction_rate_usz_layers)

        # adding layers of SZ in time
        if GP.hpipe > 0:
            O2.concentration_sz.append(O2.concentration_sz_layers_now)
            O2.reaction_rate_sz.append(O2.reaction_rate_sz_layers)

            NH4.concentration_sz.append(NH4.concentration_sz_layers_now)
            NH4.concentration_soil_sz.append(NH4.concentration_soil_sz_now)
            NH4.reaction_rate_sz.append(NH4.reaction_rate_sz_layers)

            NO3.concentration_sz.append(NO3.concentration_sz_layers_now)
            NO3.reaction_rate_sz.append(NO3.reaction_rate_sz_layers)

            DOC.concentration_sz.append(DOC.concentration_sz_layers_now)
            DOC.concentration_soil_sz.append(DOC.concentration_soil_sz_now)
            DOC.reaction_rate_sz.append(DOC.reaction_rate_sz_layers)

    # **5. Transforming in dataframe **
    # Oxygen
    data_o2 = O2.water_quality_results(USZ, SZ)
    data_nh4 = NH4.water_quality_results(USZ, SZ)
    data_no3 = NO3.water_quality_results(USZ, SZ)
    data_doc = DOC.water_quality_results(USZ, SZ)

    return data_o2, data_nh4, data_no3, data_doc


def run(interaction):
    SETUP_FILE = "parameters.ini"
    WATER_FLOW_INPUT_FILE = "water_inflow.csv"
    WATER_QUALITY_INPUT_FILE = "concentration_inflow.csv"
    GP = parameters.GeneralParameters(SETUP_FILE, WATER_FLOW_INPUT_FILE)
    GP.hpipe = hpipes[interaction]
    USZ = parameters.UnsaturatedZone(SETUP_FILE, GP.L, GP.hpipe, GP.dz)
    PZ = parameters.PondingZone(SETUP_FILE)
    SZ = parameters.SaturatedZone(SETUP_FILE, GP, USZ)
    SOIL_PLANT = parameters.SoilPlant(SETUP_FILE, USZ.nusz_ini)
    NH4 = parameters.Ammonia(USZ.m_usz, SZ.m_sz, SETUP_FILE, WATER_QUALITY_INPUT_FILE)
    NO3 = parameters.Nitrate(USZ.m_usz, SZ.m_sz, SETUP_FILE, WATER_QUALITY_INPUT_FILE)
    O2 = parameters.Oxygen(USZ.m_usz, SZ.m_sz, SETUP_FILE, WATER_QUALITY_INPUT_FILE)
    DOC = parameters.DissolvedOrganicCarbon(USZ.m_usz, SZ.m_sz, SETUP_FILE, WATER_QUALITY_INPUT_FILE)
    water_flow_module(GP, USZ, PZ, SZ)
    print("hpipe:", GP.hpipe)
    data_o2, data_nh4, data_no3, data_doc = water_quality_module(GP, USZ, PZ, SZ, SOIL_PLANT, NH4, NO3, O2, DOC)

    wf_path = results_paths[interaction] + "water_flow_results.csv"
    nh4_path = results_paths[interaction] + "results_Kin_pf_nh4_2.csv"
    o2_path = results_paths[interaction] + "results_Kin_pf_o2.csv"
    no3_path = results_paths[interaction] + "results_Kin_pf_no3.csv"
    doc_path = results_paths[interaction] + "results_Kin_pf_doc.csv"

    wf_test = results_tests.water_flow_comparison_test(wf_path, GP, PZ, USZ, SZ)
    nh4_test = results_tests.water_quality_comparison_test(nh4_path, data_nh4)
    o2_test = results_tests.water_quality_comparison_test(o2_path, data_o2)
    no3_test = results_tests.water_quality_comparison_test(no3_path, data_no3)
    doc_test = results_tests.water_quality_comparison_test(doc_path, data_doc)

    if len(wf_test) == 0 and len(nh4_test) == 0 and len(o2_test) == 0 and len(no3_test) == 0 and len(doc_test) == 0:
        print("Passed at this test. The results seems to be the same as the original code.")
    else:
        print("Don't panic. You can always Rollback")
        print(len(wf_test), len(nh4_test), len(o2_test), len(no3_test), len(doc_test))


if __name__ == '__main__':
    hpipes = [0, 0.1, 0.2, 0.3, 0.4]
    results_paths = ["results/results_00/", "results/results_10/",
                     "results/results_20/", "results/results_30/",
                     "results/results_40/"]
    start = datetime.datetime.now()

    for test in range(len(hpipes)):
        run(test)

    # data_nh4.to_csv('results_Kin_pf_nh4_2.csv', index = False)
    # data_o2.to_csv('results_Kin_pf_o2.csv', index = False)
    # data_no3.to_csv('results_Kin_pf_no3.csv', index = False)
    # data_doc.to_csv('results_Kin_pf_doc.csv', index = False)

    end = datetime.datetime.now()
    print('Elapsed time: ', end - start)
    print('Done!')
