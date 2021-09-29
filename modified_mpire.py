#!/usr/bin/env python

import parameters
import tests

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
































__author__ = 'thalita'

# !/usr/bin/env python
# coding: utf-8

import parameters
from modified_mpire import water_flow_module
import pandas as pd


###FUNCTIONS####
def falfa_beta_sz(j):
    if GP.hpipe <= 0.03:
        alfa2 = 0.5
        beta2 = 0.5
    else:
        alfa2 = (SZ.m_sz - 1 - j) / (SZ.m_sz - 1)
        beta2 = j / (SZ.m_sz - 1)

    return alfa2, beta2

def fV_sz(alfa2, beta2, Qfs, Qhc, Qet2, Qpipe, Qinf_sz, teta_i_sz):
    V_sz = (alfa2 * (Qfs - Qhc - Qet2) + beta2 * (Qpipe + Qinf_sz)) / (PZ.Ab * teta_i_sz)

    return V_sz

def fD_x(V_x):
    D_x = EC.lamta1 * V_x
    return D_x

def fPe_usz(V_usz, D_usz):
    if D_usz > 0:
        Pe_usz = V_usz * GP.dz / D_usz
    else:
        Pe_usz = 100

    return Pe_usz

def fPe_sz(V_sz, D_sz):
    D_sz = EC.lamta2 * V_sz
    if D_sz > 0:
        Pe_sz = V_sz * GP.dz / D_sz
    else:
        Pe_sz = 100

    return Pe_sz


##### Transport equations #####
def ftransport(teta_i, teta_iplus1, ci, csoil_i, dc, dc_dz, kads, kdes, D_x, V_x, Rx):
    if teta_i == 0:
        delta_c_i1 = 0
    elif teta_iplus1 == 0:
        delta_c_i1 = 0
    else:
        delta_c_i1 = ((1 / teta_iplus1) * GP.dt * (-teta_i * kads * ci
                                                + SOIL_PLANT.ro * kdes * csoil_i
                                                + teta_i * (D_x * (dc / GP.dz ** 2)
                                                            - V_x * dc_dz)
                                                + Rx))
    return delta_c_i1

def fcsoil(csoil_a, teta, kads, ci, kdes, Rxs):
    csoil_abs = csoil_a + ((teta / SOIL_PLANT.ro) * kads * ci - kdes * csoil_a + Rxs) * GP.dt

    if csoil_abs <= 0:
        csoil = 0
    else:
        csoil = csoil_abs

    return csoil

def fR_dieoff_l(teta_i, muel, ci):
    dieoff_l = -teta_i * muel * ci

    return dieoff_l

def fR_dieoff_s(mues, csoil):
    dieoff_s = -mues * csoil
    return dieoff_s

###ENDFUNCTIONS####


####### **4. Model routine**
def ecoli_module(GP, PZ, USZ, SZ, SOIL_PLANT, EC):
    ##### Ponding Zone
    for time in range(len(GP.rain_inflow) - 1):
        PZ.height_now = PZ.height_after[time]
        if PZ.height_now < 0.001:
            PZ.height_now = 0
        else:
            PZ.height_now = PZ.height_now

        if time == 0:
            PZ.height_before = 0
        else:
            PZ.height_before = PZ.height_after[time - 1]

        if time == 0:
            SZ.pipe_outflow_now = 0
        else:
            SZ.pipe_outflow_now = SZ.pipe_outflow[time]

        ### Reacao na ponding zone ###
        EC.muel = EC.mue1 * EC.sita ** (EC.temperature[time] - 20)
        EC.reaction_rate_pz_now = EC.f_reaction_pz()

        if time < (len(GP.rain_inflow) - 1):
            USZ.theta_after = USZ.theta[time + 1]
            SZ.theta_after = SZ.theta[time + 1]
        else:
            USZ.theta_after = USZ.theta[time]
            SZ.theta_after = SZ.theta[time]

        if PZ.height_now == 0:
            EC.concentration_pz_now = 0
        else:
            EC.concentration_pz_now = PZ.f_concentration(time, GP, EC)

        EC.concentration_pz.append(EC.concentration_pz_now)
        EC.concentration_pz_before = EC.concentration_pz[-1]

        #### SOIL MIX
        ## die-off in USZ, according to s[time]
        EC.mues1 = EC.mue1 * EC.sita ** (EC.temperature[time] - 20)
        EC.muel1 = EC.mue1 * EC.sita ** (EC.temperature[time] - 20)
        EC.muel2 = EC.mue2 * EC.sita ** (EC.temperature[time] - 20)
        EC.mues2 = EC.mue2 * EC.sita ** (EC.temperature[time] - 20)

        # USZ
        EC.concentration_usz_layers = EC.concentration_usz[time].copy()

        # SZ
        if GP.hpipe >= 0.03:
            EC.concentration_sz_layers = EC.concentration_sz[time].copy()  # copiar lista

        if USZ.m_usz != 0:
            EC.concentration_usz_layers_now = []
            EC.concentration_soil_usz_now = []
            EC.reaction_rate_usz_layers = []

            # Predictive step
            ###   USZ   ####
            for usz_layer in range(USZ.m_usz):
                if usz_layer < (USZ.m_usz - 1):
                    EC.concentration_usz_next_layer = EC.concentration_usz_layers[usz_layer + 1]
                else:
                    EC.concentration_usz_next_layer = 0

                EC.concentration_soil_usz_before = EC.concentration_soil_usz[time][usz_layer]
                EC.concentration_soil_usz_layer = EC.f_straining()
                EC.concentration_soil_usz_now.append(EC.f_concentration_soil_usz(time, usz_layer, GP, USZ, SOIL_PLANT))

                USZ.unit_flux_now = USZ.f_unit_flux(time, usz_layer, GP, PZ, SZ)
                USZ.unit_flux.append(USZ.unit_flux_now)

                ## Reacoes ## dieoff_l (parte liquida) + straining + dieoff (parte solida)
                EC.dieoff_liquid_phase = EC.f_dieoff_liquid_phase_usz(time, usz_layer, USZ)
                EC.straining = EC.f_straining()
                EC.dieoff_solid_phase = EC.f_dieoff_solid_phase_usz(usz_layer)
                EC.reaction_rate_usz_now = EC.dieoff_liquid_phase + EC.straining + EC.dieoff_solid_phase
                EC.reaction_rate_usz_layers.append(EC.reaction_rate_usz_now * (1 / USZ.theta_after) * GP.dt)

                ## Coeficiente de difusao ##
                EC.D = EC.f_diffusion_coefficient_usz(USZ)
                EC.kads = EC.kads1
                EC.kdes = EC.kdes1
                EC.peclet = USZ.f_peclet(GP, EC)

                dc, dc_dz = EC.f_delta_concentration_layer_usz(usz_layer, GP.dz, USZ.m_usz)

                ### Transport usz ###
                delta_c_usz = ftransport(USZ.theta[time], USZ.theta_after, EC.concentration_usz_layers[usz_layer], EC.concentration_soil_usz_now[usz_layer], dc, dc_dz, EC.kads1, EC.kdes1,
                                         EC.D, USZ.unit_flux_now, EC.reaction_rate_usz_now)


                if USZ.theta_after > 0:
                    ci1 = EC.concentration_usz_layers[usz_layer] + delta_c_usz
                else:
                    ci1 = 0

                if ci1 < 0.01:
                    ci1 = 0
                else:
                    ci1 = ci1

                EC.concentration_usz_layers_now.append(ci1)

            ### SZ ###

            cj_i1 = []
            csoil_i_sz = []
            Rxj = []

            for j in range(SZ.m_sz):
                cj_sz = EC.concentration_sz_layers[j]
                cjminus1 = EC.concentration_sz_layers[j - 1]
                if j < (SZ.m_sz - 1):
                    cjplus1 = EC.concentration_sz_layers[j + 1]
                else:
                    cjplus1 = 0
                cmminus1 = EC.concentration_usz_layers[USZ.m_usz - 1]

                ## Concentracao no solo - sz ##
                EC.csoil_sz_a = EC.csoil_sz_list[time][j]  # chamar a lista mais ampla

                Rxi_soil_sz = 0

                csoil_sz = fcsoil(EC.csoil_sz_a, SZ.theta[time], EC.kads2, cj_sz, EC.kdes2, Rxi_soil_sz)

                if csoil_sz < 0.00000000000001:
                    csoil_sz = 0

                csoil_i_sz.append(csoil_sz)

                V_sz = []
                alfa2, beta2 = falfa_beta_sz(j)

                Vi_sz = fV_sz(alfa2, beta2, USZ.infiltration_to_sz[time], USZ.capillary_rise[time], USZ.evapotranspiration[time], SZ.pipe_outflow_now, SZ.infiltration_to_surround[time], SZ.theta[time])
                V_sz.append(Vi_sz)

                ## Reacoes ## dieoff_l (parte liquida) + dieoff_s (parte solida)
                Rxi_sz = fR_dieoff_l(SZ.theta[time], EC.muel2, cj_sz) + fR_dieoff_s(EC.mues2, csoil_sz)


                Rxj.append(Rxi_sz * (1 / SZ.theta_after) * GP.dt)

                ## Coeficiente de difusao ##
                Di_sz = fD_x(Vi_sz)
                Pe_sz = fPe_sz(Vi_sz, Di_sz)

                if USZ.m_usz < (GP.n - 1):
                    if Pe_sz <= 2:
                        if j == 0:  # first cell
                            dc = cjplus1 - 2 * cj_sz + cmminus1
                            dc_dz = (cjplus1 - cmminus1) / (2 * GP.dz)

                        elif j == (SZ.m_sz - 1):  # last cell
                            dc = cj_sz - 2 * cj_sz + cjminus1
                            dc_dz = (cj_sz - cjminus1) / GP.dz

                        else:
                            dc = cjplus1 - 2 * cj_sz + cjminus1
                            dc_dz = (cjplus1 - cjminus1) / (2 * GP.dz)

                    else:  # Pesz > 2
                        if j == 0:  # first cell
                            dc = cjplus1 - 2 * cj_sz + cmminus1
                            dc_dz = (cj_sz - cmminus1) / GP.dz

                        elif j == (SZ.m_sz - 1):  # last cell
                            dc = cj_sz - 2 * cj_sz + cjminus1
                            dc_dz = (cj_sz - cjminus1) / GP.dz

                        else:
                            dc = cjplus1 - 2 * cj_sz + cjminus1
                            dc_dz = (cj_sz - cjminus1) / GP.dz

                if USZ.m_usz == (GP.n - 1):
                    dc = cj_sz - 2 * cj_sz + cjminus1
                    dc_dz = (cj_sz - cjminus1) / (2 * GP.dz)

                if USZ.m_usz == 0:
                    if Pe_sz <= 2:
                        if j == 0:  # first cell
                            dc = cjplus1 - 2 * cj_sz + EC.concentration_pz_now
                            dc_dz = (cjplus1 - EC.concentration_pz_now) / (2 * GP.dz)

                        elif j == (SZ.m_sz - 1):  # last cell
                            dc = cj_sz - 2 * cj_sz + cjminus1
                            dc_dz = (cj_sz - cjminus1) / GP.dz

                        else:
                            dc = cjplus1 - 2 * cj_sz + cjminus1
                            dc_dz = (cjplus1 - cjminus1) / (2 * GP.dz)

                    else:  # EC.peclet > 2
                        if j == 0:  # first cell
                            dc = cjplus1 - 2 * cj_sz + EC.concentration_pz_now
                            dc_dz = (cj_sz - EC.concentration_pz_now) / GP.dz

                        elif j == (SZ.m_sz - 1):  # last cell
                            dc = cj_sz - 2 * cj_sz + cjminus1
                            dc_dz = (cj_sz - cjminus1) / GP.dz

                        else:
                            dc = cjplus1 - 2 * cj_sz + cjminus1
                            dc_dz = (cj_sz - cjminus1) / GP.dz


                delta_c_sz = ftransport(SZ.theta[time], SZ.theta_after, cj_sz, csoil_sz, dc, dc_dz, EC.kads2, EC.kdes2,
                                        Di_sz, Vi_sz, Rxi_sz)

                if SZ.theta_after > 0:
                    ci1 = cj_sz + delta_c_sz
                else:
                    ci1 = 0

                if ci1 < 0.01:
                    ci1 = 0
                else:
                    ci1 = ci1

                cj_i1.append(ci1)

        # Corrective step

        if GP.hpipe > 0.03:

            Mstor_ast = sum(EC.concentration_usz[time]) * PZ.Ab * USZ.theta[time] * USZ.height[time] * 1 / USZ.m_usz + sum(EC.concentration_sz[time]) * PZ.Ab * \
                        SZ.theta[time] * SZ.height[time] * 1 / SZ.m_sz
            EC.Mstor_ast_list.append(Mstor_ast)

            Msoil_a = EC.concentration_soil_usz_before * SOIL_PLANT.ro * PZ.Ab * USZ.height[time] * 1 + EC.csoil_sz_a * SOIL_PLANT.ro * PZ.Ab * SZ.height[time] * 1
            Msoil = sum(EC.concentration_soil_usz[time]) * SOIL_PLANT.ro * PZ.Ab * USZ.height[time] * 1 / USZ.m_usz + sum(EC.concentration_soil_usz[time]) * SOIL_PLANT.ro * PZ.Ab * SZ.height[
                time] * 1 / SZ.m_sz
            EC.Msoil_acum.append(Msoil - Msoil_a)

            MRx = - (sum(EC.Rx_usz_list[time]) * PZ.Ab * USZ.theta[time] * USZ.height[time] * 1 / USZ.m_usz + sum(EC.Rx_sz_list[time]) * PZ.Ab *
                     SZ.theta[time] * SZ.height[time] * 1 / SZ.m_sz)

            if time == 0:
                EC.MRx_acum.append(MRx)
            else:
                EC.MRx_acum.append(MRx + EC.MRx_acum[-1])

            M_in = GP.inflow[time] * EC.concentration_inflow[time] * GP.dt * 1
            if time == 0:
                EC.M_in_acum.append(M_in)
            else:
                EC.M_in_acum.append(M_in + EC.M_in_acum[-1])

            M_over = PZ.overflow[time] * EC.concentration_pz[time] * GP.dt * 1
            if time == 0:
                EC.M_over_acum.append(M_over)
            else:
                EC.M_over_acum.append(M_over + EC.M_over_acum[-1])

            Mpipe = SZ.pipe_outflow[time] * EC.concentration_sz[time][SZ.m_sz - 1] * GP.dt * 1
            if time == 0:
                EC.Mpipe_acum.append(Mpipe)
            else:
                EC.Mpipe_acum.append(Mpipe + EC.Mpipe_acum[-1])

            Minf_sz = SZ.infiltration_to_surround[time] * EC.concentration_sz[time][SZ.m_sz - 1] * GP.dt * 1
            if time == 0:
                EC.Minf_sz_acum.append(Minf_sz)
            else:
                EC.Minf_sz_acum.append(Minf_sz + EC.Minf_sz_acum[-1])

            M_et = PZ.evapotranspiration_overall[time] * EC.concentration_usz_layers_now[0] * GP.dt * 1
            if time == 0:
                EC.M_et_acum.append(M_et)
            else:
                EC.M_et_acum.append(M_et + EC.M_et_acum[-1])

            M_pz = PZ.height_after[time] * PZ.Ab * 1 * EC.concentration_pz[time]
            EC.M_pz_list.append(M_pz)

            Mstor_mb = EC.M_in_acum[-1] - EC.M_pz_list[-1] - EC.M_over_acum[-1] - EC.Mpipe_acum[-1] - \
                       EC.Minf_sz_acum[-1] - EC.M_et_acum[-1] - EC.Msoil_acum[-1] - EC.MRx_acum[-1]
            EC.Mstor_mb_list.append(Mstor_mb)

            delta_ast_usz = (Mstor_mb - Mstor_ast) / (PZ.Ab * USZ.theta[time] * USZ.height[time] * 1 / USZ.m_usz)
            EC.concentration_usz_layers_now[1] = EC.concentration_usz_layers_now[1] + delta_ast_usz
            if EC.concentration_usz_layers_now[1] > 0:
                EC.concentration_usz_layers_now[1] = EC.concentration_usz_layers_now[1]
            else:
                EC.concentration_usz_layers_now[1] = 0

        if GP.hpipe == 0:
            Mstor_ast = sum(EC.concentration_usz[time]) * PZ.Ab * USZ.theta[time] * USZ.height[time] * 1 / USZ.m_usz
            EC.Mstor_ast_list.append(Mstor_ast)

            Msoil_a = EC.concentration_soil_usz_before * SOIL_PLANT.ro * PZ.Ab * USZ.height[time] * 1
            Msoil = sum(EC.concentration_soil_usz[time]) * SOIL_PLANT.ro * PZ.Ab * USZ.height[time] * 1 / USZ.m_usz
            EC.Msoil_acum.append(Msoil - Msoil_a)

            MRx = - (sum(EC.Rx_usz_list[time]) * PZ.Ab * USZ.theta[time] * USZ.height[time] * 1 / USZ.m_usz)
            if time == 0:
                EC.MRx_acum.append(MRx)
            else:
                EC.MRx_acum.append(MRx + EC.MRx_acum[-1])

            M_in = GP.inflow[time] * EC.concentration_inflow[time] * GP.dt * 1
            if time == 0:
                EC.M_in_acum.append(M_in)
            else:
                EC.M_in_acum.append(M_in + EC.M_in_acum[-1])

            M_over = PZ.overflow[time] * EC.concentration_pz[time] * GP.dt * 1
            if time == 0:
                EC.M_over_acum.append(M_over)
            else:
                EC.M_over_acum.append(M_over + EC.M_over_acum[-1])

            Mpipe = SZ.pipe_outflow[time] * EC.concentration_usz[time][USZ.m_usz - 1] * GP.dt * 1
            if time == 0:
                EC.Mpipe_acum.append(Mpipe)
            else:
                EC.Mpipe_acum.append(Mpipe + EC.Mpipe_acum[-1])

            Minf_sz = SZ.infiltration_to_surround[time] * EC.concentration_usz[time][USZ.m_usz - 1] * GP.dt * 1
            if time == 0:
                EC.Minf_sz_acum.append(Minf_sz)
            else:
                EC.Minf_sz_acum.append(Minf_sz + EC.Minf_sz_acum[-1])

            M_et = PZ.evapotranspiration_overall[time] * EC.concentration_usz_layers_now[0] * GP.dt * 1
            if time == 0:
                EC.M_et_acum.append(M_et)
            else:
                EC.M_et_acum.append(M_et + EC.M_et_acum[-1])

            M_pz = PZ.height_after[time] * PZ.Ab * 1 * EC.concentration_pz[time]
            EC.M_pz_list.append(M_pz)

            Mstor_mb = EC.M_in_acum[-1] - EC.M_pz_list[-1] - EC.M_over_acum[-1] - EC.Mpipe_acum[-1] - EC.Minf_sz_acum[-1] - EC.M_et_acum[
                -1] - EC.Msoil_acum[-1] - EC.MRx_acum[-1]
            EC.Mstor_mb_list.append(Mstor_mb)

            delta_ast_usz = (Mstor_mb - Mstor_ast) / (PZ.Ab * USZ.theta[time] * USZ.height[time] * 1 / USZ.m_usz)
            EC.concentration_usz_layers_now[1] = EC.concentration_usz_layers_now[1] + delta_ast_usz
            if EC.concentration_usz_layers_now[1] > 0:
                EC.concentration_usz_layers_now[1] = EC.concentration_usz_layers_now[1]
            else:
                EC.concentration_usz_layers_now[1] = 0

        EC.concentration_soil_usz.append(EC.concentration_soil_usz_now)
        EC.concentration_usz.append(EC.concentration_usz_layers_now)
        EC.Rx_usz_list.append(EC.reaction_rate_usz_layers)

        if GP.hpipe >= 0.03:
            EC.csoil_sz_list.append(csoil_i_sz)
            EC.concentration_sz.append(cj_i1)
            EC.Rx_sz_list.append(Rxj)

    # **5. Transforming in dataframe **
    ## E Coli
    data_usz = pd.DataFrame(EC.concentration_usz)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz' + str(a)
        column_name.append(name)
        a = a + 1
    data_usz.set_axis(column_name, axis='columns', inplace=True)

    data_sz = pd.DataFrame(EC.concentration_sz)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz' + str(a)
        column_name.append(name)
        a = a + 1
    data_sz.set_axis(column_name, axis='columns', inplace=True)


    data_soil_usz = pd.DataFrame(EC.concentration_soil_usz)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz_soil' + str(a)
        column_name.append(name)
        a = a + 1
    data_soil_usz.set_axis(column_name, axis='columns', inplace=True)

    data_soil_sz = pd.DataFrame(EC.csoil_sz_list)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz_soil' + str(a)
        column_name.append(name)
        a = a + 1
    data_soil_sz.set_axis(column_name, axis='columns', inplace=True)

    data_rx_usz = pd.DataFrame(EC.Rx_usz_list)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_usz.set_axis(column_name, axis='columns', inplace=True)

    data_rx_sz = pd.DataFrame(EC.Rx_sz_list)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_sz.set_axis(column_name, axis='columns', inplace=True)

    frames = [data_usz, data_sz, data_soil_usz, data_soil_sz, data_rx_usz, data_rx_sz]
    data_EColi = pd.concat((frames), axis=1)

    # ----- teste
    cpz_list_calib = EC.concentration_pz[:]
    cpz_list_calib.append(0)
    data_EColi['cpz'] = cpz_list_calib
    EC.concentration_inflow[time] = EC.concentration_inflow[:]
    EC.concentration_inflow[time].append(0)
    c_in = EC.concentration_inflow[time][:len(GP.rain_inflow)]
    data_EColi['c_in'] = c_in

    data_EColi['t'] = range(len(GP.rain_inflow))

    Qorif_list = SZ.pipe_outflow[:len(GP.rain_inflow)]

    data_EColi['Qorif'] = Qorif_list

    if GP.hpipe > 0:
        data_EColi['Morif'] = data_EColi['Qorif'] * data_EColi['sz3'] * 1000  # m3/s
    else:
        data_EColi['Morif'] = data_EColi['Qorif'] * data_EColi['usz10'] * 1000  # m3/s

    return data_EColi


if __name__ == '__main__':
    SETUP_FILE = "parameters.ini"
    WATER_FLOW_INPUT_FILE = "input_files/ecoli/test_event/water_inflow.csv"
    WATER_QUALITY_INPUT_FILE = "input_files/ecoli/test_event/concentration_inflow.csv"
    GP = parameters.GeneralParameters(SETUP_FILE, WATER_FLOW_INPUT_FILE)
    USZ = parameters.UnsaturatedZone(SETUP_FILE, GP.L, GP.hpipe, GP.dz)
    PZ = parameters.PondingZone(SETUP_FILE)
    SZ = parameters.SaturatedZone(SETUP_FILE, GP, USZ)
    SOIL_PLANT = parameters.SoilPlant(SETUP_FILE, USZ.nusz_ini)
    EC = parameters.Ecoli(USZ.nusz_ini, USZ.m_usz, SZ.m_sz, SETUP_FILE, WATER_QUALITY_INPUT_FILE)
    water_flow_module(GP, USZ, PZ, SZ)
    data_EColi = ecoli_module(GP, PZ, USZ, SZ, SOIL_PLANT, EC)

    #data_EColi.to_csv('EColi_results_calib.csv', index=False, sep=';', decimal=',')
    print('Done!')


from tests import water_quality_comparison_test

errors = water_quality_comparison_test("test_files/ecoli_tested_thalita_parameters.csv", data_EColi)
print(errors)


