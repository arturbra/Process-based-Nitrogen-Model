#!/usr/bin/env python

def water_flow_module(GP, USZ, PZ, SZ, RTC=False):
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

        if str(type(RTC)) == "<class 'parameters.RealTimeControl'>":
            if RTC.strategy == 1:
                SZ.pipe_outflow.append(RTC.control(time, GP, USZ, SZ))
            if RTC.strategy == 2:
                SZ.pipe_outflow.append(RTC.control(time, GP, USZ, SZ))

        else:
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


def nitrogen_module(GP, PZ, USZ, SZ, SOIL_PLANT, NH4, NO3, O2, DOC):
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


def ecoli_module(GP, PZ, USZ, SZ, SOIL_PLANT, EC, log_transformation=True):
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
                EC.concentration_soil_usz_layer = EC.f_concentration_soil_usz(time, usz_layer, GP, USZ, SOIL_PLANT)
                EC.concentration_soil_usz_now.append(EC.concentration_soil_usz_layer)

                USZ.unit_flux_now = USZ.f_unit_flux(time, usz_layer, GP, PZ, SZ)
                USZ.unit_flux.append(USZ.unit_flux_now)

                # Reactions
                EC.dieoff_liquid_phase = EC.f_dieoff_liquid_phase_usz(time, usz_layer, USZ)
                EC.straining = EC.f_straining()
                EC.dieoff_solid_phase = EC.f_dieoff_solid_phase_usz(usz_layer)
                EC.reaction_rate_usz_now = EC.dieoff_liquid_phase + EC.straining + EC.dieoff_solid_phase
                EC.reaction_rate_usz_layers.append(EC.reaction_rate_usz_now * (1 / USZ.theta_after) * GP.dt)

                EC.D = EC.f_diffusion_coefficient_usz(USZ)
                EC.kads = EC.kads1
                EC.kdes = EC.kdes1
                EC.peclet = USZ.f_peclet(GP, EC)
                EC.delta_concentration = EC.f_delta_concentration_usz(time, usz_layer, GP, USZ, SOIL_PLANT, threshold=0.01)
                EC.concentration_usz_layers_now.append(EC.delta_concentration)

            ### SZ ###

            EC.concentration_sz_layers_now = []
            EC.concentration_soil_sz_now = []
            EC.reaction_rate_sz_layers = []

            for sz_layer in range(SZ.m_sz):
                if sz_layer < (SZ.m_sz - 1):
                    EC.concentration_sz_next_layer = EC.concentration_sz_layers[sz_layer + 1]
                else:
                    EC.concentration_sz_next_layer = 0
                ## Concentracao no solo - sz ##
                EC.concentration_soil_sz_before = EC.concentration_soil_sz[time][sz_layer]
                EC.concentration_soil_sz_layer = EC.f_concentration_soil_sz(time, sz_layer, GP, SZ, SOIL_PLANT)
                EC.concentration_soil_sz_now.append(EC.concentration_soil_sz_layer)

                SZ.unit_flux_now = SZ.f_unit_flux(time, sz_layer, PZ, USZ)
                SZ.unit_flux.append(SZ.unit_flux_now)

                ## Reacoes ## dieoff_l (parte liquida) + dieoff_s (parte solida)
                EC.dieoff_liquid_phase = EC.f_dieoff_liquid_phase_sz(time, sz_layer, SZ)
                EC.straining = EC.f_straining()
                EC.dieoff_solid_phase = EC.f_dieoff_solid_phase_sz(sz_layer)
                EC.reaction_rate_sz_now = EC.dieoff_liquid_phase + EC.straining + EC.dieoff_solid_phase
                EC.reaction_rate_sz_layers.append(EC.reaction_rate_sz_now * (1 / SZ.theta_after) * GP.dt)

                ## Coeficiente de difusao ##
                EC.D = EC.f_diffusion_coefficient_sz(SZ)
                EC.peclet = SZ.f_peclet(GP, EC)
                EC.delta_concentration = EC.f_delta_concentration_sz(time, sz_layer, GP, USZ, SZ, SOIL_PLANT, threshold=0.01)
                EC.concentration_sz_layers_now.append(EC.delta_concentration)


        # Corrective step
        EC.mass_balance.append(EC.f_mass_balance_stormwater(time, USZ, SZ, PZ, GP, SOIL_PLANT))
        EC.mass_stormwater_now = EC.f_mass_stormwater(time, GP, PZ, USZ, SZ)
        EC.f_mass_balance_concentration(time, PZ, USZ)

        # adding layers of USZ in time
        EC.concentration_soil_usz.append(EC.concentration_soil_usz_now)
        EC.concentration_usz.append(EC.concentration_usz_layers_now)
        EC.reaction_rate_usz.append(EC.reaction_rate_usz_layers)

        # adding layers of SZ in time
        if GP.hpipe >= 0.03:
            EC.concentration_soil_sz.append(EC.concentration_soil_sz_now)
            EC.concentration_sz.append(EC.concentration_sz_layers_now)
            EC.reaction_rate_sz.append(EC.reaction_rate_sz_layers)

    if log_transformation:
        EC.concentration_inflow_log = EC.f_log_transformation(EC.concentration_inflow)
        EC.concentration_pz_log = EC.f_log_transformation(EC.concentration_pz)
        EC.concentration_soil_usz_log = EC.f_log_transformation(EC.concentration_soil_usz)
        EC.concentration_usz_log = EC.f_log_transformation(EC.concentration_usz)
        EC.reaction_rate_usz_log = EC.f_log_transformation(EC.reaction_rate_usz)
        EC.concentration_soil_sz_log = EC.f_log_transformation(EC.concentration_soil_sz)
        EC.concentration_sz_log = EC.f_log_transformation(EC.concentration_sz)
        EC.reaction_rate_sz_log = EC.f_log_transformation(EC.reaction_rate_sz)

    # **5. Transforming in dataframe **
    data_ecoli = EC.water_quality_results(GP, USZ, SZ, log_transformation)

    return data_ecoli
