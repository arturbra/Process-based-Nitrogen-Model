#!/usr/bin/env python

import parameters
import results_tests
import datetime


def water_flow_module(WFR, GP, USZ, PZ, SZ):
    for time in range(len(WFR.tQrain)):
        # PZ
        PZ.overflow.append(PZ.f_weir_overflow(PZ.height_after[time], GP.dt, GP.inflow[time], GP.rain_inflow[time]))
        PZ.height.append(PZ.f_height(time, GP))
        PZ.infiltration_to_surrounding.append(PZ.f_infiltration_to_surrounding(USZ.Kf, USZ.A, PZ.height[time]))
        PZ.infiltration_to_filter_material.append(PZ.f_infiltration_to_filter_material(USZ.Ks, PZ.height[time], USZ.height[time], USZ.A, GP.dt, USZ.wilting_point_moisture[time], USZ.porosity[time], PZ.infiltration_to_surrounding[time]))
        PZ.height_after.append(PZ.f_height_after(time, GP))

        # USZ
        USZ.wilting_point_estimated.append(USZ.f_wilting_point_moisture_estimation(time, PZ, GP))
        USZ.capillary_rise.append(USZ.f_capillary_rise(GP.evapotranspiration_max[time], USZ.wilting_point_estimated[time], GP.Kc))
        USZ.infiltration_to_sz.append(USZ.f_infiltration_to_sz(PZ.height_after[-1], USZ.height[time], USZ.porosity[time], GP.dt, USZ.wilting_point_estimated[time]))
        USZ.wilting_point_estimated_after.append(USZ.f_wilting_point_moisture_next_interaction_estimation(time, SZ))

        # EVT
        PZ.evapotranspiration_overall.append(PZ.f_evapotranspiration_overall(USZ.sw, USZ.sh, USZ.ss, GP.Kc, GP.evapotranspiration_max[time], USZ.A, USZ.wilting_point_estimated_after[time],GP.dt))
        PZ.evapotranspiration.append(PZ.f_evapotranspiration(time, USZ, SZ))
        USZ.evapotranspiration.append(USZ.f_evapotranspiration(time, PZ))

        # SZ
        SZ.height_estimated.append(SZ.f_height_estimation(time, GP, USZ))
        SZ.infiltration_to_surround.append(SZ.f_infiltration_to_surround(USZ.Kf, USZ.A, PZ.Cs, SZ.height_estimated[time]))
        SZ.pipe_outflow.append(SZ.f_underdrain_flow(GP.hpipe, USZ.A, SZ.porosity[time], GP.dt, SZ.infiltration_to_surround[time], GP.Apipe, SZ.height_estimated[time], GP.Cd))
        SZ.height.append(SZ.f_height(time, GP, USZ))
        USZ.height.append(USZ.f_height(GP, SZ))

        # Porosity#
        SZ.porosity.append(SZ.f_porosity(SZ.height[-1], GP.L, GP.Dt, GP.Dg, GP.nf, GP.nt, GP.ng))
        USZ.porosity.append(USZ.f_porosity(USZ.height[-1], SZ.height[-1], GP.ng, GP.Dg, GP.Df))

        if time == 0:
            husz_a = USZ.husz_ini
            nusz_a = USZ.nusz_ini
            s_a = USZ.sw
        else:
            husz_a = WFR.thusz[time - 1]
            nusz_a = WFR.tnusz[time - 1]
            s_a = WFR.ts[time - 1]

        USZ.wilting_point_moisture.append(max(min(1.0, (s_a * husz_a * nusz_a * USZ.A + GP.dt * (PZ.infiltration_to_filter_material[time] + USZ.capillary_rise[time] - USZ.infiltration_to_sz[time] - PZ.evapotranspiration[time])) / (USZ.A * USZ.height[-1] * USZ.porosity[-1])), USZ.sh))

        # save all results to WFR
        WFR.tt.append(time)
        WFR.tQover.append(PZ.overflow[time])
        WFR.tQpf.append(PZ.infiltration_to_filter_material[time])
        WFR.tQinfp.append(PZ.infiltration_to_surrounding[time])
        WFR.tQfs.append(USZ.infiltration_to_sz[time])
        WFR.tQhc.append(USZ.capillary_rise[time])
        WFR.tQet.append(PZ.evapotranspiration_overall[time])
        WFR.tQinfsz.append(SZ.infiltration_to_surround[time])
        WFR.tQpipe.append(SZ.pipe_outflow[time])
        WFR.tQet1.append(PZ.evapotranspiration[time])
        WFR.tQet2.append(USZ.evapotranspiration[time])

        WFR.thp.append(PZ.height[time])
        WFR.ts.append(USZ.wilting_point_moisture[-1])
        WFR.thusz.append(USZ.height[-1])
        WFR.thsz.append(SZ.height[-1])
        WFR.thszEST.append(SZ.height_estimated[time])
        WFR.tnsz.append(SZ.porosity[-1])
        WFR.tnusz.append(USZ.porosity[-1])
        WFR.thpEND.append(PZ.height_after[-1])
        WFR.tteta_usz.append(USZ.wilting_point_moisture[-1] * USZ.porosity[-1])
        WFR.tteta_sz.append(SZ.porosity[-1])


# **4. Model routine**
def water_quality_module(WFR, GP, USZ, PZ, SZ, SOIL_PLANT, NH4, NO3, O2, DOC):
    for time in range(len(WFR.indice) - 1):
        # PZ
        height_pz = WFR.thpEND[time]
        if height_pz < 0.001:
            height_pz = 0

        if time == 0:
            height_pz_before = 0
            Qorif = 0
        else:
            height_pz_before = WFR.thpEND[time - 1]
            Qorif = WFR.tQpipe[time]

        if time < (len(WFR.indice) - 1):
            theta_usz_after = WFR.tteta_usz[time + 1]
            theta_sz_after = WFR.tteta_sz[time + 1]
        else:
            theta_usz_after = WFR.tteta_usz[time]
            theta_sz_after = WFR.tteta_sz[time]

        O2.reaction_rate_pz_now = O2.f_reaction_pz()
        NH4.reaction_rate_pz_now = NH4.f_reaction_pz()
        NO3.reaction_rate_pz_now = NO3.f_reaction_pz()
        DOC.reaction_rate_pz_now = DOC.f_reaction_pz()

        if height_pz == 0:
            O2.concentration_pz_now = 0
            NH4.concentration_pz_now = 0
            NO3.concentration_pz_now = 0
            DOC.concentration_pz_now = 0

        else:
            O2.concentration_pz_now = PZ.f_concentration(O2.concentration_inflow[time], WFR.tQin[time], O2.concentration_pz_before, WFR.tQpf[time], WFR.tQover[time], O2.reaction_rate_pz_now, height_pz, height_pz_before, GP.dt)
            NH4.concentration_pz_now = PZ.f_concentration(NH4.concentration_inflow[time], WFR.tQin[time], NH4.concentration_pz_before, WFR.tQpf[time], WFR.tQover[time], NH4.reaction_rate_pz_now, height_pz, height_pz_before, GP.dt)
            NO3.concentration_pz_now = PZ.f_concentration(NO3.concentration_inflow[time], WFR.tQin[time], NO3.concentration_pz_before, WFR.tQpf[time], WFR.tQover[time], NO3.reaction_rate_pz_now, height_pz, height_pz_before, GP.dt)
            DOC.concentration_pz_now = PZ.f_concentration(DOC.concentration_inflow[time], WFR.tQin[time], DOC.concentration_pz_before, WFR.tQpf[time], WFR.tQover[time], DOC.reaction_rate_pz_now, height_pz, height_pz_before, GP.dt)

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
                NH4.concentration_soil_usz_now.append(NH4.f_concentration_soil(NH4.concentration_soil_usz_layer, WFR.tteta_usz[time], NH4.kads, NH4.concentration_usz_layers[usz_layer], NH4.kdes, NH4.k_mb, SOIL_PLANT.ro, GP.dt))

                DOC.concentration_soil_usz_before = DOC.concentration_soil_usz[time][usz_layer]
                DOC.concentration_soil_usz_layer = DOC.f_concentration_soil(DOC.concentration_soil_usz_before, WFR.tteta_usz[time], DOC.kads, DOC.concentration_usz_layers[usz_layer], DOC.kdes, DOC.k_mb, SOIL_PLANT.ro, GP.dt)
                DOC.concentration_soil_usz_now.append(DOC.concentration_soil_usz_layer)

                USZ.unit_flux_now = USZ.f_unit_flux(usz_layer, WFR.tQpf[time], WFR.tQet1[time], WFR.tQfs[time], WFR.tQhc[time], Qorif, WFR.tQinfsz[time], WFR.tteta_usz[time], GP.hpipe, PZ.Ab)
                USZ.unit_flux.append(USZ.unit_flux_now)

                O2.reaction_rate_usz_now = O2.f_reaction_usz(O2.concentration_usz_layers[usz_layer], NH4.concentration_usz_layers[usz_layer], GP.k_nit) + O2.f_plant_uptake_usz(O2.concentration_usz_layers[usz_layer], SOIL_PLANT.c_o2_root, WFR.tteta_usz[time], SOIL_PLANT.root_fraction, SOIL_PLANT.lamda)
                NH4.reaction_rate_usz_now = NH4.f_reaction_usz(NH4.concentration_usz_layers[usz_layer], GP.k_nit) + NH4.f_plant_uptake_usz(NH4.concentration_usz_layers[usz_layer], WFR.tteta_usz[time], SOIL_PLANT.root_fraction)
                NO3.reaction_rate_usz_now = NO3.f_reaction_usz(NH4.concentration_usz_layers[usz_layer], GP.k_nit) + NO3.f_plant_uptake_usz(NO3.concentration_usz_layers[usz_layer], WFR.tteta_usz[time], SOIL_PLANT.root_fraction)
                DOC.reaction_rate_usz_now = DOC.f_reaction_usz(DOC.concentration_usz_layers[usz_layer])

                O2.reaction_rate_usz_layers.append(O2.reaction_rate_usz_now * (1 / theta_usz_after) * GP.dt)
                NH4.reaction_rate_usz_layers.append(NH4.reaction_rate_usz_now * (1 / theta_usz_after) * GP.dt)
                NO3.reaction_rate_usz_layers.append(NO3.reaction_rate_usz_now * (1 / theta_usz_after) * GP.dt)
                DOC.reaction_rate_usz_layers.append(DOC.reaction_rate_usz_now * (1 / theta_usz_after) * GP.dt)

                O2.peclet = USZ.f_peclet(USZ.unit_flux_now, O2.D, GP.dz)
                O2.delta_concentration = O2.f_delta_concentration_usz(time, usz_layer, USZ.m_usz, WFR.tteta_usz, theta_usz_after, USZ.unit_flux_now, SOIL_PLANT.ro, SOIL_PLANT.f, GP.dt, GP.dz)
                O2.concentration_usz_layers_now.append(O2.delta_concentration)

                NH4.peclet = USZ.f_peclet(USZ.unit_flux_now, NH4.D, GP.dz)
                NH4.delta_concentration = NH4.f_delta_concentration_usz(time, usz_layer, USZ.m_usz, WFR.tteta_usz, theta_usz_after, USZ.unit_flux_now, SOIL_PLANT.ro, SOIL_PLANT.f, GP.dt, GP.dz)
                NH4.concentration_usz_layers_now.append(NH4.delta_concentration)

                NO3.peclet = USZ.f_peclet(USZ.unit_flux_now, NO3.D, GP.dz)
                NO3.delta_concentration = NO3.f_delta_concentration_usz(time, usz_layer, USZ.m_usz, WFR.tteta_usz, theta_usz_after, USZ.unit_flux_now, SOIL_PLANT.ro, SOIL_PLANT.f, GP.dt, GP.dz)
                NO3.concentration_usz_layers_now.append(NO3.delta_concentration)

                DOC.peclet = USZ.f_peclet(USZ.unit_flux_now, DOC.D, GP.dz)
                DOC.delta_concentration = DOC.f_delta_concentration_usz(time, usz_layer, USZ.m_usz, WFR.tteta_usz, theta_usz_after, USZ.unit_flux_now, SOIL_PLANT.ro, SOIL_PLANT.f, GP.dt, GP.dz)
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
            
            NH4.concentration_soil_sz_before = NH4.concentration_soil_sz[time][sz_layer]
            NH4.concentration_soil_sz_layer = NH4.f_concentration_soil(NH4.concentration_soil_sz_before, WFR.tteta_sz[time], NH4.kads2, NH4.concentration_sz_layers[sz_layer], NH4.kdes2, NH4.k_mb, SOIL_PLANT.ro, GP.dt)
            NH4.concentration_soil_sz_now.append(NH4.concentration_soil_sz_layer)

            DOC.concentration_soil_sz_before = DOC.concentration_soil_sz[time][sz_layer]
            DOC.concentration_soil_sz_layer = DOC.f_concentration_soil(DOC.concentration_soil_sz_before, WFR.tteta_sz[time], DOC.kads2, DOC.concentration_sz_layers[sz_layer], DOC.kdes2, DOC.k_mb, SOIL_PLANT.ro, GP.dt)
            DOC.concentration_soil_sz_now.append(DOC.concentration_soil_sz_layer)

            SZ.unit_flux_now = SZ.f_unit_flux(sz_layer, WFR.tQfs[time], WFR.tQhc[time], WFR.tQet2[time], Qorif, WFR.tQinfsz[time], WFR.tteta_sz[time], PZ.Ab)
            SZ.unit_flux.append(SZ.unit_flux_now)

            O2.reaction_rate_sz_now = O2.f_reaction_sz(O2.concentration_sz_layers[sz_layer], NH4.concentration_sz_layers[sz_layer], GP.k_nit) + O2.f_plant_uptake_sz(O2.concentration_sz_layers[sz_layer], SOIL_PLANT.c_o2_root, WFR.tteta_sz[time], SOIL_PLANT.root_fraction, SOIL_PLANT.lamda)
            NH4.reaction_rate_sz_now = NH4.f_reaction_sz() + NH4.f_plant_uptake_sz(NH4.concentration_sz_layers[sz_layer], WFR.tteta_sz[time], SOIL_PLANT.root_fraction)
            NO3.reaction_rate_sz_now = NO3.f_reaction_sz(NO3.concentration_sz_layers[sz_layer], O2.concentration_sz_layers[sz_layer], DOC.concentration_sz_layers[sz_layer], GP.k_denit) + NO3.f_plant_uptake_sz(NO3.concentration_sz_layers[sz_layer], WFR.tteta_sz[time], SOIL_PLANT.root_fraction)
            DOC.reaction_rate_sz_now = DOC.f_reaction_sz(DOC.concentration_sz_layers[sz_layer])

            O2.reaction_rate_sz_layers.append(O2.reaction_rate_sz_now * (1 / theta_sz_after) * GP.dt)
            NH4.reaction_rate_sz_layers.append(NH4.reaction_rate_sz_now * (1 / theta_sz_after) * GP.dt)
            NO3.reaction_rate_sz_layers.append(NO3.reaction_rate_sz_now * (1 / theta_sz_after) * GP.dt)
            DOC.reaction_rate_sz_layers.append(DOC.reaction_rate_sz_now * (1 / theta_sz_after) * GP.dt)

            ### Oxygen
            O2.peclet = SZ.f_peclet(SZ.unit_flux_now, O2.D, GP.dz)
            O2.delta_concentration = O2.f_delta_concentration_sz(time, sz_layer, USZ.m_usz, SZ.m_sz, WFR.tteta_sz, theta_sz_after, SZ.unit_flux_now, SOIL_PLANT.ro, SOIL_PLANT.f, GP.dt, GP.dz, GP.n)
            O2.concentration_sz_layers_now.append(O2.delta_concentration)

            ### Amonia
            NH4.peclet = SZ.f_peclet(SZ.unit_flux_now, NH4.D, GP.dz)
            NH4.delta_concentration = NH4.f_delta_concentration_sz(time, sz_layer, USZ.m_usz, SZ.m_sz, WFR.tteta_sz, theta_sz_after, SZ.unit_flux_now, SOIL_PLANT.ro, SOIL_PLANT.f, GP.dt, GP.dz, GP.n)
            NH4.concentration_sz_layers_now.append(NH4.delta_concentration)

            ### Nitrate
            NO3.peclet = SZ.f_peclet(SZ.unit_flux_now, NO3.D, GP.dz)
            NO3.delta_concentration = NO3.f_delta_concentration_sz(time, sz_layer, USZ.m_usz, SZ.m_sz, WFR.tteta_sz, theta_sz_after, SZ.unit_flux_now, SOIL_PLANT.ro, SOIL_PLANT.f, GP.dt, GP.dz, GP.n)
            NO3.concentration_sz_layers_now.append(NO3.delta_concentration)

            ### DOC
            DOC.peclet = SZ.f_peclet(SZ.unit_flux_now, DOC.D, GP.dz)
            DOC.delta_concentration = DOC.f_delta_concentration_sz(time, sz_layer, USZ.m_usz, SZ.m_sz, WFR.tteta_sz, theta_sz_after, SZ.unit_flux_now, SOIL_PLANT.ro, SOIL_PLANT.f, GP.dt, GP.dz, GP.n)
            DOC.concentration_sz_layers_now.append(DOC.delta_concentration)

        # Corrective step
        O2.mass_balance.append(O2.f_mass_balance_stormwater(time, USZ, SZ, WFR, PZ, GP))
        O2.mass_stormwater_now = O2.f_mass_stormwater(time, USZ.m_usz, SZ.m_sz, PZ.Ab, WFR.tteta_usz, WFR.thusz, WFR.tteta_sz, WFR.thsz, GP.hpipe)
        O2.f_mass_balance_concentration(time, USZ.m_usz, PZ.Ab, WFR.tteta_usz, WFR.thusz)
        
        NH4.mass_balance.append(NH4.f_mass_balance_stormwater(time, USZ, SZ, WFR, PZ, GP, SOIL_PLANT))
        NH4.mass_stormwater_now = NH4.f_mass_stormwater(time, USZ.m_usz, SZ.m_sz, PZ.Ab, WFR.tteta_usz, WFR.thusz, WFR.tteta_sz, WFR.thsz, GP.hpipe)
        NH4.f_mass_balance_concentration(time, USZ.m_usz, PZ.Ab, WFR.tteta_usz, WFR.thusz)

        NO3.mass_balance.append(NO3.f_mass_balance_stormwater(time, USZ, SZ, WFR, PZ, GP))
        NO3.mass_stormwater_now = NO3.f_mass_stormwater(time, USZ.m_usz, SZ.m_sz, PZ.Ab, WFR.tteta_usz, WFR.thusz, WFR.tteta_sz, WFR.thsz, GP.hpipe)
        NO3.f_mass_balance_concentration(time, USZ.m_usz, PZ.Ab, WFR.tteta_usz, WFR.thusz)

        DOC.mass_balance.append(DOC.f_mass_balance_stormwater(time, USZ, SZ, WFR, PZ, GP, SOIL_PLANT))
        DOC.mass_stormwater_now = DOC.f_mass_stormwater(time, USZ.m_usz, SZ.m_sz, PZ.Ab, WFR.tteta_usz, WFR.thusz, WFR.tteta_sz, WFR.thsz, GP.hpipe)
        DOC.f_mass_balance_concentration(time, USZ.m_usz, PZ.Ab, WFR.tteta_usz, WFR.thusz)

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

        ## adding layers of SZ in time
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
    ## Oxygen
    data_o2 = O2.water_quality_results(USZ.m_usz, SZ.m_sz)
    data_nh4 = NH4.water_quality_results(USZ.m_usz, SZ.m_sz)
    data_no3 = NO3.water_quality_results(USZ.m_usz, SZ.m_sz)
    data_doc = DOC.water_quality_results(USZ.m_usz, SZ.m_sz)

    return data_o2, data_nh4, data_no3, data_doc

def run(interaction):
    SETUP_FILE = "parameters.ini"
    WATER_FLOW_INPUT_FILE = "water_inflow.csv"
    WATER_QUALITY_INPUT_FILE = "concentration_inflow.csv"
    GP = parameters.GeneralParameters(SETUP_FILE, WATER_FLOW_INPUT_FILE)
    WFR = parameters.WaterFlowResults(GP.rain_inflow, GP.inflow, GP.evapotranspiration_max)
    GP.hpipe = hpipes[interaction]
    USZ = parameters.UnsaturatedZone(SETUP_FILE, GP.L, GP.hpipe, GP.dz)
    PZ = parameters.PondingZone(SETUP_FILE)
    SZ = parameters.SaturatedZone(SETUP_FILE, GP, USZ)
    SOIL_PLANT = parameters.SoilPlant(SETUP_FILE, USZ.nusz_ini)
    NH4 = parameters.Ammonia(USZ.m_usz, SZ.m_sz, SETUP_FILE, WATER_QUALITY_INPUT_FILE)
    NO3 = parameters.Nitrate(USZ.m_usz, SZ.m_sz, SETUP_FILE, WATER_QUALITY_INPUT_FILE)
    O2 = parameters.Oxygen(USZ.m_usz, SZ.m_sz, SETUP_FILE, WATER_QUALITY_INPUT_FILE)
    DOC = parameters.DissolvedOrganicCarbon(USZ.m_usz, SZ.m_sz, SETUP_FILE, WATER_QUALITY_INPUT_FILE)
    water_flow_module(WFR, GP, USZ, PZ, SZ)
    #WFR.water_balance(GP.dt)
    print("hpipe:", GP.hpipe)
    data_o2, data_nh4, data_no3, data_doc = water_quality_module(WFR, GP, USZ, PZ,
                                                                 SZ, SOIL_PLANT, NH4, NO3, O2, DOC)

    wf_path = results_paths[interaction] + "water_flow_results.csv"
    nh4_path = results_paths[interaction] + "results_Kin_pf_nh4_2.csv"
    o2_path = results_paths[interaction] + "results_Kin_pf_o2.csv"
    no3_path = results_paths[interaction] + "results_Kin_pf_no3.csv"
    doc_path = results_paths[interaction] + "results_Kin_pf_doc.csv"

    wf_test = results_tests.water_flow_comparison_test(wf_path, WFR)
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
    inicio = datetime.datetime.now()

    for test in range(len(hpipes)):
        run(test)


    #data_nh4.to_csv('results_Kin_pf_nh4_2.csv', index = False)
    #data_o2.to_csv('results_Kin_pf_o2.csv', index = False)
    #data_no3.to_csv('results_Kin_pf_no3.csv', index = False)
    #data_doc.to_csv('results_Kin_pf_doc.csv', index = False)

    fim = datetime.datetime.now()
    print('Elapsed time: ', fim - inicio)
    print('Done!')

