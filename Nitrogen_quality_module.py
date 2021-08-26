#!/usr/bin/env python

import parameters
import results_tests
import datetime
import pandas as pd

def water_flow_module(WFR, GENERAL_PARAMETERS, USZ, PZ, SZ):
    hpEND = 0
    husz = USZ.husz_ini

    if GENERAL_PARAMETERS.hpipe > 0:
        hsz = GENERAL_PARAMETERS.hpipe
    else:
        hsz = GENERAL_PARAMETERS.dpipe / 1000 + 0.03

    nusz = USZ.nusz_ini
    nsz = GENERAL_PARAMETERS.ng
    s = USZ.sw
    for time in range(len(WFR.tQrain)):
        Qin = WFR.tQin[time]
        Qrain = WFR.tQrain[time]
        Emax = WFR.tEmax[time]

        # PZ#
        Qover = PZ.f_weir_overflow(hpEND, GENERAL_PARAMETERS.dt, Qin, Qrain)
        height_pz = max(hpEND + GENERAL_PARAMETERS.dt / PZ.Ab * (WFR.tQrain[time] + Qin - Qover), 0)  # beginning
        Qinfp = PZ.f_infiltration_to_surrounding(USZ.Kf, USZ.A, height_pz)
        Qpf = PZ.f_infiltration_to_filter_material(USZ.Ks, height_pz, husz, USZ.A, GENERAL_PARAMETERS.dt, s, nusz, Qinfp)
        hpEND = max(height_pz - GENERAL_PARAMETERS.dt / PZ.Ab * (Qpf + Qinfp), 0)  # end

        # USZ#
        sEST = max(min(s + Qpf * GENERAL_PARAMETERS.dt / (nusz * USZ.A * husz), 1), 0)
        Qhc = USZ.f_capillary_rise(Emax, sEST, GENERAL_PARAMETERS.Kc)
        Qfs = USZ.f_infiltration_to_sz(hpEND, husz, nusz, GENERAL_PARAMETERS.dt, sEST)
        sEST2 = (sEST * nusz * husz + nsz * hsz) / (nusz * husz + nsz * hsz)
        Qet = PZ.f_evapotranspiration(USZ.sw, USZ.sh, USZ.ss, GENERAL_PARAMETERS.Kc, Emax, USZ.A, sEST2,
                                      GENERAL_PARAMETERS.dt)
        Qet1 = Qet * (sEST * nusz * husz) / (sEST * nusz * husz + nsz * hsz)
        Qet2 = Qet - Qet1

        # SZ#
        hszEST = hsz + GENERAL_PARAMETERS.dt * (Qfs - Qhc - Qet2) / USZ.A / nsz
        Qinfsz = SZ.f_infiltration_to_surround(USZ.Kf, USZ.A, PZ.Cs, hszEST)
        Qpipe = SZ.f_underdrain_flow(GENERAL_PARAMETERS.hpipe, USZ.A, nsz, GENERAL_PARAMETERS.dt, Qinfsz,
                                     GENERAL_PARAMETERS.Apipe, hszEST, GENERAL_PARAMETERS.Cd)
        hsz = hsz + GENERAL_PARAMETERS.dt * (Qfs - Qhc - Qinfsz - Qpipe - Qet2) / USZ.A / nsz
        husz = GENERAL_PARAMETERS.L - hsz

        # Porosity#
        nsz = SZ.f_porosity(hsz, GENERAL_PARAMETERS.L, GENERAL_PARAMETERS.Dt, GENERAL_PARAMETERS.Dg,
                            GENERAL_PARAMETERS.nf, GENERAL_PARAMETERS.nt, GENERAL_PARAMETERS.ng)
        nusz = USZ.f_porosity(husz, hsz, GENERAL_PARAMETERS.ng, GENERAL_PARAMETERS.Dg, GENERAL_PARAMETERS.Df)

        if time == 0:
            husz_a = USZ.husz_ini
            nusz_a = USZ.nusz_ini
            s_a = USZ.sw
        else:
            husz_a = WFR.thusz[time - 1]
            nusz_a = WFR.tnusz[time - 1]
            s_a = WFR.ts[time - 1]

        s = max(min(1.0, (s_a * husz_a * nusz_a * USZ.A + GENERAL_PARAMETERS.dt * (Qpf + Qhc - Qfs - Qet1)) / (
                    USZ.A * husz * nusz)), USZ.sh)

        # save all results to WFR
        WFR.tt.append(time)
        WFR.tQover.append(Qover)
        WFR.tQpf.append(Qpf)
        WFR.tQinfp.append(Qinfp)
        WFR.tQfs.append(Qfs)
        WFR.tQhc.append(Qhc)
        WFR.tQet.append(Qet)
        WFR.tQinfsz.append(Qinfsz)
        WFR.tQpipe.append(Qpipe)
        WFR.tQet1.append(Qet1)
        WFR.tQet2.append(Qet2)

        WFR.thp.append(height_pz)
        WFR.ts.append(s)
        WFR.thusz.append(husz)
        WFR.thsz.append(hsz)
        WFR.thszEST.append(hszEST)
        WFR.tnsz.append(nsz)
        WFR.tnusz.append(nusz)
        WFR.thpEND.append(hpEND)
        WFR.tteta_usz.append(s * nusz)
        WFR.tteta_sz.append(nsz)


# **4. Model routine**
def water_quality_module(WFR, GENERAL_PARAMETERS, USZ, PZ, SZ, SOIL_PLANT, NH4, NO3, O2, DOC):
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
            O2.concentration_pz_now = PZ.f_concentration(O2.concentration_inflow[time], WFR.tQin[time], O2.concentration_pz_before, WFR.tQpf[time], WFR.tQover[time], O2.reaction_rate_pz_now, height_pz, height_pz_before, GENERAL_PARAMETERS.dt)
            NH4.concentration_pz_now = PZ.f_concentration(NH4.concentration_inflow[time], WFR.tQin[time], NH4.concentration_pz_before, WFR.tQpf[time], WFR.tQover[time], NH4.reaction_rate_pz_now, height_pz, height_pz_before, GENERAL_PARAMETERS.dt)
            NO3.concentration_pz_now = PZ.f_concentration(NO3.concentration_inflow[time], WFR.tQin[time], NO3.concentration_pz_before, WFR.tQpf[time], WFR.tQover[time], NO3.reaction_rate_pz_now, height_pz, height_pz_before, GENERAL_PARAMETERS.dt)
            DOC.concentration_pz_now = PZ.f_concentration(DOC.concentration_inflow[time], WFR.tQin[time], DOC.concentration_pz_before, WFR.tQpf[time], WFR.tQover[time], DOC.reaction_rate_pz_now, height_pz, height_pz_before, GENERAL_PARAMETERS.dt)

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
        if GENERAL_PARAMETERS.hpipe > 0:
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

                O2.concentration_soil_usz_now = O2.f_concentration_soil()
                O2.concentration_soil_usz_layer = O2.f_concentration_soil()
                NO3.concentration_soil_usz_now = NO3.f_concentration_soil()
                NO3.concentration_soil_usz_layer = NO3.f_concentration_soil()

                NH4.concentration_soil_usz_layer = NH4.concentration_soil_usz[time][usz_layer]
                NH4.concentration_soil_usz_now.append(NH4.f_concentration_soil(NH4.concentration_soil_usz_layer, WFR.tteta_usz[time], NH4.kads, NH4.concentration_usz_layers[usz_layer], NH4.kdes, NH4.k_mb, SOIL_PLANT.ro, GENERAL_PARAMETERS.dt))
                NO3.concentration_soil_usz_now = NO3.f_concentration_soil()

                DOC.concentration_soil_usz_before = DOC.concentration_soil_usz[time][usz_layer]
                DOC.concentration_soil_usz_layer = DOC.f_concentration_soil(DOC.concentration_soil_usz_before, WFR.tteta_usz[time], DOC.kads, DOC.concentration_usz_layers[usz_layer], DOC.kdes, DOC.k_mb, SOIL_PLANT.ro, GENERAL_PARAMETERS.dt)
                DOC.concentration_soil_usz_now.append(DOC.concentration_soil_usz_layer)

                USZ.unit_flux_now = USZ.f_unit_flux(usz_layer, WFR.tQpf[time], WFR.tQet1[time], WFR.tQfs[time], WFR.tQhc[time], Qorif, WFR.tQinfsz[time], WFR.tteta_usz[time], GENERAL_PARAMETERS.hpipe, PZ.Ab)
                USZ.unit_flux.append(USZ.unit_flux_now)

                O2.reaction_rate_usz_now = O2.f_reaction_usz(O2.concentration_usz_layers[usz_layer], NH4.concentration_usz_layers[usz_layer], GENERAL_PARAMETERS.k_nit) + O2.f_plant_uptake_usz(O2.concentration_usz_layers[usz_layer], SOIL_PLANT.c_o2_root, WFR.tteta_usz[time], SOIL_PLANT.root_fraction, SOIL_PLANT.lamda)
                NH4.reaction_rate_usz_now = NH4.f_reaction_usz(NH4.concentration_usz_layers[usz_layer], GENERAL_PARAMETERS.k_nit) + NH4.f_plant_uptake_usz(NH4.concentration_usz_layers[usz_layer], WFR.tteta_usz[time], SOIL_PLANT.root_fraction)
                NO3.reaction_rate_usz_now = NO3.f_reaction_usz(NH4.concentration_usz_layers[usz_layer], GENERAL_PARAMETERS.k_nit) + NO3.f_plant_uptake_usz(NO3.concentration_usz_layers[usz_layer], WFR.tteta_usz[time], SOIL_PLANT.root_fraction)
                DOC.reaction_rate_usz_now = DOC.f_reaction_usz(DOC.concentration_usz_layers[usz_layer])

                O2.reaction_rate_usz_layers.append(O2.reaction_rate_usz_now * (1 / theta_usz_after) * GENERAL_PARAMETERS.dt)
                NH4.reaction_rate_usz_layers.append(NH4.reaction_rate_usz_now * (1 / theta_usz_after) * GENERAL_PARAMETERS.dt)
                NO3.reaction_rate_usz_layers.append(NO3.reaction_rate_usz_now * (1 / theta_usz_after) * GENERAL_PARAMETERS.dt)
                DOC.reaction_rate_usz_layers.append(DOC.reaction_rate_usz_now * (1 / theta_usz_after) * GENERAL_PARAMETERS.dt)

                O2.peclet = USZ.f_peclet(USZ.unit_flux_now, O2.D, GENERAL_PARAMETERS.dz)
                O2.delta_concentration = O2.f_delta_concentration_usz(time, usz_layer, USZ.m_usz, WFR.tteta_usz, theta_usz_after, USZ.unit_flux_now, SOIL_PLANT.ro, SOIL_PLANT.f, GENERAL_PARAMETERS.dt, GENERAL_PARAMETERS.dz)
                O2.concentration_usz_layers_now.append(O2.delta_concentration)

                NH4.peclet = USZ.f_peclet(USZ.unit_flux_now, NH4.D, GENERAL_PARAMETERS.dz)
                NH4.delta_concentration = NH4.f_delta_concentration_usz(time, usz_layer, USZ.m_usz, WFR.tteta_usz, theta_usz_after, USZ.unit_flux_now, SOIL_PLANT.ro, SOIL_PLANT.f, GENERAL_PARAMETERS.dt, GENERAL_PARAMETERS.dz)
                NH4.concentration_usz_layers_now.append(NH4.delta_concentration)

                NO3.peclet = USZ.f_peclet(USZ.unit_flux_now, NO3.D, GENERAL_PARAMETERS.dz)
                NO3.delta_concentration = NO3.f_delta_concentration_usz(time, usz_layer, USZ.m_usz, WFR.tteta_usz, theta_usz_after, USZ.unit_flux_now, SOIL_PLANT.ro, SOIL_PLANT.f, GENERAL_PARAMETERS.dt, GENERAL_PARAMETERS.dz)
                NO3.concentration_usz_layers_now.append(NO3.delta_concentration)

                ### Nitrate
                DOC.peclet = USZ.f_peclet(USZ.unit_flux_now, DOC.D, GENERAL_PARAMETERS.dz)
                DOC.delta_concentration = DOC.f_delta_concentration_usz(time, usz_layer, USZ.m_usz, WFR.tteta_usz, theta_usz_after, USZ.unit_flux_now, SOIL_PLANT.ro, SOIL_PLANT.f, GENERAL_PARAMETERS.dt, GENERAL_PARAMETERS.dz)
                DOC.concentration_usz_layers_now.append(DOC.delta_concentration)

        #####   SZ   #####

        cj_i1_o2 = []
        Rxj_o2 = []

        cj_i1_nh4 = []
        csi_sz_nh4 = []
        Rxj_nh4 = []

        cj_i1_no3 = []
        Rxj_no3 = []

        cj_i1_doc = []
        csi_sz_doc = []
        Rxj_doc = []

        for sz_layer in range(SZ.m_sz):
            cj_o2 = O2.concentration_sz_layers[sz_layer]
            cjminus1_o2 = O2.concentration_sz_layers[sz_layer - 1]
            if sz_layer < (SZ.m_sz - 1):
                cjplus1_o2 = O2.concentration_sz_layers[sz_layer + 1]
            else:
                cjplus1_o2 = 0
            cmminus1_o2 = O2.concentration_usz_layers[USZ.m_usz - 1]

            cj_nh4 = NH4.concentration_sz_layers[sz_layer]
            cjminus1_nh4 = NH4.concentration_sz_layers[sz_layer - 1]
            if sz_layer < (SZ.m_sz - 1):
                cjplus1_nh4 = NH4.concentration_sz_layers[sz_layer + 1]
            else:
                cjplus1_nh4 = 0
            cmminus1_nh4 = NH4.concentration_usz_layers[USZ.m_usz - 1]

            cj_no3 = NO3.concentration_sz_layers[sz_layer]
            cjminus1_no3 = NO3.concentration_sz_layers[sz_layer - 1]
            if sz_layer < (SZ.m_sz - 1):
                cjplus1_no3 = NO3.concentration_sz_layers[sz_layer + 1]
            else:
                cjplus1_no3 = 0
            cmminus1_no3 = NO3.concentration_usz_layers[USZ.m_usz - 1]

            cj_doc = DOC.concentration_sz_layers[sz_layer]
            cjminus1_doc = DOC.concentration_sz_layers[sz_layer - 1]
            if sz_layer < (SZ.m_sz - 1):
                cjplus1_doc = DOC.concentration_sz_layers[sz_layer + 1]
            else:
                cjplus1_doc = 0
            cmminus1_doc = DOC.concentration_usz_layers[USZ.m_usz - 1]

            O2.concentration_soil_sz_now = 0
            NH4.concentration_soil_sz_before = NH4.concentration_soil_sz[time][sz_layer]
            NH4.concentration_soil_sz_layer = NH4.f_concentration_soil(NH4.concentration_soil_sz_before, WFR.tteta_sz[time], NH4.kads2, cj_nh4, NH4.kdes2, NH4.k_mb, SOIL_PLANT.ro, GENERAL_PARAMETERS.dt)
            csi_sz_nh4.append(NH4.concentration_soil_sz_layer)

            NO3.concentration_soil_sz_now = 0

            DOC.concentration_soil_sz_before = DOC.concentration_soil_sz[time][sz_layer]
            DOC.concentration_soil_sz_layer = DOC.f_concentration_soil(DOC.concentration_soil_sz_before, WFR.tteta_sz[time], DOC.kads2, cj_doc, DOC.kdes2, DOC.k_mb, SOIL_PLANT.ro, GENERAL_PARAMETERS.dt)
            csi_sz_doc.append(DOC.concentration_soil_sz_layer)

            UF_sz = []
            alfa2, beta2 = SZ.f_alfa_beta(sz_layer)
            # WFR.tQet2[time] = 0
            UFi_sz = SZ.f_unit_flux(alfa2, beta2, WFR.tQfs[time], WFR.tQhc[time], WFR.tQet2[time], Qorif, WFR.tQinfsz[time], WFR.tteta_sz[time], PZ.Ab)
            #                 print('WFR.tQfs[time], WFR.tQhc[time], WFR.tQet2[time], Qorif, WFR.tQinfsz[time]:', WFR.tQfs[time], WFR.tQhc[time], WFR.tQet2[time], Qorif, WFR.tQinfsz[time])
            #                 print('UF: ', UFi_sz, ', alfa2: ', alfa2, ', beta2: ' ,beta2)

            UF_sz.append(UFi_sz)

            Rxi_3_o2 = O2.f_reaction_sz(cj_o2, cj_nh4, GENERAL_PARAMETERS.k_nit) + O2.f_plant_uptake_sz(cj_o2, SOIL_PLANT.c_o2_root, WFR.tteta_sz[time], SOIL_PLANT.root_fraction, SOIL_PLANT.lamda)
            Rxi_3_nh4 = NH4.f_reaction_sz() + NH4.f_plant_uptake_sz(cj_nh4, WFR.tteta_sz[time], SOIL_PLANT.root_fraction)
            Rxi_3_no3 = NO3.f_reaction_sz(cj_no3, cj_o2, cj_doc, GENERAL_PARAMETERS.k_denit) + NO3.f_plant_uptake_sz(cj_no3, WFR.tteta_sz[time], SOIL_PLANT.root_fraction)
            Rxi_3_doc = DOC.f_reaction_sz(cj_doc)

            Rxj_o2.append(Rxi_3_o2 * (1 / theta_sz_after) * GENERAL_PARAMETERS.dt)
            Rxj_nh4.append(Rxi_3_nh4 * (1 / theta_sz_after) * GENERAL_PARAMETERS.dt)
            Rxj_no3.append(Rxi_3_no3 * (1 / theta_sz_after) * GENERAL_PARAMETERS.dt)
            Rxj_doc.append(Rxi_3_doc * (1 / theta_sz_after) * GENERAL_PARAMETERS.dt)
            # if time < 20:
            # print('time: ', time, 'sz_layer: ', sz_layer, 'O2.reaction_rate_sz: ', Rxi_3_o2, 'NH4.reaction_rate_sz: ', Rxi_3_nh4, 'NO3.reaction_rate_sz: ', Rxi_3_no3)

            ### Oxygen
            Pesz_o2 = SZ.f_peclet(UFi_sz, O2.D, GENERAL_PARAMETERS.dz)

            if USZ.m_usz < (GENERAL_PARAMETERS.n - 1):
                if Pesz_o2 <= 2:
                    if sz_layer == 0:  # first cell
                        dc_o2 = cjplus1_o2 - 2 * cj_o2 + cmminus1_o2
                        dc_dz_o2 = (cjplus1_o2 - cmminus1_o2) / (2 * GENERAL_PARAMETERS.dz)

                    elif sz_layer == (SZ.m_sz - 1):  # last cell
                        dc_o2 = cj_o2 - 2 * cj_o2 + cjminus1_o2
                        dc_dz_o2 = (cj_o2 - cjminus1_o2) / GENERAL_PARAMETERS.dz

                    else:
                        dc_o2 = cjplus1_o2 - 2 * cj_o2 + cjminus1_o2
                        dc_dz_o2 = (cjplus1_o2 - cjminus1_o2) / (2 * GENERAL_PARAMETERS.dz)

                else:  # Pesz > 2
                    if sz_layer == 0:  # first cell
                        dc_o2 = cjplus1_o2 - 2 * cj_o2 + cmminus1_o2
                        dc_dz_o2 = (cj_o2 - cmminus1_o2) / GENERAL_PARAMETERS.dz

                    elif sz_layer == (SZ.m_sz - 1):  # last cell
                        dc_o2 = cj_o2 - 2 * cj_o2 + cjminus1_o2
                        dc_dz_o2 = (cj_o2 - cjminus1_o2) / GENERAL_PARAMETERS.dz

                    else:
                        dc_o2 = cjplus1_o2 - 2 * cj_o2 + cjminus1_o2
                        dc_dz_o2 = (cj_o2 - cjminus1_o2) / GENERAL_PARAMETERS.dz

            if USZ.m_usz == (GENERAL_PARAMETERS.n - 1):
                dc_o2 = cj_o2 - 2 * cj_o2 + cjminus1_o2
                dc_dz_o2 = (cj_o2 - cjminus1_o2) / (2 * GENERAL_PARAMETERS.dz)

            if USZ.m_usz == 0:
                if Pesz_o2 <= 2:
                    if sz_layer == 0:  # first cell
                        dc_o2 = cjplus1_o2 - 2 * cj_o2 + O2.concentration_pz_now
                        dc_dz_o2 = (cjplus1_o2 - O2.concentration_pz_now) / (2 * GENERAL_PARAMETERS.dz)

                    elif sz_layer == (SZ.m_sz - 1):  # last cell
                        dc_o2 = cj_o2 - 2 * cj_o2 + cjminus1_o2
                        dc_dz_o2 = (cj_o2 - cjminus1_o2) / GENERAL_PARAMETERS.dz

                    else:
                        dc_o2 = cjplus1_o2 - 2 * cj_o2 + cjminus1_o2
                        dc_dz_o2 = (cjplus1_o2 - cjminus1_o2) / (2 * GENERAL_PARAMETERS.dz)

                else:  # Pusz > 2
                    if sz_layer == 0:  # first cell
                        dc_o2 = cjplus1_o2 - 2 * cj_o2 + O2.concentration_pz_now
                        dc_dz_o2 = (cj_o2 - O2.concentration_pz_now) / GENERAL_PARAMETERS.dz

                    elif sz_layer == (SZ.m_sz - 1):  # last cell
                        dc_o2 = cj_o2 - 2 * cj_o2 + cjminus1_o2
                        dc_dz_o2 = (cj_o2 - cjminus1_o2) / GENERAL_PARAMETERS.dz

                    else:
                        dc_o2 = cjplus1_o2 - 2 * cj_o2 + cjminus1_o2
                        dc_dz_o2 = (cj_o2 - cjminus1_o2) / GENERAL_PARAMETERS.dz

            delta_c_o2 = O2.f_transport(WFR.tteta_sz[time], theta_sz_after, cj_o2, O2.concentration_soil_sz_now, dc_o2, dc_dz_o2, 0, 0, O2.D, UFi_sz, Rxi_3_o2, SOIL_PLANT.ro, SOIL_PLANT.f, GENERAL_PARAMETERS.dt, GENERAL_PARAMETERS.dz)
            if theta_sz_after > 0:
                ci1_o2 = cj_o2 + delta_c_o2
            else:
                ci1_o2 = 0

            if ci1_o2 <= 0.0000000000000001:
                ci1_o2 = 0
            else:
                ci1_o2 = ci1_o2

            cj_i1_o2.append(ci1_o2)

            ### Amonia
            Pesz_nh4 = SZ.f_peclet(UFi_sz, NH4.D, GENERAL_PARAMETERS.dz)

            if USZ.m_usz < (GENERAL_PARAMETERS.n - 1):
                if Pesz_nh4 <= 2:
                    if sz_layer == 0:  # first cell
                        dc_nh4 = cjplus1_nh4 - 2 * cj_nh4 + cmminus1_nh4
                        dc_dz_nh4 = (cjplus1_nh4 - cmminus1_nh4) / (2 * GENERAL_PARAMETERS.dz)

                    elif sz_layer == (SZ.m_sz - 1):  # last cell
                        dc_nh4 = cj_nh4 - 2 * cj_nh4 + cjminus1_nh4
                        dc_dz_nh4 = (cj_nh4 - cjminus1_nh4) / GENERAL_PARAMETERS.dz

                    else:
                        dc_nh4 = cjplus1_nh4 - 2 * cj_nh4 + cjminus1_nh4
                        dc_dz_nh4 = (cjplus1_nh4 - cjminus1_nh4) / (2 * GENERAL_PARAMETERS.dz)

                else:  # Pesz > 2
                    if sz_layer == 0:  # first cell
                        dc_nh4 = cjplus1_nh4 - 2 * cj_nh4 + cmminus1_nh4
                        dc_dz_nh4 = (cj_nh4 - cmminus1_nh4) / GENERAL_PARAMETERS.dz

                    elif sz_layer == (SZ.m_sz - 1):  # last cell
                        dc_nh4 = cj_nh4 - 2 * cj_nh4 + cjminus1_nh4
                        dc_dz_nh4 = (cj_nh4 - cjminus1_nh4) / GENERAL_PARAMETERS.dz

                    else:
                        dc_nh4 = cjplus1_nh4 - 2 * cj_nh4 + cjminus1_nh4
                        dc_dz_nh4 = (cj_nh4 - cjminus1_nh4) / GENERAL_PARAMETERS.dz

            if USZ.m_usz == (GENERAL_PARAMETERS.n - 1):
                dc_nh4 = cj_nh4 - 2 * cj_nh4 + cjminus1_nh4
                dc_dz_nh4 = (cj_nh4 - cjminus1_nh4) / (2 * GENERAL_PARAMETERS.dz)

            if USZ.m_usz == 0:
                if Pesz_nh4 <= 2:
                    if sz_layer == 0:  # first cell
                        dc_nh4 = cjplus1_nh4 - 2 * cj_nh4 + NH4.concentration_pz_now
                        dc_dz_nh4 = (cjplus1_nh4 - NH4.concentration_pz_now) / (2 * GENERAL_PARAMETERS.dz)

                    elif sz_layer == (SZ.m_sz - 1):  # last cell
                        dc_nh4 = cj_nh4 - 2 * cj_nh4 + cjminus1_nh4
                        dc_dz_nh4 = (cj_nh4 - cjminus1_nh4) / GENERAL_PARAMETERS.dz

                    else:
                        dc_nh4 = cjplus1_nh4 - 2 * cj_nh4 + cjminus1_nh4
                        dc_dz_nh4 = (cjplus1_nh4 - cjminus1_nh4) / (2 * GENERAL_PARAMETERS.dz)

                else:  # Pusz > 2
                    if sz_layer == 0:  # first cell
                        dc_nh4 = cjplus1_nh4 - 2 * cj_nh4 + NH4.concentration_pz_now
                        dc_dz_nh4 = (cj_nh4 - NH4.concentration_pz_now) / GENERAL_PARAMETERS.dz

                    elif sz_layer == (SZ.m_sz - 1):  # last cell
                        dc_nh4 = cj_nh4 - 2 * cj_nh4 + cjminus1_nh4
                        dc_dz_nh4 = (cj_nh4 - cjminus1_nh4) / GENERAL_PARAMETERS.dz

                    else:
                        dc_nh4 = cjplus1_nh4 - 2 * cj_nh4 + cjminus1_nh4
                        dc_dz_nh4 = (cj_nh4 - cjminus1_nh4) / GENERAL_PARAMETERS.dz

            delta_c_nh4 = NH4.f_transport(WFR.tteta_sz[time], theta_sz_after, cj_nh4, NH4.concentration_soil_sz_layer, dc_nh4, dc_dz_nh4, NH4.kads2, NH4.kdes2,
                                  NH4.D, UFi_sz, Rxi_3_nh4, SOIL_PLANT.ro, SOIL_PLANT.f, GENERAL_PARAMETERS.dt, GENERAL_PARAMETERS.dz)
            #                 print('WFR.tteta_sz[time]: ', WFR.tteta_sz[time], ', theta_sz_after: ', theta_sz_after, ', cj_nh4: ', cj_nh4, ', NH4.concentration_soil_sz_layer: ', NH4.concentration_soil_sz_layer, ', dc_nh4: ', dc_nh4, ', dc_dz_nh4: ', dc_dz_nh4, ', NH4.kads2: ', NH4.kads2, ', NH4.kdes2: ', NH4.kdes2, ', NH4.D: ', NH4.D, ', UFi_sz: ', UFi_sz, ', Rxi_3_nh4: ', Rxi_3_nh4)
            #                 print('delta_c_nh4: ', delta_c_nh4)

            if theta_sz_after > 0:
                ci1_nh4 = cj_nh4 + delta_c_nh4
                # print('calculo final - cj_nh4, WFR.tteta_sz[time], theta_sz_after, delta_c_nh4: ', cj_nh4, ', ', WFR.tteta_sz[time], ', ', theta_sz_after, ', ', delta_c_nh4)
            else:
                ci1_nh4 = 0
            # print('ci1_nh4: ', ci1_nh4)

            if ci1_nh4 <= 0.0000000000000001:
                ci1_nh4 = 0
            else:
                ci1_nh4 = ci1_nh4

            cj_i1_nh4.append(ci1_nh4)

            ### Nitrate
            Pesz_no3 = SZ.f_peclet(UFi_sz, NO3.D, GENERAL_PARAMETERS.dz)

            if USZ.m_usz < (GENERAL_PARAMETERS.n - 1):
                if Pesz_no3 <= 2:
                    if sz_layer == 0:  # first cell
                        dc_no3 = cjplus1_no3 - 2 * cj_no3 + cmminus1_no3
                        dc_dz_no3 = (cjplus1_no3 - cmminus1_no3) / (2 * GENERAL_PARAMETERS.dz)

                    elif sz_layer == (SZ.m_sz - 1):  # last cell
                        dc_no3 = cj_no3 - 2 * cj_no3 + cjminus1_no3
                        dc_dz_no3 = (cj_no3 - cjminus1_no3) / GENERAL_PARAMETERS.dz

                    else:
                        dc_no3 = cjplus1_no3 - 2 * cj_no3 + cjminus1_no3
                        dc_dz_no3 = (cjplus1_no3 - cjminus1_no3) / (2 * GENERAL_PARAMETERS.dz)

                else:  # Pesz > 2
                    if sz_layer == 0:  # first cell
                        dc_no3 = cjplus1_no3 - 2 * cj_no3 + cmminus1_no3
                        dc_dz_no3 = (cj_no3 - cmminus1_no3) / GENERAL_PARAMETERS.dz

                    elif sz_layer == (SZ.m_sz - 1):  # last cell
                        dc_no3 = cj_no3 - 2 * cj_no3 + cjminus1_no3
                        dc_dz_no3 = (cj_no3 - cjminus1_no3) / GENERAL_PARAMETERS.dz

                    else:
                        dc_no3 = cjplus1_no3 - 2 * cj_no3 + cjminus1_no3
                        dc_dz_no3 = (cj_no3 - cjminus1_no3) / GENERAL_PARAMETERS.dz

            if USZ.m_usz == (GENERAL_PARAMETERS.n - 1):
                dc_no3 = cj_no3 - 2 * cj_no3 + cjminus1_no3
                dc_dz_no3 = (cj_no3 - cjminus1_no3) / (2 * GENERAL_PARAMETERS.dz)

            if USZ.m_usz == 0:
                if Pesz_no3 <= 2:
                    if sz_layer == 0:  # first cell
                        dc_no3 = cjplus1_no3 - 2 * cj_no3 + NO3.concentration_pz_now
                        dc_dz_no3 = (cjplus1_no3 - NO3.concentration_pz_now) / (2 * GENERAL_PARAMETERS.dz)

                    elif sz_layer == (SZ.m_sz - 1):  # last cell
                        dc_no3 = cj_no3 - 2 * cj_no3 + cjminus1_no3
                        dc_dz_no3 = (cj_no3 - cjminus1_no3) / GENERAL_PARAMETERS.dz

                    else:
                        dc_no3 = cjplus1_no3 - 2 * cj_no3 + cjminus1_no3
                        dc_dz_no3 = (cjplus1_no3 - cjminus1_no3) / (2 * GENERAL_PARAMETERS.dz)

                else:  # Pusz > 2
                    if sz_layer == 0:  # first cell
                        dc_no3 = cjplus1_no3 - 2 * cj_no3 + NO3.concentration_pz_now
                        dc_dz_no3 = (cj_no3 - NO3.concentration_pz_now) / GENERAL_PARAMETERS.dz

                    elif sz_layer == (SZ.m_sz - 1):  # last cell
                        dc_no3 = cj_no3 - 2 * cj_no3 + cjminus1_no3
                        dc_dz_no3 = (cj_no3 - cjminus1_no3) / GENERAL_PARAMETERS.dz

                    else:
                        dc_no3 = cjplus1_no3 - 2 * cj_no3 + cjminus1_no3
                        dc_dz_no3 = (cj_no3 - cjminus1_no3) / GENERAL_PARAMETERS.dz

            delta_c_no3 = NO3.f_transport(WFR.tteta_sz[time], theta_sz_after, cj_no3, NO3.concentration_soil_sz_now, dc_no3, dc_dz_no3, 0, 0, NO3.D, UFi_sz,
                                  Rxi_3_no3, SOIL_PLANT.ro, SOIL_PLANT.f, GENERAL_PARAMETERS.dt, GENERAL_PARAMETERS.dz)
            if theta_sz_after > 0:
                ci1_no3 = cj_no3 + delta_c_no3
            else:
                ci1_no3 = 0

            if ci1_no3 <= 0.0000000000000001:
                ci1_no3 = 0
            else:
                ci1_no3 = ci1_no3

            cj_i1_no3.append(ci1_no3)

            ### DOC
            Pesz_doc = SZ.f_peclet(UFi_sz, DOC.D,GENERAL_PARAMETERS.dz)

            if USZ.m_usz < (GENERAL_PARAMETERS.n - 1):
                if Pesz_doc <= 2:
                    if sz_layer == 0:  # first cell
                        dc_doc = cjplus1_doc - 2 * cj_doc + cmminus1_doc
                        dc_dz_doc = (cjplus1_doc - cmminus1_doc) / (2 * GENERAL_PARAMETERS.dz)

                    elif sz_layer == (SZ.m_sz - 1):  # last cell
                        dc_doc = cj_doc - 2 * cj_doc + cjminus1_doc
                        dc_dz_doc = (cj_doc - cjminus1_doc) / GENERAL_PARAMETERS.dz

                    else:
                        dc_doc = cjplus1_doc - 2 * cj_doc + cjminus1_doc
                        dc_dz_doc = (cjplus1_doc - cjminus1_doc) / (2 * GENERAL_PARAMETERS.dz)

                else:  # Pesz > 2
                    if sz_layer == 0:  # first cell
                        dc_doc = cjplus1_doc - 2 * cj_doc + cmminus1_doc
                        dc_dz_doc = (cj_doc - cmminus1_doc) / GENERAL_PARAMETERS.dz

                    elif sz_layer == (SZ.m_sz - 1):  # last cell
                        dc_doc = cj_doc - 2 * cj_doc + cjminus1_doc
                        dc_dz_doc = (cj_doc - cjminus1_doc) / GENERAL_PARAMETERS.dz

                    else:
                        dc_doc = cjplus1_doc - 2 * cj_doc + cjminus1_doc
                        dc_dz_doc = (cj_doc - cjminus1_doc) / GENERAL_PARAMETERS.dz

            if USZ.m_usz == (GENERAL_PARAMETERS.n - 1):
                dc_doc = cj_doc - 2 * cj_doc + cjminus1_doc
                dc_dz_doc = (cj_doc - cjminus1_doc) / (2 * GENERAL_PARAMETERS.dz)

            if USZ.m_usz == 0:
                if Pesz_doc <= 2:
                    if sz_layer == 0:  # first cell
                        dc_doc = cjplus1_doc - 2 * cj_doc + DOC.concentration_pz_now
                        dc_dz_doc = (cjplus1_doc - DOC.concentration_pz_now) / (2 * GENERAL_PARAMETERS.dz)

                    elif sz_layer == (SZ.m_sz - 1):  # last cell
                        dc_doc = cj_doc - 2 * cj_doc + cjminus1_doc
                        dc_dz_doc = (cj_doc - cjminus1_doc) / GENERAL_PARAMETERS.dz

                    else:
                        dc_doc = cjplus1_doc - 2 * cj_doc + cjminus1_doc
                        dc_dz_doc = (cjplus1_doc - cjminus1_doc) / (2 * GENERAL_PARAMETERS.dz)

                else:  # Pusz > 2
                    if sz_layer == 0:  # first cell
                        dc_doc = cjplus1_doc - 2 * cj_doc + DOC.concentration_pz_now
                        dc_dz_doc = (cj_doc - DOC.concentration_pz_now) / GENERAL_PARAMETERS.dz

                    elif sz_layer == (SZ.m_sz - 1):  # last cell
                        dc_doc = cj_doc - 2 * cj_doc + cjminus1_doc
                        dc_dz_doc = (cj_doc - cjminus1_doc) / GENERAL_PARAMETERS.dz

                    else:
                        dc_doc = cjplus1_doc - 2 * cj_doc + cjminus1_doc
                        dc_dz_doc = (cj_doc - cjminus1_doc) / GENERAL_PARAMETERS.dz

            delta_c_doc = DOC.f_transport(WFR.tteta_sz[time], theta_sz_after, cj_doc, DOC.concentration_soil_sz_layer, dc_doc, dc_dz_doc, DOC.kads2, DOC.kdes2,
                                  DOC.D, UFi_sz, Rxi_3_doc, SOIL_PLANT.ro, SOIL_PLANT.f, GENERAL_PARAMETERS.dt, GENERAL_PARAMETERS.dz)
            if theta_sz_after > 0:
                ci1_doc = cj_doc + delta_c_doc
            else:
                ci1_doc = 0

            if ci1_doc <= 0.0000000000000001:
                ci1_doc = 0
            else:
                ci1_doc = ci1_doc

            cj_i1_doc.append(ci1_doc)

        # Corrective step

        if GENERAL_PARAMETERS.hpipe > 0:
            ### Oxygen
            Mstor_o2_ast = sum(O2.concentration_usz[time]) * PZ.Ab * WFR.tteta_usz[time] * WFR.thusz[time] * 1000 / USZ.m_usz + sum(O2.concentration_sz[time]) * PZ.Ab * \
                           WFR.tteta_sz[time] * WFR.thsz[time] * 1000 / SZ.m_sz
            O2.mass_stormwater.append(Mstor_o2_ast)

            Msoil_o2_a = 0
            Msoil_o2 = 0
            O2.mass_soil.append(0)

            MRx_o2 = - (sum(O2.reaction_rate_usz[time]) * PZ.Ab * WFR.tteta_usz[time] * WFR.thusz[time] * 1000 / USZ.m_usz + sum(O2.reaction_rate_sz[time]) * PZ.Ab *
                        WFR.tteta_sz[time] * WFR.thsz[time] * 1000 / SZ.m_sz)
            if time == 0:
                O2.mass_accumulated_reaction.append(MRx_o2)
            else:
                O2.mass_accumulated_reaction.append(MRx_o2 + O2.mass_accumulated_reaction[-1])

            Min_o2 = WFR.tQin[time] * O2.concentration_inflow[time] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                O2.mass_accumulated_inflow.append(Min_o2)
            else:
                O2.mass_accumulated_inflow.append(Min_o2 + O2.mass_accumulated_inflow[-1])

            Mover_o2 = WFR.tQover[time] * O2.concentration_pz[time] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                O2.mass_accumulated_overflow.append(Mover_o2)
            else:
                O2.mass_accumulated_overflow.append(Mover_o2 + O2.mass_accumulated_overflow[-1])

            Mpipe_o2 = WFR.tQpipe[time] * O2.concentration_sz[time][SZ.m_sz - 1] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                O2.mass_accumulated_pipe_outflow.append(Mpipe_o2)
            else:
                O2.mass_accumulated_pipe_outflow.append(Mpipe_o2 + O2.mass_accumulated_pipe_outflow[-1])

            Minfsz_o2 = WFR.tQinfsz[time] * O2.concentration_sz[time][SZ.m_sz - 1] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                O2.mass_accumulated_infiltrated_to_SZ.append(Minfsz_o2)
            else:
                O2.mass_accumulated_infiltrated_to_SZ.append(Minfsz_o2 + O2.mass_accumulated_infiltrated_to_SZ[-1])

            Met_o2 = WFR.tQet[time] * O2.concentration_usz_layers_now[0] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                O2.mass_accumulated_evapotranspiration.append(Met_o2)
            else:
                O2.mass_accumulated_evapotranspiration.append(Met_o2 + O2.mass_accumulated_evapotranspiration[-1])

            Mpz_o2 = WFR.thpEND[time] * PZ.Ab * 1000 * O2.concentration_pz[time]
            O2.mass_PZ.append(Mpz_o2)

            Mstor_o2_mb = O2.mass_accumulated_inflow[-1] - O2.mass_PZ[-1] - O2.mass_accumulated_overflow[-1] - O2.mass_accumulated_pipe_outflow[-1] - O2.mass_accumulated_infiltrated_to_SZ[
                -1] - O2.mass_accumulated_evapotranspiration[-1] - O2.mass_soil[-1] - O2.mass_accumulated_reaction[-1]
            O2.mass_balance.append(Mstor_o2_mb)

            delta_ast_o2_usz = (Mstor_o2_mb - Mstor_o2_ast) / (PZ.Ab * WFR.tteta_usz[time] * WFR.thusz[time] * 1000 / USZ.m_usz)
            O2.concentration_usz_layers_now[1] = O2.concentration_usz_layers_now[1] + delta_ast_o2_usz
            if O2.concentration_usz_layers_now[1] > 0:
                O2.concentration_usz_layers_now[1] = O2.concentration_usz_layers_now[1]
            else:
                O2.concentration_usz_layers_now[1] = 0

            ### Amonia

            Mstor_nh4_ast = sum(NH4.concentration_usz[time]) * PZ.Ab * WFR.tteta_usz[time] * WFR.thusz[time] * 1000 / USZ.m_usz + sum(NH4.concentration_sz[time]) * PZ.Ab * \
                            WFR.tteta_sz[time] * WFR.thsz[time] * 1000 / SZ.m_sz
            NH4.mass_stormwater.append(Mstor_nh4_ast)

            Msoil_nh4_a = NH4.concentration_soil_usz_before * SOIL_PLANT.ro * PZ.Ab * WFR.thusz[time] * 1000 + NH4.concentration_soil_sz_before * SOIL_PLANT.ro * PZ.Ab * WFR.thsz[time] * 1000
            Msoil_nh4 = sum(NH4.concentration_soil_usz[time]) * SOIL_PLANT.ro * PZ.Ab * WFR.thusz[time] * 1000 / USZ.m_usz + sum(NH4.concentration_soil_sz[time]) * SOIL_PLANT.ro * PZ.Ab * WFR.thsz[
                time] * 1000 / SZ.m_sz
            NH4.mass_soil.append(Msoil_nh4 - Msoil_nh4_a)

            MRx_nh4 = - (sum(NH4.reaction_rate_usz[time]) * PZ.Ab * WFR.tteta_usz[time] * WFR.thusz[time] * 1000 / USZ.m_usz + sum(NH4.reaction_rate_sz[time]) * PZ.Ab *
                         WFR.tteta_sz[time] * WFR.thsz[time] * 1000 / SZ.m_sz)
            if time == 0:
                NH4.mass_accumulated_reaction.append(MRx_nh4)
            else:
                NH4.mass_accumulated_reaction.append(MRx_nh4 + NH4.mass_accumulated_reaction[-1])

            Min_nh4 = WFR.tQin[time] * NH4.concentration_inflow[time] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                NH4.mass_accumulated_inflow.append(Min_nh4)
            else:
                NH4.mass_accumulated_inflow.append(Min_nh4 + NH4.mass_accumulated_inflow[-1])

            Mover_nh4 = WFR.tQover[time] * NH4.concentration_pz[time] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                NH4.mass_accumulated_overflow.append(Mover_nh4)
            else:
                NH4.mass_accumulated_overflow.append(Mover_nh4 + NH4.mass_accumulated_overflow[-1])

            Mpipe_nh4 = WFR.tQpipe[time] * NH4.concentration_sz[time][SZ.m_sz - 1] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                NH4.mass_accumulated_pipe_outflow.append(Mpipe_nh4)
            else:
                NH4.mass_accumulated_pipe_outflow.append(Mpipe_nh4 + NH4.mass_accumulated_pipe_outflow[-1])

            Minfsz_nh4 = WFR.tQinfsz[time] * NH4.concentration_sz[time][SZ.m_sz - 1] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                NH4.mass_accumulated_infiltrated_to_SZ.append(Minfsz_nh4)
            else:
                NH4.mass_accumulated_infiltrated_to_SZ.append(Minfsz_nh4 + NH4.mass_accumulated_infiltrated_to_SZ[-1])

            Met_nh4 = WFR.tQet[time] * NH4.concentration_usz_layers_now[0] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                NH4.mass_accumulated_evapotranspiration.append(Met_nh4)
            else:
                NH4.mass_accumulated_evapotranspiration.append(Met_nh4 + NH4.mass_accumulated_evapotranspiration[-1])

            Mpz_nh4 = WFR.thpEND[time] * PZ.Ab * 1000 * NH4.concentration_pz[time]
            NH4.mass_PZ.append(Mpz_nh4)

            Mstor_nh4_mb = NH4.mass_accumulated_inflow[-1] - NH4.mass_PZ[-1] - NH4.mass_accumulated_overflow[-1] - NH4.mass_accumulated_pipe_outflow[-1] - \
                           NH4.mass_accumulated_infiltrated_to_SZ[-1] - NH4.mass_accumulated_evapotranspiration[-1] - NH4.mass_soil[-1] - NH4.mass_accumulated_reaction[-1]
            NH4.mass_balance.append(Mstor_nh4_mb)

            #         delta_ast_nh4_usz = 0
            #         delta_ast_nh4_sz = 0

            #         delta_ast_nh4_usz = (Mstor_nh4_mb - Mstor_nh4_ast)/(GENERAL_PARAMETERS.n*(PZ.Ab*WFR.tteta_usz[time]*WFR.thusz[time]*1000/USZ.m_usz))
            #         delta_ast_nh4_sz = (Mstor_nh4_mb - Mstor_nh4_ast)/(GENERAL_PARAMETERS.n*(PZ.Ab*WFR.tteta_sz[time]*WFR.thsz[time]*1000/SZ.m_sz))

            delta_ast_nh4_usz = (Mstor_nh4_mb - Mstor_nh4_ast) / (PZ.Ab * WFR.tteta_usz[time] * WFR.thusz[time] * 1000 / USZ.m_usz)
            NH4.concentration_usz_layers_now[1] = NH4.concentration_usz_layers_now[1] + delta_ast_nh4_usz
            if NH4.concentration_usz_layers_now[1] > 0:
                NH4.concentration_usz_layers_now[1] = NH4.concentration_usz_layers_now[1]
            else:
                NH4.concentration_usz_layers_now[1] = 0
            #         sum_cusz_nh4 = sum(NH4.concentration_usz[time])
            #         sum_csz_nh4 = sum(NH4.concentration_sz[time])
            #         print('time: ', time, 'Mstor_trace: ', Mstor_trace, ', Mstor_mb: ', Mstor_mb, ', dif: ', (Mstor_mb - Mstor_trace))
            #         print('time: ', time, ', delta_usz: ', delta_trace_nh4_usz, ', delta_sz: ', delta_trace_nh4_sz)
            #         print('teta_sm: ', WFR.tteta_usz[time], ', husz: ', WFR.thusz[time], ', teta_b: ', WFR.tteta_sz[time], ', hsz: ', WFR.thsz[time])
            #         print('sum_cusz_nh4: ', sum_cusz_nh4, ', sum_csz_nh4: ', sum_csz_nh4)
            #         print('Qin: ', WFR.tQin[time], ',hpEND_i: ', WFR.thpEND[time], ', hpEND_i-1: ', WFR.thpEND[time-1], ',Qpipe: ', WFR.tQpipe[time], ', Qover: ', WFR.tQover[time], ', Qinfsz: ', WFR.tQinfsz[time], ', Qet: ', WFR.tQet[time])
            #         print('Cin: ', NH4.concentration_inflow[time], ', Cp_i: ', NH4.concentration_pz[time], ', Cp_i-1: ', NH4.concentration_pz[time-1], ', Csz_ultimo: ', c_pipe)
            #         input()
            #
            #         for usz_layer in range(USZ.m_usz):
            #             NH4.concentration_usz_layers_now[usz_layer] = NH4.concentration_usz_layers_now[usz_layer] + delta_ast_nh4_usz
            #         for sz_layer in range(SZ.m_sz):
            #             cj_i1_nh4[sz_layer] = cj_i1_nh4[sz_layer] + delta_ast_nh4_sz

            ### Nitrate

            Mstor_no3_ast = sum(NO3.concentration_usz[time]) * PZ.Ab * WFR.tteta_usz[time] * WFR.thusz[time] * 1000 / USZ.m_usz + sum(NO3.concentration_sz[time]) * PZ.Ab * \
                            WFR.tteta_sz[time] * WFR.thsz[time] * 1000 / SZ.m_sz
            NO3.mass_stormwater.append(Mstor_no3_ast)

            Msoil_no3_a = 0
            Msoil_no3 = 0
            NO3.mass_soil.append(0)

            MRx_no3 = - (sum(NO3.reaction_rate_usz[time]) * PZ.Ab * WFR.tteta_usz[time] * WFR.thusz[time] * 1000 / USZ.m_usz + sum(NO3.reaction_rate_sz[time]) * PZ.Ab *
                         WFR.tteta_sz[time] * WFR.thsz[time] * 1000 / SZ.m_sz)
            if time == 0:
                NO3.mass_accumulated_reaction.append(MRx_no3)
            else:
                NO3.mass_accumulated_reaction.append(MRx_no3 + NO3.mass_accumulated_reaction[-1])

            Min_no3 = WFR.tQin[time] * NO3.concentration_inflow[time] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                NO3.mass_accumulated_inflow.append(Min_no3)
            else:
                NO3.mass_accumulated_inflow.append(Min_no3 + NO3.mass_accumulated_inflow[-1])

            Mover_no3 = WFR.tQover[time] * NO3.concentration_pz[time] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                NO3.mass_accumulated_overflow.append(Mover_no3)
            else:
                NO3.mass_accumulated_overflow.append(Mover_no3 + NO3.mass_accumulated_overflow[-1])

            Mpipe_no3 = WFR.tQpipe[time] * NO3.concentration_sz[time][SZ.m_sz - 1] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                NO3.mass_accumulated_pipe_outflow.append(Mpipe_no3)
            else:
                NO3.mass_accumulated_pipe_outflow.append(Mpipe_no3 + NO3.mass_accumulated_pipe_outflow[-1])

            Minfsz_no3 = WFR.tQinfsz[time] * NO3.concentration_sz[time][SZ.m_sz - 1] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                NO3.mass_accumulated_infiltrated_to_SZ.append(Minfsz_no3)
            else:
                NO3.mass_accumulated_infiltrated_to_SZ.append(Minfsz_no3 + NO3.mass_accumulated_infiltrated_to_SZ[-1])

            Met_no3 = WFR.tQet[time] * NO3.concentration_usz_layers_now[0] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                NO3.mass_accumulated_evapotranspiration.append(Met_no3)
            else:
                NO3.mass_accumulated_evapotranspiration.append(Met_no3 + NO3.mass_accumulated_evapotranspiration[-1])

            Mpz_no3 = WFR.thpEND[time] * PZ.Ab * 1000 * NO3.concentration_pz[time]
            NO3.mass_PZ.append(Mpz_no3)

            Mstor_no3_mb = NO3.mass_accumulated_inflow[-1] - NO3.mass_PZ[-1] - NO3.mass_accumulated_overflow[-1] - NO3.mass_accumulated_pipe_outflow[-1] - \
                           NO3.mass_accumulated_infiltrated_to_SZ[-1] - NO3.mass_accumulated_evapotranspiration[-1] - NO3.mass_soil[-1] - NO3.mass_accumulated_reaction[-1]
            NO3.mass_balance.append(Mstor_no3_mb)

            delta_ast_no3_usz = (Mstor_no3_mb - Mstor_no3_ast) / (PZ.Ab * WFR.tteta_usz[time] * WFR.thusz[time] * 1000 / USZ.m_usz)
            NO3.concentration_usz_layers_now[1] = NO3.concentration_usz_layers_now[1] + delta_ast_no3_usz
            if NO3.concentration_usz_layers_now[1] > 0:
                NO3.concentration_usz_layers_now[1] = NO3.concentration_usz_layers_now[1]
            else:
                NO3.concentration_usz_layers_now[1] = 0

                ### DOC

            Mstor_doc_ast = sum(DOC.concentration_usz[time]) * PZ.Ab * WFR.tteta_usz[time] * WFR.thusz[time] * 1000 / USZ.m_usz + sum(DOC.concentration_sz[time]) * PZ.Ab * \
                            WFR.tteta_sz[time] * WFR.thsz[time] * 1000 / SZ.m_sz
            DOC.mass_stormwater.append(Mstor_doc_ast)

            Msoil_doc_a = DOC.concentration_soil_usz_before * SOIL_PLANT.ro * PZ.Ab * WFR.thusz[time] * 1000 + DOC.concentration_soil_sz_before * SOIL_PLANT.ro * PZ.Ab * WFR.thsz[time] * 1000
            Msoil_doc = sum(DOC.concentration_soil_usz[time]) * SOIL_PLANT.ro * PZ.Ab * WFR.thusz[time] * 1000 / USZ.m_usz + sum(DOC.concentration_soil_sz[time]) * SOIL_PLANT.ro * PZ.Ab * WFR.thsz[
                time] * 1000 / SZ.m_sz
            DOC.mass_soil.append(Msoil_doc - Msoil_doc_a)

            MRx_doc = - (sum(DOC.reaction_rate_usz[time]) * PZ.Ab * WFR.tteta_usz[time] * WFR.thusz[time] * 1000 / USZ.m_usz + sum(DOC.reaction_rate_sz[time]) * PZ.Ab *
                         WFR.tteta_sz[time] * WFR.thsz[time] * 1000 / SZ.m_sz)
            if time == 0:
                DOC.mass_accumulated_reaction.append(MRx_doc)
            else:
                DOC.mass_accumulated_reaction.append(MRx_doc + DOC.mass_accumulated_reaction[-1])

            Min_doc = WFR.tQin[time] * DOC.concentration_inflow[time] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                DOC.mass_accumulated_inflow.append(Min_doc)
            else:
                DOC.mass_accumulated_inflow.append(Min_doc + DOC.mass_accumulated_inflow[-1])

            Mover_doc = WFR.tQover[time] * DOC.concentration_pz[time] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                DOC.mass_accumulated_overflow.append(Mover_doc)
            else:
                DOC.mass_accumulated_overflow.append(Mover_doc + DOC.mass_accumulated_overflow[-1])

            Mpipe_doc = WFR.tQpipe[time] * DOC.concentration_sz[time][SZ.m_sz - 1] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                DOC.mass_accumulated_pipe_outflow.append(Mpipe_doc)
            else:
                DOC.mass_accumulated_pipe_outflow.append(Mpipe_doc + DOC.mass_accumulated_pipe_outflow[-1])

            Minfsz_doc = WFR.tQinfsz[time] * DOC.concentration_sz[time][SZ.m_sz - 1] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                DOC.mass_accumulated_infiltrated_to_SZ.append(Minfsz_doc)
            else:
                DOC.mass_accumulated_infiltrated_to_SZ.append(Minfsz_doc + DOC.mass_accumulated_infiltrated_to_SZ[-1])

            Met_doc = WFR.tQet[time] * DOC.concentration_usz_layers_now[0] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                DOC.mass_accumulated_evapotranspiration.append(Met_doc)
            else:
                DOC.mass_accumulated_evapotranspiration.append(Met_doc + DOC.mass_accumulated_evapotranspiration[-1])

            Mpz_doc = WFR.thpEND[time] * PZ.Ab * 1000 * DOC.concentration_pz[time]
            DOC.mass_PZ.append(Mpz_doc)

            Mstor_doc_mb = DOC.mass_accumulated_inflow[-1] - DOC.mass_PZ[-1] - DOC.mass_accumulated_overflow[-1] - DOC.mass_accumulated_pipe_outflow[-1] - \
                           DOC.mass_accumulated_infiltrated_to_SZ[-1] - DOC.mass_accumulated_evapotranspiration[-1] - DOC.mass_soil[-1] - DOC.mass_accumulated_reaction[-1]
            DOC.mass_balance.append(Mstor_doc_mb)

            delta_ast_doc_usz = (Mstor_doc_mb - Mstor_doc_ast) / (PZ.Ab * WFR.tteta_usz[time] * WFR.thusz[time] * 1000 / USZ.m_usz)
            DOC.concentration_usz_layers_now[1] = DOC.concentration_usz_layers_now[1] + delta_ast_doc_usz
            if DOC.concentration_usz_layers_now[1] > 0:
                DOC.concentration_usz_layers_now[1] = DOC.concentration_usz_layers_now[1]
            else:
                DOC.concentration_usz_layers_now[1] = 0


        else:  # if GENERAL_PARAMETERS.hpipe == 0

            ### Oxygen

            Mstor_o2_ast = sum(O2.concentration_usz[time]) * PZ.Ab * WFR.tteta_usz[time] * WFR.thusz[time] * 1000 / USZ.m_usz
            O2.mass_stormwater.append(Mstor_o2_ast)

            Msoil_o2_a = 0
            Msoil_o2 = 0
            O2.mass_soil.append(0)

            MRx_o2 = - (sum(O2.reaction_rate_usz[time]) * PZ.Ab * WFR.tteta_usz[time] * WFR.thusz[time] * 1000 / USZ.m_usz)
            if time == 0:
                O2.mass_accumulated_reaction.append(MRx_o2)
            else:
                O2.mass_accumulated_reaction.append(MRx_o2 + O2.mass_accumulated_reaction[-1])

            Min_o2 = WFR.tQin[time] * O2.concentration_inflow[time] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                O2.mass_accumulated_inflow.append(Min_o2)
            else:
                O2.mass_accumulated_inflow.append(Min_o2 + O2.mass_accumulated_inflow[-1])

            Mover_o2 = WFR.tQover[time] * O2.concentration_pz[time] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                O2.mass_accumulated_overflow.append(Mover_o2)
            else:
                O2.mass_accumulated_overflow.append(Mover_o2 + O2.mass_accumulated_overflow[-1])

            Mpipe_o2 = WFR.tQpipe[time] * O2.concentration_usz[time][USZ.m_usz - 1] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                O2.mass_accumulated_pipe_outflow.append(Mpipe_o2)
            else:
                O2.mass_accumulated_pipe_outflow.append(Mpipe_o2 + O2.mass_accumulated_pipe_outflow[-1])

            Minfsz_o2 = WFR.tQinfsz[time] * O2.concentration_usz[time][USZ.m_usz - 1] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                O2.mass_accumulated_infiltrated_to_SZ.append(Minfsz_o2)
            else:
                O2.mass_accumulated_infiltrated_to_SZ.append(Minfsz_o2 + O2.mass_accumulated_infiltrated_to_SZ[-1])

            Met_o2 = WFR.tQet[time] * O2.concentration_usz_layers_now[0] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                O2.mass_accumulated_evapotranspiration.append(Met_o2)
            else:
                O2.mass_accumulated_evapotranspiration.append(Met_o2 + O2.mass_accumulated_evapotranspiration[-1])

            Mpz_o2 = WFR.thpEND[time] * PZ.Ab * 1000 * O2.concentration_pz[time]
            O2.mass_PZ.append(Mpz_o2)

            Mstor_o2_mb = O2.mass_accumulated_inflow[-1] - O2.mass_PZ[-1] - O2.mass_accumulated_overflow[-1] - O2.mass_accumulated_pipe_outflow[-1] - O2.mass_accumulated_infiltrated_to_SZ[
                -1] - O2.mass_accumulated_evapotranspiration[-1] - O2.mass_soil[-1] - O2.mass_accumulated_reaction[-1]
            O2.mass_balance.append(Mstor_o2_mb)

            delta_ast_o2_usz = (Mstor_o2_mb - Mstor_o2_ast) / (PZ.Ab * WFR.tteta_usz[time] * WFR.thusz[time] * 1000 / USZ.m_usz)
            O2.concentration_usz_layers_now[1] = O2.concentration_usz_layers_now[1] + delta_ast_o2_usz
            if O2.concentration_usz_layers_now[1] > 0:
                O2.concentration_usz_layers_now[1] = O2.concentration_usz_layers_now[1]
            else:
                O2.concentration_usz_layers_now[1] = 0

            ### Amonia

            Mstor_nh4_ast = sum(NH4.concentration_usz[time]) * PZ.Ab * WFR.tteta_usz[time] * WFR.thusz[time] * 1000 / USZ.m_usz
            NH4.mass_stormwater.append(Mstor_nh4_ast)

            Msoil_nh4_a = NH4.concentration_soil_usz_before * SOIL_PLANT.ro * PZ.Ab * WFR.thusz[time] * 1000
            Msoil_nh4 = sum(NH4.concentration_soil_usz[time]) * SOIL_PLANT.ro * PZ.Ab * WFR.thusz[time] * 1000 / USZ.m_usz
            NH4.mass_soil.append(Msoil_nh4 - Msoil_nh4_a)

            MRx_nh4 = - (sum(NH4.reaction_rate_usz[time]) * PZ.Ab * WFR.tteta_usz[time] * WFR.thusz[time] * 1000 / USZ.m_usz)
            if time == 0:
                NH4.mass_accumulated_reaction.append(MRx_nh4)
            else:
                NH4.mass_accumulated_reaction.append(MRx_nh4 + NH4.mass_accumulated_reaction[-1])

            Min_nh4 = WFR.tQin[time] * NH4.concentration_inflow[time] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                NH4.mass_accumulated_inflow.append(Min_nh4)
            else:
                NH4.mass_accumulated_inflow.append(Min_nh4 + NH4.mass_accumulated_inflow[-1])

            Mover_nh4 = WFR.tQover[time] * NH4.concentration_pz[time] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                NH4.mass_accumulated_overflow.append(Mover_nh4)
            else:
                NH4.mass_accumulated_overflow.append(Mover_nh4 + NH4.mass_accumulated_overflow[-1])

            Mpipe_nh4 = WFR.tQpipe[time] * NH4.concentration_usz[time][USZ.m_usz - 1] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                NH4.mass_accumulated_pipe_outflow.append(Mpipe_nh4)
            else:
                NH4.mass_accumulated_pipe_outflow.append(Mpipe_nh4 + NH4.mass_accumulated_pipe_outflow[-1])

            Minfsz_nh4 = WFR.tQinfsz[time] * NH4.concentration_usz[time][USZ.m_usz - 1] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                NH4.mass_accumulated_infiltrated_to_SZ.append(Minfsz_nh4)
            else:
                NH4.mass_accumulated_infiltrated_to_SZ.append(Minfsz_nh4 + NH4.mass_accumulated_infiltrated_to_SZ[-1])

            Met_nh4 = WFR.tQet[time] * NH4.concentration_usz_layers_now[0] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                NH4.mass_accumulated_evapotranspiration.append(Met_nh4)
            else:
                NH4.mass_accumulated_evapotranspiration.append(Met_nh4 + NH4.mass_accumulated_evapotranspiration[-1])

            Mpz_nh4 = WFR.thpEND[time] * PZ.Ab * 1000 * NH4.concentration_pz[time]
            NH4.mass_PZ.append(Mpz_nh4)

            Mstor_nh4_mb = NH4.mass_accumulated_inflow[-1] - NH4.mass_PZ[-1] - NH4.mass_accumulated_overflow[-1] - NH4.mass_accumulated_pipe_outflow[-1] - \
                           NH4.mass_accumulated_infiltrated_to_SZ[-1] - NH4.mass_accumulated_evapotranspiration[-1] - NH4.mass_soil[-1] - NH4.mass_accumulated_reaction[-1]
            NH4.mass_balance.append(Mstor_nh4_mb)

            delta_ast_nh4_usz = (Mstor_nh4_mb - Mstor_nh4_ast) / (PZ.Ab * WFR.tteta_usz[time] * WFR.thusz[time] * 1000 / USZ.m_usz)
            NH4.concentration_usz_layers_now[1] = NH4.concentration_usz_layers_now[1] + delta_ast_nh4_usz
            if NH4.concentration_usz_layers_now[1] > 0:
                NH4.concentration_usz_layers_now[1] = NH4.concentration_usz_layers_now[1]
            else:
                NH4.concentration_usz_layers_now[1] = 0

            ### Nitrate

            Mstor_no3_ast = sum(NO3.concentration_usz[time]) * PZ.Ab * WFR.tteta_usz[time] * WFR.thusz[time] * 1000 / USZ.m_usz
            NO3.mass_stormwater.append(Mstor_no3_ast)

            Msoil_no3_a = 0
            Msoil_no3 = 0
            NO3.mass_soil.append(0)

            MRx_no3 = - (sum(NO3.reaction_rate_usz[time]) * PZ.Ab * WFR.tteta_usz[time] * WFR.thusz[time] * 1000 / USZ.m_usz)
            if time == 0:
                NO3.mass_accumulated_reaction.append(MRx_no3)
            else:
                NO3.mass_accumulated_reaction.append(MRx_no3 + NO3.mass_accumulated_reaction[-1])

            Min_no3 = WFR.tQin[time] * NO3.concentration_inflow[time] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                NO3.mass_accumulated_inflow.append(Min_no3)
            else:
                NO3.mass_accumulated_inflow.append(Min_no3 + NO3.mass_accumulated_inflow[-1])

            Mover_no3 = WFR.tQover[time] * NO3.concentration_pz[time] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                NO3.mass_accumulated_overflow.append(Mover_no3)
            else:
                NO3.mass_accumulated_overflow.append(Mover_no3 + NO3.mass_accumulated_overflow[-1])

            Mpipe_no3 = WFR.tQpipe[time] * NO3.concentration_usz[time][USZ.m_usz - 1] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                NO3.mass_accumulated_pipe_outflow.append(Mpipe_no3)
            else:
                NO3.mass_accumulated_pipe_outflow.append(Mpipe_no3 + NO3.mass_accumulated_pipe_outflow[-1])

            Minfsz_no3 = WFR.tQinfsz[time] * NO3.concentration_usz[time][USZ.m_usz - 1] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                NO3.mass_accumulated_infiltrated_to_SZ.append(Minfsz_no3)
            else:
                NO3.mass_accumulated_infiltrated_to_SZ.append(Minfsz_no3 + NO3.mass_accumulated_infiltrated_to_SZ[-1])

            Met_no3 = WFR.tQet[time] * NO3.concentration_usz_layers_now[0] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                NO3.mass_accumulated_evapotranspiration.append(Met_no3)
            else:
                NO3.mass_accumulated_evapotranspiration.append(Met_no3 + NO3.mass_accumulated_evapotranspiration[-1])

            Mpz_no3 = WFR.thpEND[time] * PZ.Ab * 1000 * NO3.concentration_pz[time]
            NO3.mass_PZ.append(Mpz_no3)

            Mstor_no3_mb = NO3.mass_accumulated_inflow[-1] - NO3.mass_PZ[-1] - NO3.mass_accumulated_overflow[-1] - NO3.mass_accumulated_pipe_outflow[-1] - \
                           NO3.mass_accumulated_infiltrated_to_SZ[-1] - NO3.mass_accumulated_evapotranspiration[-1] - NO3.mass_soil[-1] - NO3.mass_accumulated_reaction[-1]
            NO3.mass_balance.append(Mstor_no3_mb)

            delta_ast_no3_usz = (Mstor_no3_mb - Mstor_no3_ast) / (PZ.Ab * WFR.tteta_usz[time] * WFR.thusz[time] * 1000 / USZ.m_usz)
            NO3.concentration_usz_layers_now[1] = NO3.concentration_usz_layers_now[1] + delta_ast_no3_usz
            if NO3.concentration_usz_layers_now[1] > 0:
                NO3.concentration_usz_layers_now[1] = NO3.concentration_usz_layers_now[1]
            else:
                NO3.concentration_usz_layers_now[1] = 0

                ### DOC

            Mstor_doc_ast = sum(DOC.concentration_usz[time]) * PZ.Ab * WFR.tteta_usz[time] * WFR.thusz[time] * 1000 / USZ.m_usz
            DOC.mass_stormwater.append(Mstor_doc_ast)

            Msoil_doc_a = DOC.concentration_soil_usz_before * SOIL_PLANT.ro * PZ.Ab * WFR.thusz[time] * 1000
            Msoil_doc = sum(DOC.concentration_soil_usz[time]) * SOIL_PLANT.ro * PZ.Ab * WFR.thusz[time] * 1000 / USZ.m_usz
            DOC.mass_soil.append(Msoil_doc - Msoil_doc_a)

            MRx_doc = - (sum(DOC.reaction_rate_usz[time]) * PZ.Ab * WFR.tteta_usz[time] * WFR.thusz[time] * 1000 / USZ.m_usz)
            if time == 0:
                DOC.mass_accumulated_reaction.append(MRx_doc)
            else:
                DOC.mass_accumulated_reaction.append(MRx_doc + DOC.mass_accumulated_reaction[-1])

            Min_doc = WFR.tQin[time] * DOC.concentration_inflow[time] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                DOC.mass_accumulated_inflow.append(Min_doc)
            else:
                DOC.mass_accumulated_inflow.append(Min_doc + DOC.mass_accumulated_inflow[-1])

            Mover_doc = WFR.tQover[time] * DOC.concentration_pz[time] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                DOC.mass_accumulated_overflow.append(Mover_doc)
            else:
                DOC.mass_accumulated_overflow.append(Mover_doc + DOC.mass_accumulated_overflow[-1])

            Mpipe_doc = WFR.tQpipe[time] * DOC.concentration_usz[time][USZ.m_usz - 1] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                DOC.mass_accumulated_pipe_outflow.append(Mpipe_doc)
            else:
                DOC.mass_accumulated_pipe_outflow.append(Mpipe_doc + DOC.mass_accumulated_pipe_outflow[-1])

            Minfsz_doc = WFR.tQinfsz[time] * DOC.concentration_usz[time][USZ.m_usz - 1] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                DOC.mass_accumulated_infiltrated_to_SZ.append(Minfsz_doc)
            else:
                DOC.mass_accumulated_infiltrated_to_SZ.append(Minfsz_doc + DOC.mass_accumulated_infiltrated_to_SZ[-1])

            Met_doc = WFR.tQet[time] * DOC.concentration_usz_layers_now[0] * GENERAL_PARAMETERS.dt * 1000
            if time == 0:
                DOC.mass_accumulated_evapotranspiration.append(Met_doc)
            else:
                DOC.mass_accumulated_evapotranspiration.append(Met_doc + DOC.mass_accumulated_evapotranspiration[-1])

            Mpz_doc = WFR.thpEND[time] * PZ.Ab * 1000 * DOC.concentration_pz[time]
            DOC.mass_PZ.append(Mpz_doc)

            Mstor_doc_mb = DOC.mass_accumulated_inflow[-1] - DOC.mass_PZ[-1] - DOC.mass_accumulated_overflow[-1] - DOC.mass_accumulated_pipe_outflow[-1] - \
                           DOC.mass_accumulated_infiltrated_to_SZ[-1] - DOC.mass_accumulated_evapotranspiration[-1] - DOC.mass_soil[-1] - DOC.mass_accumulated_reaction[-1]
            DOC.mass_balance.append(Mstor_doc_mb)

            delta_ast_doc_usz = (Mstor_doc_mb - Mstor_doc_ast) / (PZ.Ab * WFR.tteta_usz[time] * WFR.thusz[time] * 1000 / USZ.m_usz)
            DOC.concentration_usz_layers_now[1] = DOC.concentration_usz_layers_now[1] + delta_ast_doc_usz
            if DOC.concentration_usz_layers_now[1] > 0:
                DOC.concentration_usz_layers_now[1] = DOC.concentration_usz_layers_now[1]
            else:
                DOC.concentration_usz_layers_now[1] = 0

        ## adding layers of USZ in time
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
        if GENERAL_PARAMETERS.hpipe > 0:
            O2.concentration_sz.append(cj_i1_o2)
            O2.reaction_rate_sz.append(Rxj_o2)

            NH4.concentration_sz.append(cj_i1_nh4)
            NH4.concentration_soil_sz.append(csi_sz_nh4)
            NH4.reaction_rate_sz.append(Rxj_nh4)

            NO3.concentration_sz.append(cj_i1_no3)
            NO3.reaction_rate_sz.append(Rxj_no3)

            DOC.concentration_sz.append(cj_i1_doc)
            DOC.concentration_soil_sz.append(csi_sz_doc)
            DOC.reaction_rate_sz.append(Rxj_doc)

    # **5. Transforming in dataframe **
    ## Oxygen
    data_usz_o2 = pd.DataFrame(O2.concentration_usz)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz' + str(a)
        column_name.append(name)
        a = a + 1
    data_usz_o2.set_axis(column_name, axis='columns', inplace=True)

    # print(data_usz_o2)

    data_sz_o2 = pd.DataFrame(O2.concentration_sz)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz' + str(a)
        column_name.append(name)
        a = a + 1
    data_sz_o2.set_axis(column_name, axis='columns', inplace=True)

    data_rx_usz_o2 = pd.DataFrame(O2.reaction_rate_usz)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_usz_o2.set_axis(column_name, axis='columns', inplace=True)

    data_rx_sz_o2 = pd.DataFrame(O2.reaction_rate_sz)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_sz_o2.set_axis(column_name, axis='columns', inplace=True)

    frames = [data_usz_o2, data_sz_o2, data_rx_usz_o2, data_rx_sz_o2]
    data_o2 = pd.concat((frames), axis=1)

    O2.concentration_pz.append(0)
    data_o2['pz'] = O2.concentration_pz

    indice_n = list(range(len(O2.concentration_pz)))
    O2.concentration_inflow.append(0)
    c_in = O2.concentration_inflow[:len(indice_n)]
    # print('len c_in:', len(c_in), ', len_indice:', len(WFR.indice))
    data_o2['c_in'] = c_in

    data_o2['t'] = indice_n

    ## Amonia
    data_usz_nh4 = pd.DataFrame(NH4.concentration_usz)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz' + str(a)
        column_name.append(name)
        a = a + 1
    data_usz_nh4.set_axis(column_name, axis='columns', inplace=True)

    data_sz_nh4 = pd.DataFrame(NH4.concentration_sz)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz' + str(a)
        column_name.append(name)
        a = a + 1
    data_sz_nh4.set_axis(column_name, axis='columns', inplace=True)

    data_s_usz_nh4 = pd.DataFrame(NH4.concentration_soil_usz)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz_soil' + str(a)
        column_name.append(name)
        a = a + 1
    data_s_usz_nh4.set_axis(column_name, axis='columns', inplace=True)

    data_s_sz_nh4 = pd.DataFrame(NH4.concentration_soil_sz)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz_soil' + str(a)
        column_name.append(name)
        a = a + 1
    data_s_sz_nh4.set_axis(column_name, axis='columns', inplace=True)

    data_rx_usz_nh4 = pd.DataFrame(NH4.reaction_rate_usz)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_usz_nh4.set_axis(column_name, axis='columns', inplace=True)

    data_rx_sz_nh4 = pd.DataFrame(NH4.reaction_rate_sz)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_sz_nh4.set_axis(column_name, axis='columns', inplace=True)

    frames = [data_usz_nh4, data_sz_nh4, data_s_usz_nh4, data_s_sz_nh4, data_rx_usz_nh4, data_rx_sz_nh4]
    data_nh4 = pd.concat((frames), axis=1)
    NH4.concentration_pz.append(0)
    NH4.concentration_inflow.append(0)
    c_in = NH4.concentration_inflow[:len(indice_n)]
    # print('len c_in:', len(c_in), ', len_indice:', len(WFR.indice))
    data_nh4['c_in'] = c_in

    data_nh4['pz'] = NH4.concentration_pz
    # WFR.indice.append(len(WFR.indice)+1)
    data_nh4['t'] = indice_n

    ## Nitrate
    data_usz_no3 = pd.DataFrame(NO3.concentration_usz)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz' + str(a)
        column_name.append(name)
        a = a + 1
    data_usz_no3.set_axis(column_name, axis='columns', inplace=True)

    data_sz_no3 = pd.DataFrame(NO3.concentration_sz)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz' + str(a)
        column_name.append(name)
        a = a + 1
    data_sz_no3.set_axis(column_name, axis='columns', inplace=True)

    data_rx_usz_no3 = pd.DataFrame(NO3.reaction_rate_usz)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_usz_no3.set_axis(column_name, axis='columns', inplace=True)

    data_rx_sz_no3 = pd.DataFrame(NO3.reaction_rate_sz)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_sz_no3.set_axis(column_name, axis='columns', inplace=True)

    frames = [data_usz_no3, data_sz_no3, data_rx_usz_no3, data_rx_sz_no3]
    data_no3 = pd.concat((frames), axis=1)
    NO3.concentration_pz.append(0)
    data_no3['pz'] = NO3.concentration_pz

    NO3.concentration_inflow.append(0)
    c_in = NO3.concentration_inflow[:len(indice_n)]
    data_no3['c_in'] = c_in

    data_no3['t'] = indice_n

    ## DOC
    data_usz_doc = pd.DataFrame(DOC.concentration_usz)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz' + str(a)
        column_name.append(name)
        a = a + 1
    data_usz_doc.set_axis(column_name, axis='columns', inplace=True)

    data_sz_doc = pd.DataFrame(DOC.concentration_sz)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz' + str(a)
        column_name.append(name)
        a = a + 1
    data_sz_doc.set_axis(column_name, axis='columns', inplace=True)

    data_rx_usz_doc = pd.DataFrame(DOC.reaction_rate_usz)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_usz_doc.set_axis(column_name, axis='columns', inplace=True)

    data_rx_sz_doc = pd.DataFrame(DOC.reaction_rate_sz)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_sz_doc.set_axis(column_name, axis='columns', inplace=True)

    frames = [data_usz_doc, data_sz_doc, data_rx_usz_doc, data_rx_sz_doc]
    data_doc = pd.concat((frames), axis=1)
    DOC.concentration_pz.append(0)
    data_doc['pz'] = DOC.concentration_pz

    DOC.concentration_inflow.append(0)
    c_in = DOC.concentration_inflow[:len(indice_n)]
    data_doc['c_in'] = c_in

    data_doc['t'] = indice_n
    return data_o2, data_nh4, data_no3, data_doc

def run(interaction):
    SETUP_FILE = "parameters.ini"
    INFLOW_FILE = "concentration_inflow.csv"
    WF_INFLOW = parameters.WaterInflow("water_inflow.csv")
    WFR = parameters.WaterFlowResults(WF_INFLOW.tQrain, WF_INFLOW.tQin, WF_INFLOW.tEmax)
    GENERAL_PARAMETERS = parameters.GeneralParameters(SETUP_FILE)
    GENERAL_PARAMETERS.hpipe = hpipes[interaction]
    USZ = parameters.UnsaturatedZone(SETUP_FILE, GENERAL_PARAMETERS.L, GENERAL_PARAMETERS.hpipe,
                                     GENERAL_PARAMETERS.dz)
    PZ = parameters.PondingZone(SETUP_FILE)
    SZ = parameters.SaturatedZone(SETUP_FILE, GENERAL_PARAMETERS.n, USZ.m_usz)
    SOIL_PLANT = parameters.SoilPlant(SETUP_FILE, USZ.nusz_ini)
    NH4 = parameters.Ammonia(USZ.m_usz, SZ.m_sz, SETUP_FILE, INFLOW_FILE)
    NO3 = parameters.Nitrate(USZ.m_usz, SZ.m_sz, SETUP_FILE, INFLOW_FILE)
    O2 = parameters.Oxygen(USZ.m_usz, SZ.m_sz, SETUP_FILE, INFLOW_FILE)
    DOC = parameters.DissolvedOrganicCarbon(USZ.m_usz, SZ.m_sz, SETUP_FILE, INFLOW_FILE)
    water_flow_module(WFR, GENERAL_PARAMETERS, USZ, PZ, SZ)
    #WFR.water_balance(GENERAL_PARAMETERS.dt)
    print("hpipe:", GENERAL_PARAMETERS.hpipe)
    data_o2, data_nh4, data_no3, data_doc = water_quality_module(WFR, GENERAL_PARAMETERS, USZ, PZ,
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

if __name__ == '__main__':
    hpipes = [0, 0.1, 0.2, 0.3, 0.4]
    results_paths = ["results/results_00/", "results/results_10/", "results/results_20/", "results/results_30/",
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

