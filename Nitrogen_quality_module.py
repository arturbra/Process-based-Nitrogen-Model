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
    for t in range(len(WFR.tQrain)):
        Qin = WFR.tQin[t]
        Qrain = WFR.tQrain[t]
        Emax = WFR.tEmax[t]

        # PZ#
        Qover = PZ.f_weir_overflow(hpEND, GENERAL_PARAMETERS.dt, Qin, Qrain)
        height_pz = max(hpEND + GENERAL_PARAMETERS.dt / PZ.Ab * (WFR.tQrain[t] + Qin - Qover), 0)  # beginning
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

        if t == 0:
            husz_a = USZ.husz_ini
            nusz_a = USZ.nusz_ini
            s_a = USZ.sw
        else:
            husz_a = WFR.thusz[t - 1]
            nusz_a = WFR.tnusz[t - 1]
            s_a = WFR.ts[t - 1]

        s = max(min(1.0, (s_a * husz_a * nusz_a * USZ.A + GENERAL_PARAMETERS.dt * (Qpf + Qhc - Qfs - Qet1)) / (
                    USZ.A * husz * nusz)), USZ.sh)

        # save all results to WFR
        WFR.tt.append(t)
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
    for t in range(len(WFR.indice) - 1):
        height_pz = WFR.thpEND[t]
        if height_pz < 0.001:
            height_pz = 0

        if t == 0:
            height_pz_before = 0
            Qorif = 0
        else:
            height_pz_before = WFR.thpEND[t - 1]
            Qorif = WFR.tQpipe[t]

        Qinf_sz = WFR.tQinfsz[t]
        teta_sm_i = WFR.tteta_usz[t]
        teta_b_i = WFR.tteta_sz[t]

        if t < (len(WFR.indice) - 1):
            teta_sm_iplus1 = WFR.tteta_usz[t + 1]
            teta_b_iplus1 = WFR.tteta_sz[t + 1]
        else:
            teta_sm_iplus1 = WFR.tteta_usz[t]
            teta_b_iplus1 = WFR.tteta_sz[t]

        Rxi_p_o2 = O2.f_reaction_pz()
        Rxi_p_nh4 = NH4.f_reaction_pz()
        Rxi_p_no3 = NO3.f_reaction_pz()
        Rxi_p_doc = DOC.f_reaction_pz()

        # if t < 20:
        # print('t: ', t, 'Rx_p_o2: ', Rxi_p_o2, 'Rx_p_nh4: ', Rxi_p_nh4, 'Rx_p_no3: ', Rxi_p_no3)

        if height_pz == 0:
            cpi_o2 = 0
            cpi_nh4 = 0
            cpi_no3 = 0
            cpi_doc = 0

        else:
            cpi_o2 = PZ.f_concentration(O2.concentration_inflow[t], WFR.tQin[t], O2.concentration_PZ_before, WFR.tQpf[t], WFR.tQover[t], Rxi_p_o2, height_pz, height_pz_before, PZ.Ab, GENERAL_PARAMETERS.dt)
            cpi_nh4 = PZ.f_concentration(NH4.concentration_inflow[t], WFR.tQin[t], NH4.concentration_PZ_before, WFR.tQpf[t], WFR.tQover[t], Rxi_p_nh4, height_pz, height_pz_before, PZ.Ab, GENERAL_PARAMETERS.dt)
            cpi_no3 = PZ.f_concentration(NO3.concentration_inflow[t], WFR.tQin[t], NO3.concentration_PZ_before, WFR.tQpf[t], WFR.tQover[t], Rxi_p_no3, height_pz, height_pz_before, PZ.Ab, GENERAL_PARAMETERS.dt)
            cpi_doc = PZ.f_concentration(DOC.concentration_inflow[t], WFR.tQin[t], DOC.concentration_PZ_before, WFR.tQpf[t], WFR.tQover[t], Rxi_p_doc, height_pz, height_pz_before, PZ.Ab, GENERAL_PARAMETERS.dt)
        #             print('NH4.concentration_inflow[t], WFR.tQin[t], NH4.concentration_PZ_before, WFR.tQpf[t], WFR.tQover[t], Rxi_p_nh4, height_pz, height_pz_before: ', NH4.concentration_inflow[t], WFR.tQin[t], NH4.concentration_PZ_before, WFR.tQpf[t], WFR.tQover[t], Rxi_p_nh4, height_pz, height_pz_before)
        #             print('cpi_nh4: ', cpi_nh4)

        if cpi_o2 < 0.00002:
            cpi_o2 = 0

        if cpi_nh4 < 0.00002:
            cpi_nh4 = 0

        if cpi_no3 < 0.00002:
            cpi_no3 = 0

        if cpi_doc < 0.00002:
            cpi_doc = 0


        O2.concentration_PZ.append(cpi_o2)
        O2.concentration_PZ_before = O2.concentration_PZ[-1]

        NH4.concentration_PZ.append(cpi_nh4)
        NH4.concentration_PZ_before = NH4.concentration_PZ[-1]

        NO3.concentration_PZ.append(cpi_no3)
        NO3.concentration_PZ_before = NO3.concentration_PZ[-1]

        DOC.concentration_PZ.append(cpi_doc)
        DOC.concentration_PZ_before = DOC.concentration_PZ[-1]

        # USZ
        cli_o2_list = O2.concentration_USZ[t].copy()
        # print('1_o2', cli_o2_list)
        cli_nh4_list = NH4.concentration_USZ[t].copy()
        # print('1_nh4', cli_nh4_list)
        cli_no3_list = NO3.concentration_USZ[t].copy()
        cli_doc_list = DOC.concentration_USZ[t].copy()

        # SZ
        if GENERAL_PARAMETERS.hpipe > 0:
            cji_o2_list = O2.concentration_SZ[t].copy()  # copiar lista
            # cji_o2_list = O2.concentration_SZ[t][:] #outro jeito de copiar a lista
            cji_nh4_list = NH4.concentration_SZ[t].copy()
            cji_no3_list = NO3.concentration_SZ[t].copy()
            cji_doc_list = DOC.concentration_SZ[t].copy()

        ####   USZ   ####
        if USZ.m_usz != 0:

            cl_i1_o2 = []
            Rxl_o2 = []

            cl_i1_nh4 = []
            csi_usz_nh4 = []
            Rxl_nh4 = []

            cl_i1_no3 = []
            Rxl_no3 = []

            cl_i1_doc = []
            csi_usz_doc = []
            Rxl_doc = []

            # Predictive step

            for l in range(USZ.m_usz):
                # print('l', l)

                cl_o2 = cli_o2_list[l]
                clminus1_o2 = cli_o2_list[l - 1]
                if l < (USZ.m_usz - 1):
                    clplus1_o2 = cli_o2_list[l + 1]
                else:
                    clplus1_o2 = 0

                cl_nh4 = cli_nh4_list[l]
                clminus1_nh4 = cli_nh4_list[l - 1]
                if l < (USZ.m_usz - 1):
                    clplus1_nh4 = cli_nh4_list[l + 1]
                else:
                    clplus1_nh4 = 0

                cl_no3 = cli_no3_list[l]
                clminus1_no3 = cli_no3_list[l - 1]
                if l < (USZ.m_usz - 1):
                    clplus1_no3 = cli_no3_list[l + 1]
                else:
                    clplus1_no3 = 0

                len(cli_doc_list)
                cl_doc = cli_doc_list[l]
                clminus1_doc = cli_doc_list[l - 1]

                if l < (USZ.m_usz - 1):
                    clplus1_doc = cli_doc_list[l + 1]
                else:
                    clplus1_doc = 0

                cs_o2 = 0

                # since we have added an initial value of cs = 0 in the NH4.concentration_soil_USZ list, when calling the index equal = 't' we are actually calling the value corresponding to t-1
                cs_nh4 = NH4.concentration_soil_USZ[t][l]
                cs_nh4_iplus1 = NH4.f_concentration_soil(cs_nh4, teta_sm_i, NH4.kads, cl_nh4, NH4.kdes, NH4.k_mb, SOIL_PLANT.ro, GENERAL_PARAMETERS.dt)
                if cs_nh4_iplus1 < 0.00000000000001:
                    cs_nh4_iplus1 = 0

                #                 sor = (teta_sm_i/SOIL_PLANT.ro)*NH4.kads*cl_nh4*GENERAL_PARAMETERS.dt
                #                 des =  NH4.kdes*NH4.concentration_soil_USZ_before*GENERAL_PARAMETERS.dt

                #                 print('t: ', t, ', l: ', l , ', cs_nh4_a: ', NH4.concentration_soil_USZ_before, 'teta_sm_i: ', teta_sm_i, 'NH4.kads: ', NH4.kads, 'cl_nh4: ', cl_nh4, 'NH4.kdes: ', NH4.kdes, 'NH4.k_mb: ', NH4.k_mb)
                #                 print('cs_nh4: ', cs_nh4)
                #                 input()

                #                 logger.debug(SOIL_PLANT.f' t: {t} \t l: {l} \t sor: {sor} \t des: {des} \t cs_nh4: {cs_nh4}')
                #                 logger.debug('----SOIL----')
                #                 logger.debug('cs_nh4_a: ', NH4.concentration_soil_USZ_before)
                #                 logger.debug('cs_nh4_a: ', NH4.concentration_soil_USZ_before, 'teta_sm_i: ', teta_sm_i, 'NH4.kads: ', NH4.kads, 'cl_nh4: ', cl_nh4, 'NH4.kdes: ', NH4.kdes, 'NH4.k_mb: ', NH4.k_mb)
                #                 logger.debug('cs_nh4: ', cs_nh4)
                #                 logger.debug('----WATER----')
                csi_usz_nh4.append(cs_nh4_iplus1)

                cs_no3 = 0

                DOC.concentration_soil_USZ_before = DOC.concentration_soil_USZ[t][l]
                cs_doc = DOC.f_concentration_soil(DOC.concentration_soil_USZ_before, teta_sm_i, DOC.kads, cl_doc, DOC.kdes, DOC.k_mb, SOIL_PLANT.ro, GENERAL_PARAMETERS.dt)
                csi_usz_doc.append(cs_doc)
                if cs_doc < 0.00000000000001:
                    cs_doc = 0

                UF_usz = []
                # WFR.tQet1[t] = 0
                alfa, beta = USZ.f_alfa_beta(l)
                UFi_usz = USZ.f_unit_flux(alfa, beta, WFR.tQpf[t], WFR.tQet1[t], WFR.tQfs[t], WFR.tQhc[t], Qorif, Qinf_sz, teta_sm_i, GENERAL_PARAMETERS.hpipe, PZ.Ab)
                #                 print(t, '  ', l)
                #                 print('UFi_usz: ', UFi_usz, ', alfa: ', alfa, ', beta: ', beta,  ', WFR.tQpf[t]: ', WFR.tQpf[t], ', WFR.tQet1[t]:', WFR.tQet1[t], ', WFR.tQfs[t]: ', WFR.tQfs[t], ', WFR.tQhc[t]: ', WFR.tQhc[t], ', Qorif: ', Qorif, ', Qinf_sz: ', Qinf_sz)
                #                 input()

                UF_usz.append(UFi_usz)

                Rxi_2_o2 = O2.f_reaction_usz(cl_o2, cl_nh4, GENERAL_PARAMETERS.k_nit) + O2.f_plant_uptake_usz(cl_o2, SOIL_PLANT.c_o2_root, teta_sm_i, SOIL_PLANT.root_fraction, SOIL_PLANT.lamda)
                Rxi_2_nh4 = NH4.f_reaction_usz(cl_nh4, GENERAL_PARAMETERS.k_nit) + NH4.f_plant_uptake_usz(cl_nh4, teta_sm_i, SOIL_PLANT.root_fraction)
                Rxi_2_no3 = NO3.f_reaction_usz(cl_nh4, GENERAL_PARAMETERS.k_nit) + NO3.f_plant_uptake_usz(cl_no3, teta_sm_i, SOIL_PLANT.root_fraction)
                Rxi_2_doc = DOC.f_reaction_usz(cl_doc)

                Rxl_o2.append(Rxi_2_o2 * (1 / teta_sm_iplus1) * GENERAL_PARAMETERS.dt)
                Rxl_nh4.append(Rxi_2_nh4 * (1 / teta_sm_iplus1) * GENERAL_PARAMETERS.dt)
                Rxl_no3.append(Rxi_2_no3 * (1 / teta_sm_iplus1) * GENERAL_PARAMETERS.dt)
                Rxl_doc.append(Rxi_2_doc * (1 / teta_sm_iplus1) * GENERAL_PARAMETERS.dt)

                #                 Rxl_o2.append(Rxi_2_o2*(GENERAL_PARAMETERS.dt/teta_sm_i))
                #                 Rxl_nh4.append(Rxi_2_nh4*(GENERAL_PARAMETERS.dt/teta_sm_i))
                #                 Rxl_no3.append(Rxi_2_no3*(GENERAL_PARAMETERS.dt/teta_sm_i))
                #                 Rxl_doc.append(Rxi_2_doc*(GENERAL_PARAMETERS.dt/teta_sm_i))
                #                 print('t: ', t, ', l: ', l, ', Rx: ', Rxi_2_nh4*(GENERAL_PARAMETERS.dt/teta_sm_i))
                #                 input()
                # if t < 20:
                # print('t: ', t, 'l: ', l, 'O2.reaction_rate_USZ: ', Rxi_2_o2, 'NH4.reaction_rate_USZ: ', Rxi_2_nh4, 'NO3.reaction_rate_USZ: ', Rxi_2_no3)

                ### Oxygen
                Peusz_o2 = USZ.f_peclet(UFi_usz, O2.D, GENERAL_PARAMETERS.dz)
                # print(Peusz_o2)

                if Peusz_o2 <= 2:
                    if l == 0:  # first cell
                        dc_o2 = clplus1_o2 - 2 * cl_o2 + cpi_o2
                        dc_dz_o2 = (clplus1_o2 - cpi_o2) / (2 * GENERAL_PARAMETERS.dz)


                    elif l == (USZ.m_usz - 1):  # last cell
                        dc_o2 = cl_o2 - 2 * cl_o2 + clminus1_o2
                        dc_dz_o2 = (cl_o2 - clminus1_o2) / GENERAL_PARAMETERS.dz

                    else:
                        dc_o2 = clplus1_o2 - 2 * cl_o2 + clminus1_o2
                        dc_dz_o2 = (clplus1_o2 - clminus1_o2) / (2 * GENERAL_PARAMETERS.dz)

                else:  # Peusz > 2
                    if l == 0:  # first cell
                        dc_o2 = clplus1_o2 - 2 * cl_o2 + cpi_o2
                        dc_dz_o2 = (cl_o2 - cpi_o2) / GENERAL_PARAMETERS.dz

                    elif l == (USZ.m_usz - 1):  # last cell
                        dc_o2 = cl_o2 - 2 * cl_o2 + clminus1_o2
                        dc_dz_o2 = (cl_o2 - clminus1_o2) / GENERAL_PARAMETERS.dz

                    else:
                        dc_o2 = clplus1_o2 - 2 * cl_o2 + clminus1_o2
                        dc_dz_o2 = (cl_o2 - clminus1_o2) / GENERAL_PARAMETERS.dz

                delta_c_o2 = O2.f_transport(teta_sm_i, teta_sm_iplus1, cl_o2, cs_o2, dc_o2, dc_dz_o2, 0, 0, O2.D, UFi_usz,
                                     Rxi_2_o2, SOIL_PLANT.ro, SOIL_PLANT.f, GENERAL_PARAMETERS.dt, GENERAL_PARAMETERS.dz)
                if teta_sm_iplus1 > 0:
                    ci1_o2 = cl_o2 + delta_c_o2
                else:
                    ci1_o2 = 0

                if ci1_o2 <= 0.0000000000000001:
                    ci1_o2 = 0
                else:
                    ci1_o2 = ci1_o2

                    # print('ci1_o2', ci1_o2)
                cl_i1_o2.append(ci1_o2)
                # print('2_o2', cl_i1_o2, l)

                ### Amonia
                Peusz_nh4 = USZ.f_peclet(UFi_usz, NH4.D, GENERAL_PARAMETERS.dz)
                # print('t:', t, 'l:', l, 'Peusz_nh4:', Peusz_nh4)

                if Peusz_nh4 <= 2:
                    if l == 0:  # first cell
                        dc_nh4 = clplus1_nh4 - 2 * cl_nh4 + cpi_nh4
                        dc_dz_nh4 = (clplus1_nh4 - cpi_nh4) / (2 * GENERAL_PARAMETERS.dz)

                    elif l == (USZ.m_usz - 1):  # last cell
                        dc_nh4 = cl_nh4 - 2 * cl_nh4 + clminus1_nh4
                        dc_dz_nh4 = (cl_nh4 - clminus1_nh4) / GENERAL_PARAMETERS.dz

                    else:
                        dc_nh4 = clplus1_nh4 - 2 * cl_nh4 + clminus1_nh4
                        dc_dz_nh4 = (clplus1_nh4 - clminus1_nh4) / (2 * GENERAL_PARAMETERS.dz)

                else:  # Peusz > 2
                    if l == 0:  # first cell
                        dc_nh4 = clplus1_nh4 - 2 * cl_nh4 + cpi_nh4
                        dc_dz_nh4 = (cl_nh4 - cpi_nh4) / GENERAL_PARAMETERS.dz

                    elif l == (USZ.m_usz - 1):  # last cell
                        dc_nh4 = cl_nh4 - 2 * cl_nh4 + clminus1_nh4
                        dc_dz_nh4 = (cl_nh4 - clminus1_nh4) / GENERAL_PARAMETERS.dz

                    else:
                        dc_nh4 = clplus1_nh4 - 2 * cl_nh4 + clminus1_nh4
                        dc_dz_nh4 = (cl_nh4 - clminus1_nh4) / GENERAL_PARAMETERS.dz

                delta_c_nh4 = NH4.f_transport(teta_sm_i, teta_sm_iplus1, cl_nh4, cs_nh4, dc_nh4, dc_dz_nh4, NH4.kads, NH4.kdes,
                                      NH4.D, UFi_usz, Rxi_2_nh4, SOIL_PLANT.ro, SOIL_PLANT.f, GENERAL_PARAMETERS.dt, GENERAL_PARAMETERS.dz)
                #                 print('t: ', t, ', l: ', l)
                #                 print('clplus1_nh4: ',  clplus1_nh4, ', cl_nh4: ' ,cl_nh4 , ', cpi_nh4:' , cpi_nh4, ', clminus1_nh4: ' , clminus1_nh4)
                #                 print('GENERAL_PARAMETERS.dt:', GENERAL_PARAMETERS.dt, 'teta_sm_i:', teta_sm_i, 'teta_sm_iplus1:', teta_sm_iplus1, 'cl_nh4:', cl_nh4, 'cs_nh4:', cs_nh4, 'dc_nh4:', dc_nh4, 'dc_dz_nh4:', dc_dz_nh4, 'NH4.kads:', NH4.kads, 'NH4.kdes:', NH4.kdes, 'NH4.D:', NH4.D, 'UFi_usz:', UFi_usz, 'Rxi_2_nh4:', Rxi_2_nh4)
                #                 print('delta_c_nh4: ', delta_c_nh4)
                #                 input()

                if teta_sm_iplus1 > 0:
                    ci1_nh4 = cl_nh4 + delta_c_nh4
                else:
                    ci1_nh4 = 0

                if ci1_nh4 <= 0.0000000000000001:
                    ci1_nh4 = 0
                else:
                    ci1_nh4 = ci1_nh4
                    # print('ci1_nh4: ', ci1_nh4)
                cl_i1_nh4.append(ci1_nh4)
                # print('2_nh4', cl_i1_nh4)

                ### Nitrate
                Peusz_no3 = USZ.f_peclet(UFi_usz, NO3.D, GENERAL_PARAMETERS.dz)

                if Peusz_no3 <= 2:
                    if l == 0:  # first cell
                        dc_no3 = clplus1_no3 - 2 * cl_no3 + cpi_no3
                        dc_dz_no3 = (clplus1_no3 - cpi_no3) / (2 * GENERAL_PARAMETERS.dz)

                    elif l == (USZ.m_usz - 1):  # last cell
                        dc_no3 = cl_no3 - 2 * cl_no3 + clminus1_no3
                        dc_dz_no3 = (cl_no3 - clminus1_no3) / GENERAL_PARAMETERS.dz

                    else:
                        dc_no3 = clplus1_no3 - 2 * cl_no3 + clminus1_no3
                        dc_dz_no3 = (clplus1_no3 - clminus1_no3) / (2 * GENERAL_PARAMETERS.dz)

                else:  # Peusz > 2
                    if l == 0:  # first cell
                        dc_no3 = clplus1_no3 - 2 * cl_no3 + cpi_no3
                        dc_dz_no3 = (cl_no3 - cpi_no3) / GENERAL_PARAMETERS.dz

                    elif l == (USZ.m_usz - 1):  # last cell
                        dc_no3 = cl_no3 - 2 * cl_no3 + clminus1_no3
                        dc_dz_no3 = (cl_no3 - clminus1_no3) / GENERAL_PARAMETERS.dz

                    else:
                        dc_no3 = clplus1_no3 - 2 * cl_no3 + clminus1_no3
                        dc_dz_no3 = (cl_no3 - clminus1_no3) / GENERAL_PARAMETERS.dz

                delta_c_no3 = NO3.f_transport(teta_sm_i, teta_sm_iplus1, cl_no3, cs_no3, dc_no3, dc_dz_no3, 0, 0, NO3.D,
                                      UFi_usz, Rxi_2_no3, SOIL_PLANT.ro, SOIL_PLANT.f, GENERAL_PARAMETERS.dt, GENERAL_PARAMETERS.dz)
                if teta_sm_iplus1 > 0:
                    ci1_no3 = cl_no3 + delta_c_no3
                else:
                    ci1_no3 = 0

                if ci1_no3 <= 0.0000000000000001:
                    ci1_no3 = 0
                else:
                    ci1_no3 = ci1_no3

                cl_i1_no3.append(ci1_no3)

                ### DOC
                Peusz_doc = USZ.f_peclet(UFi_usz, DOC.D, GENERAL_PARAMETERS.dz)

                if Peusz_doc <= 2:
                    if l == 0:  # first cell
                        dc_doc = clplus1_doc - 2 * cl_doc + cpi_doc
                        dc_dz_doc = (clplus1_doc - cpi_doc) / (2 * GENERAL_PARAMETERS.dz)

                    elif l == (USZ.m_usz - 1):  # last cell
                        dc_doc = cl_doc - 2 * cl_doc + clminus1_doc
                        dc_dz_doc = (cl_doc - clminus1_doc) / GENERAL_PARAMETERS.dz

                    else:
                        dc_doc = clplus1_doc - 2 * cl_doc + clminus1_doc
                        dc_dz_doc = (clplus1_doc - clminus1_doc) / (2 * GENERAL_PARAMETERS.dz)

                else:  # Peusz > 2
                    if l == 0:  # first cell
                        dc_doc = clplus1_doc - 2 * cl_doc + cpi_doc
                        dc_dz_doc = (cl_doc - cpi_doc) / GENERAL_PARAMETERS.dz

                    elif l == (USZ.m_usz - 1):  # last cell
                        dc_doc = cl_doc - 2 * cl_doc + clminus1_doc
                        dc_dz_doc = (cl_doc - clminus1_doc) / GENERAL_PARAMETERS.dz

                    else:
                        dc_doc = clplus1_doc - 2 * cl_doc + clminus1_doc
                        dc_dz_doc = (cl_doc - clminus1_doc) / GENERAL_PARAMETERS.dz

                delta_c_doc = DOC.f_transport(teta_sm_i, teta_sm_iplus1, cl_doc, cs_doc, dc_doc, dc_dz_doc, DOC.kads, DOC.kdes,
                                      DOC.D, UFi_usz, Rxi_2_doc, SOIL_PLANT.ro, SOIL_PLANT.f, GENERAL_PARAMETERS.dt, GENERAL_PARAMETERS.dz)
                if teta_sm_iplus1 > 0:
                    ci1_doc = cl_doc + delta_c_doc
                else:
                    ci1_doc = 0

                if ci1_doc <= 0.0000000000000001:
                    ci1_doc = 0
                else:
                    ci1_doc = ci1_doc

                cl_i1_doc.append(ci1_doc)
                # print('2_nh4', cl_i1_nh4)

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

        for j in range(SZ.m_sz):
            #                 if teta_b_iplus1 > 0:
            #                     den_teta = 1/teta_b_iplus1
            #                 else:
            #                     den_teta = 0

            cj_o2 = cji_o2_list[j]
            cjminus1_o2 = cji_o2_list[j - 1]
            if j < (SZ.m_sz - 1):
                cjplus1_o2 = cji_o2_list[j + 1]
            else:
                cjplus1_o2 = 0
            cmminus1_o2 = cli_o2_list[USZ.m_usz - 1]

            cj_nh4 = cji_nh4_list[j]
            cjminus1_nh4 = cji_nh4_list[j - 1]
            if j < (SZ.m_sz - 1):
                cjplus1_nh4 = cji_nh4_list[j + 1]
            else:
                cjplus1_nh4 = 0
            cmminus1_nh4 = cli_nh4_list[USZ.m_usz - 1]

            cj_no3 = cji_no3_list[j]
            cjminus1_no3 = cji_no3_list[j - 1]
            if j < (SZ.m_sz - 1):
                cjplus1_no3 = cji_no3_list[j + 1]
            else:
                cjplus1_no3 = 0
            cmminus1_no3 = cli_no3_list[USZ.m_usz - 1]

            cj_doc = cji_doc_list[j]
            cjminus1_doc = cji_doc_list[j - 1]
            if j < (SZ.m_sz - 1):
                cjplus1_doc = cji_doc_list[j + 1]
            else:
                cjplus1_doc = 0
            cmminus1_doc = cli_doc_list[USZ.m_usz - 1]

            cs_o2 = 0

            NH4.concentration_soil_SZ_before = NH4.concentration_soil_SZ[t][j]
            cs_nh4 = NH4.f_concentration_soil(NH4.concentration_soil_SZ_before, teta_b_i, NH4.kads2, cj_nh4, NH4.kdes2, NH4.k_mb, SOIL_PLANT.ro, GENERAL_PARAMETERS.dt)
            csi_sz_nh4.append(cs_nh4)
            if cs_nh4 < 0.00000000000001:
                cs_nh4 = 0

            #                 print('====>> t: ', t, ', j: ', j)
            #                 print('---- SOIL ----')
            #                 print('cs_nh4_a: ', NH4.concentration_soil_SZ_before, 'teta_b_i: ', teta_b_i, 'NH4.kads: ', NH4.kads, 'cj_nh4: ', cj_nh4, 'NH4.kdes: ', NH4.kdes, 'NH4.k_mb: ', NH4.k_mb)
            #                 print('cs_nh4: ', cs_nh4)
            #                 print('---- WATER ----')
            #                 print('cmminus1_nh4: ', cmminus1_nh4)
            cs_no3 = 0

            DOC.concentration_soil_SZ_before = DOC.concentration_soil_SZ[t][j]
            cs_doc = DOC.f_concentration_soil(DOC.concentration_soil_SZ_before, teta_b_i, DOC.kads2, cj_doc, DOC.kdes2, DOC.k_mb, SOIL_PLANT.ro, GENERAL_PARAMETERS.dt)
            csi_sz_doc.append(cs_doc)

            UF_sz = []
            alfa2, beta2 = SZ.f_alfa_beta(j)
            # WFR.tQet2[t] = 0
            UFi_sz = SZ.f_unit_flux(alfa2, beta2, WFR.tQfs[t], WFR.tQhc[t], WFR.tQet2[t], Qorif, Qinf_sz, teta_b_i, PZ.Ab)
            #                 print('WFR.tQfs[t], WFR.tQhc[t], WFR.tQet2[t], Qorif, Qinf_sz:', WFR.tQfs[t], WFR.tQhc[t], WFR.tQet2[t], Qorif, Qinf_sz)
            #                 print('UF: ', UFi_sz, ', alfa2: ', alfa2, ', beta2: ' ,beta2)

            UF_sz.append(UFi_sz)

            Rxi_3_o2 = O2.f_reaction_sz(cj_o2, cj_nh4, GENERAL_PARAMETERS.k_nit) + O2.f_plant_uptake_sz(cj_o2, SOIL_PLANT.c_o2_root, teta_b_i, SOIL_PLANT.root_fraction, SOIL_PLANT.lamda)
            Rxi_3_nh4 = NH4.f_reaction_sz() + NH4.f_plant_uptake_sz(cj_nh4, teta_b_i, SOIL_PLANT.root_fraction)
            Rxi_3_no3 = NO3.f_reaction_sz(cj_no3, cj_o2, cj_doc, GENERAL_PARAMETERS.k_denit) + NO3.f_plant_uptake_sz(cj_no3, teta_b_i, SOIL_PLANT.root_fraction)
            Rxi_3_doc = DOC.f_reaction_sz(cj_doc)

            Rxj_o2.append(Rxi_3_o2 * (1 / teta_b_iplus1) * GENERAL_PARAMETERS.dt)
            Rxj_nh4.append(Rxi_3_nh4 * (1 / teta_b_iplus1) * GENERAL_PARAMETERS.dt)
            Rxj_no3.append(Rxi_3_no3 * (1 / teta_b_iplus1) * GENERAL_PARAMETERS.dt)
            Rxj_doc.append(Rxi_3_doc * (1 / teta_b_iplus1) * GENERAL_PARAMETERS.dt)
            # if t < 20:
            # print('t: ', t, 'j: ', j, 'O2.reaction_rate_SZ: ', Rxi_3_o2, 'NH4.reaction_rate_SZ: ', Rxi_3_nh4, 'NO3.reaction_rate_SZ: ', Rxi_3_no3)

            ### Oxygen
            Pesz_o2 = SZ.f_peclet(UFi_sz, O2.D, GENERAL_PARAMETERS.dz)

            if USZ.m_usz < (GENERAL_PARAMETERS.n - 1):
                if Pesz_o2 <= 2:
                    if j == 0:  # first cell
                        dc_o2 = cjplus1_o2 - 2 * cj_o2 + cmminus1_o2
                        dc_dz_o2 = (cjplus1_o2 - cmminus1_o2) / (2 * GENERAL_PARAMETERS.dz)

                    elif j == (SZ.m_sz - 1):  # last cell
                        dc_o2 = cj_o2 - 2 * cj_o2 + cjminus1_o2
                        dc_dz_o2 = (cj_o2 - cjminus1_o2) / GENERAL_PARAMETERS.dz

                    else:
                        dc_o2 = cjplus1_o2 - 2 * cj_o2 + cjminus1_o2
                        dc_dz_o2 = (cjplus1_o2 - cjminus1_o2) / (2 * GENERAL_PARAMETERS.dz)

                else:  # Pesz > 2
                    if j == 0:  # first cell
                        dc_o2 = cjplus1_o2 - 2 * cj_o2 + cmminus1_o2
                        dc_dz_o2 = (cj_o2 - cmminus1_o2) / GENERAL_PARAMETERS.dz

                    elif j == (SZ.m_sz - 1):  # last cell
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
                    if j == 0:  # first cell
                        dc_o2 = cjplus1_o2 - 2 * cj_o2 + cpi_o2
                        dc_dz_o2 = (cjplus1_o2 - cpi_o2) / (2 * GENERAL_PARAMETERS.dz)

                    elif j == (SZ.m_sz - 1):  # last cell
                        dc_o2 = cj_o2 - 2 * cj_o2 + cjminus1_o2
                        dc_dz_o2 = (cj_o2 - cjminus1_o2) / GENERAL_PARAMETERS.dz

                    else:
                        dc_o2 = cjplus1_o2 - 2 * cj_o2 + cjminus1_o2
                        dc_dz_o2 = (cjplus1_o2 - cjminus1_o2) / (2 * GENERAL_PARAMETERS.dz)

                else:  # Pusz > 2
                    if j == 0:  # first cell
                        dc_o2 = cjplus1_o2 - 2 * cj_o2 + cpi_o2
                        dc_dz_o2 = (cj_o2 - cpi_o2) / GENERAL_PARAMETERS.dz

                    elif j == (SZ.m_sz - 1):  # last cell
                        dc_o2 = cj_o2 - 2 * cj_o2 + cjminus1_o2
                        dc_dz_o2 = (cj_o2 - cjminus1_o2) / GENERAL_PARAMETERS.dz

                    else:
                        dc_o2 = cjplus1_o2 - 2 * cj_o2 + cjminus1_o2
                        dc_dz_o2 = (cj_o2 - cjminus1_o2) / GENERAL_PARAMETERS.dz

            delta_c_o2 = O2.f_transport(teta_b_i, teta_b_iplus1, cj_o2, cs_o2, dc_o2, dc_dz_o2, 0, 0, O2.D, UFi_sz, Rxi_3_o2, SOIL_PLANT.ro, SOIL_PLANT.f, GENERAL_PARAMETERS.dt, GENERAL_PARAMETERS.dz)
            if teta_b_iplus1 > 0:
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
                    if j == 0:  # first cell
                        dc_nh4 = cjplus1_nh4 - 2 * cj_nh4 + cmminus1_nh4
                        dc_dz_nh4 = (cjplus1_nh4 - cmminus1_nh4) / (2 * GENERAL_PARAMETERS.dz)

                    elif j == (SZ.m_sz - 1):  # last cell
                        dc_nh4 = cj_nh4 - 2 * cj_nh4 + cjminus1_nh4
                        dc_dz_nh4 = (cj_nh4 - cjminus1_nh4) / GENERAL_PARAMETERS.dz

                    else:
                        dc_nh4 = cjplus1_nh4 - 2 * cj_nh4 + cjminus1_nh4
                        dc_dz_nh4 = (cjplus1_nh4 - cjminus1_nh4) / (2 * GENERAL_PARAMETERS.dz)

                else:  # Pesz > 2
                    if j == 0:  # first cell
                        dc_nh4 = cjplus1_nh4 - 2 * cj_nh4 + cmminus1_nh4
                        dc_dz_nh4 = (cj_nh4 - cmminus1_nh4) / GENERAL_PARAMETERS.dz

                    elif j == (SZ.m_sz - 1):  # last cell
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
                    if j == 0:  # first cell
                        dc_nh4 = cjplus1_nh4 - 2 * cj_nh4 + cpi_nh4
                        dc_dz_nh4 = (cjplus1_nh4 - cpi_nh4) / (2 * GENERAL_PARAMETERS.dz)

                    elif j == (SZ.m_sz - 1):  # last cell
                        dc_nh4 = cj_nh4 - 2 * cj_nh4 + cjminus1_nh4
                        dc_dz_nh4 = (cj_nh4 - cjminus1_nh4) / GENERAL_PARAMETERS.dz

                    else:
                        dc_nh4 = cjplus1_nh4 - 2 * cj_nh4 + cjminus1_nh4
                        dc_dz_nh4 = (cjplus1_nh4 - cjminus1_nh4) / (2 * GENERAL_PARAMETERS.dz)

                else:  # Pusz > 2
                    if j == 0:  # first cell
                        dc_nh4 = cjplus1_nh4 - 2 * cj_nh4 + cpi_nh4
                        dc_dz_nh4 = (cj_nh4 - cpi_nh4) / GENERAL_PARAMETERS.dz

                    elif j == (SZ.m_sz - 1):  # last cell
                        dc_nh4 = cj_nh4 - 2 * cj_nh4 + cjminus1_nh4
                        dc_dz_nh4 = (cj_nh4 - cjminus1_nh4) / GENERAL_PARAMETERS.dz

                    else:
                        dc_nh4 = cjplus1_nh4 - 2 * cj_nh4 + cjminus1_nh4
                        dc_dz_nh4 = (cj_nh4 - cjminus1_nh4) / GENERAL_PARAMETERS.dz

            delta_c_nh4 = NH4.f_transport(teta_b_i, teta_b_iplus1, cj_nh4, cs_nh4, dc_nh4, dc_dz_nh4, NH4.kads2, NH4.kdes2,
                                  NH4.D, UFi_sz, Rxi_3_nh4, SOIL_PLANT.ro, SOIL_PLANT.f, GENERAL_PARAMETERS.dt, GENERAL_PARAMETERS.dz)
            #                 print('teta_b_i: ', teta_b_i, ', teta_b_iplus1: ', teta_b_iplus1, ', cj_nh4: ', cj_nh4, ', cs_nh4: ', cs_nh4, ', dc_nh4: ', dc_nh4, ', dc_dz_nh4: ', dc_dz_nh4, ', NH4.kads2: ', NH4.kads2, ', NH4.kdes2: ', NH4.kdes2, ', NH4.D: ', NH4.D, ', UFi_sz: ', UFi_sz, ', Rxi_3_nh4: ', Rxi_3_nh4)
            #                 print('delta_c_nh4: ', delta_c_nh4)

            if teta_b_iplus1 > 0:
                ci1_nh4 = cj_nh4 + delta_c_nh4
                # print('calculo final - cj_nh4, teta_b_i, teta_b_iplus1, delta_c_nh4: ', cj_nh4, ', ', teta_b_i, ', ', teta_b_iplus1, ', ', delta_c_nh4)
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
                    if j == 0:  # first cell
                        dc_no3 = cjplus1_no3 - 2 * cj_no3 + cmminus1_no3
                        dc_dz_no3 = (cjplus1_no3 - cmminus1_no3) / (2 * GENERAL_PARAMETERS.dz)

                    elif j == (SZ.m_sz - 1):  # last cell
                        dc_no3 = cj_no3 - 2 * cj_no3 + cjminus1_no3
                        dc_dz_no3 = (cj_no3 - cjminus1_no3) / GENERAL_PARAMETERS.dz

                    else:
                        dc_no3 = cjplus1_no3 - 2 * cj_no3 + cjminus1_no3
                        dc_dz_no3 = (cjplus1_no3 - cjminus1_no3) / (2 * GENERAL_PARAMETERS.dz)

                else:  # Pesz > 2
                    if j == 0:  # first cell
                        dc_no3 = cjplus1_no3 - 2 * cj_no3 + cmminus1_no3
                        dc_dz_no3 = (cj_no3 - cmminus1_no3) / GENERAL_PARAMETERS.dz

                    elif j == (SZ.m_sz - 1):  # last cell
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
                    if j == 0:  # first cell
                        dc_no3 = cjplus1_no3 - 2 * cj_no3 + cpi_no3
                        dc_dz_no3 = (cjplus1_no3 - cpi_no3) / (2 * GENERAL_PARAMETERS.dz)

                    elif j == (SZ.m_sz - 1):  # last cell
                        dc_no3 = cj_no3 - 2 * cj_no3 + cjminus1_no3
                        dc_dz_no3 = (cj_no3 - cjminus1_no3) / GENERAL_PARAMETERS.dz

                    else:
                        dc_no3 = cjplus1_no3 - 2 * cj_no3 + cjminus1_no3
                        dc_dz_no3 = (cjplus1_no3 - cjminus1_no3) / (2 * GENERAL_PARAMETERS.dz)

                else:  # Pusz > 2
                    if j == 0:  # first cell
                        dc_no3 = cjplus1_no3 - 2 * cj_no3 + cpi_no3
                        dc_dz_no3 = (cj_no3 - cpi_no3) / GENERAL_PARAMETERS.dz

                    elif j == (SZ.m_sz - 1):  # last cell
                        dc_no3 = cj_no3 - 2 * cj_no3 + cjminus1_no3
                        dc_dz_no3 = (cj_no3 - cjminus1_no3) / GENERAL_PARAMETERS.dz

                    else:
                        dc_no3 = cjplus1_no3 - 2 * cj_no3 + cjminus1_no3
                        dc_dz_no3 = (cj_no3 - cjminus1_no3) / GENERAL_PARAMETERS.dz

            delta_c_no3 = NO3.f_transport(teta_b_i, teta_b_iplus1, cj_no3, cs_no3, dc_no3, dc_dz_no3, 0, 0, NO3.D, UFi_sz,
                                  Rxi_3_no3, SOIL_PLANT.ro, SOIL_PLANT.f, GENERAL_PARAMETERS.dt, GENERAL_PARAMETERS.dz)
            if teta_b_iplus1 > 0:
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
                    if j == 0:  # first cell
                        dc_doc = cjplus1_doc - 2 * cj_doc + cmminus1_doc
                        dc_dz_doc = (cjplus1_doc - cmminus1_doc) / (2 * GENERAL_PARAMETERS.dz)

                    elif j == (SZ.m_sz - 1):  # last cell
                        dc_doc = cj_doc - 2 * cj_doc + cjminus1_doc
                        dc_dz_doc = (cj_doc - cjminus1_doc) / GENERAL_PARAMETERS.dz

                    else:
                        dc_doc = cjplus1_doc - 2 * cj_doc + cjminus1_doc
                        dc_dz_doc = (cjplus1_doc - cjminus1_doc) / (2 * GENERAL_PARAMETERS.dz)

                else:  # Pesz > 2
                    if j == 0:  # first cell
                        dc_doc = cjplus1_doc - 2 * cj_doc + cmminus1_doc
                        dc_dz_doc = (cj_doc - cmminus1_doc) / GENERAL_PARAMETERS.dz

                    elif j == (SZ.m_sz - 1):  # last cell
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
                    if j == 0:  # first cell
                        dc_doc = cjplus1_doc - 2 * cj_doc + cpi_doc
                        dc_dz_doc = (cjplus1_doc - cpi_doc) / (2 * GENERAL_PARAMETERS.dz)

                    elif j == (SZ.m_sz - 1):  # last cell
                        dc_doc = cj_doc - 2 * cj_doc + cjminus1_doc
                        dc_dz_doc = (cj_doc - cjminus1_doc) / GENERAL_PARAMETERS.dz

                    else:
                        dc_doc = cjplus1_doc - 2 * cj_doc + cjminus1_doc
                        dc_dz_doc = (cjplus1_doc - cjminus1_doc) / (2 * GENERAL_PARAMETERS.dz)

                else:  # Pusz > 2
                    if j == 0:  # first cell
                        dc_doc = cjplus1_doc - 2 * cj_doc + cpi_doc
                        dc_dz_doc = (cj_doc - cpi_doc) / GENERAL_PARAMETERS.dz

                    elif j == (SZ.m_sz - 1):  # last cell
                        dc_doc = cj_doc - 2 * cj_doc + cjminus1_doc
                        dc_dz_doc = (cj_doc - cjminus1_doc) / GENERAL_PARAMETERS.dz

                    else:
                        dc_doc = cjplus1_doc - 2 * cj_doc + cjminus1_doc
                        dc_dz_doc = (cj_doc - cjminus1_doc) / GENERAL_PARAMETERS.dz

            delta_c_doc = DOC.f_transport(teta_b_i, teta_b_iplus1, cj_doc, cs_doc, dc_doc, dc_dz_doc, DOC.kads2, DOC.kdes2,
                                  DOC.D, UFi_sz, Rxi_3_doc, SOIL_PLANT.ro, SOIL_PLANT.f, GENERAL_PARAMETERS.dt, GENERAL_PARAMETERS.dz)
            if teta_b_iplus1 > 0:
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
            Mstor_o2_ast = sum(O2.concentration_USZ[t]) * PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz + sum(O2.concentration_SZ[t]) * PZ.Ab * \
                           WFR.tteta_sz[t] * WFR.thsz[t] * 1000 / SZ.m_sz
            O2.mass_stormwater.append(Mstor_o2_ast)

            Msoil_o2_a = 0
            Msoil_o2 = 0
            O2.mass_soil.append(0)

            MRx_o2 = - (sum(O2.reaction_rate_USZ[t]) * PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz + sum(O2.reaction_rate_SZ[t]) * PZ.Ab *
                        WFR.tteta_sz[t] * WFR.thsz[t] * 1000 / SZ.m_sz)
            if t == 0:
                O2.mass_accumulated_reaction.append(MRx_o2)
            else:
                O2.mass_accumulated_reaction.append(MRx_o2 + O2.mass_accumulated_reaction[-1])

            Min_o2 = WFR.tQin[t] * O2.concentration_inflow[t] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                O2.mass_accumulated_inflow.append(Min_o2)
            else:
                O2.mass_accumulated_inflow.append(Min_o2 + O2.mass_accumulated_inflow[-1])

            Mover_o2 = WFR.tQover[t] * O2.concentration_PZ[t] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                O2.mass_accumulated_overflow.append(Mover_o2)
            else:
                O2.mass_accumulated_overflow.append(Mover_o2 + O2.mass_accumulated_overflow[-1])

            Mpipe_o2 = WFR.tQpipe[t] * O2.concentration_SZ[t][SZ.m_sz - 1] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                O2.mass_accumulated_pipe_outflow.append(Mpipe_o2)
            else:
                O2.mass_accumulated_pipe_outflow.append(Mpipe_o2 + O2.mass_accumulated_pipe_outflow[-1])

            Minfsz_o2 = WFR.tQinfsz[t] * O2.concentration_SZ[t][SZ.m_sz - 1] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                O2.mass_accumulated_infiltrated_to_SZ.append(Minfsz_o2)
            else:
                O2.mass_accumulated_infiltrated_to_SZ.append(Minfsz_o2 + O2.mass_accumulated_infiltrated_to_SZ[-1])

            Met_o2 = WFR.tQet[t] * cl_i1_o2[0] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                O2.mass_accumulated_evapotranspiration.append(Met_o2)
            else:
                O2.mass_accumulated_evapotranspiration.append(Met_o2 + O2.mass_accumulated_evapotranspiration[-1])

            Mpz_o2 = WFR.thpEND[t] * PZ.Ab * 1000 * O2.concentration_PZ[t]
            O2.mass_PZ.append(Mpz_o2)

            Mstor_o2_mb = O2.mass_accumulated_inflow[-1] - O2.mass_PZ[-1] - O2.mass_accumulated_overflow[-1] - O2.mass_accumulated_pipe_outflow[-1] - O2.mass_accumulated_infiltrated_to_SZ[
                -1] - O2.mass_accumulated_evapotranspiration[-1] - O2.mass_soil[-1] - O2.mass_accumulated_reaction[-1]
            O2.mass_balance.append(Mstor_o2_mb)

            delta_ast_o2_usz = (Mstor_o2_mb - Mstor_o2_ast) / (PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz)
            cl_i1_o2[1] = cl_i1_o2[1] + delta_ast_o2_usz
            if cl_i1_o2[1] > 0:
                cl_i1_o2[1] = cl_i1_o2[1]
            else:
                cl_i1_o2[1] = 0

            ### Amonia

            Mstor_nh4_ast = sum(NH4.concentration_USZ[t]) * PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz + sum(NH4.concentration_SZ[t]) * PZ.Ab * \
                            WFR.tteta_sz[t] * WFR.thsz[t] * 1000 / SZ.m_sz
            NH4.mass_stormwater.append(Mstor_nh4_ast)

            Msoil_nh4_a = NH4.concentration_soil_USZ_before * SOIL_PLANT.ro * PZ.Ab * WFR.thusz[t] * 1000 + NH4.concentration_soil_SZ_before * SOIL_PLANT.ro * PZ.Ab * WFR.thsz[t] * 1000
            Msoil_nh4 = sum(NH4.concentration_soil_USZ[t]) * SOIL_PLANT.ro * PZ.Ab * WFR.thusz[t] * 1000 / USZ.m_usz + sum(NH4.concentration_soil_SZ[t]) * SOIL_PLANT.ro * PZ.Ab * WFR.thsz[
                t] * 1000 / SZ.m_sz
            NH4.mass_soil.append(Msoil_nh4 - Msoil_nh4_a)

            MRx_nh4 = - (sum(NH4.reaction_rate_USZ[t]) * PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz + sum(NH4.reaction_rate_SZ[t]) * PZ.Ab *
                         WFR.tteta_sz[t] * WFR.thsz[t] * 1000 / SZ.m_sz)
            if t == 0:
                NH4.mass_accumulated_reaction.append(MRx_nh4)
            else:
                NH4.mass_accumulated_reaction.append(MRx_nh4 + NH4.mass_accumulated_reaction[-1])

            Min_nh4 = WFR.tQin[t] * NH4.concentration_inflow[t] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                NH4.mass_accumulated_inflow.append(Min_nh4)
            else:
                NH4.mass_accumulated_inflow.append(Min_nh4 + NH4.mass_accumulated_inflow[-1])

            Mover_nh4 = WFR.tQover[t] * NH4.concentration_PZ[t] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                NH4.mass_accumulated_overflow.append(Mover_nh4)
            else:
                NH4.mass_accumulated_overflow.append(Mover_nh4 + NH4.mass_accumulated_overflow[-1])

            Mpipe_nh4 = WFR.tQpipe[t] * NH4.concentration_SZ[t][SZ.m_sz - 1] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                NH4.mass_accumulated_pipe_outflow.append(Mpipe_nh4)
            else:
                NH4.mass_accumulated_pipe_outflow.append(Mpipe_nh4 + NH4.mass_accumulated_pipe_outflow[-1])

            Minfsz_nh4 = WFR.tQinfsz[t] * NH4.concentration_SZ[t][SZ.m_sz - 1] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                NH4.mass_accumulated_infiltrated_to_SZ.append(Minfsz_nh4)
            else:
                NH4.mass_accumulated_infiltrated_to_SZ.append(Minfsz_nh4 + NH4.mass_accumulated_infiltrated_to_SZ[-1])

            Met_nh4 = WFR.tQet[t] * cl_i1_nh4[0] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                NH4.mass_accumulated_evapotranspiration.append(Met_nh4)
            else:
                NH4.mass_accumulated_evapotranspiration.append(Met_nh4 + NH4.mass_accumulated_evapotranspiration[-1])

            Mpz_nh4 = WFR.thpEND[t] * PZ.Ab * 1000 * NH4.concentration_PZ[t]
            NH4.mass_PZ.append(Mpz_nh4)

            Mstor_nh4_mb = NH4.mass_accumulated_inflow[-1] - NH4.mass_PZ[-1] - NH4.mass_accumulated_overflow[-1] - NH4.mass_accumulated_pipe_outflow[-1] - \
                           NH4.mass_accumulated_infiltrated_to_SZ[-1] - NH4.mass_accumulated_evapotranspiration[-1] - NH4.mass_soil[-1] - NH4.mass_accumulated_reaction[-1]
            NH4.mass_balance.append(Mstor_nh4_mb)

            #         delta_ast_nh4_usz = 0
            #         delta_ast_nh4_sz = 0

            #         delta_ast_nh4_usz = (Mstor_nh4_mb - Mstor_nh4_ast)/(GENERAL_PARAMETERS.n*(PZ.Ab*WFR.tteta_usz[t]*WFR.thusz[t]*1000/USZ.m_usz))
            #         delta_ast_nh4_sz = (Mstor_nh4_mb - Mstor_nh4_ast)/(GENERAL_PARAMETERS.n*(PZ.Ab*WFR.tteta_sz[t]*WFR.thsz[t]*1000/SZ.m_sz))

            delta_ast_nh4_usz = (Mstor_nh4_mb - Mstor_nh4_ast) / (PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz)
            cl_i1_nh4[1] = cl_i1_nh4[1] + delta_ast_nh4_usz
            if cl_i1_nh4[1] > 0:
                cl_i1_nh4[1] = cl_i1_nh4[1]
            else:
                cl_i1_nh4[1] = 0
            #         sum_cusz_nh4 = sum(NH4.concentration_USZ[t])
            #         sum_csz_nh4 = sum(NH4.concentration_SZ[t])
            #         print('t: ', t, 'Mstor_trace: ', Mstor_trace, ', Mstor_mb: ', Mstor_mb, ', dif: ', (Mstor_mb - Mstor_trace))
            #         print('t: ', t, ', delta_usz: ', delta_trace_nh4_usz, ', delta_sz: ', delta_trace_nh4_sz)
            #         print('teta_sm: ', WFR.tteta_usz[t], ', husz: ', WFR.thusz[t], ', teta_b: ', WFR.tteta_sz[t], ', hsz: ', WFR.thsz[t])
            #         print('sum_cusz_nh4: ', sum_cusz_nh4, ', sum_csz_nh4: ', sum_csz_nh4)
            #         print('Qin: ', WFR.tQin[t], ',hpEND_i: ', WFR.thpEND[t], ', hpEND_i-1: ', WFR.thpEND[t-1], ',Qpipe: ', WFR.tQpipe[t], ', Qover: ', WFR.tQover[t], ', Qinfsz: ', WFR.tQinfsz[t], ', Qet: ', WFR.tQet[t])
            #         print('Cin: ', NH4.concentration_inflow[t], ', Cp_i: ', NH4.concentration_PZ[t], ', Cp_i-1: ', NH4.concentration_PZ[t-1], ', Csz_ultimo: ', c_pipe)
            #         input()
            #
            #         for l in range(USZ.m_usz):
            #             cl_i1_nh4[l] = cl_i1_nh4[l] + delta_ast_nh4_usz
            #         for j in range(SZ.m_sz):
            #             cj_i1_nh4[j] = cj_i1_nh4[j] + delta_ast_nh4_sz

            ### Nitrate

            Mstor_no3_ast = sum(NO3.concentration_USZ[t]) * PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz + sum(NO3.concentration_SZ[t]) * PZ.Ab * \
                            WFR.tteta_sz[t] * WFR.thsz[t] * 1000 / SZ.m_sz
            NO3.mass_stormwater.append(Mstor_no3_ast)

            Msoil_no3_a = 0
            Msoil_no3 = 0
            NO3.mass_soil.append(0)

            MRx_no3 = - (sum(NO3.reaction_rate_USZ[t]) * PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz + sum(NO3.reaction_rate_SZ[t]) * PZ.Ab *
                         WFR.tteta_sz[t] * WFR.thsz[t] * 1000 / SZ.m_sz)
            if t == 0:
                NO3.mass_accumulated_reaction.append(MRx_no3)
            else:
                NO3.mass_accumulated_reaction.append(MRx_no3 + NO3.mass_accumulated_reaction[-1])

            Min_no3 = WFR.tQin[t] * NO3.concentration_inflow[t] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                NO3.mass_accumulated_inflow.append(Min_no3)
            else:
                NO3.mass_accumulated_inflow.append(Min_no3 + NO3.mass_accumulated_inflow[-1])

            Mover_no3 = WFR.tQover[t] * NO3.concentration_PZ[t] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                NO3.mass_accumulated_overflow.append(Mover_no3)
            else:
                NO3.mass_accumulated_overflow.append(Mover_no3 + NO3.mass_accumulated_overflow[-1])

            Mpipe_no3 = WFR.tQpipe[t] * NO3.concentration_SZ[t][SZ.m_sz - 1] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                NO3.mass_accumulated_pipe_outflow.append(Mpipe_no3)
            else:
                NO3.mass_accumulated_pipe_outflow.append(Mpipe_no3 + NO3.mass_accumulated_pipe_outflow[-1])

            Minfsz_no3 = WFR.tQinfsz[t] * NO3.concentration_SZ[t][SZ.m_sz - 1] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                NO3.mass_accumulated_infiltrated_to_SZ.append(Minfsz_no3)
            else:
                NO3.mass_accumulated_infiltrated_to_SZ.append(Minfsz_no3 + NO3.mass_accumulated_infiltrated_to_SZ[-1])

            Met_no3 = WFR.tQet[t] * cl_i1_no3[0] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                NO3.mass_accumulated_evapotranspiration.append(Met_no3)
            else:
                NO3.mass_accumulated_evapotranspiration.append(Met_no3 + NO3.mass_accumulated_evapotranspiration[-1])

            Mpz_no3 = WFR.thpEND[t] * PZ.Ab * 1000 * NO3.concentration_PZ[t]
            NO3.mass_PZ.append(Mpz_no3)

            Mstor_no3_mb = NO3.mass_accumulated_inflow[-1] - NO3.mass_PZ[-1] - NO3.mass_accumulated_overflow[-1] - NO3.mass_accumulated_pipe_outflow[-1] - \
                           NO3.mass_accumulated_infiltrated_to_SZ[-1] - NO3.mass_accumulated_evapotranspiration[-1] - NO3.mass_soil[-1] - NO3.mass_accumulated_reaction[-1]
            NO3.mass_balance.append(Mstor_no3_mb)

            delta_ast_no3_usz = (Mstor_no3_mb - Mstor_no3_ast) / (PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz)
            cl_i1_no3[1] = cl_i1_no3[1] + delta_ast_no3_usz
            if cl_i1_no3[1] > 0:
                cl_i1_no3[1] = cl_i1_no3[1]
            else:
                cl_i1_no3[1] = 0

                ### DOC

            Mstor_doc_ast = sum(DOC.concentration_USZ[t]) * PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz + sum(DOC.concentration_SZ[t]) * PZ.Ab * \
                            WFR.tteta_sz[t] * WFR.thsz[t] * 1000 / SZ.m_sz
            DOC.mass_stormwater.append(Mstor_doc_ast)

            Msoil_doc_a = DOC.concentration_soil_USZ_before * SOIL_PLANT.ro * PZ.Ab * WFR.thusz[t] * 1000 + DOC.concentration_soil_SZ_before * SOIL_PLANT.ro * PZ.Ab * WFR.thsz[t] * 1000
            Msoil_doc = sum(DOC.concentration_soil_USZ[t]) * SOIL_PLANT.ro * PZ.Ab * WFR.thusz[t] * 1000 / USZ.m_usz + sum(DOC.concentration_soil_SZ[t]) * SOIL_PLANT.ro * PZ.Ab * WFR.thsz[
                t] * 1000 / SZ.m_sz
            DOC.mass_soil.append(Msoil_doc - Msoil_doc_a)

            MRx_doc = - (sum(DOC.reaction_rate_USZ[t]) * PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz + sum(DOC.reaction_rate_SZ[t]) * PZ.Ab *
                         WFR.tteta_sz[t] * WFR.thsz[t] * 1000 / SZ.m_sz)
            if t == 0:
                DOC.mass_accumulated_reaction.append(MRx_doc)
            else:
                DOC.mass_accumulated_reaction.append(MRx_doc + DOC.mass_accumulated_reaction[-1])

            Min_doc = WFR.tQin[t] * DOC.concentration_inflow[t] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                DOC.mass_accumulated_inflow.append(Min_doc)
            else:
                DOC.mass_accumulated_inflow.append(Min_doc + DOC.mass_accumulated_inflow[-1])

            Mover_doc = WFR.tQover[t] * DOC.concentration_PZ[t] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                DOC.mass_accumulated_overflow.append(Mover_doc)
            else:
                DOC.mass_accumulated_overflow.append(Mover_doc + DOC.mass_accumulated_overflow[-1])

            Mpipe_doc = WFR.tQpipe[t] * DOC.concentration_SZ[t][SZ.m_sz - 1] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                DOC.mass_accumulated_pipe_outflow.append(Mpipe_doc)
            else:
                DOC.mass_accumulated_pipe_outflow.append(Mpipe_doc + DOC.mass_accumulated_pipe_outflow[-1])

            Minfsz_doc = WFR.tQinfsz[t] * DOC.concentration_SZ[t][SZ.m_sz - 1] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                DOC.mass_accumulated_infiltrated_to_SZ.append(Minfsz_doc)
            else:
                DOC.mass_accumulated_infiltrated_to_SZ.append(Minfsz_doc + DOC.mass_accumulated_infiltrated_to_SZ[-1])

            Met_doc = WFR.tQet[t] * cl_i1_doc[0] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                DOC.mass_accumulated_evapotranspiration.append(Met_doc)
            else:
                DOC.mass_accumulated_evapotranspiration.append(Met_doc + DOC.mass_accumulated_evapotranspiration[-1])

            Mpz_doc = WFR.thpEND[t] * PZ.Ab * 1000 * DOC.concentration_PZ[t]
            DOC.mass_PZ.append(Mpz_doc)

            Mstor_doc_mb = DOC.mass_accumulated_inflow[-1] - DOC.mass_PZ[-1] - DOC.mass_accumulated_overflow[-1] - DOC.mass_accumulated_pipe_outflow[-1] - \
                           DOC.mass_accumulated_infiltrated_to_SZ[-1] - DOC.mass_accumulated_evapotranspiration[-1] - DOC.mass_soil[-1] - DOC.mass_accumulated_reaction[-1]
            DOC.mass_balance.append(Mstor_doc_mb)

            delta_ast_doc_usz = (Mstor_doc_mb - Mstor_doc_ast) / (PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz)
            cl_i1_doc[1] = cl_i1_doc[1] + delta_ast_doc_usz
            if cl_i1_doc[1] > 0:
                cl_i1_doc[1] = cl_i1_doc[1]
            else:
                cl_i1_doc[1] = 0


        else:  # if GENERAL_PARAMETERS.hpipe == 0

            ### Oxygen

            Mstor_o2_ast = sum(O2.concentration_USZ[t]) * PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz
            O2.mass_stormwater.append(Mstor_o2_ast)

            Msoil_o2_a = 0
            Msoil_o2 = 0
            O2.mass_soil.append(0)

            MRx_o2 = - (sum(O2.reaction_rate_USZ[t]) * PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz)
            if t == 0:
                O2.mass_accumulated_reaction.append(MRx_o2)
            else:
                O2.mass_accumulated_reaction.append(MRx_o2 + O2.mass_accumulated_reaction[-1])

            Min_o2 = WFR.tQin[t] * O2.concentration_inflow[t] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                O2.mass_accumulated_inflow.append(Min_o2)
            else:
                O2.mass_accumulated_inflow.append(Min_o2 + O2.mass_accumulated_inflow[-1])

            Mover_o2 = WFR.tQover[t] * O2.concentration_PZ[t] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                O2.mass_accumulated_overflow.append(Mover_o2)
            else:
                O2.mass_accumulated_overflow.append(Mover_o2 + O2.mass_accumulated_overflow[-1])

            Mpipe_o2 = WFR.tQpipe[t] * O2.concentration_USZ[t][USZ.m_usz - 1] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                O2.mass_accumulated_pipe_outflow.append(Mpipe_o2)
            else:
                O2.mass_accumulated_pipe_outflow.append(Mpipe_o2 + O2.mass_accumulated_pipe_outflow[-1])

            Minfsz_o2 = WFR.tQinfsz[t] * O2.concentration_USZ[t][USZ.m_usz - 1] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                O2.mass_accumulated_infiltrated_to_SZ.append(Minfsz_o2)
            else:
                O2.mass_accumulated_infiltrated_to_SZ.append(Minfsz_o2 + O2.mass_accumulated_infiltrated_to_SZ[-1])

            Met_o2 = WFR.tQet[t] * cl_i1_o2[0] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                O2.mass_accumulated_evapotranspiration.append(Met_o2)
            else:
                O2.mass_accumulated_evapotranspiration.append(Met_o2 + O2.mass_accumulated_evapotranspiration[-1])

            Mpz_o2 = WFR.thpEND[t] * PZ.Ab * 1000 * O2.concentration_PZ[t]
            O2.mass_PZ.append(Mpz_o2)

            Mstor_o2_mb = O2.mass_accumulated_inflow[-1] - O2.mass_PZ[-1] - O2.mass_accumulated_overflow[-1] - O2.mass_accumulated_pipe_outflow[-1] - O2.mass_accumulated_infiltrated_to_SZ[
                -1] - O2.mass_accumulated_evapotranspiration[-1] - O2.mass_soil[-1] - O2.mass_accumulated_reaction[-1]
            O2.mass_balance.append(Mstor_o2_mb)

            delta_ast_o2_usz = (Mstor_o2_mb - Mstor_o2_ast) / (PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz)
            cl_i1_o2[1] = cl_i1_o2[1] + delta_ast_o2_usz
            if cl_i1_o2[1] > 0:
                cl_i1_o2[1] = cl_i1_o2[1]
            else:
                cl_i1_o2[1] = 0

            ### Amonia

            Mstor_nh4_ast = sum(NH4.concentration_USZ[t]) * PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz
            NH4.mass_stormwater.append(Mstor_nh4_ast)

            Msoil_nh4_a = NH4.concentration_soil_USZ_before * SOIL_PLANT.ro * PZ.Ab * WFR.thusz[t] * 1000
            Msoil_nh4 = sum(NH4.concentration_soil_USZ[t]) * SOIL_PLANT.ro * PZ.Ab * WFR.thusz[t] * 1000 / USZ.m_usz
            NH4.mass_soil.append(Msoil_nh4 - Msoil_nh4_a)

            MRx_nh4 = - (sum(NH4.reaction_rate_USZ[t]) * PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz)
            if t == 0:
                NH4.mass_accumulated_reaction.append(MRx_nh4)
            else:
                NH4.mass_accumulated_reaction.append(MRx_nh4 + NH4.mass_accumulated_reaction[-1])

            Min_nh4 = WFR.tQin[t] * NH4.concentration_inflow[t] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                NH4.mass_accumulated_inflow.append(Min_nh4)
            else:
                NH4.mass_accumulated_inflow.append(Min_nh4 + NH4.mass_accumulated_inflow[-1])

            Mover_nh4 = WFR.tQover[t] * NH4.concentration_PZ[t] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                NH4.mass_accumulated_overflow.append(Mover_nh4)
            else:
                NH4.mass_accumulated_overflow.append(Mover_nh4 + NH4.mass_accumulated_overflow[-1])

            Mpipe_nh4 = WFR.tQpipe[t] * NH4.concentration_USZ[t][USZ.m_usz - 1] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                NH4.mass_accumulated_pipe_outflow.append(Mpipe_nh4)
            else:
                NH4.mass_accumulated_pipe_outflow.append(Mpipe_nh4 + NH4.mass_accumulated_pipe_outflow[-1])

            Minfsz_nh4 = WFR.tQinfsz[t] * NH4.concentration_USZ[t][USZ.m_usz - 1] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                NH4.mass_accumulated_infiltrated_to_SZ.append(Minfsz_nh4)
            else:
                NH4.mass_accumulated_infiltrated_to_SZ.append(Minfsz_nh4 + NH4.mass_accumulated_infiltrated_to_SZ[-1])

            Met_nh4 = WFR.tQet[t] * cl_i1_nh4[0] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                NH4.mass_accumulated_evapotranspiration.append(Met_nh4)
            else:
                NH4.mass_accumulated_evapotranspiration.append(Met_nh4 + NH4.mass_accumulated_evapotranspiration[-1])

            Mpz_nh4 = WFR.thpEND[t] * PZ.Ab * 1000 * NH4.concentration_PZ[t]
            NH4.mass_PZ.append(Mpz_nh4)

            Mstor_nh4_mb = NH4.mass_accumulated_inflow[-1] - NH4.mass_PZ[-1] - NH4.mass_accumulated_overflow[-1] - NH4.mass_accumulated_pipe_outflow[-1] - \
                           NH4.mass_accumulated_infiltrated_to_SZ[-1] - NH4.mass_accumulated_evapotranspiration[-1] - NH4.mass_soil[-1] - NH4.mass_accumulated_reaction[-1]
            NH4.mass_balance.append(Mstor_nh4_mb)

            delta_ast_nh4_usz = (Mstor_nh4_mb - Mstor_nh4_ast) / (PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz)
            cl_i1_nh4[1] = cl_i1_nh4[1] + delta_ast_nh4_usz
            if cl_i1_nh4[1] > 0:
                cl_i1_nh4[1] = cl_i1_nh4[1]
            else:
                cl_i1_nh4[1] = 0

            ### Nitrate

            Mstor_no3_ast = sum(NO3.concentration_USZ[t]) * PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz
            NO3.mass_stormwater.append(Mstor_no3_ast)

            Msoil_no3_a = 0
            Msoil_no3 = 0
            NO3.mass_soil.append(0)

            MRx_no3 = - (sum(NO3.reaction_rate_USZ[t]) * PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz)
            if t == 0:
                NO3.mass_accumulated_reaction.append(MRx_no3)
            else:
                NO3.mass_accumulated_reaction.append(MRx_no3 + NO3.mass_accumulated_reaction[-1])

            Min_no3 = WFR.tQin[t] * NO3.concentration_inflow[t] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                NO3.mass_accumulated_inflow.append(Min_no3)
            else:
                NO3.mass_accumulated_inflow.append(Min_no3 + NO3.mass_accumulated_inflow[-1])

            Mover_no3 = WFR.tQover[t] * NO3.concentration_PZ[t] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                NO3.mass_accumulated_overflow.append(Mover_no3)
            else:
                NO3.mass_accumulated_overflow.append(Mover_no3 + NO3.mass_accumulated_overflow[-1])

            Mpipe_no3 = WFR.tQpipe[t] * NO3.concentration_USZ[t][USZ.m_usz - 1] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                NO3.mass_accumulated_pipe_outflow.append(Mpipe_no3)
            else:
                NO3.mass_accumulated_pipe_outflow.append(Mpipe_no3 + NO3.mass_accumulated_pipe_outflow[-1])

            Minfsz_no3 = WFR.tQinfsz[t] * NO3.concentration_USZ[t][USZ.m_usz - 1] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                NO3.mass_accumulated_infiltrated_to_SZ.append(Minfsz_no3)
            else:
                NO3.mass_accumulated_infiltrated_to_SZ.append(Minfsz_no3 + NO3.mass_accumulated_infiltrated_to_SZ[-1])

            Met_no3 = WFR.tQet[t] * cl_i1_no3[0] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                NO3.mass_accumulated_evapotranspiration.append(Met_no3)
            else:
                NO3.mass_accumulated_evapotranspiration.append(Met_no3 + NO3.mass_accumulated_evapotranspiration[-1])

            Mpz_no3 = WFR.thpEND[t] * PZ.Ab * 1000 * NO3.concentration_PZ[t]
            NO3.mass_PZ.append(Mpz_no3)

            Mstor_no3_mb = NO3.mass_accumulated_inflow[-1] - NO3.mass_PZ[-1] - NO3.mass_accumulated_overflow[-1] - NO3.mass_accumulated_pipe_outflow[-1] - \
                           NO3.mass_accumulated_infiltrated_to_SZ[-1] - NO3.mass_accumulated_evapotranspiration[-1] - NO3.mass_soil[-1] - NO3.mass_accumulated_reaction[-1]
            NO3.mass_balance.append(Mstor_no3_mb)

            delta_ast_no3_usz = (Mstor_no3_mb - Mstor_no3_ast) / (PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz)
            cl_i1_no3[1] = cl_i1_no3[1] + delta_ast_no3_usz
            if cl_i1_no3[1] > 0:
                cl_i1_no3[1] = cl_i1_no3[1]
            else:
                cl_i1_no3[1] = 0

                ### DOC

            Mstor_doc_ast = sum(DOC.concentration_USZ[t]) * PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz
            DOC.mass_stormwater.append(Mstor_doc_ast)

            Msoil_doc_a = DOC.concentration_soil_USZ_before * SOIL_PLANT.ro * PZ.Ab * WFR.thusz[t] * 1000
            Msoil_doc = sum(DOC.concentration_soil_USZ[t]) * SOIL_PLANT.ro * PZ.Ab * WFR.thusz[t] * 1000 / USZ.m_usz
            DOC.mass_soil.append(Msoil_doc - Msoil_doc_a)

            MRx_doc = - (sum(DOC.reaction_rate_USZ[t]) * PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz)
            if t == 0:
                DOC.mass_accumulated_reaction.append(MRx_doc)
            else:
                DOC.mass_accumulated_reaction.append(MRx_doc + DOC.mass_accumulated_reaction[-1])

            Min_doc = WFR.tQin[t] * DOC.concentration_inflow[t] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                DOC.mass_accumulated_inflow.append(Min_doc)
            else:
                DOC.mass_accumulated_inflow.append(Min_doc + DOC.mass_accumulated_inflow[-1])

            Mover_doc = WFR.tQover[t] * DOC.concentration_PZ[t] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                DOC.mass_accumulated_overflow.append(Mover_doc)
            else:
                DOC.mass_accumulated_overflow.append(Mover_doc + DOC.mass_accumulated_overflow[-1])

            Mpipe_doc = WFR.tQpipe[t] * DOC.concentration_USZ[t][USZ.m_usz - 1] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                DOC.mass_accumulated_pipe_outflow.append(Mpipe_doc)
            else:
                DOC.mass_accumulated_pipe_outflow.append(Mpipe_doc + DOC.mass_accumulated_pipe_outflow[-1])

            Minfsz_doc = WFR.tQinfsz[t] * DOC.concentration_USZ[t][USZ.m_usz - 1] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                DOC.mass_accumulated_infiltrated_to_SZ.append(Minfsz_doc)
            else:
                DOC.mass_accumulated_infiltrated_to_SZ.append(Minfsz_doc + DOC.mass_accumulated_infiltrated_to_SZ[-1])

            Met_doc = WFR.tQet[t] * cl_i1_doc[0] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                DOC.mass_accumulated_evapotranspiration.append(Met_doc)
            else:
                DOC.mass_accumulated_evapotranspiration.append(Met_doc + DOC.mass_accumulated_evapotranspiration[-1])

            Mpz_doc = WFR.thpEND[t] * PZ.Ab * 1000 * DOC.concentration_PZ[t]
            DOC.mass_PZ.append(Mpz_doc)

            Mstor_doc_mb = DOC.mass_accumulated_inflow[-1] - DOC.mass_PZ[-1] - DOC.mass_accumulated_overflow[-1] - DOC.mass_accumulated_pipe_outflow[-1] - \
                           DOC.mass_accumulated_infiltrated_to_SZ[-1] - DOC.mass_accumulated_evapotranspiration[-1] - DOC.mass_soil[-1] - DOC.mass_accumulated_reaction[-1]
            DOC.mass_balance.append(Mstor_doc_mb)

            delta_ast_doc_usz = (Mstor_doc_mb - Mstor_doc_ast) / (PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz)
            cl_i1_doc[1] = cl_i1_doc[1] + delta_ast_doc_usz
            if cl_i1_doc[1] > 0:
                cl_i1_doc[1] = cl_i1_doc[1]
            else:
                cl_i1_doc[1] = 0

        ## adding layers of USZ in time
        O2.concentration_USZ.append(cl_i1_o2)
        O2.reaction_rate_USZ.append(Rxl_o2)

        NH4.concentration_USZ.append(cl_i1_nh4)
        NH4.concentration_soil_USZ.append(csi_usz_nh4)
        NH4.reaction_rate_USZ.append(Rxl_nh4)

        NO3.concentration_USZ.append(cl_i1_no3)
        NO3.reaction_rate_USZ.append(Rxl_no3)

        DOC.concentration_USZ.append(cl_i1_doc)
        DOC.concentration_soil_USZ.append(csi_usz_doc)
        DOC.reaction_rate_USZ.append(Rxl_doc)

        ## adding layers of SZ in time
        if GENERAL_PARAMETERS.hpipe > 0:
            O2.concentration_SZ.append(cj_i1_o2)
            O2.reaction_rate_SZ.append(Rxj_o2)

            NH4.concentration_SZ.append(cj_i1_nh4)
            NH4.concentration_soil_SZ.append(csi_sz_nh4)
            NH4.reaction_rate_SZ.append(Rxj_nh4)

            NO3.concentration_SZ.append(cj_i1_no3)
            NO3.reaction_rate_SZ.append(Rxj_no3)

            DOC.concentration_SZ.append(cj_i1_doc)
            DOC.concentration_soil_SZ.append(csi_sz_doc)
            DOC.reaction_rate_SZ.append(Rxj_doc)

    # **5. Transforming in dataframe **
    ## Oxygen
    data_usz_o2 = pd.DataFrame(O2.concentration_USZ)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz' + str(a)
        column_name.append(name)
        a = a + 1
    data_usz_o2.set_axis(column_name, axis='columns', inplace=True)

    # print(data_usz_o2)

    data_sz_o2 = pd.DataFrame(O2.concentration_SZ)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz' + str(a)
        column_name.append(name)
        a = a + 1
    data_sz_o2.set_axis(column_name, axis='columns', inplace=True)

    data_rx_usz_o2 = pd.DataFrame(O2.reaction_rate_USZ)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_usz_o2.set_axis(column_name, axis='columns', inplace=True)

    data_rx_sz_o2 = pd.DataFrame(O2.reaction_rate_SZ)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_sz_o2.set_axis(column_name, axis='columns', inplace=True)

    frames = [data_usz_o2, data_sz_o2, data_rx_usz_o2, data_rx_sz_o2]
    data_o2 = pd.concat((frames), axis=1)

    O2.concentration_PZ.append(0)
    data_o2['pz'] = O2.concentration_PZ

    indice_n = list(range(len(O2.concentration_PZ)))
    O2.concentration_inflow.append(0)
    c_in = O2.concentration_inflow[:len(indice_n)]
    # print('len c_in:', len(c_in), ', len_indice:', len(WFR.indice))
    data_o2['c_in'] = c_in

    data_o2['t'] = indice_n

    ## Amonia
    data_usz_nh4 = pd.DataFrame(NH4.concentration_USZ)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz' + str(a)
        column_name.append(name)
        a = a + 1
    data_usz_nh4.set_axis(column_name, axis='columns', inplace=True)

    data_sz_nh4 = pd.DataFrame(NH4.concentration_SZ)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz' + str(a)
        column_name.append(name)
        a = a + 1
    data_sz_nh4.set_axis(column_name, axis='columns', inplace=True)

    data_s_usz_nh4 = pd.DataFrame(NH4.concentration_soil_USZ)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz_soil' + str(a)
        column_name.append(name)
        a = a + 1
    data_s_usz_nh4.set_axis(column_name, axis='columns', inplace=True)

    data_s_sz_nh4 = pd.DataFrame(NH4.concentration_soil_SZ)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz_soil' + str(a)
        column_name.append(name)
        a = a + 1
    data_s_sz_nh4.set_axis(column_name, axis='columns', inplace=True)

    data_rx_usz_nh4 = pd.DataFrame(NH4.reaction_rate_USZ)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_usz_nh4.set_axis(column_name, axis='columns', inplace=True)

    data_rx_sz_nh4 = pd.DataFrame(NH4.reaction_rate_SZ)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_sz_nh4.set_axis(column_name, axis='columns', inplace=True)

    frames = [data_usz_nh4, data_sz_nh4, data_s_usz_nh4, data_s_sz_nh4, data_rx_usz_nh4, data_rx_sz_nh4]
    data_nh4 = pd.concat((frames), axis=1)
    NH4.concentration_PZ.append(0)
    NH4.concentration_inflow.append(0)
    c_in = NH4.concentration_inflow[:len(indice_n)]
    # print('len c_in:', len(c_in), ', len_indice:', len(WFR.indice))
    data_nh4['c_in'] = c_in

    data_nh4['pz'] = NH4.concentration_PZ
    # WFR.indice.append(len(WFR.indice)+1)
    data_nh4['t'] = indice_n

    ## Nitrate
    data_usz_no3 = pd.DataFrame(NO3.concentration_USZ)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz' + str(a)
        column_name.append(name)
        a = a + 1
    data_usz_no3.set_axis(column_name, axis='columns', inplace=True)

    data_sz_no3 = pd.DataFrame(NO3.concentration_SZ)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz' + str(a)
        column_name.append(name)
        a = a + 1
    data_sz_no3.set_axis(column_name, axis='columns', inplace=True)

    data_rx_usz_no3 = pd.DataFrame(NO3.reaction_rate_USZ)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_usz_no3.set_axis(column_name, axis='columns', inplace=True)

    data_rx_sz_no3 = pd.DataFrame(NO3.reaction_rate_SZ)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_sz_no3.set_axis(column_name, axis='columns', inplace=True)

    frames = [data_usz_no3, data_sz_no3, data_rx_usz_no3, data_rx_sz_no3]
    data_no3 = pd.concat((frames), axis=1)
    NO3.concentration_PZ.append(0)
    data_no3['pz'] = NO3.concentration_PZ

    NO3.concentration_inflow.append(0)
    c_in = NO3.concentration_inflow[:len(indice_n)]
    data_no3['c_in'] = c_in

    data_no3['t'] = indice_n

    ## DOC
    data_usz_doc = pd.DataFrame(DOC.concentration_USZ)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz' + str(a)
        column_name.append(name)
        a = a + 1
    data_usz_doc.set_axis(column_name, axis='columns', inplace=True)

    data_sz_doc = pd.DataFrame(DOC.concentration_SZ)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz' + str(a)
        column_name.append(name)
        a = a + 1
    data_sz_doc.set_axis(column_name, axis='columns', inplace=True)

    data_rx_usz_doc = pd.DataFrame(DOC.reaction_rate_USZ)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_usz_doc.set_axis(column_name, axis='columns', inplace=True)

    data_rx_sz_doc = pd.DataFrame(DOC.reaction_rate_SZ)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_sz_doc.set_axis(column_name, axis='columns', inplace=True)

    frames = [data_usz_doc, data_sz_doc, data_rx_usz_doc, data_rx_sz_doc]
    data_doc = pd.concat((frames), axis=1)
    DOC.concentration_PZ.append(0)
    data_doc['pz'] = DOC.concentration_PZ

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

