#!/usr/bin/env python

import parameters
import results_tests
import datetime
import pandas as pd 


SETUP_FILE = "parameters.ini"
INFLOW_PARAMETERS = parameters.ConcentrationInflow("concentration_inflow.csv")
WF_INFLOW = parameters.WaterInflow("water_inflow.csv")
WFR = parameters.WaterFlowResults(WF_INFLOW.tQrain, WF_INFLOW.tQin, WF_INFLOW.tEmax)
GENERAL_PARAMETERS = parameters.GeneralParameters(SETUP_FILE)
USZ = parameters.UnsaturatedZone(SETUP_FILE, GENERAL_PARAMETERS.L, GENERAL_PARAMETERS.hpipe, GENERAL_PARAMETERS.dz)
PZ = parameters.PondingZone(SETUP_FILE)
SZ = parameters.SaturatedZone(SETUP_FILE, GENERAL_PARAMETERS.n, USZ.m_usz)
SOIL_PLANT = parameters.SoilPlant(SETUP_FILE, USZ.nusz_ini)
NH4 = parameters.Ammonia(USZ.m_usz, SZ.m_sz, SETUP_FILE)
NO3 = parameters.Nitrate(USZ.m_usz, SZ.m_sz, SETUP_FILE)
O2 = parameters.Oxygen(USZ.m_usz, SZ.m_sz, SETUP_FILE)
DOC = parameters.DissolvedOrganicCarbon(USZ.m_usz, SZ.m_sz, SETUP_FILE)


def run_W():
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
        hp = max(hpEND + GENERAL_PARAMETERS.dt / PZ.Ab * (WFR.tQrain[t] + Qin - Qover), 0)  # beginning
        Qinfp = PZ.f_infiltration_to_surrounding(USZ.Kf, USZ.A, hp)
        Qpf = PZ.f_infiltration_to_filter_material(USZ.Ks, hp, husz, USZ.A, GENERAL_PARAMETERS.dt, s, nusz, Qinfp)
        hpEND = max(hp - GENERAL_PARAMETERS.dt / PZ.Ab * (Qpf + Qinfp), 0)  # end

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

        WFR.thp.append(hp)
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
def run_Kin():
    # O2
    cp_o2 = []
    c_usz_o2 = []
    cs_usz_o2 = []
    c_sz_o2 = []
    cs_sz_o2 = []

    Rx_usz_o2 = []
    Rx_sz_o2 = []
    cp_a_o2 = 0

    c0_usz_o2 = [0]
    c0_usz_o2 = c0_usz_o2 * USZ.m_usz
    c_usz_o2.append(c0_usz_o2)
    cs_usz_o2.append(c0_usz_o2)
    c0_sz_o2 = [0]
    c0_sz_o2 = c0_sz_o2 * SZ.m_sz
    c_sz_o2.append(c0_sz_o2)
    cs_sz_o2.append(c0_sz_o2)
    Rx_usz_o2.append(c0_usz_o2)
    Rx_sz_o2.append(c0_sz_o2)

    Mstor_o2_ast_list = []
    Msoil_o2_acum = []
    MRx_o2_acum = []
    Min_o2_acum = []
    Mover_o2_acum = []
    Mpipe_o2_acum = []
    Minfsz_o2_acum = []
    Met_o2_acum = []
    Mpz_o2_list = []
    Mstor_o2_mb_list = []

    # NH4
    cp_nh4 = []
    c_usz_nh4 = []
    cs_usz_nh4 = []
    c_sz_nh4 = []
    cs_sz_nh4 = []

    Rx_usz_nh4 = []
    Rx_sz_nh4 = []

    cp_a_nh4 = 0

    cs_usz_nh4_a = 0.0
    cs_sz_nh4_a = 0

    c0_usz_nh4 = [0]
    c0_usz_nh4 = c0_usz_nh4 * USZ.m_usz
    cs0_usz_nh4 = [cs_usz_nh4_a]
    cs0_usz_nh4 = cs0_usz_nh4 * USZ.m_usz
    c_usz_nh4.append(c0_usz_nh4)
    cs_usz_nh4.append(cs0_usz_nh4)
    c0_sz_nh4 = [0]
    c0_sz_nh4 = c0_sz_nh4 * SZ.m_sz
    c_sz_nh4.append(c0_sz_nh4)
    cs_sz_nh4.append(c0_sz_nh4)
    Rx_usz_nh4.append(c0_usz_nh4)
    Rx_sz_nh4.append(c0_sz_nh4)

    Mstor_nh4_ast_list = []
    Msoil_nh4_acum = []
    MRx_nh4_acum = []
    Min_nh4_acum = []
    Mover_nh4_acum = []
    Mpipe_nh4_acum = []
    Minfsz_nh4_acum = []
    Met_nh4_acum = []
    Mpz_nh4_list = []
    Mstor_nh4_mb_list = []

    # NO3
    cp_no3 = []
    c_usz_no3 = []
    cs_usz_no3 = []
    c_sz_no3 = []
    cs_sz_no3 = []

    Rx_usz_no3 = []
    Rx_sz_no3 = []

    cp_a_no3 = 0

    c0_usz_no3 = [0]
    c0_usz_no3 = c0_usz_no3 * USZ.m_usz
    c_usz_no3.append(c0_usz_no3)
    cs_usz_no3.append(c0_usz_no3)
    c0_sz_no3 = [0]
    c0_sz_no3 = c0_sz_no3 * SZ.m_sz
    c_sz_no3.append(c0_sz_no3)
    cs_sz_no3.append(c0_sz_no3)
    Rx_usz_no3.append(c0_usz_no3)
    Rx_sz_no3.append(c0_sz_no3)

    Mstor_no3_ast_list = []
    Msoil_no3_acum = []
    MRx_no3_acum = []
    Min_no3_acum = []
    Mover_no3_acum = []
    Mpipe_no3_acum = []
    Minfsz_no3_acum = []
    Met_no3_acum = []
    Mpz_no3_list = []
    Mstor_no3_mb_list = []

    # DOC
    cp_doc = []
    c_usz_doc = []
    cs_usz_doc = []
    c_sz_doc = []
    cs_sz_doc = []

    Rx_usz_doc = []
    Rx_sz_doc = []

    cp_a_doc = 0

    c0_usz_doc = [0]
    c0_usz_doc = c0_usz_doc * USZ.m_usz
    c_usz_doc.append(c0_usz_doc)
    cs_usz_doc.append(c0_usz_doc)
    c0_sz_doc = [0]
    c0_sz_doc = c0_sz_doc * SZ.m_sz
    c_sz_doc.append(c0_sz_doc)
    cs_sz_doc.append(c0_sz_doc)
    Rx_usz_doc.append(c0_usz_doc)
    Rx_sz_doc.append(c0_sz_doc)

    cs_usz_doc_a = 0
    cs_sz_doc_a = 0

    Mstor_doc_ast_list = []
    Msoil_doc_acum = []
    MRx_doc_acum = []
    Min_doc_acum = []
    Mover_doc_acum = []
    Mpipe_doc_acum = []
    Minfsz_doc_acum = []
    Met_doc_acum = []
    Mpz_doc_list = []
    Mstor_doc_mb_list = []

    for t in range(len(WFR.indice) - 1):
        # Ponding zone
        cin_o2 = INFLOW_PARAMETERS.cin_o2[t]
        cin_nh4 = INFLOW_PARAMETERS.cin_nh4[t]
        cin_no3 = INFLOW_PARAMETERS.cin_no3[t]
        cin_doc = INFLOW_PARAMETERS.cin_doc[t]

        # print('t', t)
        hp = WFR.thpEND[t]
        if hp < 0.001:
            hp = 0
        else:
            hp = hp

        if t == 0:
            hp_a = 0
        else:
            hp_a = WFR.thpEND[t - 1]
        # print('hp', hp)
        Qin_p = WFR.tQin[t]
        I1 = WFR.tQpf[t]
        Qv = WFR.tQover[t]
        Qet_1 = WFR.tQet1[t]
        I2 = WFR.tQfs[t]
        Qhc = WFR.tQhc[t]
        Qet_2 = WFR.tQet2[t]
        if t == 0:
            Qorif = 0
        else:
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

        if hp == 0:
            cpi_o2 = 0
            cpi_nh4 = 0
            cpi_no3 = 0
            cpi_doc = 0

        else:
            cpi_o2 = PZ.f_concentration(cin_o2, Qin_p, cp_a_o2, I1, Qv, Rxi_p_o2, hp, hp_a, PZ.Ab, GENERAL_PARAMETERS.dt)
            cpi_nh4 = PZ.f_concentration(cin_nh4, Qin_p, cp_a_nh4, I1, Qv, Rxi_p_nh4, hp, hp_a, PZ.Ab, GENERAL_PARAMETERS.dt)
            cpi_no3 = PZ.f_concentration(cin_no3, Qin_p, cp_a_no3, I1, Qv, Rxi_p_no3, hp, hp_a, PZ.Ab, GENERAL_PARAMETERS.dt)
            cpi_doc = PZ.f_concentration(cin_doc, Qin_p, cp_a_doc, I1, Qv, Rxi_p_doc, hp, hp_a, PZ.Ab, GENERAL_PARAMETERS.dt)
        #             print('cin_nh4, Qin_p, cp_a_nh4, I1, Qv, Rxi_p_nh4, hp, hp_a: ', cin_nh4, Qin_p, cp_a_nh4, I1, Qv, Rxi_p_nh4, hp, hp_a)
        #             print('cpi_nh4: ', cpi_nh4)

        if cpi_o2 < 0.00002:
            cpi_o2 = 0
        else:
            cpi_o2 = cpi_o2
        if cpi_nh4 < 0.00002:
            cpi_nh4 = 0
        else:
            cpi_nh4 = cpi_nh4
        if cpi_no3 < 0.00002:
            cpi_no3 = 0
        else:
            cpi_no3 = cpi_no3
        if cpi_doc < 0.00002:
            cpi_doc = 0
        else:
            cpi_doc = cpi_doc

        cp_o2.append(cpi_o2)
        cp_a_o2 = cp_o2[-1]

        cp_nh4.append(cpi_nh4)
        cp_a_nh4 = cp_nh4[-1]

        cp_no3.append(cpi_no3)
        cp_a_no3 = cp_no3[-1]

        cp_doc.append(cpi_doc)
        cp_a_doc = cp_doc[-1]

        # USZ
        cli_o2_list = c_usz_o2[t].copy()
        # print('1_o2', cli_o2_list)
        cli_nh4_list = c_usz_nh4[t].copy()
        # print('1_nh4', cli_nh4_list)
        cli_no3_list = c_usz_no3[t].copy()
        cli_doc_list = c_usz_doc[t].copy()

        # SZ
        if GENERAL_PARAMETERS.hpipe > 0:
            cji_o2_list = c_sz_o2[t].copy()  # copiar lista
            # cji_o2_list = c_sz_o2[t][:] #outro jeito de copiar a lista
            cji_nh4_list = c_sz_nh4[t].copy()
            cji_no3_list = c_sz_no3[t].copy()
            cji_doc_list = c_sz_doc[t].copy()

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

                # since we have added an initial value of cs = 0 in the cs_usz_nh4 list, when calling the index equal = 't' we are actually calling the value corresponding to t-1 
                cs_nh4 = cs_usz_nh4[t][l]
                cs_nh4_iplus1 = NH4.f_concentration_soil(cs_nh4, teta_sm_i, NH4.kads, cl_nh4, NH4.kdes, NH4.k_mb, SOIL_PLANT.ro, GENERAL_PARAMETERS.dt)
                if cs_nh4_iplus1 < 0.00000000000001:
                    cs_nh4_iplus1 = 0

                #                 sor = (teta_sm_i/SOIL_PLANT.ro)*NH4.kads*cl_nh4*GENERAL_PARAMETERS.dt
                #                 des =  NH4.kdes*cs_usz_nh4_a*GENERAL_PARAMETERS.dt

                #                 print('t: ', t, ', l: ', l , ', cs_nh4_a: ', cs_usz_nh4_a, 'teta_sm_i: ', teta_sm_i, 'NH4.kads: ', NH4.kads, 'cl_nh4: ', cl_nh4, 'NH4.kdes: ', NH4.kdes, 'NH4.k_mb: ', NH4.k_mb)
                #                 print('cs_nh4: ', cs_nh4)
                #                 input()

                #                 logger.debug(SOIL_PLANT.f' t: {t} \t l: {l} \t sor: {sor} \t des: {des} \t cs_nh4: {cs_nh4}')
                #                 logger.debug('----SOIL----')
                #                 logger.debug('cs_nh4_a: ', cs_usz_nh4_a)
                #                 logger.debug('cs_nh4_a: ', cs_usz_nh4_a, 'teta_sm_i: ', teta_sm_i, 'NH4.kads: ', NH4.kads, 'cl_nh4: ', cl_nh4, 'NH4.kdes: ', NH4.kdes, 'NH4.k_mb: ', NH4.k_mb)
                #                 logger.debug('cs_nh4: ', cs_nh4)
                #                 logger.debug('----WATER----')
                csi_usz_nh4.append(cs_nh4_iplus1)

                cs_no3 = 0

                cs_usz_doc_a = cs_usz_doc[t][l]
                cs_doc = DOC.f_concentration_soil(cs_usz_doc_a, teta_sm_i, DOC.kads, cl_doc, DOC.kdes, DOC.k_mb, SOIL_PLANT.ro, GENERAL_PARAMETERS.dt)
                csi_usz_doc.append(cs_doc)
                if cs_doc < 0.00000000000001:
                    cs_doc = 0

                UF_usz = []
                # Qet_1 = 0                        
                alfa, beta = USZ.f_alfa_beta(l)
                UFi_usz = USZ.f_unit_flux(alfa, beta, I1, Qet_1, I2, Qhc, Qorif, Qinf_sz, teta_sm_i, GENERAL_PARAMETERS.hpipe, PZ.Ab)
                #                 print(t, '  ', l)
                #                 print('UFi_usz: ', UFi_usz, ', alfa: ', alfa, ', beta: ', beta,  ', I1: ', I1, ', Qet_1:', Qet_1, ', I2: ', I2, ', Qhc: ', Qhc, ', Qorif: ', Qorif, ', Qinf_sz: ', Qinf_sz)            
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
                # print('t: ', t, 'l: ', l, 'Rx_usz_o2: ', Rxi_2_o2, 'Rx_usz_nh4: ', Rxi_2_nh4, 'Rx_usz_no3: ', Rxi_2_no3)   

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

            cs_sz_nh4_a = cs_sz_nh4[t][j]
            cs_nh4 = NH4.f_concentration_soil(cs_sz_nh4_a, teta_b_i, NH4.kads2, cj_nh4, NH4.kdes2, NH4.k_mb, SOIL_PLANT.ro, GENERAL_PARAMETERS.dt)
            csi_sz_nh4.append(cs_nh4)
            if cs_nh4 < 0.00000000000001:
                cs_nh4 = 0

            #                 print('====>> t: ', t, ', j: ', j)
            #                 print('---- SOIL ----')
            #                 print('cs_nh4_a: ', cs_sz_nh4_a, 'teta_b_i: ', teta_b_i, 'NH4.kads: ', NH4.kads, 'cj_nh4: ', cj_nh4, 'NH4.kdes: ', NH4.kdes, 'NH4.k_mb: ', NH4.k_mb)
            #                 print('cs_nh4: ', cs_nh4)
            #                 print('---- WATER ----')
            #                 print('cmminus1_nh4: ', cmminus1_nh4)                
            cs_no3 = 0

            cs_sz_doc_a = cs_sz_doc[t][j]
            cs_doc = DOC.f_concentration_soil(cs_sz_doc_a, teta_b_i, DOC.kads2, cj_doc, DOC.kdes2, DOC.k_mb, SOIL_PLANT.ro, GENERAL_PARAMETERS.dt)
            csi_sz_doc.append(cs_doc)

            UF_sz = []
            alfa2, beta2 = SZ.f_alfa_beta(j)
            # Qet_2 = 0
            UFi_sz = SZ.f_unit_flux(alfa2, beta2, I2, Qhc, Qet_2, Qorif, Qinf_sz, teta_b_i, PZ.Ab)
            #                 print('I2, Qhc, Qet_2, Qorif, Qinf_sz:', I2, Qhc, Qet_2, Qorif, Qinf_sz)
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
            # print('t: ', t, 'j: ', j, 'Rx_sz_o2: ', Rxi_3_o2, 'Rx_sz_nh4: ', Rxi_3_nh4, 'Rx_sz_no3: ', Rxi_3_no3)

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
            Mstor_o2_ast = sum(c_usz_o2[t]) * PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz + sum(c_sz_o2[t]) * PZ.Ab * \
                           WFR.tteta_sz[t] * WFR.thsz[t] * 1000 / SZ.m_sz
            Mstor_o2_ast_list.append(Mstor_o2_ast)

            Msoil_o2_a = 0
            Msoil_o2 = 0
            Msoil_o2_acum.append(0)

            MRx_o2 = - (sum(Rx_usz_o2[t]) * PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz + sum(Rx_sz_o2[t]) * PZ.Ab *
                        WFR.tteta_sz[t] * WFR.thsz[t] * 1000 / SZ.m_sz)
            if t == 0:
                MRx_o2_acum.append(MRx_o2)
            else:
                MRx_o2_acum.append(MRx_o2 + MRx_o2_acum[-1])

            Min_o2 = WFR.tQin[t] * INFLOW_PARAMETERS.cin_o2[t] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Min_o2_acum.append(Min_o2)
            else:
                Min_o2_acum.append(Min_o2 + Min_o2_acum[-1])

            Mover_o2 = WFR.tQover[t] * cp_o2[t] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Mover_o2_acum.append(Mover_o2)
            else:
                Mover_o2_acum.append(Mover_o2 + Mover_o2_acum[-1])

            Mpipe_o2 = WFR.tQpipe[t] * c_sz_o2[t][SZ.m_sz - 1] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Mpipe_o2_acum.append(Mpipe_o2)
            else:
                Mpipe_o2_acum.append(Mpipe_o2 + Mpipe_o2_acum[-1])

            Minfsz_o2 = WFR.tQinfsz[t] * c_sz_o2[t][SZ.m_sz - 1] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Minfsz_o2_acum.append(Minfsz_o2)
            else:
                Minfsz_o2_acum.append(Minfsz_o2 + Minfsz_o2_acum[-1])

            Met_o2 = WFR.tQet[t] * cl_i1_o2[0] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Met_o2_acum.append(Met_o2)
            else:
                Met_o2_acum.append(Met_o2 + Met_o2_acum[-1])

            Mpz_o2 = WFR.thpEND[t] * PZ.Ab * 1000 * cp_o2[t]
            Mpz_o2_list.append(Mpz_o2)

            Mstor_o2_mb = Min_o2_acum[-1] - Mpz_o2_list[-1] - Mover_o2_acum[-1] - Mpipe_o2_acum[-1] - Minfsz_o2_acum[
                -1] - Met_o2_acum[-1] - Msoil_o2_acum[-1] - MRx_o2_acum[-1]
            Mstor_o2_mb_list.append(Mstor_o2_mb)

            delta_ast_o2_usz = (Mstor_o2_mb - Mstor_o2_ast) / (PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz)
            cl_i1_o2[1] = cl_i1_o2[1] + delta_ast_o2_usz
            if cl_i1_o2[1] > 0:
                cl_i1_o2[1] = cl_i1_o2[1]
            else:
                cl_i1_o2[1] = 0

            ### Amonia

            Mstor_nh4_ast = sum(c_usz_nh4[t]) * PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz + sum(c_sz_nh4[t]) * PZ.Ab * \
                            WFR.tteta_sz[t] * WFR.thsz[t] * 1000 / SZ.m_sz
            Mstor_nh4_ast_list.append(Mstor_nh4_ast)

            Msoil_nh4_a = cs_usz_nh4_a * SOIL_PLANT.ro * PZ.Ab * WFR.thusz[t] * 1000 + cs_sz_nh4_a * SOIL_PLANT.ro * PZ.Ab * WFR.thsz[t] * 1000
            Msoil_nh4 = sum(cs_usz_nh4[t]) * SOIL_PLANT.ro * PZ.Ab * WFR.thusz[t] * 1000 / USZ.m_usz + sum(cs_sz_nh4[t]) * SOIL_PLANT.ro * PZ.Ab * WFR.thsz[
                t] * 1000 / SZ.m_sz
            Msoil_nh4_acum.append(Msoil_nh4 - Msoil_nh4_a)

            MRx_nh4 = - (sum(Rx_usz_nh4[t]) * PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz + sum(Rx_sz_nh4[t]) * PZ.Ab *
                         WFR.tteta_sz[t] * WFR.thsz[t] * 1000 / SZ.m_sz)
            if t == 0:
                MRx_nh4_acum.append(MRx_nh4)
            else:
                MRx_nh4_acum.append(MRx_nh4 + MRx_nh4_acum[-1])

            Min_nh4 = WFR.tQin[t] * INFLOW_PARAMETERS.cin_nh4[t] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Min_nh4_acum.append(Min_nh4)
            else:
                Min_nh4_acum.append(Min_nh4 + Min_nh4_acum[-1])

            Mover_nh4 = WFR.tQover[t] * cp_nh4[t] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Mover_nh4_acum.append(Mover_nh4)
            else:
                Mover_nh4_acum.append(Mover_nh4 + Mover_nh4_acum[-1])

            Mpipe_nh4 = WFR.tQpipe[t] * c_sz_nh4[t][SZ.m_sz - 1] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Mpipe_nh4_acum.append(Mpipe_nh4)
            else:
                Mpipe_nh4_acum.append(Mpipe_nh4 + Mpipe_nh4_acum[-1])

            Minfsz_nh4 = WFR.tQinfsz[t] * c_sz_nh4[t][SZ.m_sz - 1] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Minfsz_nh4_acum.append(Minfsz_nh4)
            else:
                Minfsz_nh4_acum.append(Minfsz_nh4 + Minfsz_nh4_acum[-1])

            Met_nh4 = WFR.tQet[t] * cl_i1_nh4[0] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Met_nh4_acum.append(Met_nh4)
            else:
                Met_nh4_acum.append(Met_nh4 + Met_nh4_acum[-1])

            Mpz_nh4 = WFR.thpEND[t] * PZ.Ab * 1000 * cp_nh4[t]
            Mpz_nh4_list.append(Mpz_nh4)

            Mstor_nh4_mb = Min_nh4_acum[-1] - Mpz_nh4_list[-1] - Mover_nh4_acum[-1] - Mpipe_nh4_acum[-1] - \
                           Minfsz_nh4_acum[-1] - Met_nh4_acum[-1] - Msoil_nh4_acum[-1] - MRx_nh4_acum[-1]
            Mstor_nh4_mb_list.append(Mstor_nh4_mb)

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
            #         sum_cusz_nh4 = sum(c_usz_nh4[t])
            #         sum_csz_nh4 = sum(c_sz_nh4[t])        
            #         print('t: ', t, 'Mstor_trace: ', Mstor_trace, ', Mstor_mb: ', Mstor_mb, ', dif: ', (Mstor_mb - Mstor_trace))  
            #         print('t: ', t, ', delta_usz: ', delta_trace_nh4_usz, ', delta_sz: ', delta_trace_nh4_sz)
            #         print('teta_sm: ', WFR.tteta_usz[t], ', husz: ', WFR.thusz[t], ', teta_b: ', WFR.tteta_sz[t], ', hsz: ', WFR.thsz[t])
            #         print('sum_cusz_nh4: ', sum_cusz_nh4, ', sum_csz_nh4: ', sum_csz_nh4)             
            #         print('Qin: ', WFR.tQin[t], ',hpEND_i: ', WFR.thpEND[t], ', hpEND_i-1: ', WFR.thpEND[t-1], ',Qpipe: ', WFR.tQpipe[t], ', Qover: ', WFR.tQover[t], ', Qinfsz: ', WFR.tQinfsz[t], ', Qet: ', WFR.tQet[t])
            #         print('Cin: ', INFLOW_PARAMETERS.cin_nh4[t], ', Cp_i: ', cp_nh4[t], ', Cp_i-1: ', cp_nh4[t-1], ', Csz_ultimo: ', c_pipe)
            #         input()
            #                      
            #         for l in range(USZ.m_usz):
            #             cl_i1_nh4[l] = cl_i1_nh4[l] + delta_ast_nh4_usz
            #         for j in range(SZ.m_sz):
            #             cj_i1_nh4[j] = cj_i1_nh4[j] + delta_ast_nh4_sz 

            ### Nitrate

            Mstor_no3_ast = sum(c_usz_no3[t]) * PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz + sum(c_sz_no3[t]) * PZ.Ab * \
                            WFR.tteta_sz[t] * WFR.thsz[t] * 1000 / SZ.m_sz
            Mstor_no3_ast_list.append(Mstor_no3_ast)

            Msoil_no3_a = 0
            Msoil_no3 = 0
            Msoil_no3_acum.append(0)

            MRx_no3 = - (sum(Rx_usz_no3[t]) * PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz + sum(Rx_sz_no3[t]) * PZ.Ab *
                         WFR.tteta_sz[t] * WFR.thsz[t] * 1000 / SZ.m_sz)
            if t == 0:
                MRx_no3_acum.append(MRx_no3)
            else:
                MRx_no3_acum.append(MRx_no3 + MRx_no3_acum[-1])

            Min_no3 = WFR.tQin[t] * INFLOW_PARAMETERS.cin_no3[t] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Min_no3_acum.append(Min_no3)
            else:
                Min_no3_acum.append(Min_no3 + Min_no3_acum[-1])

            Mover_no3 = WFR.tQover[t] * cp_no3[t] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Mover_no3_acum.append(Mover_no3)
            else:
                Mover_no3_acum.append(Mover_no3 + Mover_no3_acum[-1])

            Mpipe_no3 = WFR.tQpipe[t] * c_sz_no3[t][SZ.m_sz - 1] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Mpipe_no3_acum.append(Mpipe_no3)
            else:
                Mpipe_no3_acum.append(Mpipe_no3 + Mpipe_no3_acum[-1])

            Minfsz_no3 = WFR.tQinfsz[t] * c_sz_no3[t][SZ.m_sz - 1] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Minfsz_no3_acum.append(Minfsz_no3)
            else:
                Minfsz_no3_acum.append(Minfsz_no3 + Minfsz_no3_acum[-1])

            Met_no3 = WFR.tQet[t] * cl_i1_no3[0] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Met_no3_acum.append(Met_no3)
            else:
                Met_no3_acum.append(Met_no3 + Met_no3_acum[-1])

            Mpz_no3 = WFR.thpEND[t] * PZ.Ab * 1000 * cp_no3[t]
            Mpz_no3_list.append(Mpz_no3)

            Mstor_no3_mb = Min_no3_acum[-1] - Mpz_no3_list[-1] - Mover_no3_acum[-1] - Mpipe_no3_acum[-1] - \
                           Minfsz_no3_acum[-1] - Met_no3_acum[-1] - Msoil_no3_acum[-1] - MRx_no3_acum[-1]
            Mstor_no3_mb_list.append(Mstor_no3_mb)

            delta_ast_no3_usz = (Mstor_no3_mb - Mstor_no3_ast) / (PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz)
            cl_i1_no3[1] = cl_i1_no3[1] + delta_ast_no3_usz
            if cl_i1_no3[1] > 0:
                cl_i1_no3[1] = cl_i1_no3[1]
            else:
                cl_i1_no3[1] = 0

                ### DOC

            Mstor_doc_ast = sum(c_usz_doc[t]) * PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz + sum(c_sz_doc[t]) * PZ.Ab * \
                            WFR.tteta_sz[t] * WFR.thsz[t] * 1000 / SZ.m_sz
            Mstor_doc_ast_list.append(Mstor_doc_ast)

            Msoil_doc_a = cs_usz_doc_a * SOIL_PLANT.ro * PZ.Ab * WFR.thusz[t] * 1000 + cs_sz_doc_a * SOIL_PLANT.ro * PZ.Ab * WFR.thsz[t] * 1000
            Msoil_doc = sum(cs_usz_doc[t]) * SOIL_PLANT.ro * PZ.Ab * WFR.thusz[t] * 1000 / USZ.m_usz + sum(cs_sz_doc[t]) * SOIL_PLANT.ro * PZ.Ab * WFR.thsz[
                t] * 1000 / SZ.m_sz
            Msoil_doc_acum.append(Msoil_doc - Msoil_doc_a)

            MRx_doc = - (sum(Rx_usz_doc[t]) * PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz + sum(Rx_sz_doc[t]) * PZ.Ab *
                         WFR.tteta_sz[t] * WFR.thsz[t] * 1000 / SZ.m_sz)
            if t == 0:
                MRx_doc_acum.append(MRx_doc)
            else:
                MRx_doc_acum.append(MRx_doc + MRx_doc_acum[-1])

            Min_doc = WFR.tQin[t] * INFLOW_PARAMETERS.cin_doc[t] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Min_doc_acum.append(Min_doc)
            else:
                Min_doc_acum.append(Min_doc + Min_doc_acum[-1])

            Mover_doc = WFR.tQover[t] * cp_doc[t] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Mover_doc_acum.append(Mover_doc)
            else:
                Mover_doc_acum.append(Mover_doc + Mover_doc_acum[-1])

            Mpipe_doc = WFR.tQpipe[t] * c_sz_doc[t][SZ.m_sz - 1] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Mpipe_doc_acum.append(Mpipe_doc)
            else:
                Mpipe_doc_acum.append(Mpipe_doc + Mpipe_doc_acum[-1])

            Minfsz_doc = WFR.tQinfsz[t] * c_sz_doc[t][SZ.m_sz - 1] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Minfsz_doc_acum.append(Minfsz_doc)
            else:
                Minfsz_doc_acum.append(Minfsz_doc + Minfsz_doc_acum[-1])

            Met_doc = WFR.tQet[t] * cl_i1_doc[0] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Met_doc_acum.append(Met_doc)
            else:
                Met_doc_acum.append(Met_doc + Met_doc_acum[-1])

            Mpz_doc = WFR.thpEND[t] * PZ.Ab * 1000 * cp_doc[t]
            Mpz_doc_list.append(Mpz_doc)

            Mstor_doc_mb = Min_doc_acum[-1] - Mpz_doc_list[-1] - Mover_doc_acum[-1] - Mpipe_doc_acum[-1] - \
                           Minfsz_doc_acum[-1] - Met_doc_acum[-1] - Msoil_doc_acum[-1] - MRx_doc_acum[-1]
            Mstor_doc_mb_list.append(Mstor_doc_mb)

            delta_ast_doc_usz = (Mstor_doc_mb - Mstor_doc_ast) / (PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz)
            cl_i1_doc[1] = cl_i1_doc[1] + delta_ast_doc_usz
            if cl_i1_doc[1] > 0:
                cl_i1_doc[1] = cl_i1_doc[1]
            else:
                cl_i1_doc[1] = 0


        else:  # if GENERAL_PARAMETERS.hpipe == 0

            ### Oxygen

            Mstor_o2_ast = sum(c_usz_o2[t]) * PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz
            Mstor_o2_ast_list.append(Mstor_o2_ast)

            Msoil_o2_a = 0
            Msoil_o2 = 0
            Msoil_o2_acum.append(0)

            MRx_o2 = - (sum(Rx_usz_o2[t]) * PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz)
            if t == 0:
                MRx_o2_acum.append(MRx_o2)
            else:
                MRx_o2_acum.append(MRx_o2 + MRx_o2_acum[-1])

            Min_o2 = WFR.tQin[t] * INFLOW_PARAMETERS.cin_o2[t] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Min_o2_acum.append(Min_o2)
            else:
                Min_o2_acum.append(Min_o2 + Min_o2_acum[-1])

            Mover_o2 = WFR.tQover[t] * cp_o2[t] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Mover_o2_acum.append(Mover_o2)
            else:
                Mover_o2_acum.append(Mover_o2 + Mover_o2_acum[-1])

            Mpipe_o2 = WFR.tQpipe[t] * c_usz_o2[t][USZ.m_usz - 1] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Mpipe_o2_acum.append(Mpipe_o2)
            else:
                Mpipe_o2_acum.append(Mpipe_o2 + Mpipe_o2_acum[-1])

            Minfsz_o2 = WFR.tQinfsz[t] * c_usz_o2[t][USZ.m_usz - 1] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Minfsz_o2_acum.append(Minfsz_o2)
            else:
                Minfsz_o2_acum.append(Minfsz_o2 + Minfsz_o2_acum[-1])

            Met_o2 = WFR.tQet[t] * cl_i1_o2[0] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Met_o2_acum.append(Met_o2)
            else:
                Met_o2_acum.append(Met_o2 + Met_o2_acum[-1])

            Mpz_o2 = WFR.thpEND[t] * PZ.Ab * 1000 * cp_o2[t]
            Mpz_o2_list.append(Mpz_o2)

            Mstor_o2_mb = Min_o2_acum[-1] - Mpz_o2_list[-1] - Mover_o2_acum[-1] - Mpipe_o2_acum[-1] - Minfsz_o2_acum[
                -1] - Met_o2_acum[-1] - Msoil_o2_acum[-1] - MRx_o2_acum[-1]
            Mstor_o2_mb_list.append(Mstor_o2_mb)

            delta_ast_o2_usz = (Mstor_o2_mb - Mstor_o2_ast) / (PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz)
            cl_i1_o2[1] = cl_i1_o2[1] + delta_ast_o2_usz
            if cl_i1_o2[1] > 0:
                cl_i1_o2[1] = cl_i1_o2[1]
            else:
                cl_i1_o2[1] = 0

            ### Amonia

            Mstor_nh4_ast = sum(c_usz_nh4[t]) * PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz
            Mstor_nh4_ast_list.append(Mstor_nh4_ast)

            Msoil_nh4_a = cs_usz_nh4_a * SOIL_PLANT.ro * PZ.Ab * WFR.thusz[t] * 1000
            Msoil_nh4 = sum(cs_usz_nh4[t]) * SOIL_PLANT.ro * PZ.Ab * WFR.thusz[t] * 1000 / USZ.m_usz
            Msoil_nh4_acum.append(Msoil_nh4 - Msoil_nh4_a)

            MRx_nh4 = - (sum(Rx_usz_nh4[t]) * PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz)
            if t == 0:
                MRx_nh4_acum.append(MRx_nh4)
            else:
                MRx_nh4_acum.append(MRx_nh4 + MRx_nh4_acum[-1])

            Min_nh4 = WFR.tQin[t] * INFLOW_PARAMETERS.cin_nh4[t] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Min_nh4_acum.append(Min_nh4)
            else:
                Min_nh4_acum.append(Min_nh4 + Min_nh4_acum[-1])

            Mover_nh4 = WFR.tQover[t] * cp_nh4[t] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Mover_nh4_acum.append(Mover_nh4)
            else:
                Mover_nh4_acum.append(Mover_nh4 + Mover_nh4_acum[-1])

            Mpipe_nh4 = WFR.tQpipe[t] * c_usz_nh4[t][USZ.m_usz - 1] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Mpipe_nh4_acum.append(Mpipe_nh4)
            else:
                Mpipe_nh4_acum.append(Mpipe_nh4 + Mpipe_nh4_acum[-1])

            Minfsz_nh4 = WFR.tQinfsz[t] * c_usz_nh4[t][USZ.m_usz - 1] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Minfsz_nh4_acum.append(Minfsz_nh4)
            else:
                Minfsz_nh4_acum.append(Minfsz_nh4 + Minfsz_nh4_acum[-1])

            Met_nh4 = WFR.tQet[t] * cl_i1_nh4[0] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Met_nh4_acum.append(Met_nh4)
            else:
                Met_nh4_acum.append(Met_nh4 + Met_nh4_acum[-1])

            Mpz_nh4 = WFR.thpEND[t] * PZ.Ab * 1000 * cp_nh4[t]
            Mpz_nh4_list.append(Mpz_nh4)

            Mstor_nh4_mb = Min_nh4_acum[-1] - Mpz_nh4_list[-1] - Mover_nh4_acum[-1] - Mpipe_nh4_acum[-1] - \
                           Minfsz_nh4_acum[-1] - Met_nh4_acum[-1] - Msoil_nh4_acum[-1] - MRx_nh4_acum[-1]
            Mstor_nh4_mb_list.append(Mstor_nh4_mb)

            delta_ast_nh4_usz = (Mstor_nh4_mb - Mstor_nh4_ast) / (PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz)
            cl_i1_nh4[1] = cl_i1_nh4[1] + delta_ast_nh4_usz
            if cl_i1_nh4[1] > 0:
                cl_i1_nh4[1] = cl_i1_nh4[1]
            else:
                cl_i1_nh4[1] = 0

            ### Nitrate

            Mstor_no3_ast = sum(c_usz_no3[t]) * PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz
            Mstor_no3_ast_list.append(Mstor_no3_ast)

            Msoil_no3_a = 0
            Msoil_no3 = 0
            Msoil_no3_acum.append(0)

            MRx_no3 = - (sum(Rx_usz_no3[t]) * PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz)
            if t == 0:
                MRx_no3_acum.append(MRx_no3)
            else:
                MRx_no3_acum.append(MRx_no3 + MRx_no3_acum[-1])

            Min_no3 = WFR.tQin[t] * INFLOW_PARAMETERS.cin_no3[t] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Min_no3_acum.append(Min_no3)
            else:
                Min_no3_acum.append(Min_no3 + Min_no3_acum[-1])

            Mover_no3 = WFR.tQover[t] * cp_no3[t] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Mover_no3_acum.append(Mover_no3)
            else:
                Mover_no3_acum.append(Mover_no3 + Mover_no3_acum[-1])

            Mpipe_no3 = WFR.tQpipe[t] * c_usz_no3[t][USZ.m_usz - 1] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Mpipe_no3_acum.append(Mpipe_no3)
            else:
                Mpipe_no3_acum.append(Mpipe_no3 + Mpipe_no3_acum[-1])

            Minfsz_no3 = WFR.tQinfsz[t] * c_usz_no3[t][USZ.m_usz - 1] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Minfsz_no3_acum.append(Minfsz_no3)
            else:
                Minfsz_no3_acum.append(Minfsz_no3 + Minfsz_no3_acum[-1])

            Met_no3 = WFR.tQet[t] * cl_i1_no3[0] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Met_no3_acum.append(Met_no3)
            else:
                Met_no3_acum.append(Met_no3 + Met_no3_acum[-1])

            Mpz_no3 = WFR.thpEND[t] * PZ.Ab * 1000 * cp_no3[t]
            Mpz_no3_list.append(Mpz_no3)

            Mstor_no3_mb = Min_no3_acum[-1] - Mpz_no3_list[-1] - Mover_no3_acum[-1] - Mpipe_no3_acum[-1] - \
                           Minfsz_no3_acum[-1] - Met_no3_acum[-1] - Msoil_no3_acum[-1] - MRx_no3_acum[-1]
            Mstor_no3_mb_list.append(Mstor_no3_mb)

            delta_ast_no3_usz = (Mstor_no3_mb - Mstor_no3_ast) / (PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz)
            cl_i1_no3[1] = cl_i1_no3[1] + delta_ast_no3_usz
            if cl_i1_no3[1] > 0:
                cl_i1_no3[1] = cl_i1_no3[1]
            else:
                cl_i1_no3[1] = 0

                ### DOC

            Mstor_doc_ast = sum(c_usz_doc[t]) * PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz
            Mstor_doc_ast_list.append(Mstor_doc_ast)

            Msoil_doc_a = cs_usz_doc_a * SOIL_PLANT.ro * PZ.Ab * WFR.thusz[t] * 1000
            Msoil_doc = sum(cs_usz_doc[t]) * SOIL_PLANT.ro * PZ.Ab * WFR.thusz[t] * 1000 / USZ.m_usz
            Msoil_doc_acum.append(Msoil_doc - Msoil_doc_a)

            MRx_doc = - (sum(Rx_usz_doc[t]) * PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz)
            if t == 0:
                MRx_doc_acum.append(MRx_doc)
            else:
                MRx_doc_acum.append(MRx_doc + MRx_doc_acum[-1])

            Min_doc = WFR.tQin[t] * INFLOW_PARAMETERS.cin_doc[t] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Min_doc_acum.append(Min_doc)
            else:
                Min_doc_acum.append(Min_doc + Min_doc_acum[-1])

            Mover_doc = WFR.tQover[t] * cp_doc[t] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Mover_doc_acum.append(Mover_doc)
            else:
                Mover_doc_acum.append(Mover_doc + Mover_doc_acum[-1])

            Mpipe_doc = WFR.tQpipe[t] * c_usz_doc[t][USZ.m_usz - 1] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Mpipe_doc_acum.append(Mpipe_doc)
            else:
                Mpipe_doc_acum.append(Mpipe_doc + Mpipe_doc_acum[-1])

            Minfsz_doc = WFR.tQinfsz[t] * c_usz_doc[t][USZ.m_usz - 1] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Minfsz_doc_acum.append(Minfsz_doc)
            else:
                Minfsz_doc_acum.append(Minfsz_doc + Minfsz_doc_acum[-1])

            Met_doc = WFR.tQet[t] * cl_i1_doc[0] * GENERAL_PARAMETERS.dt * 1000
            if t == 0:
                Met_doc_acum.append(Met_doc)
            else:
                Met_doc_acum.append(Met_doc + Met_doc_acum[-1])

            Mpz_doc = WFR.thpEND[t] * PZ.Ab * 1000 * cp_doc[t]
            Mpz_doc_list.append(Mpz_doc)

            Mstor_doc_mb = Min_doc_acum[-1] - Mpz_doc_list[-1] - Mover_doc_acum[-1] - Mpipe_doc_acum[-1] - \
                           Minfsz_doc_acum[-1] - Met_doc_acum[-1] - Msoil_doc_acum[-1] - MRx_doc_acum[-1]
            Mstor_doc_mb_list.append(Mstor_doc_mb)

            delta_ast_doc_usz = (Mstor_doc_mb - Mstor_doc_ast) / (PZ.Ab * WFR.tteta_usz[t] * WFR.thusz[t] * 1000 / USZ.m_usz)
            cl_i1_doc[1] = cl_i1_doc[1] + delta_ast_doc_usz
            if cl_i1_doc[1] > 0:
                cl_i1_doc[1] = cl_i1_doc[1]
            else:
                cl_i1_doc[1] = 0

        ## adding layers of USZ in time 
        c_usz_o2.append(cl_i1_o2)
        Rx_usz_o2.append(Rxl_o2)

        c_usz_nh4.append(cl_i1_nh4)
        cs_usz_nh4.append(csi_usz_nh4)
        Rx_usz_nh4.append(Rxl_nh4)

        c_usz_no3.append(cl_i1_no3)
        Rx_usz_no3.append(Rxl_no3)

        c_usz_doc.append(cl_i1_doc)
        cs_usz_doc.append(csi_usz_doc)
        Rx_usz_doc.append(Rxl_doc)

        ## adding layers of SZ in time
        if GENERAL_PARAMETERS.hpipe > 0:
            c_sz_o2.append(cj_i1_o2)
            Rx_sz_o2.append(Rxj_o2)

            c_sz_nh4.append(cj_i1_nh4)
            cs_sz_nh4.append(csi_sz_nh4)
            Rx_sz_nh4.append(Rxj_nh4)

            c_sz_no3.append(cj_i1_no3)
            Rx_sz_no3.append(Rxj_no3)

            c_sz_doc.append(cj_i1_doc)
            cs_sz_doc.append(csi_sz_doc)
            Rx_sz_doc.append(Rxj_doc)

    # **5. Transforming in dataframe **
    ## Oxygen
    data_usz_o2 = pd.DataFrame(c_usz_o2)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz' + str(a)
        column_name.append(name)
        a = a + 1
    data_usz_o2.set_axis(column_name, axis='columns', inplace=True)

    # print(data_usz_o2)

    data_sz_o2 = pd.DataFrame(c_sz_o2)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz' + str(a)
        column_name.append(name)
        a = a + 1
    data_sz_o2.set_axis(column_name, axis='columns', inplace=True)

    data_rx_usz_o2 = pd.DataFrame(Rx_usz_o2)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_usz_o2.set_axis(column_name, axis='columns', inplace=True)

    data_rx_sz_o2 = pd.DataFrame(Rx_sz_o2)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_sz_o2.set_axis(column_name, axis='columns', inplace=True)

    frames = [data_usz_o2, data_sz_o2, data_rx_usz_o2, data_rx_sz_o2]
    data_o2 = pd.concat((frames), axis=1)

    cp_o2.append(0)
    data_o2['pz'] = cp_o2

    indice_n = list(range(len(cp_o2)))
    INFLOW_PARAMETERS.cin_o2.append(0)
    c_in = INFLOW_PARAMETERS.cin_o2[:len(indice_n)]
    # print('len c_in:', len(c_in), ', len_indice:', len(WFR.indice))
    data_o2['c_in'] = c_in

    data_o2['t'] = indice_n

    ## Amonia
    data_usz_nh4 = pd.DataFrame(c_usz_nh4)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz' + str(a)
        column_name.append(name)
        a = a + 1
    data_usz_nh4.set_axis(column_name, axis='columns', inplace=True)

    data_sz_nh4 = pd.DataFrame(c_sz_nh4)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz' + str(a)
        column_name.append(name)
        a = a + 1
    data_sz_nh4.set_axis(column_name, axis='columns', inplace=True)

    data_s_usz_nh4 = pd.DataFrame(cs_usz_nh4)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz_soil' + str(a)
        column_name.append(name)
        a = a + 1
    data_s_usz_nh4.set_axis(column_name, axis='columns', inplace=True)

    data_s_sz_nh4 = pd.DataFrame(cs_sz_nh4)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz_soil' + str(a)
        column_name.append(name)
        a = a + 1
    data_s_sz_nh4.set_axis(column_name, axis='columns', inplace=True)

    data_rx_usz_nh4 = pd.DataFrame(Rx_usz_nh4)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_usz_nh4.set_axis(column_name, axis='columns', inplace=True)

    data_rx_sz_nh4 = pd.DataFrame(Rx_sz_nh4)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_sz_nh4.set_axis(column_name, axis='columns', inplace=True)

    frames = [data_usz_nh4, data_sz_nh4, data_s_usz_nh4, data_s_sz_nh4, data_rx_usz_nh4, data_rx_sz_nh4]
    data_nh4 = pd.concat((frames), axis=1)
    cp_nh4.append(0)
    INFLOW_PARAMETERS.cin_nh4.append(0)
    c_in = INFLOW_PARAMETERS.cin_nh4[:len(indice_n)]
    # print('len c_in:', len(c_in), ', len_indice:', len(WFR.indice))
    data_nh4['c_in'] = c_in

    data_nh4['pz'] = cp_nh4
    # WFR.indice.append(len(WFR.indice)+1)
    data_nh4['t'] = indice_n

    ## Nitrate
    data_usz_no3 = pd.DataFrame(c_usz_no3)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz' + str(a)
        column_name.append(name)
        a = a + 1
    data_usz_no3.set_axis(column_name, axis='columns', inplace=True)

    data_sz_no3 = pd.DataFrame(c_sz_no3)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz' + str(a)
        column_name.append(name)
        a = a + 1
    data_sz_no3.set_axis(column_name, axis='columns', inplace=True)

    data_rx_usz_no3 = pd.DataFrame(Rx_usz_no3)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_usz_no3.set_axis(column_name, axis='columns', inplace=True)

    data_rx_sz_no3 = pd.DataFrame(Rx_sz_no3)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_sz_no3.set_axis(column_name, axis='columns', inplace=True)

    frames = [data_usz_no3, data_sz_no3, data_rx_usz_no3, data_rx_sz_no3]
    data_no3 = pd.concat((frames), axis=1)
    cp_no3.append(0)
    data_no3['pz'] = cp_no3

    INFLOW_PARAMETERS.cin_no3.append(0)
    c_in = INFLOW_PARAMETERS.cin_no3[:len(indice_n)]
    data_no3['c_in'] = c_in

    data_no3['t'] = indice_n

    ## DOC
    data_usz_doc = pd.DataFrame(c_usz_doc)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz' + str(a)
        column_name.append(name)
        a = a + 1
    data_usz_doc.set_axis(column_name, axis='columns', inplace=True)

    data_sz_doc = pd.DataFrame(c_sz_doc)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz' + str(a)
        column_name.append(name)
        a = a + 1
    data_sz_doc.set_axis(column_name, axis='columns', inplace=True)

    data_rx_usz_doc = pd.DataFrame(Rx_usz_doc)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_usz_doc.set_axis(column_name, axis='columns', inplace=True)

    data_rx_sz_doc = pd.DataFrame(Rx_sz_doc)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_sz_doc.set_axis(column_name, axis='columns', inplace=True)

    frames = [data_usz_doc, data_sz_doc, data_rx_usz_doc, data_rx_sz_doc]
    data_doc = pd.concat((frames), axis=1)
    cp_doc.append(0)
    data_doc['pz'] = cp_doc

    INFLOW_PARAMETERS.cin_doc.append(0)
    c_in = INFLOW_PARAMETERS.cin_doc[:len(indice_n)]
    data_doc['c_in'] = c_in

    data_doc['t'] = indice_n
    return data_o2, data_nh4, data_no3, data_doc


if __name__ == '__main__':
    inicio = datetime.datetime.now()
    run_W()
    WFR.water_balance(GENERAL_PARAMETERS.dt)
    data_o2, data_nh4, data_no3, data_doc = run_Kin()

    #data_nh4.to_csv('results_Kin_pf_nh4_2.csv', index = False)
    #data_o2.to_csv('results_Kin_pf_o2.csv', index = False)
    #data_no3.to_csv('results_Kin_pf_no3.csv', index = False)
    #data_doc.to_csv('results_Kin_pf_doc.csv', index = False)

    fim = datetime.datetime.now()
    print('Elapsed time: ', fim - inicio)
    print('Done!')
    wf_test = results_tests.water_flow_comparison_test("results/results_00/water_flow_results.csv", WFR)
    nh4_test = results_tests.water_quality_comparison_test("results/results_00/results_Kin_pf_nh4_2.csv", data_nh4)
    o2_test = results_tests.water_quality_comparison_test("results/results_00/results_Kin_pf_o2.csv", data_o2)
    no3_test = results_tests.water_quality_comparison_test("results/results_00/results_Kin_pf_no3.csv", data_no3)
    doc_test = results_tests.water_quality_comparison_test("results/results_00/results_Kin_pf_doc.csv", data_doc)
    print(len(wf_test), len(nh4_test), len(o2_test), len(no3_test))
