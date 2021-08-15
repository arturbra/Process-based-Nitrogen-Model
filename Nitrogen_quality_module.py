import parameters
import datetime
import results_tests

SETUP_FILE = "parameters.ini"
INFLOW_PARAMETERS = parameters.ConcentrationInflow("concentration_inflow.csv")
GENERAL_PARAMETERS = parameters.GeneralParameters(SETUP_FILE)
USZ = parameters.UnsaturatedZone(SETUP_FILE, GENERAL_PARAMETERS.L, GENERAL_PARAMETERS.hpipe, GENERAL_PARAMETERS.dz)
PZ = parameters.PondingZone(SETUP_FILE)
SZ = parameters.SaturatedZone(SETUP_FILE, GENERAL_PARAMETERS.n, USZ.m_usz)
SOIL_PLANT = parameters.SoilPlant(SETUP_FILE, USZ.nusz_ini)
NH4 = parameters.Ammonia(USZ.m_usz, SZ.m_sz, SETUP_FILE)
NO3 = parameters.Nitrate(USZ.m_usz, SZ.m_sz, SETUP_FILE)
O2 = parameters.Oxygen(USZ.m_usz, SZ.m_sz, SETUP_FILE)
DOC = parameters.DissolvedOrganicCarbon(USZ.m_usz, SZ.m_sz, SETUP_FILE)


# **4. Model routine**
def run_Kin():
    for t in range(len(WFR.indice)-1):
        #Ponding zone
        cin_o2 = INFLOW_PARAMETERS.cin_o2[t]
        cin_nh4 = INFLOW_PARAMETERS.cin_nh4[t]
        cin_no3 = INFLOW_PARAMETERS.cin_no3[t]
        cin_doc = INFLOW_PARAMETERS.cin_doc[t]
        
        hp = WFR.thpEND[t]
        if hp < 0.001:
            hp = 0
            
        if t == 0:
            hp_a = 0
            Qorif = 0
        else:
            hp_a = WFR.thpEND[t-1]
            Qorif = WFR.tQpipe[t]

        if t < (len(WFR.indice) - 1):
            teta_sm_iplus1 = WFR.tteta_usz[t + 1]
            teta_b_iplus1 = WFR.tteta_sz[t + 1]
        else:
            teta_sm_iplus1 = WFR.tteta_usz[t]
            teta_b_iplus1 = WFR.tteta_sz[t]        

        O2.Rx_pz = O2.f_reaction_pz()
        NH4.Rx_pz = NH4.f_reaction_pz()
        NO3.Rx_pz = NO3.f_reaction_pz()
        DOC.Rx_pz = DOC.f_reaction_pz()
           
        if hp == 0:
            O2.cpi = 0
            NH4.cpi = 0
            NO3.cpi = 0
            DOC.cpi = 0
        else:        
            O2.cpi = PZ.f_concentration(cin_o2, WFR.tQin[t], O2.cp_a, WFR.tQpf[t], WFR.tQover[t], O2.Rx_pz, hp, hp_a, PZ.Ab, GENERAL_PARAMETERS.dt, 0.00002)
            NH4.cpi = PZ.f_concentration(cin_nh4, WFR.tQin[t], NH4.cp_a, WFR.tQpf[t], WFR.tQover[t], NH4.Rx_pz, hp, hp_a, PZ.Ab, GENERAL_PARAMETERS.dt, 0.00002)
            NO3.cpi = PZ.f_concentration(cin_no3, WFR.tQin[t], NO3.cp_a, WFR.tQpf[t], WFR.tQover[t], NO3.Rx_pz, hp, hp_a, PZ.Ab, GENERAL_PARAMETERS.dt, 0.00002)
            DOC.cpi = PZ.f_concentration(cin_doc, WFR.tQin[t], DOC.cp_a, WFR.tQpf[t], WFR.tQover[t], DOC.Rx_pz, hp, hp_a, PZ.Ab, GENERAL_PARAMETERS.dt, 0.00002)
            
        O2.cp.append(O2.cpi)
        O2.cp_a = O2.cp[-1]

        NH4.cp.append(NH4.cpi)
        NH4.cp_a = NH4.cp[-1]

        NO3.cp.append(NO3.cpi)
        NO3.cp_a = NO3.cp[-1]

        DOC.cp.append(DOC.cpi)
        DOC.cp_a = DOC.cp[-1]
        
        #USZ
        O2.cli = O2.c_usz[t].copy()
        NH4.cli = NH4.c_usz[t].copy()
        NO3.cli = NO3.c_usz[t].copy()
        DOC.cli = DOC.c_usz[t].copy()
        
        #SZ
        if GENERAL_PARAMETERS.hpipe > 0:
            O2.cji = O2.c_sz[t].copy()
            NH4.cji = NH4.c_sz[t].copy()
            NO3.cji = NO3.c_sz[t].copy()
            DOC.cji = DOC.c_sz[t].copy()
            
        if USZ.m_usz != 0:
            O2.cl_i1 = []
            NH4.cl_i1 = []
            NO3.cl_i1 = []
            DOC.cl_i1 = []

            O2.Rxl = []
            NH4.Rxl = []
            NO3.Rxl = []
            DOC.Rxl = []

            NH4.csi_usz = []
            DOC.csi_usz = []

    #Predictive step

    ####   USZ   ####
            for l in range(USZ.m_usz):
                O2.cliplus1 = O2.concentration_water_phase_usz(USZ.m_usz, l)
                NH4.cliplus1 = NH4.concentration_water_phase_usz(USZ.m_usz, l)
                NO3.cliplus1 = NO3.concentration_water_phase_usz(USZ.m_usz, l)
                DOC.cliplus1 = DOC.concentration_water_phase_usz(USZ.m_usz, l)
                #since we have added an initial value of cs = 0 in the NH4.cs_usz list, when calling the index equal = 't' we are actually calling the value corresponding to t-1
                O2.cs = O2.concentration_soil_phase_usz()
                NH4.csi_usz.append(NH4.concentration_soil_phase_usz(l, t, WFR.tteta_usz[t], NH4.cli[l], SOIL_PLANT.ro, GENERAL_PARAMETERS.dt, 0.00000000000001, kmicro=NH4.k_nh4_mb))
                NO3.cs = NO3.concentration_soil_phase_usz()
                DOC.csi_usz.append(DOC.concentration_soil_phase_usz(l, t, WFR.tteta_usz[t], DOC.cli[l], SOIL_PLANT.ro, GENERAL_PARAMETERS.dt, 0.00000000000001, kmicro=DOC.k_doc_mb))
                DOC.cs_usz_a = DOC.cs_usz[t][l]

                USZ.unit_flux = USZ.f_unit_flux(l, WFR.tQpf[t], WFR.tQet1[t], WFR.tQfs[t], WFR.tQhc[t], Qorif, WFR.tQinfsz[t], WFR.tteta_usz[t], PZ.Ab, GENERAL_PARAMETERS.hpipe)

                O2.Rxi_usz = O2.f_reaction_usz(O2.cli[l], NH4.cli[l], GENERAL_PARAMETERS.k_nit) + SOIL_PLANT.f_plant_uptake_usz(O2.cli[l], WFR.tteta_usz[t], "O2")
                NH4.Rxi_usz = NH4.f_reaction_usz(NH4.cli[l], GENERAL_PARAMETERS.k_nit) + SOIL_PLANT.f_plant_uptake_usz(NH4.cli[l], WFR.tteta_usz[t], "NH4", NH4.Fm_nh4, NH4.Km_nh4)
                NO3.Rxi_usz = NO3.f_reaction_usz(NH4.cli[l], GENERAL_PARAMETERS.k_nit) + SOIL_PLANT.f_plant_uptake_usz(NO3.cli[l], WFR.tteta_usz[t], "NO3", NO3.Fm_no3, NO3.Km_no3)
                DOC.Rxi_usz = DOC.f_reaction_usz(DOC.cli[l])

                O2.Rxl.append(O2.Rxi_usz * (1 / teta_sm_iplus1) * GENERAL_PARAMETERS.dt)
                NH4.Rxl.append(NH4.Rxi_usz * (1 / teta_sm_iplus1) * GENERAL_PARAMETERS.dt)
                NO3.Rxl.append(NO3.Rxi_usz * (1 / teta_sm_iplus1) * GENERAL_PARAMETERS.dt)
                DOC.Rxl.append(DOC.Rxi_usz * (1 / teta_sm_iplus1) * GENERAL_PARAMETERS.dt)

    ### Oxygen
                O2.peclet_usz = USZ.f_peclet(USZ.unit_flux, O2.D, GENERAL_PARAMETERS.dz)
                concentration = O2.concentration_delta_usz(O2.peclet_usz, l, USZ.m_usz, GENERAL_PARAMETERS.dz, WFR.tteta_usz[t], teta_sm_iplus1, USZ.unit_flux, GENERAL_PARAMETERS.dt, SOIL_PLANT.ro, SOIL_PLANT.f, 0.0000000000000001)
                O2.cl_i1.append(concentration)

    ### Amonia
                NH4.peclet_usz = USZ.f_peclet(USZ.unit_flux, NH4.D, GENERAL_PARAMETERS.dz)
                concentration = NH4.concentration_delta_usz(NH4.peclet_usz, l, USZ.m_usz, GENERAL_PARAMETERS.dz, WFR.tteta_usz[t], teta_sm_iplus1, USZ.unit_flux, GENERAL_PARAMETERS.dt, SOIL_PLANT.ro, SOIL_PLANT.f, 0.0000000000000001)
                NH4.cl_i1.append(concentration)

    ### Nitrate
                NO3.peclet_usz = USZ.f_peclet(USZ.unit_flux, NO3.D, GENERAL_PARAMETERS.dz)
                concentration = NO3.concentration_delta_usz(NO3.peclet_usz, l, USZ.m_usz, GENERAL_PARAMETERS.dz, WFR.tteta_usz[t], teta_sm_iplus1, USZ.unit_flux, GENERAL_PARAMETERS.dt, SOIL_PLANT.ro, SOIL_PLANT.f, 0.0000000000000001)
                NO3.cl_i1.append(concentration)

    ### DOC
                DOC.peclet_usz = USZ.f_peclet(USZ.unit_flux, DOC.D, GENERAL_PARAMETERS.dz)
                concentration = DOC.concentration_delta_usz(DOC.peclet_usz, l, USZ.m_usz, GENERAL_PARAMETERS.dz, WFR.tteta_usz[t], teta_sm_iplus1, USZ.unit_flux, GENERAL_PARAMETERS.dt, SOIL_PLANT.ro, SOIL_PLANT.f, 0.0000000000000001)
                DOC.cl_i1.append(concentration)

    #####   SZ   #####

        O2.cj_i1 = []
        NH4.cj_i1 = []
        NO3.cj_i1 = []
        DOC.cj_i1 = []
        
        O2.Rxj = []
        NH4.Rxj = []
        NO3.Rxj = []
        DOC.Rxj = []
        
        NH4.csi_sz = []
        DOC.csi_sz = []

        for j in range(SZ.m_sz):
            O2.cjplus1 = O2.concentration_water_phase_sz(SZ.m_sz, j)
            NH4.cjplus1 = NH4.concentration_water_phase_sz(SZ.m_sz, j)
            NO3.cjplus1 = NO3.concentration_water_phase_sz(SZ.m_sz, j)
            DOC.cjplus1 = DOC.concentration_water_phase_sz(SZ.m_sz, j)

            O2.cs = O2.concentration_soil_phase_sz()
            NO3.cs = NO3.concentration_soil_phase_sz()
            NH4.csi_sz.append(NH4.concentration_soil_phase_sz(j, t, WFR.tteta_sz[t], NH4.cji[l], SOIL_PLANT.ro, GENERAL_PARAMETERS.dt, 0.0000000000000001, kmicro=NH4.k_nh4_mb))
            NH4.cs_sz_a = NH4.cs_sz[t][j]
            DOC.csi_sz.append(DOC.concentration_soil_phase(j, t,  WFR.tteta_sz[t], DOC.cji[l], SOIL_PLANT.ro, GENERAL_PARAMETERS.dt, 0.0000000000000001, kmicro=DOC.k_doc_mb))
            DOC.cs_sz_a = DOC.cs_sz[t][j]
            
            SZ.unit_flux = SZ.f_unit_flux(j, WFR.tQfs[t], WFR.tQhc[t], WFR.tQet2[t], Qorif, WFR.tQinfsz[t], WFR.tteta_sz[t], PZ.Ab)
            
            O2.Rxi_sz = O2.f_reaction_sz(O2.cji[j], NH4.cji[j], GENERAL_PARAMETERS.k_nit) + SOIL_PLANT.f_plant_uptake_sz(O2.cji[j], WFR.tteta_sz[t], "O2")
            NH4.Rxi_sz = NH4.f_reaction_sz() + SOIL_PLANT.f_plant_uptake_sz(NH4.cji[j], WFR.tteta_sz[t], "NH4", NH4.Fm_nh4, NH4.Km_nh4)
            NO3.Rxi_sz = NO3.f_reaction_sz(NO3.cji[j], O2.cji[j], DOC.cji[j], GENERAL_PARAMETERS.k_denit) + SOIL_PLANT.f_plant_uptake_sz(NO3.cji[j], WFR.tteta_sz[t], "NO3", NO3.Fm_no3, NO3.Km_no3)
            DOC.Rxi_sz = DOC.f_reaction_sz(DOC.cji[j])

            O2.Rxj.append(O2.Rxi_sz * (1 / teta_b_iplus1) * GENERAL_PARAMETERS.dt)
            NH4.Rxj.append(NH4.Rxi_sz * (1 / teta_b_iplus1) * GENERAL_PARAMETERS.dt)
            NO3.Rxj.append(NO3.Rxi_sz * (1 / teta_b_iplus1) * GENERAL_PARAMETERS.dt)
            DOC.Rxj.append(DOC.Rxi_sz * (1 / teta_b_iplus1) * GENERAL_PARAMETERS.dt)

### Oxygen
            O2.peclet_sz = SZ.f_peclet(SZ.unit_flux, O2.D, GENERAL_PARAMETERS.dz)
            concentration = O2.concentration_delta_sz(USZ.m_usz, SZ.m_sz, GENERAL_PARAMETERS.n, O2.peclet_sz, j, GENERAL_PARAMETERS.dz, WFR.tteta_sz[t], teta_b_iplus1, SZ.unit_flux, GENERAL_PARAMETERS.dt, SOIL_PLANT.ro, SOIL_PLANT.f, 0.0000000000000001)
            O2.cj_i1.append(concentration)

### Amonia
            NH4.peclet_sz = SZ.f_peclet(SZ.unit_flux, NH4.D, GENERAL_PARAMETERS.dz)
            concentration = NH4.concentration_delta_sz(USZ.m_usz, SZ.m_sz, GENERAL_PARAMETERS.n, NH4.peclet_sz, j, GENERAL_PARAMETERS.dz, WFR.tteta_sz[t], teta_b_iplus1, SZ.unit_flux, GENERAL_PARAMETERS.dt, SOIL_PLANT.ro, SOIL_PLANT.f, 0.0000000000000001)
            NH4.cj_i1.append(concentration)

### Nitrate
            NO3.peclet_sz = SZ.f_peclet(SZ.unit_flux, NO3.D, GENERAL_PARAMETERS.dz)
            concentration = NO3.concentration_delta_sz(USZ.m_usz, SZ.m_sz, GENERAL_PARAMETERS.n, NO3.peclet_sz, j, GENERAL_PARAMETERS.dz, WFR.tteta_sz[t], teta_b_iplus1, SZ.unit_flux, GENERAL_PARAMETERS.dt, SOIL_PLANT.ro, SOIL_PLANT.f, 0.0000000000000001)
            NO3.cj_i1.append(concentration)

### DOC
            DOC.peclet_sz = SZ.f_peclet(SZ.unit_flux, DOC.D, GENERAL_PARAMETERS.dz)
            concentration = DOC.concentration_delta_sz(USZ.m_usz, SZ.m_sz, GENERAL_PARAMETERS.n, DOC.peclet_sz, j, GENERAL_PARAMETERS.dz, WFR.tteta_sz[t], teta_b_iplus1, SZ.unit_flux, GENERAL_PARAMETERS.dt, SOIL_PLANT.ro, SOIL_PLANT.f, 0.0000000000000001)
            DOC.cj_i1.append(concentration)

        #Corrective step
            
        ### Oxygen
        O2.mass_storm_asterisk = O2.f_mass_storm_asterisk(t, PZ.Ab, WFR.tteta_usz, WFR.tteta_sz, USZ.m_usz, SZ.m_sz, WFR.thusz, WFR.thsz, GENERAL_PARAMETERS.hpipe)
        O2.mass_soil = O2.f_mass_soil()
        O2.mass_reaction_rate = O2.f_mass_reaction_rate(t, PZ.Ab, WFR.tteta_usz, WFR.tteta_sz, USZ.m_usz, SZ.m_sz, WFR.thusz, WFR.thsz, GENERAL_PARAMETERS.hpipe)
        O2.mass_inflow = O2.f_mass_inflow(t, WFR.tQin, INFLOW_PARAMETERS.cin_o2, GENERAL_PARAMETERS.dt)
        O2.mass_overflow = O2.f_mass_overflow(t, WFR.tQover, GENERAL_PARAMETERS.dt)
        O2.mass_pipe_outflow = O2.f_mass_pipe_outflow(t, WFR.tQpipe, USZ.m_usz, SZ.m_sz, GENERAL_PARAMETERS.dt, GENERAL_PARAMETERS.hpipe)
        O2.mass_infiltration_sz = O2.f_mass_infiltration_sz(t, WFR.tQinfsz, USZ.m_usz, SZ.m_sz, GENERAL_PARAMETERS.dt, GENERAL_PARAMETERS.hpipe)
        O2.mass_evapotranspiration = O2.f_mass_evapotranspiration(t, WFR.tQet, GENERAL_PARAMETERS.dt)
        O2.mass_pz = O2.f_mass_pz(t, WFR.thpEND, PZ.Ab)

        O2.mass_balance.append(O2.f_mass_balance(t, O2.mass_storm_asterisk, O2.mass_soil, O2.mass_reaction_rate,
                                                 O2.mass_inflow, O2.mass_overflow, O2.mass_pipe_outflow,
                                                 O2.mass_infiltration_sz, O2.mass_evapotranspiration, O2.mass_pz))

        O2.cl_i1[1] = O2.f_delta_mass_balance_usz(t, PZ.Ab, WFR.tteta_usz, WFR.thusz, USZ.m_usz)

        ### Amonia
        NH4.mass_storm_asterisk = NH4.f_mass_storm_asterisk(t, PZ.Ab, WFR.tteta_usz, WFR.tteta_sz, USZ.m_usz, SZ.m_sz, WFR.thusz, WFR.thsz, GENERAL_PARAMETERS.hpipe)
        NH4.mass_soil = NH4.f_mass_soil(t, SOIL_PLANT.ro, PZ.Ab, WFR.thusz, WFR.thsz, USZ.m_usz, SZ.m_sz, GENERAL_PARAMETERS.hpipe)
        NH4.mass_reaction_rate = NH4.f_mass_reaction_rate(t, PZ.Ab, WFR.tteta_usz, WFR.tteta_sz, USZ.m_usz, SZ.m_sz, WFR.thusz, WFR.thsz, GENERAL_PARAMETERS.hpipe)
        NH4.mass_inflow = NH4.f_mass_inflow(t, WFR.tQin, INFLOW_PARAMETERS.cin_nh4, GENERAL_PARAMETERS.dt)
        NH4.mass_overflow = NH4.f_mass_overflow(t, WFR.tQover, GENERAL_PARAMETERS.dt)
        NH4.mass_pipe_outflow = NH4.f_mass_pipe_outflow(t, WFR.tQpipe, USZ.m_usz, SZ.m_sz, GENERAL_PARAMETERS.dt, GENERAL_PARAMETERS.hpipe)
        NH4.mass_infiltration_sz = NH4.f_mass_infiltration_sz(t, WFR.tQinfsz, USZ.m_usz, SZ.m_sz, GENERAL_PARAMETERS.dt, GENERAL_PARAMETERS.hpipe)
        NH4.mass_evapotranspiration = NH4.f_mass_evapotranspiration(t, WFR.tQet, GENERAL_PARAMETERS.dt)
        NH4.mass_pz = NH4.f_mass_pz(t, WFR.thpEND, PZ.Ab)

        NH4.mass_balance.append(NH4.f_mass_balance(t, NH4.mass_storm_asterisk, NH4.mass_soil, NH4.mass_reaction_rate,
                                                 NH4.mass_inflow, NH4.mass_overflow, NH4.mass_pipe_outflow,
                                                 NH4.mass_infiltration_sz, NH4.mass_evapotranspiration, NH4.mass_pz))

        NH4.cl_i1[1] = NH4.f_delta_mass_balance_usz(t, PZ.Ab, WFR.tteta_usz, WFR.thusz, USZ.m_usz)


        ### Nitrate
        NO3.mass_storm_asterisk = NO3.f_mass_storm_asterisk(t, PZ.Ab, WFR.tteta_usz, WFR.tteta_sz, USZ.m_usz, SZ.m_sz, WFR.thusz, WFR.thsz, GENERAL_PARAMETERS.hpipe)
        NO3.mass_soil = NO3.f_mass_soil()
        NO3.mass_reaction_rate = NO3.f_mass_reaction_rate(t, PZ.Ab, WFR.tteta_usz, WFR.tteta_sz, USZ.m_usz, SZ.m_sz, WFR.thusz, WFR.thsz, GENERAL_PARAMETERS.hpipe)
        NO3.mass_inflow = NO3.f_mass_inflow(t, WFR.tQin, INFLOW_PARAMETERS.cin_no3, GENERAL_PARAMETERS.dt)
        NO3.mass_overflow = NO3.f_mass_overflow(t, WFR.tQover, GENERAL_PARAMETERS.dt)
        NO3.mass_pipe_outflow = NO3.f_mass_pipe_outflow(t, WFR.tQpipe, USZ.m_usz, SZ.m_sz, GENERAL_PARAMETERS.dt, GENERAL_PARAMETERS.hpipe)
        NO3.mass_infiltration_sz = NO3.f_mass_infiltration_sz(t, WFR.tQinfsz, USZ.m_usz, SZ.m_sz, GENERAL_PARAMETERS.dt, GENERAL_PARAMETERS.hpipe)
        NO3.mass_evapotranspiration = NO3.f_mass_evapotranspiration(t, WFR.tQet, GENERAL_PARAMETERS.dt)
        NO3.mass_pz = NO3.f_mass_pz(t, WFR.thpEND, PZ.Ab)

        NO3.mass_balance.append(NO3.f_mass_balance(t, NO3.mass_storm_asterisk, NO3.mass_soil, NO3.mass_reaction_rate,
                                                 NO3.mass_inflow, NO3.mass_overflow, NO3.mass_pipe_outflow,
                                                 NO3.mass_infiltration_sz, NO3.mass_evapotranspiration, NO3.mass_pz))

        NO3.cl_i1[1] = NO3.f_delta_mass_balance_usz(t, PZ.Ab, WFR.tteta_usz, WFR.thusz, USZ.m_usz)


        ### DOC

        DOC.mass_storm_asterisk = DOC.f_mass_storm_asterisk(t, PZ.Ab, WFR.tteta_usz, WFR.tteta_sz, USZ.m_usz, SZ.m_sz, WFR.thusz, WFR.thsz, GENERAL_PARAMETERS.hpipe)
        DOC.mass_soil = DOC.f_mass_soil(t, SOIL_PLANT.ro, PZ.Ab, WFR.thusz, WFR.thsz, USZ.m_usz, SZ.m_sz, GENERAL_PARAMETERS.hpipe)
        DOC.mass_reaction_rate = DOC.f_mass_reaction_rate(t, PZ.Ab, WFR.tteta_usz, WFR.tteta_sz, USZ.m_usz, SZ.m_sz, WFR.thusz, WFR.thsz, GENERAL_PARAMETERS.hpipe)
        DOC.mass_inflow = DOC.f_mass_inflow(t, WFR.tQin, INFLOW_PARAMETERS.cin_doc, GENERAL_PARAMETERS.dt)
        DOC.mass_overflow = DOC.f_mass_overflow(t, WFR.tQover, GENERAL_PARAMETERS.dt)
        DOC.mass_pipe_outflow = DOC.f_mass_pipe_outflow(t, WFR.tQpipe, USZ.m_usz, SZ.m_sz, GENERAL_PARAMETERS.dt, GENERAL_PARAMETERS.hpipe)
        DOC.mass_infiltration_sz = DOC.f_mass_infiltration_sz(t, WFR.tQinfsz, USZ.m_usz, SZ.m_sz, GENERAL_PARAMETERS.dt, GENERAL_PARAMETERS.hpipe)
        DOC.mass_evapotranspiration = DOC.f_mass_evapotranspiration(t, WFR.tQet, GENERAL_PARAMETERS.dt)
        DOC.mass_pz = DOC.f_mass_pz(t, WFR.thpEND, PZ.Ab)

        DOC.mass_balance.append(DOC.f_mass_balance(t, DOC.mass_storm_asterisk, DOC.mass_soil, DOC.mass_reaction_rate,
                                                 DOC.mass_inflow, DOC.mass_overflow, DOC.mass_pipe_outflow,
                                                 DOC.mass_infiltration_sz, DOC.mass_evapotranspiration, DOC.mass_pz))

        DOC.cl_i1[1] = DOC.f_delta_mass_balance_usz(t, PZ.Ab, WFR.tteta_usz, WFR.thusz, USZ.m_usz)


    ## adding layers of USZ in time
        O2.c_usz.append(O2.cl_i1)
        O2.Rx_usz.append(O2.Rxl)

        NH4.c_usz.append(NH4.cl_i1 )
        NH4.cs_usz.append(NH4.csi_usz )
        NH4.Rx_usz.append(NH4.Rxl)
        
        NO3.c_usz.append(NO3.cl_i1)
        NO3.Rx_usz.append(NO3.Rxl)
        
        DOC.c_usz.append(DOC.cl_i1)
        DOC.cs_usz.append(DOC.csi_usz)
        DOC.Rx_usz.append(DOC.Rxl)
        
    ## adding layers of SZ in time
        if GENERAL_PARAMETERS.hpipe > 0:
            O2.c_sz.append(O2.cj_i1)
            O2.Rx_sz.append(O2.Rxj)
            
            NH4.c_sz.append(NH4.cj_i1)
            NH4.cs_sz.append(NH4.csi_sz)
            NH4.Rx_sz.append(NH4.Rxj)
            
            NO3.c_sz.append(NO3.cj_i1)
            NO3.Rx_sz.append(NO3.Rxj)
            
            DOC.c_sz.append(DOC.cj_i1)
            DOC.cs_sz.append(DOC.csi_sz)
            DOC.Rx_sz.append(DOC.Rxj)
            
    # **5. Transforming in dataframe **
    ## Oxygen
    data_o2 = O2.water_quality_results(USZ.m_usz, SZ.m_sz, INFLOW_PARAMETERS.cin_o2)
    data_nh4 = NH4.water_quality_results(USZ.m_usz, SZ.m_sz, INFLOW_PARAMETERS.cin_nh4)
    data_no3 = NO3.water_quality_results(USZ.m_usz, SZ.m_sz, INFLOW_PARAMETERS.cin_no3)
    data_doc = DOC.water_quality_results(USZ.m_usz, SZ.m_sz, INFLOW_PARAMETERS.cin_doc)

    return data_o2, data_nh4, data_no3, data_doc
    
if __name__ == '__main__':
    inicio = datetime.datetime.now()

    data_o2, data_nh4, data_no3, data_doc = run_Kin()

    #data_nh4.to_csv('results_Kin_pf_nh4_2.csv', index = False)
    #data_o2.to_csv('results_Kin_pf_o2.csv', index = False)
    #data_no3.to_csv('results_Kin_pf_no3.csv', index = False)
    #data_doc.to_csv('results_Kin_pf_doc.csv', index = False)
    
    fim = datetime.datetime.now()
    print ('Elapsed time: ', fim - inicio)
    print('Done!')
    wf_test = results_tests.water_flow_comparison_test("results/water_flow_results.csv", WFR)
    nh4_test = results_tests.water_quality_comparison_test("results/results_Kin_pf_nh4_2.csv", data_nh4)
    o2_test = results_tests.water_quality_comparison_test("results/results_Kin_pf_o2.csv", data_o2)
    no3_test = results_tests.water_quality_comparison_test("results/results_Kin_pf_no3.csv", data_no3)
    doc_test = results_tests.water_quality_comparison_test("results/results_Kin_pf_doc.csv", data_doc)
    print(len(wf_test), len(nh4_test), len(o2_test), len(no3_test), len(doc_test))