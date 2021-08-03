import parameters
import datetime
import pandas as pd

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
        else:
            hp = hp
            
        if t == 0:
            hp_a = 0
            Qorif = 0
        else:
            hp_a = WFR.thpEND[t-1]
            Qorif = WFR.tQpipe[t]


        if t < (len(WFR.indice)-1):
            teta_sm_iplus1 = WFR.tteta_usz[t+1]
            teta_b_iplus1 = WFR.tteta_sz[t+1]
        else:
            teta_sm_iplus1 = WFR.tteta_usz[t]
            teta_b_iplus1 = WFR.tteta_sz[t]        

        
        Rxi_p_o2 = O2.f_reaction_pz()
        Rxi_p_nh4 = NH4.f_reaction_pz()
        Rxi_p_no3 = NO3.f_reaction_pz()
        Rxi_p_doc = DOC.f_reaction_pz()
        
        #if t < 20:
            #print('t: ', t, 'Rx_p_o2: ', Rxi_p_o2, 'Rx_p_nh4: ', Rxi_p_nh4, 'Rx_p_no3: ', Rxi_p_no3)   
           
        if hp == 0:
            cpi_o2 = 0
            cpi_nh4 = 0
            cpi_no3 = 0
            cpi_doc = 0
       
        else:        
            cpi_o2 = PZ.f_concentration(cin_o2, WFR.tQin[t], O2.cp_a, WFR.tQpf[t], WFR.tQover[t], Rxi_p_o2, hp, hp_a, PZ.Ab, GENERAL_PARAMETERS.dt)
            cpi_nh4 = PZ.f_concentration(cin_nh4, WFR.tQin[t], NH4.cp_a, WFR.tQpf[t], WFR.tQover[t], Rxi_p_nh4, hp, hp_a, PZ.Ab, GENERAL_PARAMETERS.dt)
            cpi_no3 = PZ.f_concentration(cin_no3, WFR.tQin[t], NO3.cp_a, WFR.tQpf[t], WFR.tQover[t], Rxi_p_no3, hp, hp_a, PZ.Ab, GENERAL_PARAMETERS.dt)
            cpi_doc = PZ.f_concentration(cin_doc, WFR.tQin[t], DOC.cp_a, WFR.tQpf[t], WFR.tQover[t], Rxi_p_doc, hp, hp_a, PZ.Ab, GENERAL_PARAMETERS.dt)
#             print('cin_nh4, WFR.tQin[t], NH4.cp_a, WFR.tQpf[t], WFR.tQover[t], Rxi_p_nh4, hp, hp_a: ', cin_nh4, WFR.tQin[t], NH4.cp_a, WFR.tQpf[t], WFR.tQover[t], Rxi_p_nh4, hp, hp_a)
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
            
        O2.cp.append(cpi_o2)
        O2.cp_a = O2.cp[-1]
    
        NH4.cp.append(cpi_nh4)
        NH4.cp_a = NH4.cp[-1]
    
        NO3.cp.append(cpi_no3)
        NO3.cp_a = NO3.cp[-1]
    
        DOC.cp.append(cpi_doc)
        DOC.cp_a = DOC.cp[-1]
        
        #USZ
        cli_o2_list = O2.c_usz[t].copy()
        #print('1_o2', cli_o2_list)
        cli_nh4_list = NH4.c_usz[t].copy()
        #print('1_nh4', cli_nh4_list)
        cli_no3_list = NO3.c_usz[t].copy()
        cli_doc_list = DOC.c_usz[t].copy()
        
        #SZ
        if hpipe > 0:
            cji_o2_list = O2.c_sz[t].copy() #copiar lista
            #cji_o2_list = O2.c_sz[t][:] #outro jeito de copiar a lista
            cji_nh4_list = NH4.c_sz[t].copy()
            cji_no3_list = NO3.c_sz[t].copy()
            cji_doc_list = DOC.c_sz[t].copy()
            
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
                    
    #Predictive step
    
    ####   USZ   ####
            for l in range(USZ.m_usz): 
                #print('l', l)
                
                cl_o2 = cli_o2_list[l] 
                clminus1_o2 = cli_o2_list[l-1]
                if l < (USZ.m_usz - 1):  
                    clplus1_o2 = cli_o2_list[l+1]
                else:
                    clplus1_o2 = 0           
    
                cl_nh4 = cli_nh4_list[l]
                clminus1_nh4 = cli_nh4_list[l-1]
                if l < (USZ.m_usz - 1):
                    clplus1_nh4 = cli_nh4_list[l+1]
                else:
                    clplus1_nh4 = 0
    
                cl_no3 = cli_no3_list[l]
                clminus1_no3 = cli_no3_list[l-1]
                if l < (USZ.m_usz - 1):
                    clplus1_no3 = cli_no3_list[l+1]
                else:
                    clplus1_no3 = 0
                    
                len(cli_doc_list)
                cl_doc = cli_doc_list[l]
                clminus1_doc = cli_doc_list[l-1]

                if l < (USZ.m_usz - 1):
                    clplus1_doc = cli_doc_list[l+1]
                else:
                    clplus1_doc = 0 
    
                            
                cs_o2 = 0
                
                #since we have added an initial value of cs = 0 in the NH4.cs_usz list, when calling the index equal = 't' we are actually calling the value corresponding to t-1
                cs_nh4 = NH4.cs_usz[t][l]
                cs_nh4_iplus1 = NH4.f_concentration_soil(cs_nh4, WFR.tteta_usz[t], NH4.kads_nh4, cl_nh4, NH4.kdes_nh4, SOIL_PLANT.ro, GENERAL_PARAMETERS.dt, kmicro=NH4.k_nh4_mb)
                if cs_nh4_iplus1 < 0.00000000000001:
                    cs_nh4_iplus1 = 0
                    
#                 sor = (WFR.tteta_usz[t]/SOIL_PLANT.ro)*NH4.kads_nh4*cl_nh4*dt
#                 des =  NH4.kdes_nh4*NH4.cs_usz_a*dt
                
#                 print('t: ', t, ', l: ', l , ', cs_nh4_a: ', NH4.cs_usz_a, 'WFR.tteta_usz[t]: ', WFR.tteta_usz[t], 'NH4.kads_nh4: ', NH4.kads_nh4, 'cl_nh4: ', cl_nh4, 'NH4.kdes_nh4: ', NH4.kdes_nh4, 'NH4.k_nh4_mb: ', NH4.k_nh4_mb)
#                 print('cs_nh4: ', cs_nh4)
#                 input()
                
#                 logger.debug(f' t: {t} \t l: {l} \t sor: {sor} \t des: {des} \t cs_nh4: {cs_nh4}')
#                 logger.debug('----SOIL----')
#                 logger.debug('cs_nh4_a: ', NH4.cs_usz_a)
#                 logger.debug('cs_nh4_a: ', NH4.cs_usz_a, 'WFR.tteta_usz[t]: ', WFR.tteta_usz[t], 'NH4.kads_nh4: ', NH4.kads_nh4, 'cl_nh4: ', cl_nh4, 'NH4.kdes_nh4: ', NH4.kdes_nh4, 'NH4.k_nh4_mb: ', NH4.k_nh4_mb)
#                 logger.debug('cs_nh4: ', cs_nh4)
#                 logger.debug('----WATER----')
                csi_usz_nh4.append(cs_nh4_iplus1)

                
                cs_no3 = 0
    
                DOC.cs_usz_a = DOC.cs_usz[t][l]
                cs_doc = DOC.f_concentration_soil(DOC.cs_usz_a, WFR.tteta_usz[t], DOC.kads_doc, cl_doc, DOC.kdes_doc, SOIL_PLANT.ro, GENERAL_PARAMETERS.dt, kmicro=DOC.k_doc_mb)
                csi_usz_doc.append(cs_doc)
                if cs_doc < 0.00000000000001:
                    cs_doc = 0                
                                        
                UF_usz = []
                #WFR.tQet1[t] = 0                        
                alfa, beta = USZ.f_alfa_beta(l)
                UFi_usz = USZ.f_unit_flux(alfa, beta, WFR.tQpf[t], WFR.tQet1[t], WFR.tQfs[t], WFR.tQhc[t], Qorif, WFR.tQinfsz[t], WFR.tteta_usz[t], hpipe, PZ.Ab)
#                 print(t, '  ', l)
#                 print('UFi_usz: ', UFi_usz, ', alfa: ', alfa, ', beta: ', beta,  ', WFR.tQpf[t]: ', WFR.tQpf[t], ', WFR.tQet1[t]:', WFR.tQet1[t], ', WFR.tQfs[t]: ', WFR.tQfs[t], ', WFR.tQhc[t]: ', WFR.tQhc[t], ', Qorif: ', Qorif, ', WFR.tQinfsz[t]: ', WFR.tQinfsz[t])            
#                 input()
                
                UF_usz.append(UFi_usz)
                
                Rxi_2_o2 = O2.f_reaction_usz(cl_o2, cl_nh4, GENERAL_PARAMETERS.k_nit) + SOIL_PLANT.f_plant_uptake_usz(cl_o2, WFR.tteta_usz[t], "O2")
                Rxi_2_nh4 = NH4.f_reaction_usz(cl_nh4, GENERAL_PARAMETERS.k_nit) + SOIL_PLANT.f_plant_uptake_usz(cl_nh4, WFR.tteta_usz[t], "NH4", NH4.Fm_nh4, NH4.Km_nh4)
                Rxi_2_no3 = NO3.f_reaction_usz(cl_nh4, GENERAL_PARAMETERS.k_nit) + SOIL_PLANT.f_plant_uptake_usz(cl_no3, WFR.tteta_usz[t], "NO3", NO3.Fm_no3, NO3.Km_no3)
                Rxi_2_doc = DOC.f_reaction_usz(cl_doc)
                
                Rxl_o2.append(Rxi_2_o2*(1/teta_sm_iplus1)*dt)
                Rxl_nh4.append(Rxi_2_nh4*(1/teta_sm_iplus1)*dt)
                Rxl_no3.append(Rxi_2_no3*(1/teta_sm_iplus1)*dt)
                Rxl_doc.append(Rxi_2_doc*(1/teta_sm_iplus1)*dt)

#                 Rxl_o2.append(Rxi_2_o2*(dt/WFR.tteta_usz[t]))
#                 Rxl_nh4.append(Rxi_2_nh4*(dt/WFR.tteta_usz[t]))
#                 Rxl_no3.append(Rxi_2_no3*(dt/WFR.tteta_usz[t]))
#                 Rxl_doc.append(Rxi_2_doc*(dt/WFR.tteta_usz[t]))
#                 print('t: ', t, ', l: ', l, ', Rx: ', Rxi_2_nh4*(dt/WFR.tteta_usz[t]))
#                 input()
                #if t < 20:
                    #print('t: ', t, 'l: ', l, 'O2.Rx_usz: ', Rxi_2_o2, 'NH4.Rx_usz: ', Rxi_2_nh4, 'NO3.Rx_usz: ', Rxi_2_no3)
                
    
    ### Oxygen
                Peusz_o2 = USZ.f_peclet(UFi_usz, O2.D_o2, GENERAL_PARAMETERS.dz)
                #print(Peusz_o2)
                
                if Peusz_o2 <= 2:
                    if l == 0: #first cell
                        dc_o2 = clplus1_o2 - 2*cl_o2 + cpi_o2
                        dc_dz_o2 = (clplus1_o2 - cpi_o2)/(2*GENERAL_PARAMETERS.dz)
   
        
                    elif l == (USZ.m_usz - 1): #last cell
                        dc_o2 = cl_o2 - 2*cl_o2 + clminus1_o2
                        dc_dz_o2 = (cl_o2 - clminus1_o2)/GENERAL_PARAMETERS.dz
                        
                    else:                 
                        dc_o2 = clplus1_o2 - 2*cl_o2 + clminus1_o2
                        dc_dz_o2 = (clplus1_o2 - clminus1_o2)/(2*GENERAL_PARAMETERS.dz)               
                
                else: #Peusz > 2
                    if l == 0: #first cell
                        dc_o2 = clplus1_o2 - 2*cl_o2 + cpi_o2
                        dc_dz_o2 = (cl_o2 - cpi_o2)/GENERAL_PARAMETERS.dz
        
                    elif l == (USZ.m_usz - 1): #last cell
                        dc_o2 = cl_o2 - 2*cl_o2 + clminus1_o2
                        dc_dz_o2 = (cl_o2 - clminus1_o2)/GENERAL_PARAMETERS.dz
                        
                    else:
                        dc_o2 = clplus1_o2 - 2*cl_o2 + clminus1_o2
                        dc_dz_o2 = (cl_o2 - clminus1_o2)/GENERAL_PARAMETERS.dz            

                delta_c_o2 = O2.f_transport(WFR.tteta_usz[t], teta_sm_iplus1, cl_o2, cs_o2, dc_o2, dc_dz_o2, 0, 0, O2.D_o2, UFi_usz, Rxi_2_o2, GENERAL_PARAMETERS.dt, SOIL_PLANT.ro, SOIL_PLANT.f, GENERAL_PARAMETERS.dz)
                if teta_sm_iplus1 > 0:
                    ci1_o2 = cl_o2 + delta_c_o2
                else:
                    ci1_o2 = 0
                
                if ci1_o2 <= 0.0000000000000001:
                    ci1_o2 = 0
                else:
                    ci1_o2 = ci1_o2          
                     
                #print('ci1_o2', ci1_o2)
                cl_i1_o2.append(ci1_o2)
                #print('2_o2', cl_i1_o2, l) 
                 
    ### Amonia            
                Peusz_nh4 = USZ.f_peclet(UFi_usz, NH4.D_nh4, GENERAL_PARAMETERS.dz)
                #print('t:', t, 'l:', l, 'Peusz_nh4:', Peusz_nh4)
                
                if Peusz_nh4 <= 2:
                    if l == 0: #first cell
                        dc_nh4 = clplus1_nh4 - 2*cl_nh4 + cpi_nh4
                        dc_dz_nh4 = (clplus1_nh4 - cpi_nh4)/(2*GENERAL_PARAMETERS.dz)
        
                    elif l == (USZ.m_usz - 1): #last cell
                        dc_nh4 = cl_nh4 - 2*cl_nh4 + clminus1_nh4
                        dc_dz_nh4 = (cl_nh4 - clminus1_nh4)/GENERAL_PARAMETERS.dz
                        
                    else:
                        dc_nh4 = clplus1_nh4 - 2*cl_nh4 + clminus1_nh4
                        dc_dz_nh4 = (clplus1_nh4 - clminus1_nh4)/(2*GENERAL_PARAMETERS.dz)               
                
                else: #Peusz > 2
                    if l == 0: #first cell
                        dc_nh4 = clplus1_nh4 - 2*cl_nh4 + cpi_nh4
                        dc_dz_nh4 = (cl_nh4 - cpi_nh4)/GENERAL_PARAMETERS.dz
        
                    elif l == (USZ.m_usz - 1): #last cell
                        dc_nh4 = cl_nh4 - 2*cl_nh4 + clminus1_nh4
                        dc_dz_nh4 = (cl_nh4 - clminus1_nh4)/GENERAL_PARAMETERS.dz
                        
                    else:
                        dc_nh4 = clplus1_nh4 - 2*cl_nh4 + clminus1_nh4
                        dc_dz_nh4 = (cl_nh4 - clminus1_nh4)/GENERAL_PARAMETERS.dz            
                
                delta_c_nh4 = NH4.f_transport(WFR.tteta_usz[t], teta_sm_iplus1, cl_nh4, cs_nh4, dc_nh4, dc_dz_nh4, NH4.kads_nh4, NH4.kdes_nh4, NH4.D_nh4, UFi_usz, Rxi_2_nh4, GENERAL_PARAMETERS.dt, SOIL_PLANT.ro, SOIL_PLANT.f, GENERAL_PARAMETERS.dz)

#                 print('t: ', t, ', l: ', l)
#                 print('clplus1_nh4: ',  clplus1_nh4, ', cl_nh4: ' ,cl_nh4 , ', cpi_nh4:' , cpi_nh4, ', clminus1_nh4: ' , clminus1_nh4)
#                 print('dt:', dt, 'WFR.tteta_usz[t]:', WFR.tteta_usz[t], 'teta_sm_iplus1:', teta_sm_iplus1, 'cl_nh4:', cl_nh4, 'cs_nh4:', cs_nh4, 'dc_nh4:', dc_nh4, 'dc_dz_nh4:', dc_dz_nh4, 'NH4.kads_nh4:', NH4.kads_nh4, 'NH4.kdes_nh4:', NH4.kdes_nh4, 'NH4.D_nh4:', NH4.D_nh4, 'UFi_usz:', UFi_usz, 'Rxi_2_nh4:', Rxi_2_nh4)
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
                #print('ci1_nh4: ', ci1_nh4)
                cl_i1_nh4.append(ci1_nh4)
                #print('2_nh4', cl_i1_nh4)
                            
    ### Nitrate            
                Peusz_no3 = USZ.f_peclet(UFi_usz, NO3.D_no3, GENERAL_PARAMETERS.dz)
                
                if Peusz_no3 <= 2:
                    if l == 0: #first cell
                        dc_no3 = clplus1_no3 - 2*cl_no3 + cpi_no3
                        dc_dz_no3 = (clplus1_no3 - cpi_no3)/(2*GENERAL_PARAMETERS.dz)
        
                    elif l == (USZ.m_usz - 1): #last cell
                        dc_no3 = cl_no3 - 2*cl_no3 + clminus1_no3
                        dc_dz_no3 = (cl_no3 - clminus1_no3)/GENERAL_PARAMETERS.dz
                        
                    else:
                        dc_no3 = clplus1_no3 - 2*cl_no3 + clminus1_no3
                        dc_dz_no3 = (clplus1_no3 - clminus1_no3)/(2*GENERAL_PARAMETERS.dz)               
                
                else: #Peusz > 2
                    if l == 0: #first cell
                        dc_no3 = clplus1_no3 - 2*cl_no3 + cpi_no3
                        dc_dz_no3 = (cl_no3 - cpi_no3)/GENERAL_PARAMETERS.dz
        
                    elif l == (USZ.m_usz - 1): #last cell
                        dc_no3 = cl_no3 - 2*cl_no3 + clminus1_no3
                        dc_dz_no3 = (cl_no3 - clminus1_no3)/GENERAL_PARAMETERS.dz
                        
                    else:
                        dc_no3 = clplus1_no3 - 2*cl_no3 + clminus1_no3
                        dc_dz_no3 = (cl_no3 - clminus1_no3)/GENERAL_PARAMETERS.dz            
                
                delta_c_no3 = NO3.f_transport(WFR.tteta_usz[t], teta_sm_iplus1, cl_no3, cs_no3, dc_no3, dc_dz_no3, 0, 0, NO3.D_no3, UFi_usz, Rxi_2_no3, GENERAL_PARAMETERS.dt, SOIL_PLANT.ro, SOIL_PLANT.f, GENERAL_PARAMETERS.dz)
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
                Peusz_doc = USZ.f_peclet(UFi_usz, DOC.D_doc, GENERAL_PARAMETERS.dz)
                
                if Peusz_doc <= 2:
                    if l == 0: #first cell
                        dc_doc = clplus1_doc - 2*cl_doc + cpi_doc
                        dc_dz_doc = (clplus1_doc - cpi_doc)/(2*GENERAL_PARAMETERS.dz)
        
                    elif l == (USZ.m_usz - 1): #last cell
                        dc_doc = cl_doc - 2*cl_doc + clminus1_doc
                        dc_dz_doc = (cl_doc - clminus1_doc)/GENERAL_PARAMETERS.dz
                        
                    else:
                        dc_doc = clplus1_doc - 2*cl_doc + clminus1_doc
                        dc_dz_doc = (clplus1_doc - clminus1_doc)/(2*GENERAL_PARAMETERS.dz)               
                
                else: #Peusz > 2
                    if l == 0: #first cell
                        dc_doc = clplus1_doc - 2*cl_doc + cpi_doc
                        dc_dz_doc = (cl_doc - cpi_doc)/GENERAL_PARAMETERS.dz
        
                    elif l == (USZ.m_usz - 1): #last cell
                        dc_doc = cl_doc - 2*cl_doc + clminus1_doc
                        dc_dz_doc = (cl_doc - clminus1_doc)/GENERAL_PARAMETERS.dz
                        
                    else:
                        dc_doc = clplus1_doc - 2*cl_doc + clminus1_doc
                        dc_dz_doc = (cl_doc - clminus1_doc)/GENERAL_PARAMETERS.dz            
                
                delta_c_doc = DOC.f_transport(WFR.tteta_usz[t], teta_sm_iplus1, cl_doc, cs_doc, dc_doc, dc_dz_doc, DOC.kads_doc, DOC.kdes_doc, DOC.D_doc, UFi_usz, Rxi_2_doc, GENERAL_PARAMETERS.dt, SOIL_PLANT.ro, SOIL_PLANT.f, GENERAL_PARAMETERS.dz)
                if teta_sm_iplus1 > 0:
                    ci1_doc = cl_doc + delta_c_doc
                else:
                    ci1_doc = 0
                    
                if ci1_doc <= 0.0000000000000001:
                    ci1_doc = 0
                else:
                    ci1_doc = ci1_doc    
                
                cl_i1_doc.append(ci1_doc)
                #print('2_nh4', cl_i1_nh4)
                        
    
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
                cjminus1_o2 = cji_o2_list[j-1]
                if j < (SZ.m_sz - 1):
                    cjplus1_o2 = cji_o2_list[j+1]
                else:
                    cjplus1_o2 = 0
                cmminus1_o2 = cli_o2_list[USZ.m_usz-1]
                
                cj_nh4 = cji_nh4_list[j]
                cjminus1_nh4 = cji_nh4_list[j-1]
                if j < (SZ.m_sz - 1):
                    cjplus1_nh4 = cji_nh4_list[j+1]
                else:
                    cjplus1_nh4 = 0
                cmminus1_nh4 = cli_nh4_list[USZ.m_usz-1]
                
                cj_no3 = cji_no3_list[j]
                cjminus1_no3 = cji_no3_list[j-1]
                if j < (SZ.m_sz - 1):    
                    cjplus1_no3 = cji_no3_list[j+1]
                else:
                    cjplus1_no3 = 0
                cmminus1_no3 = cli_no3_list[USZ.m_usz-1] 
    
                cj_doc = cji_doc_list[j]
                cjminus1_doc = cji_doc_list[j-1]
                if j < (SZ.m_sz - 1):    
                    cjplus1_doc = cji_doc_list[j+1]
                else:
                    cjplus1_doc = 0
                cmminus1_doc = cli_doc_list[USZ.m_usz-1]            
                
                cs_o2 = 0
                
                NH4.cs_sz_a = NH4.cs_sz[t][j]
                cs_nh4 = NH4.f_concentration_soil(NH4.cs_sz_a, WFR.tteta_sz[t], kads2_nh4, cj_nh4, kdes2_nh4, SOIL_PLANT.ro, GENERAL_PARAMETERS.dt, kmicro=NH4.k_nh4_mb)
                csi_sz_nh4.append(cs_nh4)
                if cs_nh4 < 0.00000000000001:
                    cs_nh4 = 0
                
#                 print('====>> t: ', t, ', j: ', j)
#                 print('---- SOIL ----')
#                 print('cs_nh4_a: ', NH4.cs_sz_a, 'WFR.tteta_sz[t]: ', WFR.tteta_sz[t], 'NH4.kads_nh4: ', NH4.kads_nh4, 'cj_nh4: ', cj_nh4, 'NH4.kdes_nh4: ', NH4.kdes_nh4, 'NH4.k_nh4_mb: ', NH4.k_nh4_mb)
#                 print('cs_nh4: ', cs_nh4)
#                 print('---- WATER ----')
#                 print('cmminus1_nh4: ', cmminus1_nh4)                
                cs_no3 = 0
    
                DOC.cs_sz_a = DOC.cs_sz[t][j]
                cs_doc = DOC.f_concentration_soil(DOC.cs_sz_a, WFR.tteta_sz[t], kads2_doc, cj_doc, kdes2_doc, kmicro=DOC.k_doc_mb)
                csi_sz_doc.append(cs_doc)
                                            
                UF_sz = []                                 
                alfa, beta = SZ.f_alfa_beta(j)
                #WFR.tQet2[t] = 0
                UFi_sz = SZ.f_unit_flux(alfa, beta, WFR.tQfs[t], WFR.tQhc[t], WFR.tQet2[t], Qorif, WFR.tQinfsz[t], WFR.tteta_sz[t], PZ.Ab)
#                 print('WFR.tQfs[t], WFR.tQhc[t], WFR.tQet2[t], Qorif, WFR.tQinfsz[t]:', WFR.tQfs[t], WFR.tQhc[t], WFR.tQet2[t], Qorif, WFR.tQinfsz[t])
#                 print('UF: ', UFi_sz, ', alfa2: ', alfa2, ', beta2: ' ,beta2)
                            
                UF_sz.append(UFi_sz)
                
                Rxi_3_o2 = O2.f_reaction_sz(cj_o2, cj_nh4, GENERAL_PARAMETERS.k_nit) + SOIL_PLANT.f_plant_uptake_sz(cj_o2, WFR.tteta_sz[t], "O2")
                Rxi_3_nh4 = NH4.f_reaction_sz() + SOIL_PLANT.f_plant_uptake_sz(cj_nh4, WFR.tteta_sz[t], "NH4", NH4.Fm_nh4, NH4.Km_nh4)
                Rxi_3_no3 = NO3.f_reaction_sz(cj_no3, cj_o2, cj_doc, GENERAL_PARAMETERS.k_denit) + SOIL_PLANT.f_plant_uptake_sz(cj_no3, WFR.tteta_sz[t], "NO3", NO3.Fm_no3, NO3.Km_no3)
                Rxi_3_doc = DOC.f_reaction_sz(cj_doc)
                
                Rxj_o2.append(Rxi_3_o2*(1/teta_b_iplus1)*dt)
                Rxj_nh4.append(Rxi_3_nh4*(1/teta_b_iplus1)*dt)
                Rxj_no3.append(Rxi_3_no3*(1/teta_b_iplus1)*dt)
                Rxj_doc.append(Rxi_3_doc*(1/teta_b_iplus1)*dt)
                #if t < 20:
                    #print('t: ', t, 'j: ', j, 'O2.Rx_sz: ', Rxi_3_o2, 'NH4.Rx_sz: ', Rxi_3_nh4, 'NO3.Rx_sz: ', Rxi_3_no3)
                
    ### Oxygen                       
                Pesz_o2 = SZ.f_peclet(UFi_sz, O2.D_o2, GENERAL_PARAMETERS.dz)
                
                if USZ.m_usz < (n-1):
                    if Pesz_o2 <= 2:
                        if j == 0: #first cell
                            dc_o2 = cjplus1_o2 - 2*cj_o2 + cmminus1_o2
                            dc_dz_o2 = (cjplus1_o2 - cmminus1_o2)/(2*GENERAL_PARAMETERS.dz)
            
                        elif j == (SZ.m_sz - 1): #last cell
                            dc_o2 = cj_o2 - 2*cj_o2 + cjminus1_o2
                            dc_dz_o2 = (cj_o2 - cjminus1_o2)/GENERAL_PARAMETERS.dz
                            
                        else:
                            dc_o2 = cjplus1_o2 - 2*cj_o2 + cjminus1_o2
                            dc_dz_o2 = (cjplus1_o2 - cjminus1_o2)/(2*GENERAL_PARAMETERS.dz)               
                    
                    else: #Pesz > 2
                        if j == 0: #first cell
                            dc_o2 = cjplus1_o2 - 2*cj_o2 + cmminus1_o2
                            dc_dz_o2 = (cj_o2 - cmminus1_o2)/GENERAL_PARAMETERS.dz
            
                        elif j == (SZ.m_sz - 1): #last cell
                            dc_o2 = cj_o2 - 2*cj_o2 + cjminus1_o2
                            dc_dz_o2 = (cj_o2 - cjminus1_o2)/GENERAL_PARAMETERS.dz
                            
                        else:
                            dc_o2 = cjplus1_o2 - 2*cj_o2 + cjminus1_o2
                            dc_dz_o2 = (cj_o2 - cjminus1_o2)/GENERAL_PARAMETERS.dz            
                    
                
                if USZ.m_usz == (n - 1):
                    dc_o2 = cj_o2 - 2*cj_o2 + cjminus1_o2
                    dc_dz_o2 = (cj_o2 - cjminus1_o2)/(2*GENERAL_PARAMETERS.dz)
                    
                if USZ.m_usz == 0:
                    if Pesz_o2 <= 2:
                        if j == 0: #first cell
                            dc_o2 = cjplus1_o2 - 2*cj_o2 + cpi_o2
                            dc_dz_o2 = (cjplus1_o2 - cpi_o2)/(2*GENERAL_PARAMETERS.dz)
            
                        elif j == (SZ.m_sz - 1): #last cell
                            dc_o2 = cj_o2 - 2*cj_o2 + cjminus1_o2
                            dc_dz_o2 = (cj_o2 - cjminus1_o2)/GENERAL_PARAMETERS.dz
                            
                        else:
                            dc_o2 = cjplus1_o2 - 2*cj_o2 + cjminus1_o2
                            dc_dz_o2 = (cjplus1_o2 - cjminus1_o2)/(2*GENERAL_PARAMETERS.dz)               
                    
                    else: #Pusz > 2
                        if j == 0: #first cell
                            dc_o2 = cjplus1_o2 - 2*cj_o2 + cpi_o2
                            dc_dz_o2 = (cj_o2 - cpi_o2)/GENERAL_PARAMETERS.dz
            
                        elif j == (SZ.m_sz - 1): #last cell
                            dc_o2 = cj_o2 - 2*cj_o2 + cjminus1_o2
                            dc_dz_o2 = (cj_o2 - cjminus1_o2)/GENERAL_PARAMETERS.dz
                            
                        else:
                            dc_o2 = cjplus1_o2 - 2*cj_o2 + cjminus1_o2
                            dc_dz_o2 = (cj_o2 - cjminus1_o2)/GENERAL_PARAMETERS.dz                
                    
                delta_c_o2 = O2.f_transport(WFR.tteta_sz[t], teta_b_iplus1, cj_o2, cs_o2, dc_o2, dc_dz_o2, 0, 0, O2.D_o2, UFi_sz, Rxi_3_o2, GENERAL_PARAMETERS.dt, SOIL_PLANT.ro, SOIL_PLANT.f, GENERAL_PARAMETERS.dz)
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
                Pesz_nh4 = SZ.f_peclet(UFi_sz, NH4.D_nh4, GENERAL_PARAMETERS.dz)
                
                if USZ.m_usz < (n-1):
                    if Pesz_nh4 <= 2:
                        if j == 0: #first cell
                            dc_nh4 = cjplus1_nh4 - 2*cj_nh4 + cmminus1_nh4
                            dc_dz_nh4 = (cjplus1_nh4 - cmminus1_nh4)/(2*GENERAL_PARAMETERS.dz)
            
                        elif j == (SZ.m_sz - 1): #last cell
                            dc_nh4 = cj_nh4 - 2*cj_nh4 + cjminus1_nh4
                            dc_dz_nh4 = (cj_nh4 - cjminus1_nh4)/GENERAL_PARAMETERS.dz
                            
                        else:
                            dc_nh4 = cjplus1_nh4 - 2*cj_nh4 + cjminus1_nh4
                            dc_dz_nh4 = (cjplus1_nh4 - cjminus1_nh4)/(2*GENERAL_PARAMETERS.dz)               
                    
                    else: #Pesz > 2
                        if j == 0: #first cell
                            dc_nh4 = cjplus1_nh4 - 2*cj_nh4 + cmminus1_nh4
                            dc_dz_nh4 = (cj_nh4 - cmminus1_nh4)/GENERAL_PARAMETERS.dz
            
                        elif j == (SZ.m_sz - 1): #last cell
                            dc_nh4 = cj_nh4 - 2*cj_nh4 + cjminus1_nh4
                            dc_dz_nh4 = (cj_nh4 - cjminus1_nh4)/GENERAL_PARAMETERS.dz
                            
                        else:
                            dc_nh4 = cjplus1_nh4 - 2*cj_nh4 + cjminus1_nh4
                            dc_dz_nh4 = (cj_nh4 - cjminus1_nh4)/GENERAL_PARAMETERS.dz            
                    
                
                if USZ.m_usz == (n - 1):
                    dc_nh4 = cj_nh4 - 2*cj_nh4 + cjminus1_nh4
                    dc_dz_nh4 = (cj_nh4 - cjminus1_nh4)/(2*GENERAL_PARAMETERS.dz)
                    
                if USZ.m_usz == 0:
                    if Pesz_nh4 <= 2:
                        if j == 0: #first cell
                            dc_nh4 = cjplus1_nh4 - 2*cj_nh4 + cpi_nh4
                            dc_dz_nh4 = (cjplus1_nh4 - cpi_nh4)/(2*GENERAL_PARAMETERS.dz)
            
                        elif j == (SZ.m_sz -1): #last cell
                            dc_nh4 = cj_nh4 - 2*cj_nh4 + cjminus1_nh4
                            dc_dz_nh4 = (cj_nh4 - cjminus1_nh4)/GENERAL_PARAMETERS.dz
                            
                        else:
                            dc_nh4 = cjplus1_nh4 - 2*cj_nh4 + cjminus1_nh4
                            dc_dz_nh4 = (cjplus1_nh4 - cjminus1_nh4)/(2*GENERAL_PARAMETERS.dz)               
                    
                    else: #Pusz > 2
                        if j == 0: #first cell
                            dc_nh4 = cjplus1_nh4 - 2*cj_nh4 + cpi_nh4
                            dc_dz_nh4 = (cj_nh4 - cpi_nh4)/GENERAL_PARAMETERS.dz
            
                        elif j == (SZ.m_sz - 1): #last cell
                            dc_nh4 = cj_nh4 - 2*cj_nh4 + cjminus1_nh4
                            dc_dz_nh4 = (cj_nh4 - cjminus1_nh4)/GENERAL_PARAMETERS.dz
                            
                        else:
                            dc_nh4 = cjplus1_nh4 - 2*cj_nh4 + cjminus1_nh4
                            dc_dz_nh4 = (cj_nh4 - cjminus1_nh4)/GENERAL_PARAMETERS.dz                
                    
                delta_c_nh4 = NH4.f_transport(WFR.tteta_sz[t], teta_b_iplus1, cj_nh4, cs_nh4, dc_nh4, dc_dz_nh4, kads2_nh4, kdes2_nh4, NH4.D_nh4, UFi_sz, Rxi_3_nh4, GENERAL_PARAMETERS.dt, SOIL_PLANT.ro, SOIL_PLANT.f, GENERAL_PARAMETERS.dz)
#                 print('WFR.tteta_sz[t]: ', WFR.tteta_sz[t], ', teta_b_iplus1: ', teta_b_iplus1, ', cj_nh4: ', cj_nh4, ', cs_nh4: ', cs_nh4, ', dc_nh4: ', dc_nh4, ', dc_dz_nh4: ', dc_dz_nh4, ', kads2_nh4: ', kads2_nh4, ', kdes2_nh4: ', kdes2_nh4, ', NH4.D_nh4: ', NH4.D_nh4, ', UFi_sz: ', UFi_sz, ', Rxi_3_nh4: ', Rxi_3_nh4)
#                 print('delta_c_nh4: ', delta_c_nh4)
                
                if teta_b_iplus1 > 0:
                    ci1_nh4 = cj_nh4 + delta_c_nh4
                    #print('calculo final - cj_nh4, WFR.tteta_sz[t], teta_b_iplus1, delta_c_nh4: ', cj_nh4, ', ', WFR.tteta_sz[t], ', ', teta_b_iplus1, ', ', delta_c_nh4)
                else:
                    ci1_nh4 = 0
                #print('ci1_nh4: ', ci1_nh4)
                    
                if ci1_nh4 <= 0.0000000000000001:
                    ci1_nh4 = 0
                else:
                    ci1_nh4 = ci1_nh4 
                                         
                cj_i1_nh4.append(ci1_nh4)                
                    
    ### Nitrate           
                Pesz_no3 = SZ.f_peclet(UFi_sz, NO3.D_no3, GENERAL_PARAMETERS.dz)
                
                if USZ.m_usz < (n-1):
                    if Pesz_no3 <= 2:
                        if j == 0: #first cell
                            dc_no3 = cjplus1_no3 - 2*cj_no3 + cmminus1_no3
                            dc_dz_no3 = (cjplus1_no3 - cmminus1_no3)/(2*GENERAL_PARAMETERS.dz)
            
                        elif j == (SZ.m_sz - 1): #last cell
                            dc_no3 = cj_no3 - 2*cj_no3 + cjminus1_no3
                            dc_dz_no3 = (cj_no3 - cjminus1_no3)/GENERAL_PARAMETERS.dz
                            
                        else:
                            dc_no3 = cjplus1_no3 - 2*cj_no3 + cjminus1_no3
                            dc_dz_no3 = (cjplus1_no3 - cjminus1_no3)/(2*GENERAL_PARAMETERS.dz)               
                    
                    else: #Pesz > 2
                        if j == 0: #first cell
                            dc_no3 = cjplus1_no3 - 2*cj_no3 + cmminus1_no3
                            dc_dz_no3 = (cj_no3 - cmminus1_no3)/GENERAL_PARAMETERS.dz
            
                        elif j == (SZ.m_sz - 1): #last cell
                            dc_no3 = cj_no3 - 2*cj_no3 + cjminus1_no3
                            dc_dz_no3 = (cj_no3 - cjminus1_no3)/GENERAL_PARAMETERS.dz
                            
                        else:
                            dc_no3 = cjplus1_no3 - 2*cj_no3 + cjminus1_no3
                            dc_dz_no3 = (cj_no3 - cjminus1_no3)/GENERAL_PARAMETERS.dz            
                    
                
                if USZ.m_usz == (n - 1):
                    dc_no3 = cj_no3 - 2*cj_no3 + cjminus1_no3
                    dc_dz_no3 = (cj_no3 - cjminus1_no3)/(2*GENERAL_PARAMETERS.dz)
                    
                if USZ.m_usz == 0:
                    if Pesz_no3 <= 2:
                        if j == 0: #first cell
                            dc_no3 = cjplus1_no3 - 2*cj_no3 + cpi_no3
                            dc_dz_no3 = (cjplus1_no3 - cpi_no3)/(2*GENERAL_PARAMETERS.dz)
            
                        elif j == (SZ.m_sz - 1): #last cell
                            dc_no3 = cj_no3 - 2*cj_no3 + cjminus1_no3
                            dc_dz_no3 = (cj_no3 - cjminus1_no3)/GENERAL_PARAMETERS.dz
                            
                        else:
                            dc_no3 = cjplus1_no3 - 2*cj_no3 + cjminus1_no3
                            dc_dz_no3 = (cjplus1_no3 - cjminus1_no3)/(2*GENERAL_PARAMETERS.dz)               
                    
                    else: #Pusz > 2
                        if j == 0: #first cell
                            dc_no3 = cjplus1_no3 - 2*cj_no3 + cpi_no3
                            dc_dz_no3 = (cj_no3 - cpi_no3)/GENERAL_PARAMETERS.dz
            
                        elif j == (SZ.m_sz - 1): #last cell
                            dc_no3 = cj_no3 - 2*cj_no3 + cjminus1_no3
                            dc_dz_no3 = (cj_no3 - cjminus1_no3)/GENERAL_PARAMETERS.dz
                            
                        else:
                            dc_no3 = cjplus1_no3 - 2*cj_no3 + cjminus1_no3
                            dc_dz_no3 = (cj_no3 - cjminus1_no3)/GENERAL_PARAMETERS.dz                
                    
                delta_c_no3 = NO3.f_transportftransp(WFR.tteta_sz[t], teta_b_iplus1, cj_no3, cs_no3, dc_no3, dc_dz_no3, 0, 0, NO3.D_no3, UFi_sz, Rxi_3_no3, GENERAL_PARAMETERS.dt, SOIL_PLANT.ro, SOIL_PLANT.f, GENERAL_PARAMETERS.dz)
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
                Pesz_doc = SZ.f_peclet(UFi_sz, DOC.D_doc, GENERAL_PARAMETERS.dz)
                
                if USZ.m_usz < (n-1):
                    if Pesz_doc <= 2:
                        if j == 0: #first cell
                            dc_doc = cjplus1_doc - 2*cj_doc + cmminus1_doc
                            dc_dz_doc = (cjplus1_doc - cmminus1_doc)/(2*GENERAL_PARAMETERS.dz)
            
                        elif j == (SZ.m_sz - 1): #last cell
                            dc_doc = cj_doc - 2*cj_doc + cjminus1_doc
                            dc_dz_doc = (cj_doc - cjminus1_doc)/GENERAL_PARAMETERS.dz
                            
                        else:
                            dc_doc = cjplus1_doc - 2*cj_doc + cjminus1_doc
                            dc_dz_doc = (cjplus1_doc - cjminus1_doc)/(2*GENERAL_PARAMETERS.dz)               
                    
                    else: #Pesz > 2
                        if j == 0: #first cell
                            dc_doc = cjplus1_doc - 2*cj_doc + cmminus1_doc
                            dc_dz_doc = (cj_doc - cmminus1_doc)/GENERAL_PARAMETERS.dz
            
                        elif j == (SZ.m_sz - 1): #last cell
                            dc_doc = cj_doc - 2*cj_doc + cjminus1_doc
                            dc_dz_doc = (cj_doc - cjminus1_doc)/GENERAL_PARAMETERS.dz
                            
                        else:
                            dc_doc = cjplus1_doc - 2*cj_doc + cjminus1_doc
                            dc_dz_doc = (cj_doc - cjminus1_doc)/GENERAL_PARAMETERS.dz            
                    
                
                if USZ.m_usz == (n - 1):
                    dc_doc = cj_doc - 2*cj_doc + cjminus1_doc
                    dc_dz_doc = (cj_doc - cjminus1_doc)/(2*GENERAL_PARAMETERS.dz)
                    
                if USZ.m_usz == 0:
                    if Pesz_doc <= 2:
                        if j == 0: #first cell
                            dc_doc = cjplus1_doc - 2*cj_doc + cpi_doc
                            dc_dz_doc = (cjplus1_doc - cpi_doc)/(2*GENERAL_PARAMETERS.dz)
            
                        elif j == (SZ.m_sz - 1): #last cell
                            dc_doc = cj_doc - 2*cj_doc + cjminus1_doc
                            dc_dz_doc = (cj_doc - cjminus1_doc)/GENERAL_PARAMETERS.dz
                            
                        else:
                            dc_doc = cjplus1_doc - 2*cj_doc + cjminus1_doc
                            dc_dz_doc = (cjplus1_doc - cjminus1_doc)/(2*GENERAL_PARAMETERS.dz)               
                    
                    else: #Pusz > 2
                        if j == 0: #first cell
                            dc_doc = cjplus1_doc - 2*cj_doc + cpi_doc
                            dc_dz_doc = (cj_doc - cpi_doc)/GENERAL_PARAMETERS.dz
            
                        elif j == (SZ.m_sz - 1): #last cell
                            dc_doc = cj_doc - 2*cj_doc + cjminus1_doc
                            dc_dz_doc = (cj_doc - cjminus1_doc)/GENERAL_PARAMETERS.dz
                            
                        else:
                            dc_doc = cjplus1_doc - 2*cj_doc + cjminus1_doc
                            dc_dz_doc = (cj_doc - cjminus1_doc)/GENERAL_PARAMETERS.dz                
                    
                delta_c_doc = DOC.f_transport(WFR.tteta_sz[t], teta_b_iplus1, cj_doc, cs_doc, dc_doc, dc_dz_doc, kads2_doc, kdes2_doc, DOC.D_doc, UFi_sz, Rxi_3_doc, GENERAL_PARAMETERS.dt, SOIL_PLANT.ro, SOIL_PLANT.f, GENERAL_PARAMETERS.dz)
                if teta_b_iplus1 > 0:
                    ci1_doc = cj_doc + delta_c_doc
                else:
                    ci1_doc = 0
                    
                if ci1_doc <= 0.0000000000000001:
                    ci1_doc = 0
                else:
                    ci1_doc = ci1_doc 
                                         
                cj_i1_doc.append(ci1_doc)
                        
        
        
        
        #Corrective step 
        
        if hpipe > 0:
            ### Oxygen
            Mstor_o2_ast = sum(O2.c_usz[t])*PZ.Ab*WFR.tteta_usz[t]*WFR.thusz[t]*1000/USZ.m_usz + sum(O2.c_sz[t])*PZ.Ab*WFR.tteta_sz[t]*thsz[t]*1000/SZ.m_sz
            O2.Mstor_ast_list.append(Mstor_o2_ast)
            
            Msoil_o2_a = 0
            Msoil_o2 = 0
            O2.Msoil_acum.append(0)
            
            MRx_o2 = - (sum(O2.Rx_usz[t])*PZ.Ab*WFR.tteta_usz[t]*WFR.thusz[t]*1000/USZ.m_usz + sum(O2.Rx_sz[t])*PZ.Ab*WFR.tteta_sz[t]*thsz[t]*1000/SZ.m_sz)
            if t == 0:
                O2.MRx_acum.append(MRx_o2)
            else:
                O2.MRx_acum.append(MRx_o2 + O2.MRx_acum[-1])
    
    
            Min_o2 = tQin[t]*INFLOW_PARAMETERS.cin_o2[t]*dt*1000
            if t == 0:
                O2.Min_acum.append(Min_o2)
            else:
                O2.Min_acum.append(Min_o2 + O2.Min_acum[-1])
    
            Mover_o2 = WFR.tQover[t]*O2.cp[t]*dt*1000
            if t == 0:
                O2.Mover_acum.append(Mover_o2)
            else:
                O2.Mover_acum.append(Mover_o2 + O2.Mover_acum[-1])
             
            Mpipe_o2 = WFR.tQpipe[t]*O2.c_sz[t][SZ.m_sz-1]*dt*1000
            if t == 0:
                O2.Mpipe_acum.append(Mpipe_o2)
            else:
                O2.Mpipe_acum.append(Mpipe_o2 + O2.Mpipe_acum[-1])
             
            Minfsz_o2 = WFR.tQinfsz[t]*O2.c_sz[t][SZ.m_sz-1]*dt*1000
            if t == 0:
                O2.Minfsz_acum.append(Minfsz_o2)
            else:
                O2.Minfsz_acum.append(Minfsz_o2 + O2.Minfsz_acum[-1])
             
            Met_o2 = WFR.tQet[t]*cl_i1_o2[0]*dt*1000
            if t == 0:
                O2.Met_acum.append(Met_o2)
            else:
                O2.Met_acum.append(Met_o2 + O2.Met_acum[-1])
             
            Mpz_o2 = WFR.thpEND[t]*PZ.Ab*1000*O2.cp[t]
            O2.Mpz_list.append(Mpz_o2)
             
        
            Mstor_o2_mb = O2.Min_acum[-1] - O2.Mpz_list[-1] - O2.Mover_acum[-1] - O2.Mpipe_acum[-1] - O2.Minfsz_acum[-1] - O2.Met_acum[-1] - O2.Msoil_acum[-1] - O2.MRx_acum[-1]
            O2.Mstor_mb_list.append(Mstor_o2_mb)
            
            delta_ast_o2_usz = (Mstor_o2_mb - Mstor_o2_ast)/(PZ.Ab*WFR.tteta_usz[t]*WFR.thusz[t]*1000/USZ.m_usz)
            cl_i1_o2[1] = cl_i1_o2[1] + delta_ast_o2_usz
            if cl_i1_o2[1] > 0:
                cl_i1_o2[1] = cl_i1_o2[1]
            else:
                cl_i1_o2[1] = 0
                   
            ### Amonia
            
            Mstor_nh4_ast = sum(NH4.c_usz[t])*PZ.Ab*WFR.tteta_usz[t]*WFR.thusz[t]*1000/USZ.m_usz + sum(NH4.c_sz[t])*PZ.Ab*WFR.tteta_sz[t]*thsz[t]*1000/SZ.m_sz
            NH4.Mstor_ast_list.append(Mstor_nh4_ast)
            
            Msoil_nh4_a = NH4.cs_usz_a*SOIL_PLANT.ro*PZ.Ab*WFR.thusz[t]*1000 + NH4.cs_sz_a*SOIL_PLANT.ro*PZ.Ab*thsz[t]*1000
            Msoil_nh4 = sum(NH4.cs_usz[t])*SOIL_PLANT.ro*PZ.Ab*WFR.thusz[t]*1000/USZ.m_usz + sum(NH4.cs_sz[t])*ro*PZ.Ab*thsz[t]*1000/SZ.m_sz
            NH4.Msoil_acum.append(Msoil_nh4 - Msoil_nh4_a)
            
       
            MRx_nh4 = - (sum(NH4.Rx_usz[t])*PZ.Ab*WFR.tteta_usz[t]*WFR.thusz[t]*1000/USZ.m_usz + sum(NH4.Rx_sz[t])*PZ.Ab*WFR.tteta_sz[t]*thsz[t]*1000/SZ.m_sz)
            if t == 0:
                NH4.MRx_acum.append(MRx_nh4)
            else:
                NH4.MRx_acum.append(MRx_nh4 + NH4.MRx_acum[-1])
    
    
            Min_nh4 = tQin[t]*INFLOW_PARAMETERS.cin_nh4[t]*dt*1000
            if t == 0:
                NH4.Min_acum.append(Min_nh4)
            else:
                NH4.Min_acum.append(Min_nh4 + NH4.Min_acum[-1])
    
            Mover_nh4 = WFR.tQover[t]*NH4.cp[t]*dt*1000
            if t == 0:
                NH4.Mover_acum.append(Mover_nh4)
            else:
                NH4.Mover_acum.append(Mover_nh4 + NH4.Mover_acum[-1])
             
            Mpipe_nh4 = WFR.tQpipe[t]*NH4.c_sz[t][SZ.m_sz-1]*dt*1000
            if t == 0:
                NH4.Mpipe_acum.append(Mpipe_nh4)
            else:
                NH4.Mpipe_acum.append(Mpipe_nh4 + NH4.Mpipe_acum[-1])
             
            Minfsz_nh4 = WFR.tQinfsz[t]*NH4.c_sz[t][SZ.m_sz-1]*dt*1000
            if t == 0:
                NH4.Minfsz_acum.append(Minfsz_nh4)
            else:
                NH4.Minfsz_acum.append(Minfsz_nh4 + NH4.Minfsz_acum[-1])
             
            Met_nh4 = WFR.tQet[t]*cl_i1_nh4[0]*dt*1000
            if t == 0:
                NH4.Met_acum.append(Met_nh4)
            else:
                NH4.Met_acum.append(Met_nh4 + NH4.Met_acum[-1])
             
            Mpz_nh4 = WFR.thpEND[t]*PZ.Ab*1000*NH4.cp[t]
            NH4.Mpz_list.append(Mpz_nh4)
             
        
            Mstor_nh4_mb = NH4.Min_acum[-1] - NH4.Mpz_list[-1] - NH4.Mover_acum[-1] - NH4.Mpipe_acum[-1] - NH4.Minfsz_acum[-1] - NH4.Met_acum[-1] - NH4.Msoil_acum[-1] - NH4.MRx_acum[-1]
            NH4.Mstor_mb_list.append(Mstor_nh4_mb)
             
    
    #         delta_ast_nh4_usz = 0
    #         delta_ast_nh4_sz = 0
            
    #         delta_ast_nh4_usz = (Mstor_nh4_mb - Mstor_nh4_ast)/(n*(PZ.Ab*WFR.tteta_usz[t]*WFR.thusz[t]*1000/USZ.m_usz))
    #         delta_ast_nh4_sz = (Mstor_nh4_mb - Mstor_nh4_ast)/(n*(PZ.Ab*WFR.tteta_sz[t]*thsz[t]*1000/SZ.m_sz))
            
            delta_ast_nh4_usz = (Mstor_nh4_mb - Mstor_nh4_ast)/(PZ.Ab*WFR.tteta_usz[t]*WFR.thusz[t]*1000/USZ.m_usz)
            cl_i1_nh4[1] = cl_i1_nh4[1] + delta_ast_nh4_usz
            if cl_i1_nh4[1] > 0:
                cl_i1_nh4[1] = cl_i1_nh4[1]
            else:
                cl_i1_nh4[1] = 0
    #         sum_cusz_nh4 = sum(NH4.c_usz[t])
    #         sum_csz_nh4 = sum(NH4.c_sz[t])
    #         print('t: ', t, 'Mstor_trace: ', Mstor_trace, ', Mstor_mb: ', Mstor_mb, ', dif: ', (Mstor_mb - Mstor_trace))  
    #         print('t: ', t, ', delta_usz: ', delta_trace_nh4_usz, ', delta_sz: ', delta_trace_nh4_sz)
    #         print('teta_sm: ', WFR.tteta_usz[t], ', husz: ', WFR.thusz[t], ', teta_b: ', WFR.tteta_sz[t], ', hsz: ', thsz[t])
    #         print('sum_cusz_nh4: ', sum_cusz_nh4, ', sum_csz_nh4: ', sum_csz_nh4)             
    #         print('Qin: ', tQin[t], ',hpEND_i: ', WFR.thpEND[t], ', hpEND_i-1: ', WFR.thpEND[t-1], ',Qpipe: ', WFR.tQpipe[t], ', Qover: ', WFR.tQover[t], ', Qinfsz: ', WFR.tQinfsz[t], ', Qet: ', WFR.tQet[t])
    #         print('Cin: ', INFLOW_PARAMETERS.cin_nh4[t], ', Cp_i: ', NH4.cp[t], ', Cp_i-1: ', NH4.cp[t-1], ', Csz_ultimo: ', c_pipe)
    #         input()
    #                      
    #         for l in range(USZ.m_usz):
    #             cl_i1_nh4[l] = cl_i1_nh4[l] + delta_ast_nh4_usz
    #         for j in range(SZ.m_sz):
    #             cj_i1_nh4[j] = cj_i1_nh4[j] + delta_ast_nh4_sz 
        
            ### Nitrate
            
            Mstor_no3_ast = sum(NO3.c_usz[t])*PZ.Ab*WFR.tteta_usz[t]*WFR.thusz[t]*1000/USZ.m_usz + sum(NO3.c_sz[t])*PZ.Ab*WFR.tteta_sz[t]*thsz[t]*1000/SZ.m_sz
            NO3.Mstor_ast_list.append(Mstor_no3_ast)
            
            Msoil_no3_a = 0
            Msoil_no3 = 0
            NO3.Msoil_acum.append(0)
            
       
            MRx_no3 = - (sum(NO3.Rx_usz[t])*PZ.Ab*WFR.tteta_usz[t]*WFR.thusz[t]*1000/USZ.m_usz + sum(NO3.Rx_sz[t])*PZ.Ab*WFR.tteta_sz[t]*thsz[t]*1000/SZ.m_sz)
            if t == 0:
                NO3.MRx_acum.append(MRx_no3)
            else:
                NO3.MRx_acum.append(MRx_no3 + NO3.MRx_acum[-1])
    
    
            Min_no3 = tQin[t]*INFLOW_PARAMETERS.cin_no3[t]*dt*1000
            if t == 0:
                NO3.Min_acum.append(Min_no3)
            else:
                NO3.Min_acum.append(Min_no3 + NO3.Min_acum[-1])
    
            Mover_no3 = WFR.tQover[t]*NO3.cp[t]*dt*1000
            if t == 0:
                NO3.Mover_acum.append(Mover_no3)
            else:
                NO3.Mover_acum.append(Mover_no3 + NO3.Mover_acum[-1])
             
            Mpipe_no3 = WFR.tQpipe[t]*NO3.c_sz[t][SZ.m_sz-1]*dt*1000
            if t == 0:
                NO3.Mpipe_acum.append(Mpipe_no3)
            else:
                NO3.Mpipe_acum.append(Mpipe_no3 + NO3.Mpipe_acum[-1])
             
            Minfsz_no3 = WFR.tQinfsz[t]*NO3.c_sz[t][SZ.m_sz-1]*dt*1000
            if t == 0:
                NO3.Minfsz_acum.append(Minfsz_no3)
            else:
                NO3.Minfsz_acum.append(Minfsz_no3 + NO3.Minfsz_acum[-1])
             
            Met_no3 = WFR.tQet[t]*cl_i1_no3[0]*dt*1000
            if t == 0:
                NO3.Met_acum.append(Met_no3)
            else:
                NO3.Met_acum.append(Met_no3 + NO3.Met_acum[-1])
             
            Mpz_no3 = WFR.thpEND[t]*PZ.Ab*1000*NO3.cp[t]
            NO3.Mpz_list.append(Mpz_no3)
             
        
            Mstor_no3_mb = NO3.Min_acum[-1] - NO3.Mpz_list[-1] - NO3.Mover_acum[-1] - NO3.Mpipe_acum[-1] - NO3.Minfsz_acum[-1] - NO3.Met_acum[-1] - NO3.Msoil_acum[-1] - NO3.MRx_acum[-1]
            NO3.Mstor_mb_list.append(Mstor_no3_mb)
            
            delta_ast_no3_usz = (Mstor_no3_mb - Mstor_no3_ast)/(PZ.Ab*WFR.tteta_usz[t]*WFR.thusz[t]*1000/USZ.m_usz)
            cl_i1_no3[1] = cl_i1_no3[1] + delta_ast_no3_usz
            if cl_i1_no3[1] > 0:
                cl_i1_no3[1] = cl_i1_no3[1]
            else:
                cl_i1_no3[1] = 0  
            
            ### DOC
            
            Mstor_doc_ast = sum(DOC.c_usz[t])*PZ.Ab*WFR.tteta_usz[t]*WFR.thusz[t]*1000/USZ.m_usz + sum(DOC.c_sz[t])*PZ.Ab*WFR.tteta_sz[t]*thsz[t]*1000/SZ.m_sz
            DOC.Mstor_ast_list.append(Mstor_doc_ast)
            
            Msoil_doc_a = DOC.cs_usz_a*SOIL_PLANT.ro*PZ.Ab*WFR.thusz[t]*1000 + DOC.cs_sz_a*SOIL_PLANT.ro*PZ.Ab*thsz[t]*1000
            Msoil_doc = sum(DOC.cs_usz[t])*SOIL_PLANT.ro*PZ.Ab*WFR.thusz[t]*1000/USZ.m_usz + sum(DOC.cs_sz[t])*SOIL_PLANT.ro*PZ.Ab*thsz[t]*1000/SZ.m_sz
            DOC.Msoil_acum.append(Msoil_doc - Msoil_doc_a)
            
       
            MRx_doc = - (sum(DOC.Rx_usz[t])*PZ.Ab*WFR.tteta_usz[t]*WFR.thusz[t]*1000/USZ.m_usz + sum(DOC.Rx_sz[t])*PZ.Ab*WFR.tteta_sz[t]*thsz[t]*1000/SZ.m_sz)
            if t == 0:
                DOC.MRx_acum.append(MRx_doc)
            else:
                DOC.MRx_acum.append(MRx_doc + DOC.MRx_acum[-1])
    
    
            Min_doc = tQin[t]*INFLOW_PARAMETERS.cin_doc[t]*dt*1000
            if t == 0:
                DOC.Min_acum.append(Min_doc)
            else:
                DOC.Min_acum.append(Min_doc + DOC.Min_acum[-1])
    
            Mover_doc = WFR.tQover[t]*DOC.cp[t]*dt*1000
            if t == 0:
                DOC.Mover_acum.append(Mover_doc)
            else:
                DOC.Mover_acum.append(Mover_doc + DOC.Mover_acum[-1])
             
            Mpipe_doc = WFR.tQpipe[t]*DOC.c_sz[t][SZ.m_sz-1]*dt*1000
            if t == 0:
                DOC.Mpipe_acum.append(Mpipe_doc)
            else:
                DOC.Mpipe_acum.append(Mpipe_doc + DOC.Mpipe_acum[-1])
             
            Minfsz_doc = WFR.tQinfsz[t]*DOC.c_sz[t][SZ.m_sz-1]*dt*1000
            if t == 0:
                DOC.Minfsz_acum.append(Minfsz_doc)
            else:
                DOC.Minfsz_acum.append(Minfsz_doc + DOC.Minfsz_acum[-1])
             
            Met_doc = WFR.tQet[t]*cl_i1_doc[0]*dt*1000
            if t == 0:
                DOC.Met_acum.append(Met_doc)
            else:
                DOC.Met_acum.append(Met_doc + DOC.Met_acum[-1])
             
            Mpz_doc = WFR.thpEND[t]*PZ.Ab*1000*DOC.cp[t]
            DOC.Mpz_list.append(Mpz_doc)
             
        
            Mstor_doc_mb = DOC.Min_acum[-1] - DOC.Mpz_list[-1] - DOC.Mover_acum[-1] - DOC.Mpipe_acum[-1] - DOC.Minfsz_acum[-1] - DOC.Met_acum[-1] - DOC.Msoil_acum[-1] - DOC.MRx_acum[-1]
            DOC.Mstor_mb_list.append(Mstor_doc_mb)
            
            delta_ast_doc_usz = (Mstor_doc_mb - Mstor_doc_ast)/(PZ.Ab*WFR.tteta_usz[t]*WFR.thusz[t]*1000/USZ.m_usz)
            cl_i1_doc[1] = cl_i1_doc[1] + delta_ast_doc_usz
            if cl_i1_doc[1] > 0:
                cl_i1_doc[1] = cl_i1_doc[1]
            else:
                cl_i1_doc[1] = 0
        
        
        else: #if hpipe == 0
            
            ### Oxygen
            
            Mstor_o2_ast = sum(O2.c_usz[t])*PZ.Ab*WFR.tteta_usz[t]*WFR.thusz[t]*1000/USZ.m_usz
            O2.Mstor_ast_list.append(Mstor_o2_ast)
            
            Msoil_o2_a = 0
            Msoil_o2 = 0
            O2.Msoil_acum.append(0)
            
            MRx_o2 = - (sum(O2.Rx_usz[t])*PZ.Ab*WFR.tteta_usz[t]*WFR.thusz[t]*1000/USZ.m_usz)
            if t == 0:
                O2.MRx_acum.append(MRx_o2)
            else:
                O2.MRx_acum.append(MRx_o2 + O2.MRx_acum[-1])
    
    
            Min_o2 = tQin[t]*INFLOW_PARAMETERS.cin_o2[t]*dt*1000
            if t == 0:
                O2.Min_acum.append(Min_o2)
            else:
                O2.Min_acum.append(Min_o2 + O2.Min_acum[-1])
    
            Mover_o2 = WFR.tQover[t]*O2.cp[t]*dt*1000
            if t == 0:
                O2.Mover_acum.append(Mover_o2)
            else:
                O2.Mover_acum.append(Mover_o2 + O2.Mover_acum[-1])
             
            Mpipe_o2 = WFR.tQpipe[t]*O2.c_usz[t][USZ.m_usz-1]*dt*1000
            if t == 0:
                O2.Mpipe_acum.append(Mpipe_o2)
            else:
                O2.Mpipe_acum.append(Mpipe_o2 + O2.Mpipe_acum[-1])
             
            Minfsz_o2 = WFR.tQinfsz[t]*O2.c_usz[t][USZ.m_usz-1]*dt*1000
            if t == 0:
                O2.Minfsz_acum.append(Minfsz_o2)
            else:
                O2.Minfsz_acum.append(Minfsz_o2 + O2.Minfsz_acum[-1])
             
            Met_o2 = WFR.tQet[t]*cl_i1_o2[0]*dt*1000
            if t == 0:
                O2.Met_acum.append(Met_o2)
            else:
                O2.Met_acum.append(Met_o2 + O2.Met_acum[-1])
             
            Mpz_o2 = WFR.thpEND[t]*PZ.Ab*1000*O2.cp[t]
            O2.Mpz_list.append(Mpz_o2)
             
        
            Mstor_o2_mb = O2.Min_acum[-1] - O2.Mpz_list[-1] - O2.Mover_acum[-1] - O2.Mpipe_acum[-1] - O2.Minfsz_acum[-1] - O2.Met_acum[-1] - O2.Msoil_acum[-1] - O2.MRx_acum[-1]
            O2.Mstor_mb_list.append(Mstor_o2_mb)
            
            delta_ast_o2_usz = (Mstor_o2_mb - Mstor_o2_ast)/(PZ.Ab*WFR.tteta_usz[t]*WFR.thusz[t]*1000/USZ.m_usz)
            cl_i1_o2[1] = cl_i1_o2[1] + delta_ast_o2_usz
            if cl_i1_o2[1] > 0:
                cl_i1_o2[1] = cl_i1_o2[1]
            else:
                cl_i1_o2[1] = 0
            
            ### Amonia
        
            Mstor_nh4_ast = sum(NH4.c_usz[t])*PZ.Ab*WFR.tteta_usz[t]*WFR.thusz[t]*1000/USZ.m_usz
            NH4.Mstor_ast_list.append(Mstor_nh4_ast)
            
            Msoil_nh4_a = NH4.cs_usz_a*SOIL_PLANT.ro*PZ.Ab*WFR.thusz[t]*1000
            Msoil_nh4 = sum(NH4.cs_usz[t])*SOIL_PLANT.ro*PZ.Ab*WFR.thusz[t]*1000/USZ.m_usz
            NH4.Msoil_acum.append(Msoil_nh4 - Msoil_nh4_a)
            
       
            MRx_nh4 = - (sum(NH4.Rx_usz[t])*PZ.Ab*WFR.tteta_usz[t]*WFR.thusz[t]*1000/USZ.m_usz)
            if t == 0:
                NH4.MRx_acum.append(MRx_nh4)
            else:
                NH4.MRx_acum.append(MRx_nh4 + NH4.MRx_acum[-1])
    
    
            Min_nh4 = tQin[t]*INFLOW_PARAMETERS.cin_nh4[t]*dt*1000
            if t == 0:
                NH4.Min_acum.append(Min_nh4)
            else:
                NH4.Min_acum.append(Min_nh4 + NH4.Min_acum[-1])
    
            Mover_nh4 = WFR.tQover[t]*NH4.cp[t]*dt*1000
            if t == 0:
                NH4.Mover_acum.append(Mover_nh4)
            else:
                NH4.Mover_acum.append(Mover_nh4 + NH4.Mover_acum[-1])
             
            Mpipe_nh4 = WFR.tQpipe[t]*NH4.c_usz[t][USZ.m_usz-1]*dt*1000
            if t == 0:
                NH4.Mpipe_acum.append(Mpipe_nh4)
            else:
                NH4.Mpipe_acum.append(Mpipe_nh4 + NH4.Mpipe_acum[-1])
             
            Minfsz_nh4 = WFR.tQinfsz[t]*NH4.c_usz[t][USZ.m_usz-1]*dt*1000
            if t == 0:
                NH4.Minfsz_acum.append(Minfsz_nh4)
            else:
                NH4.Minfsz_acum.append(Minfsz_nh4 + NH4.Minfsz_acum[-1])
             
            Met_nh4 = WFR.tQet[t]*cl_i1_nh4[0]*dt*1000
            if t == 0:
                NH4.Met_acum.append(Met_nh4)
            else:
                NH4.Met_acum.append(Met_nh4 + NH4.Met_acum[-1])
             
            Mpz_nh4 = WFR.thpEND[t]*PZ.Ab*1000*NH4.cp[t]
            NH4.Mpz_list.append(Mpz_nh4)
             
        
            Mstor_nh4_mb = NH4.Min_acum[-1] - NH4.Mpz_list[-1] - NH4.Mover_acum[-1] - NH4.Mpipe_acum[-1] - NH4.Minfsz_acum[-1] - NH4.Met_acum[-1] - NH4.Msoil_acum[-1] - NH4.MRx_acum[-1]
            NH4.Mstor_mb_list.append(Mstor_nh4_mb)
             
            
            delta_ast_nh4_usz = (Mstor_nh4_mb - Mstor_nh4_ast)/(PZ.Ab*WFR.tteta_usz[t]*WFR.thusz[t]*1000/USZ.m_usz)
            cl_i1_nh4[1] = cl_i1_nh4[1] + delta_ast_nh4_usz
            if cl_i1_nh4[1] > 0:
                cl_i1_nh4[1] = cl_i1_nh4[1]
            else:
                cl_i1_nh4[1] = 0

            ### Nitrate
            
            Mstor_no3_ast = sum(NO3.c_usz[t])*PZ.Ab*WFR.tteta_usz[t]*WFR.thusz[t]*1000/USZ.m_usz
            NO3.Mstor_ast_list.append(Mstor_no3_ast)
            
            Msoil_no3_a = 0
            Msoil_no3 = 0
            NO3.Msoil_acum.append(0)
            
       
            MRx_no3 = - (sum(NO3.Rx_usz[t])*PZ.Ab*WFR.tteta_usz[t]*WFR.thusz[t]*1000/USZ.m_usz)
            if t == 0:
                NO3.MRx_acum.append(MRx_no3)
            else:
                NO3.MRx_acum.append(MRx_no3 + NO3.MRx_acum[-1])
    
    
            Min_no3 = tQin[t]*INFLOW_PARAMETERS.cin_no3[t]*dt*1000
            if t == 0:
                NO3.Min_acum.append(Min_no3)
            else:
                NO3.Min_acum.append(Min_no3 + NO3.Min_acum[-1])
    
            Mover_no3 = WFR.tQover[t]*NO3.cp[t]*dt*1000
            if t == 0:
                NO3.Mover_acum.append(Mover_no3)
            else:
                NO3.Mover_acum.append(Mover_no3 + NO3.Mover_acum[-1])
             
            Mpipe_no3 = WFR.tQpipe[t]*NO3.c_usz[t][USZ.m_usz-1]*dt*1000
            if t == 0:
                NO3.Mpipe_acum.append(Mpipe_no3)
            else:
                NO3.Mpipe_acum.append(Mpipe_no3 + NO3.Mpipe_acum[-1])
             
            Minfsz_no3 = WFR.tQinfsz[t]*NO3.c_usz[t][USZ.m_usz-1]*dt*1000
            if t == 0:
                NO3.Minfsz_acum.append(Minfsz_no3)
            else:
                NO3.Minfsz_acum.append(Minfsz_no3 + NO3.Minfsz_acum[-1])
             
            Met_no3 = WFR.tQet[t]*cl_i1_no3[0]*dt*1000
            if t == 0:
                NO3.Met_acum.append(Met_no3)
            else:
                NO3.Met_acum.append(Met_no3 + NO3.Met_acum[-1])
             
            Mpz_no3 = WFR.thpEND[t]*PZ.Ab*1000*NO3.cp[t]
            NO3.Mpz_list.append(Mpz_no3)
             
        
            Mstor_no3_mb = NO3.Min_acum[-1] - NO3.Mpz_list[-1] - NO3.Mover_acum[-1] - NO3.Mpipe_acum[-1] - NO3.Minfsz_acum[-1] - NO3.Met_acum[-1] - NO3.Msoil_acum[-1] - NO3.MRx_acum[-1]
            NO3.Mstor_mb_list.append(Mstor_no3_mb)
            
            delta_ast_no3_usz = (Mstor_no3_mb - Mstor_no3_ast)/(PZ.Ab*WFR.tteta_usz[t]*WFR.thusz[t]*1000/USZ.m_usz)
            cl_i1_no3[1] = cl_i1_no3[1] + delta_ast_no3_usz
            if cl_i1_no3[1] > 0:
                cl_i1_no3[1] = cl_i1_no3[1]
            else:
                cl_i1_no3[1] = 0  
            
            ### DOC
            
            Mstor_doc_ast = sum(DOC.c_usz[t])*PZ.Ab*WFR.tteta_usz[t]*WFR.thusz[t]*1000/USZ.m_usz
            DOC.Mstor_ast_list.append(Mstor_doc_ast)
            
            Msoil_doc_a = DOC.cs_usz_a*SOIL_PLANT.ro*PZ.Ab*WFR.thusz[t]*1000
            Msoil_doc = sum(DOC.cs_usz[t])*SOIL_PLANT.ro*PZ.Ab*WFR.thusz[t]*1000/USZ.m_usz
            DOC.Msoil_acum.append(Msoil_doc - Msoil_doc_a)
            
       
            MRx_doc = - (sum(DOC.Rx_usz[t])*PZ.Ab*WFR.tteta_usz[t]*WFR.thusz[t]*1000/USZ.m_usz)
            if t == 0:
                DOC.MRx_acum.append(MRx_doc)
            else:
                DOC.MRx_acum.append(MRx_doc + DOC.MRx_acum[-1])
    
    
            Min_doc = tQin[t]*INFLOW_PARAMETERS.cin_doc[t]*dt*1000
            if t == 0:
                DOC.Min_acum.append(Min_doc)
            else:
                DOC.Min_acum.append(Min_doc + DOC.Min_acum[-1])
    
            Mover_doc = WFR.tQover[t]*DOC.cp[t]*dt*1000
            if t == 0:
                DOC.Mover_acum.append(Mover_doc)
            else:
                DOC.Mover_acum.append(Mover_doc + DOC.Mover_acum[-1])
             
            Mpipe_doc = WFR.tQpipe[t]*DOC.c_usz[t][USZ.m_usz-1]*dt*1000
            if t == 0:
                DOC.Mpipe_acum.append(Mpipe_doc)
            else:
                DOC.Mpipe_acum.append(Mpipe_doc + DOC.Mpipe_acum[-1])
             
            Minfsz_doc = WFR.tQinfsz[t]*DOC.c_usz[t][USZ.m_usz-1]*dt*1000
            if t == 0:
                DOC.Minfsz_acum.append(Minfsz_doc)
            else:
                DOC.Minfsz_acum.append(Minfsz_doc + DOC.Minfsz_acum[-1])
             
            Met_doc = WFR.tQet[t]*cl_i1_doc[0]*dt*1000
            if t == 0:
                DOC.Met_acum.append(Met_doc)
            else:
                DOC.Met_acum.append(Met_doc + DOC.Met_acum[-1])
             
            Mpz_doc = WFR.thpEND[t]*PZ.Ab*1000*DOC.cp[t]
            DOC.Mpz_list.append(Mpz_doc)
             
        
            Mstor_doc_mb = DOC.Min_acum[-1] - DOC.Mpz_list[-1] - DOC.Mover_acum[-1] - DOC.Mpipe_acum[-1] - DOC.Minfsz_acum[-1] - DOC.Met_acum[-1] - DOC.Msoil_acum[-1] - DOC.MRx_acum[-1]
            DOC.Mstor_mb_list.append(Mstor_doc_mb)
            
            delta_ast_doc_usz = (Mstor_doc_mb - Mstor_doc_ast)/(PZ.Ab*WFR.tteta_usz[t]*WFR.thusz[t]*1000/USZ.m_usz)
            cl_i1_doc[1] = cl_i1_doc[1] + delta_ast_doc_usz
            if cl_i1_doc[1] > 0:
                cl_i1_doc[1] = cl_i1_doc[1]
            else:
                cl_i1_doc[1] = 0
            
       
        
    ## adding layers of USZ in time 
        O2.c_usz.append(cl_i1_o2)
        O2.Rx_usz.append(Rxl_o2)
               
        NH4.c_usz.append(cl_i1_nh4)
        NH4.cs_usz.append(csi_usz_nh4)
        NH4.Rx_usz.append(Rxl_nh4)
        
        NO3.c_usz.append(cl_i1_no3)
        NO3.Rx_usz.append(Rxl_no3)
        
        DOC.c_usz.append(cl_i1_doc)
        DOC.cs_usz.append(csi_usz_doc)
        DOC.Rx_usz.append(Rxl_doc)
        
    ## adding layers of SZ in time
        if hpipe > 0:
            O2.c_sz.append(cj_i1_o2)
            O2.Rx_sz.append(Rxj_o2)
            
            NH4.c_sz.append(cj_i1_nh4)
            NH4.cs_sz.append(csi_sz_nh4)
            NH4.Rx_sz.append(Rxj_nh4)
            
            NO3.c_sz.append(cj_i1_no3)
            NO3.Rx_sz.append(Rxj_no3)
            
            DOC.c_sz.append(cj_i1_doc)
            DOC.cs_sz.append(csi_sz_doc)
            DOC.Rx_sz.append(Rxj_doc)
            
    # **5. Transforming in dataframe **
    ## Oxygen
    data_usz_o2 = pd.DataFrame(O2.c_usz)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz' + str(a)
        column_name.append(name)
        a = a + 1
    data_usz_o2.set_axis(column_name, axis = 'columns', inplace = True)
    
    #print(data_usz_o2)
    
    data_sz_o2 = pd.DataFrame(O2.c_sz)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz' + str(a)
        column_name.append(name)
        a = a + 1
    data_sz_o2.set_axis(column_name, axis = 'columns', inplace = True)

    data_rx_usz_o2 = pd.DataFrame(O2.Rx_usz)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_usz_o2.set_axis(column_name, axis = 'columns', inplace = True)
    
    data_rx_sz_o2 = pd.DataFrame(O2.Rx_sz)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_sz_o2.set_axis(column_name, axis = 'columns', inplace = True)
    
    frames = [data_usz_o2, data_sz_o2, data_rx_usz_o2, data_rx_sz_o2]
    data_o2 = pd.concat((frames), axis = 1)

    O2.cp.append(0)
    data_o2['pz'] = O2.cp
    
    WFR.indice_n = list(range(len(O2.cp)))
    INFLOW_PARAMETERS.cin_o2.append(0)
    c_in = INFLOW_PARAMETERS.cin_o2[:len(WFR.indice_n)]
    #print('len c_in:', len(c_in), ', len_WFR.indice:', len(WFR.indice))
    data_o2['c_in'] = c_in 
    
    data_o2['t'] = WFR.indice_n

   
    ## Amonia
    data_usz_nh4 = pd.DataFrame(NH4.c_usz)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz' + str(a)
        column_name.append(name)
        a = a + 1
    data_usz_nh4.set_axis(column_name, axis = 'columns', inplace = True)
    
    data_sz_nh4 = pd.DataFrame(NH4.c_sz)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz' + str(a)
        column_name.append(name)
        a = a + 1
    data_sz_nh4.set_axis(column_name, axis = 'columns', inplace = True)
    
    data_s_usz_nh4 = pd.DataFrame(NH4.cs_usz)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz_soil' + str(a)
        column_name.append(name)
        a = a + 1
    data_s_usz_nh4.set_axis(column_name, axis = 'columns', inplace = True)
    
    data_s_sz_nh4 = pd.DataFrame(NH4.cs_sz)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz_soil' + str(a)
        column_name.append(name)
        a = a + 1
    data_s_sz_nh4.set_axis(column_name, axis = 'columns', inplace = True)
    
    data_rx_usz_nh4 = pd.DataFrame(NH4.Rx_usz)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_usz_nh4.set_axis(column_name, axis = 'columns', inplace = True)
    
    data_rx_sz_nh4 = pd.DataFrame(NH4.Rx_sz)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_sz_nh4.set_axis(column_name, axis = 'columns', inplace = True)
    
    
    frames = [data_usz_nh4, data_sz_nh4, data_s_usz_nh4, data_s_sz_nh4, data_rx_usz_nh4, data_rx_sz_nh4]
    data_nh4 = pd.concat((frames), axis = 1)
    NH4.cp.append(0)
    INFLOW_PARAMETERS.cin_nh4.append(0)
    c_in = INFLOW_PARAMETERS.cin_nh4[:len(WFR.indice_n)]
    #print('len c_in:', len(c_in), ', len_WFR.indice:', len(WFR.indice))
    data_nh4['c_in'] = c_in
    
    data_nh4['pz'] = NH4.cp
    #WFR.indice.append(len(WFR.indice)+1)
    data_nh4['t'] = WFR.indice_n
    
    
    ## Nitrate
    data_usz_no3 = pd.DataFrame(NO3.c_usz)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz' + str(a)
        column_name.append(name)
        a = a + 1
    data_usz_no3.set_axis(column_name, axis = 'columns', inplace = True)
    
    data_sz_no3 = pd.DataFrame(NO3.c_sz)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz' + str(a)
        column_name.append(name)
        a = a + 1
    data_sz_no3.set_axis(column_name, axis = 'columns', inplace = True)
    
    data_rx_usz_no3 = pd.DataFrame(NO3.Rx_usz)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_usz_no3.set_axis(column_name, axis = 'columns', inplace = True)
    
    data_rx_sz_no3 = pd.DataFrame(NO3.Rx_sz)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_sz_no3.set_axis(column_name, axis = 'columns', inplace = True)    
    
    frames = [data_usz_no3, data_sz_no3, data_rx_usz_no3, data_rx_sz_no3]
    data_no3 = pd.concat((frames), axis = 1)
    NO3.cp.append(0)
    data_no3['pz'] = NO3.cp
    
    INFLOW_PARAMETERS.cin_no3.append(0)
    c_in = INFLOW_PARAMETERS.cin_no3[:len(WFR.indice_n)]
    data_no3['c_in'] = c_in

    data_no3['t'] = WFR.indice_n
    
    ## DOC
    data_usz_doc = pd.DataFrame(DOC.c_usz)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz' + str(a)
        column_name.append(name)
        a = a + 1
    data_usz_doc.set_axis(column_name, axis = 'columns', inplace = True)
    
    data_sz_doc = pd.DataFrame(DOC.c_sz)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz' + str(a)
        column_name.append(name)
        a = a + 1
    data_sz_doc.set_axis(column_name, axis = 'columns', inplace = True)

    data_rx_usz_doc = pd.DataFrame(DOC.Rx_usz)
    a = 0
    column_name = []
    for i in range(USZ.m_usz):
        name = 'usz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_usz_doc.set_axis(column_name, axis = 'columns', inplace = True)
    
    data_rx_sz_doc = pd.DataFrame(DOC.Rx_sz)
    a = 0
    column_name = []
    for i in range(SZ.m_sz):
        name = 'sz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_sz_doc.set_axis(column_name, axis = 'columns', inplace = True)
    
    frames = [data_usz_doc, data_sz_doc, data_rx_usz_doc, data_rx_sz_doc]
    data_doc = pd.concat((frames), axis = 1)
    DOC.cp.append(0)
    data_doc['pz'] = DOC.cp
    
    INFLOW_PARAMETERS.cin_doc.append(0)
    c_in = INFLOW_PARAMETERS.cin_doc[:len(WFR.indice_n)]
    data_doc['c_in'] = c_in
    
    data_doc['t'] = WFR.indice_n



    return data_o2, data_nh4, data_no3, data_doc
    
if __name__ == '__main__':
    inicio = datetime.datetime.now()
    
    data_o2, data_nh4, data_no3, data_doc = run_Kin()

    data_nh4.to_csv('results_Kin_pf_nh4_2.csv', index = False)
    data_o2.to_csv('results_Kin_pf_o2.csv', index = False)
    data_no3.to_csv('results_Kin_pf_no3.csv', index = False)
    data_doc.to_csv('results_Kin_pf_doc.csv', index = False)
    
    fim = datetime.datetime.now()
    print ('Elapsed time: ', fim - inicio)
    print('Done!')
