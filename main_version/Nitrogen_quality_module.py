__author__ = 'marina batalini de macedo'
#!/usr/bin/env python
# coding: utf-8

import logging
import datetime
import numpy as np 
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
import configparser
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import math
#import seaborn as sns

print('Starting..')
inicio = datetime.datetime.now()
rodada = 'teste_4'

#logging.basicConfig(filename = 'testes.txt', level=logging.DEBUG, format=f'{rodada}: - %(asctime)s - %(levelname)s - %(message)s')
#logger = logging.getLogger(__name__)


setup_file = ".\\main_version\\Nitrogen_quality_module.ini"
setup = configparser.ConfigParser()
setup.read(setup_file)


# **1. Parameters and variables**

#Soil and plant characteristics
ro_pd = float(setup['SOIL_PLANT']['ro'])
f = float(setup['SOIL_PLANT']['f'])
lamda = float(setup['SOIL_PLANT']['lamda'])


root_fraction = float(setup['SOIL_PLANT']['root_fraction'])
c_o2_root = float(setup['SOIL_PLANT']['c_o2_root'])

#NH4 constants
D_nh4 = float(setup['NH4']['D_nh4'])
kads_nh4 = float(setup['NH4']['kads_nh4'])
kdes_nh4 = float(setup['NH4']['kdes_nh4'])
kads2_nh4 = float(setup['NH4']['kads2_nh4'])
kdes2_nh4 = float(setup['NH4']['kdes2_nh4'])
k_nh4_mb = float(setup['NH4']['k_nh4_mb'])  #Varia para as duas camadas ou nao?

Fm_nh4 = float(setup['NH4']['Fm_nh4'])
Km_nh4 = float(setup['NH4']['Km_nh4'])

#NO3 constants
D_no3 = float(setup['NO3']['D_no3'])

Fm_no3 = float(setup['NO3']['Fm_no3'])
Km_no3 = float(setup['NO3']['Km_no3'])

#O2 constants
D_o2 = float(setup['O2']['D_o2'])
K_o2 = float(setup['O2']['k_inib_o2'])
k_o2 = float(setup['O2']['k_o2'])

#DOC constants
D_doc = float(setup['DOC']['D_doc'])
fb_doc = float(setup['DOC']['fb_doc'])
bDOCd = float(setup['DOC']['bDOCd'])
KbDOC = float(setup['DOC']['KbDOC'])
k_doc = float(setup['DOC']['k_doc'])

kads_doc = float(setup['DOC']['kads_doc'])
kdes_doc = float(setup['DOC']['kdes_doc'])
kads2_doc = float(setup['DOC']['kads2_doc'])
kdes2_doc = float(setup['DOC']['kdes2_doc'])
k_doc_mb = float(setup['DOC']['k_doc_mb'])

#Nitrification constants
k_nit = float(setup['NITRIFICATION']['k_nit'])

#Denitrification constants
k_denit_pz = float(setup['DENITRIFICATION']['k_denit_pz'])
k_denit = float(setup['DENITRIFICATION']['k_denit'])


# **2. Input variables**
from main_version.Water_flow_module import *
run_W()

cin_nh4_file = pd.read_csv('.\\main_version\\Cinflow_NH4.csv')
cin_nh4_list = cin_nh4_file['nh4'].tolist()

cin_no3_file = pd.read_csv('.\\main_version\\Cinflow_NO3.csv')
cin_no3_list = cin_no3_file['no3'].tolist()

cin_o2_file = pd.read_csv('.\\main_version\\Cinflow_O2.csv')
cin_o2_list = cin_o2_file['o2'].tolist()

cin_doc_file = pd.read_csv('.\\main_version\\Cinflow_DOC.csv')
cin_doc_list = cin_doc_file['doc'].tolist()

# **3. Definitions of the equations**

n = 11 #number of cells to calculation
dz = L/n #0.0004 m valor bom.
d50 = 0.00028 #median grain size of the porous medium.
m_usz = round((L-hpipe)/dz) #numero de celulas na USZ
m_sz = n - m_usz #numero de celulas na SZ
ro = (1-nusz_ini)*ro_pd

##### Unitary flow to transport equations #####
def falfa_beta_usz(l):
# l: 0 - (m_usz - 1)
    alfa = (m_usz-1-l)/(m_usz-1)
    beta = l/(m_usz-1)
    
    return alfa, beta

def falfa_beta_sz(j):
# j: 0 - (m_sz-1)
    alfa2 = (m_sz -1 -j)/(m_sz - 1)
    beta2 = j/(m_sz - 1)
    
    return alfa2, beta2    
    
def fUF_usz(alfa, beta, I1, Qet_1, I2, Qhc, Qorif, Qinf_sz, teta_sm_i):
    if hpipe > 0:
        UF_usz = (alfa*(I1 - Qet_1) + beta*(I2 - Qhc))/(Ab*teta_sm_i)
    
    else:
        UF_usz = (alfa*(I1 - Qet_1) + beta*(Qorif + Qinf_sz - Qhc))/(Ab*teta_sm_i)
    
    return UF_usz

def fUF_sz(alfa2, beta2, I2, Qhc, Qet_2, Qorif, Qinf_sz, teta_b_i):
    UF_sz = (alfa2*(I2 - Qhc - Qet_2) + beta2*(Qorif + Qinf_sz))/(Ab*teta_b_i)
    
    return UF_sz

def fPeusz(UF_usz, D):
    if D > 0:
        Peusz = UF_usz*dz/D
    else:
        Peusz = 100
    
    return Peusz

def fPesz(UF_sz, D):
    if D > 0:
        Pesz = UF_sz*dz/D
    else:
        Pesz = 100
    
    return Pesz

##### Water mass balance to ponding zone #####
def fcp(cin,Qin_p, cp_a, I1, Qv, Rxi, hp, hp_a):
    #delta_cp = ((cin*Qin_p - cp_a*(I1 + Qv))*dt)/(hp*Ab) + Rxi*dt
    #cp = cp_a + delta_cp
    
    #cp = Rxi*dt + (cin*Qin_p*dt)/(hp*Ab + dt*(I1 + Qv))
    
    cp = (cp_a*hp_a*Ab + (cin*Qin_p - cp_a*(I1 + Qv) + Rxi*hp*Ab)*dt)/(hp*Ab)    
            
    return cp

##### Transport equations #####
def ftransp(teta_i, teta_iplus1, ci, cs_i, dc, dc_dz,kads, kdes, D, UF, Rx):
    if teta_i == 0:
        delta_c_i1 = 0
    elif teta_iplus1 == 0:
        delta_c_i1 = 0
    else:
        delta_c_i1 = ((1/teta_iplus1)*dt*(-teta_i*kads*ci + ro*kdes*cs_i + teta_i*(D*f*(dc/dz**2) - UF*dc_dz) + Rx)) 
        #delta_c_i1 = ((1/teta_i)*dt*(-teta_i*kads*ci + ro*kdes*cs_i + teta_i*(D*f*(dc/dz**2) - UF*dc_dz) + Rx)) 
    return delta_c_i1

def fcs(cs_a, teta, kads, ci, kdes, kmicro):
    #Rxs = kmicro*cs_a
    Rxs = 0
    
    #Rxs = Um*(teta*cs_a/(Km + teta*cs_a))
    cs_abs = cs_a + ((teta/ro)*kads*ci - kdes*cs_a - Rxs)*dt
    
    if cs_abs <= 0:
        cs = 0
    else:
        cs = cs_abs
    
    return cs
   
######  Ammonia  ####
# Ponding zone
def fR_nh4_1():

    return 0
#Unsaturated zone
def fR_nh4_2(C_nh4_iminus1):
    R_nit = - k_nit*C_nh4_iminus1 
    
    return R_nit

#Saturated zone
def fR_nh4_3():
    
    return 0

######  Nitrate  #####
#Ponding zone
def fR_no3_1():
#vou considerar por enquanto que nao ocorre nada
    return 0
    
#Unsaturated zone
def fR_no3_2(C_nh4_iminus1):
    R_nit = k_nit*C_nh4_iminus1 
    
    return R_nit

#Saturated zone
def fR_no3_3(C_no3_iminus1, C_o2_i, C_doc_iminus1):
#     Of = K_o2/(K_o2+C_o2_i)
#     bDOCf = (C_doc_iminus1 + bDOCd*dt)/(C_doc_iminus1 + bDOCd*dt + KbDOC)
#     
#     k2 = k_denit*Of*bDOCf

    k2 = k_denit  ###testando sem influencia de DOC e O2
    
    R_denit = - k2*C_no3_iminus1
    
    return R_denit

    
######  Oxygen  ######
#Ponding zone
def fR_o2_1():
#vou deixar vazio por enquanto    
    return 0

#Unsaturated zone
def fR_o2_2(C_o2_iminus1, C_nh4_iminus1):
    R_o2 = -k_o2*C_o2_iminus1 - k_nit*C_nh4_iminus1/2
    
    return R_o2

#Saturated zone
def fR_o2_3(C_o2_iminus1, C_nh4_iminus1):
    R_o2 = -k_o2*C_o2_iminus1 - k_nit*C_nh4_iminus1/2
    
    return R_o2

######  DOC  ######
#Ponding zone
def fR_doc_1():
#vou deixar vazio por enquanto    
    return 0

#Unsaturated zone
def fR_doc_2(C_doc_iminus1):
    R_doc = -k_doc*C_doc_iminus1 + bDOCd
    
    return R_doc 

#Saturated zone
def fR_doc_3(C_doc_iminus1):
    R_doc = -k_doc*C_doc_iminus1 + bDOCd
    
    return R_doc 


######  Plant Uptake  ######

### Ammonia
#Unsaturated zone
def fPU_nh4_2(C_nh4_2, teta_sm, root_fraction):
    PU_nh4_2 = - root_fraction * (Fm_nh4*teta_sm*C_nh4_2/(Km_nh4 + teta_sm*C_nh4_2))
    
    return PU_nh4_2

#Saturated zone
def fPU_nh4_3(C_nh4_3, teta_b, root_fraction):
    PU_nh4_3 = - root_fraction * (Fm_nh4*teta_b*C_nh4_3/(Km_nh4 + teta_b*C_nh4_3))
    
    return PU_nh4_3

### Nitrate
#Unsaturated zone
def fPU_no3_2(C_no3_2, teta_sm, root_fraction):
    PU_no3_2 = -root_fraction * (Fm_no3*teta_sm*C_no3_2/(Km_no3 + teta_sm*C_no3_2))
    
    return PU_no3_2

#Saturated zone
def fPU_no3_3(C_no3_3, teta_b, root_fraction):
    PU_no3_3 = -root_fraction * (Fm_no3*teta_b*C_no3_3/(Km_no3 + teta_b*C_no3_3))
    
    return PU_no3_3
    
### Oxygen
#Unsaturated zone
def fPU_o2_2(C_o2_2, C_o2_root, teta_sm, root_fraction):
    PU_o2_2 = -root_fraction * (lamda*(teta_sm*C_o2_root - teta_sm*C_o2_2))
    
    return PU_o2_2

#Saturated zone
def fPU_o2_3(C_o2_3, C_o2_root, teta_b, root_fraction):
    PU_o2_3 = -root_fraction * (lamda*(teta_b*C_o2_root - teta_b*C_o2_3))
    
    return PU_o2_3


    
# **4. Model routine**
def run_Kin():    
    #O2
    cp_o2 = []
    c_usz_o2 = []
    cs_usz_o2 = []
    c_sz_o2 = []
    cs_sz_o2 = []
    
    Rx_usz_o2 = []
    Rx_sz_o2 = []
    cp_a_o2 = 0
    
    c0_usz_o2 = [0]
    c0_usz_o2 = c0_usz_o2 * m_usz
    c_usz_o2.append(c0_usz_o2)
    cs_usz_o2.append(c0_usz_o2)
    c0_sz_o2 = [0]
    c0_sz_o2 = c0_sz_o2 * m_sz
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
    
    #NH4
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
    c0_usz_nh4 = c0_usz_nh4 * m_usz
    cs0_usz_nh4 = [cs_usz_nh4_a]
    cs0_usz_nh4 = cs0_usz_nh4* m_usz 
    c_usz_nh4.append(c0_usz_nh4)
    cs_usz_nh4.append(cs0_usz_nh4)
    c0_sz_nh4 = [0]
    c0_sz_nh4 = c0_sz_nh4 * m_sz
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
    
    #NO3
    cp_no3 = []
    c_usz_no3 = []
    cs_usz_no3 = []
    c_sz_no3 = []
    cs_sz_no3 = []
    
    Rx_usz_no3 = []
    Rx_sz_no3 = []
    
    cp_a_no3 = 0
    
    c0_usz_no3 = [0]
    c0_usz_no3 = c0_usz_no3 * m_usz
    c_usz_no3.append(c0_usz_no3)
    cs_usz_no3.append(c0_usz_no3)
    c0_sz_no3 = [0]
    c0_sz_no3 = c0_sz_no3 * m_sz
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
    
    #DOC
    cp_doc = []
    c_usz_doc = []
    cs_usz_doc = []
    c_sz_doc = []
    cs_sz_doc = []
    
    Rx_usz_doc = []
    Rx_sz_doc = []
    
    cp_a_doc = 0
    
    c0_usz_doc = [0]
    c0_usz_doc = c0_usz_doc * m_usz
    c_usz_doc.append(c0_usz_doc)
    cs_usz_doc.append(c0_usz_doc)
    c0_sz_doc = [0]
    c0_sz_doc = c0_sz_doc * m_sz
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
    
    
    for t in range(len(indice)-1):
        #Ponding zone
        cin_o2 = cin_o2_list[t]
        cin_nh4 = cin_nh4_list[t]
        cin_no3 = cin_no3_list[t]
        cin_doc = cin_doc_list[t]
        
        #print('t', t)
        hp = thpEND[t]
        if hp < 0.001:
            hp = 0
        else:
            hp = hp
            
        if t == 0:
            hp_a = 0
        else:
            hp_a = thpEND[t-1]
        #print('hp', hp)
        Qin_p = tQin[t]
        I1 = tQpf[t]
        Qv = tQover[t]
        Qet_1 = tQet1[t]
        I2 = tQfs[t]
        Qhc = tQhc[t]
        Qet_2 = tQet2[t]
        if t == 0:
            Qorif = 0
        else:
            Qorif = tQpipe[t] 
        Qinf_sz = tQinfsz[t]
        teta_sm_i = tteta_usz[t]
        teta_b_i = tteta_sz[t]
        if t < (len(indice)-1):
            teta_sm_iplus1 = tteta_usz[t+1]
            teta_b_iplus1 = tteta_sz[t+1]
        else:
            teta_sm_iplus1 = tteta_usz[t]
            teta_b_iplus1 = tteta_sz[t]        

        
        Rxi_p_o2 = fR_o2_1()
        Rxi_p_nh4 = fR_nh4_1() 
        Rxi_p_no3 = fR_no3_1()
        Rxi_p_doc = fR_doc_1()
        
        #if t < 20:
            #print('t: ', t, 'Rx_p_o2: ', Rxi_p_o2, 'Rx_p_nh4: ', Rxi_p_nh4, 'Rx_p_no3: ', Rxi_p_no3)   
           
        if hp == 0:
            cpi_o2 = 0
            cpi_nh4 = 0
            cpi_no3 = 0
            cpi_doc = 0
       
        else:        
            cpi_o2 = fcp(cin_o2, Qin_p, cp_a_o2, I1, Qv, Rxi_p_o2, hp, hp_a)
            cpi_nh4 = fcp(cin_nh4, Qin_p, cp_a_nh4, I1, Qv, Rxi_p_nh4, hp, hp_a)
            cpi_no3 = fcp(cin_no3, Qin_p, cp_a_no3, I1, Qv, Rxi_p_no3, hp, hp_a)
            cpi_doc = fcp(cin_doc, Qin_p, cp_a_doc, I1, Qv, Rxi_p_doc, hp, hp_a)
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
        
        #USZ
        cli_o2_list = c_usz_o2[t].copy()
        #print('1_o2', cli_o2_list)
        cli_nh4_list = c_usz_nh4[t].copy()   
        #print('1_nh4', cli_nh4_list)
        cli_no3_list = c_usz_no3[t].copy()
        cli_doc_list = c_usz_doc[t].copy()
        
        #SZ
        if hpipe > 0:
            cji_o2_list = c_sz_o2[t].copy() #copiar lista
            #cji_o2_list = c_sz_o2[t][:] #outro jeito de copiar a lista
            cji_nh4_list = c_sz_nh4[t].copy()
            cji_no3_list = c_sz_no3[t].copy()
            cji_doc_list = c_sz_doc[t].copy()

        ####   USZ   ####
        if m_usz != 0:
    
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

            for l in range(m_usz): 
                #print('l', l)
                
                cl_o2 = cli_o2_list[l] 
                clminus1_o2 = cli_o2_list[l-1]
                if l < (m_usz - 1):  
                    clplus1_o2 = cli_o2_list[l+1]
                else:
                    clplus1_o2 = 0           
    
                cl_nh4 = cli_nh4_list[l]
                clminus1_nh4 = cli_nh4_list[l-1]
                if l < (m_usz - 1):
                    clplus1_nh4 = cli_nh4_list[l+1]
                else:
                    clplus1_nh4 = 0
    
                cl_no3 = cli_no3_list[l]
                clminus1_no3 = cli_no3_list[l-1]
                if l < (m_usz - 1):
                    clplus1_no3 = cli_no3_list[l+1]
                else:
                    clplus1_no3 = 0
                    
                len(cli_doc_list)
                cl_doc = cli_doc_list[l]
                clminus1_doc = cli_doc_list[l-1]

                if l < (m_usz - 1):
                    clplus1_doc = cli_doc_list[l+1]
                else:
                    clplus1_doc = 0 
    
                            
                cs_o2 = 0
                
                #since we have added an initial value of cs = 0 in the cs_usz_nh4 list, when calling the index equal = 't' we are actually calling the value corresponding to t-1 
                cs_nh4 = cs_usz_nh4[t][l]
                cs_nh4_iplus1 = fcs(cs_nh4, teta_sm_i, kads_nh4, cl_nh4, kdes_nh4, k_nh4_mb)
                if cs_nh4_iplus1 < 0.00000000000001:
                    cs_nh4_iplus1 = 0
                    
#                 sor = (teta_sm_i/ro)*kads_nh4*cl_nh4*dt
#                 des =  kdes_nh4*cs_usz_nh4_a*dt
                
#                 print('t: ', t, ', l: ', l , ', cs_nh4_a: ', cs_usz_nh4_a, 'teta_sm_i: ', teta_sm_i, 'kads_nh4: ', kads_nh4, 'cl_nh4: ', cl_nh4, 'kdes_nh4: ', kdes_nh4, 'k_nh4_mb: ', k_nh4_mb)
#                 print('cs_nh4: ', cs_nh4)
#                 input()
                
#                 logger.debug(f' t: {t} \t l: {l} \t sor: {sor} \t des: {des} \t cs_nh4: {cs_nh4}')
#                 logger.debug('----SOIL----')
#                 logger.debug('cs_nh4_a: ', cs_usz_nh4_a)
#                 logger.debug('cs_nh4_a: ', cs_usz_nh4_a, 'teta_sm_i: ', teta_sm_i, 'kads_nh4: ', kads_nh4, 'cl_nh4: ', cl_nh4, 'kdes_nh4: ', kdes_nh4, 'k_nh4_mb: ', k_nh4_mb)
#                 logger.debug('cs_nh4: ', cs_nh4)
#                 logger.debug('----WATER----')
                csi_usz_nh4.append(cs_nh4_iplus1)

                
                cs_no3 = 0
    
                cs_usz_doc_a = cs_usz_doc[t][l]
                cs_doc = fcs(cs_usz_doc_a, teta_sm_i, kads_doc, cl_doc, kdes_doc, k_doc_mb)
                csi_usz_doc.append(cs_doc)
                if cs_doc < 0.00000000000001:
                    cs_doc = 0                
                                        
                UF_usz = []
                #Qet_1 = 0                        
                alfa,beta = falfa_beta_usz(l)
                UFi_usz = fUF_usz(alfa, beta, I1, Qet_1, I2, Qhc, Qorif, Qinf_sz, teta_sm_i)
#                 print(t, '  ', l)
#                 print('UFi_usz: ', UFi_usz, ', alfa: ', alfa, ', beta: ', beta,  ', I1: ', I1, ', Qet_1:', Qet_1, ', I2: ', I2, ', Qhc: ', Qhc, ', Qorif: ', Qorif, ', Qinf_sz: ', Qinf_sz)            
#                 input()
                
                UF_usz.append(UFi_usz)
                
                Rxi_2_o2 = fR_o2_2(cl_o2, cl_nh4) + fPU_o2_2(cl_o2, c_o2_root, teta_sm_i, root_fraction)
                Rxi_2_nh4 = fR_nh4_2(cl_nh4) + fPU_nh4_2(cl_nh4, teta_sm_i, root_fraction)
                Rxi_2_no3 = fR_no3_2(cl_nh4) + fPU_no3_2(cl_no3, teta_sm_i, root_fraction)
                Rxi_2_doc = fR_doc_2(cl_doc)
                
                Rxl_o2.append(Rxi_2_o2*(1/teta_sm_iplus1)*dt)
                Rxl_nh4.append(Rxi_2_nh4*(1/teta_sm_iplus1)*dt)
                Rxl_no3.append(Rxi_2_no3*(1/teta_sm_iplus1)*dt)
                Rxl_doc.append(Rxi_2_doc*(1/teta_sm_iplus1)*dt)

#                 Rxl_o2.append(Rxi_2_o2*(dt/teta_sm_i))
#                 Rxl_nh4.append(Rxi_2_nh4*(dt/teta_sm_i))
#                 Rxl_no3.append(Rxi_2_no3*(dt/teta_sm_i))
#                 Rxl_doc.append(Rxi_2_doc*(dt/teta_sm_i))
#                 print('t: ', t, ', l: ', l, ', Rx: ', Rxi_2_nh4*(dt/teta_sm_i))
#                 input()
                #if t < 20:
                    #print('t: ', t, 'l: ', l, 'Rx_usz_o2: ', Rxi_2_o2, 'Rx_usz_nh4: ', Rxi_2_nh4, 'Rx_usz_no3: ', Rxi_2_no3)   
                
    
    ### Oxygen
                Peusz_o2 = fPeusz(UFi_usz, D_o2)
                #print(Peusz_o2)
                
                if Peusz_o2 <= 2:
                    if l == 0: #first cell
                        dc_o2 = clplus1_o2 - 2*cl_o2 + cpi_o2
                        dc_dz_o2 = (clplus1_o2 - cpi_o2)/(2*dz)
   
        
                    elif l == (m_usz - 1): #last cell
                        dc_o2 = cl_o2 - 2*cl_o2 + clminus1_o2
                        dc_dz_o2 = (cl_o2 - clminus1_o2)/dz
                        
                    else:                 
                        dc_o2 = clplus1_o2 - 2*cl_o2 + clminus1_o2
                        dc_dz_o2 = (clplus1_o2 - clminus1_o2)/(2*dz)               
                
                else: #Peusz > 2
                    if l == 0: #first cell
                        dc_o2 = clplus1_o2 - 2*cl_o2 + cpi_o2
                        dc_dz_o2 = (cl_o2 - cpi_o2)/dz
        
                    elif l == (m_usz - 1): #last cell
                        dc_o2 = cl_o2 - 2*cl_o2 + clminus1_o2
                        dc_dz_o2 = (cl_o2 - clminus1_o2)/dz
                        
                    else:
                        dc_o2 = clplus1_o2 - 2*cl_o2 + clminus1_o2
                        dc_dz_o2 = (cl_o2 - clminus1_o2)/dz            
                
                delta_c_o2 = ftransp(teta_sm_i, teta_sm_iplus1, cl_o2, cs_o2, dc_o2, dc_dz_o2, 0, 0, D_o2, UFi_usz, Rxi_2_o2)
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
                Peusz_nh4 = fPeusz(UFi_usz, D_nh4)
                #print('t:', t, 'l:', l, 'Peusz_nh4:', Peusz_nh4)
                
                if Peusz_nh4 <= 2:
                    if l == 0: #first cell
                        dc_nh4 = clplus1_nh4 - 2*cl_nh4 + cpi_nh4
                        dc_dz_nh4 = (clplus1_nh4 - cpi_nh4)/(2*dz)
        
                    elif l == (m_usz - 1): #last cell
                        dc_nh4 = cl_nh4 - 2*cl_nh4 + clminus1_nh4
                        dc_dz_nh4 = (cl_nh4 - clminus1_nh4)/dz
                        
                    else:
                        dc_nh4 = clplus1_nh4 - 2*cl_nh4 + clminus1_nh4
                        dc_dz_nh4 = (clplus1_nh4 - clminus1_nh4)/(2*dz)               
                
                else: #Peusz > 2
                    if l == 0: #first cell
                        dc_nh4 = clplus1_nh4 - 2*cl_nh4 + cpi_nh4
                        dc_dz_nh4 = (cl_nh4 - cpi_nh4)/dz
        
                    elif l == (m_usz - 1): #last cell
                        dc_nh4 = cl_nh4 - 2*cl_nh4 + clminus1_nh4
                        dc_dz_nh4 = (cl_nh4 - clminus1_nh4)/dz
                        
                    else:
                        dc_nh4 = clplus1_nh4 - 2*cl_nh4 + clminus1_nh4
                        dc_dz_nh4 = (cl_nh4 - clminus1_nh4)/dz            
                
                delta_c_nh4 = ftransp(teta_sm_i, teta_sm_iplus1, cl_nh4, cs_nh4, dc_nh4, dc_dz_nh4, kads_nh4, kdes_nh4, D_nh4, UFi_usz, Rxi_2_nh4)
#                 print('t: ', t, ', l: ', l)
#                 print('clplus1_nh4: ',  clplus1_nh4, ', cl_nh4: ' ,cl_nh4 , ', cpi_nh4:' , cpi_nh4, ', clminus1_nh4: ' , clminus1_nh4)
#                 print('dt:', dt, 'teta_sm_i:', teta_sm_i, 'teta_sm_iplus1:', teta_sm_iplus1, 'cl_nh4:', cl_nh4, 'cs_nh4:', cs_nh4, 'dc_nh4:', dc_nh4, 'dc_dz_nh4:', dc_dz_nh4, 'kads_nh4:', kads_nh4, 'kdes_nh4:', kdes_nh4, 'D_nh4:', D_nh4, 'UFi_usz:', UFi_usz, 'Rxi_2_nh4:', Rxi_2_nh4)
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
                Peusz_no3 = fPeusz(UFi_usz, D_no3)
                
                if Peusz_no3 <= 2:
                    if l == 0: #first cell
                        dc_no3 = clplus1_no3 - 2*cl_no3 + cpi_no3
                        dc_dz_no3 = (clplus1_no3 - cpi_no3)/(2*dz)
        
                    elif l == (m_usz - 1): #last cell
                        dc_no3 = cl_no3 - 2*cl_no3 + clminus1_no3
                        dc_dz_no3 = (cl_no3 - clminus1_no3)/dz
                        
                    else:
                        dc_no3 = clplus1_no3 - 2*cl_no3 + clminus1_no3
                        dc_dz_no3 = (clplus1_no3 - clminus1_no3)/(2*dz)               
                
                else: #Peusz > 2
                    if l == 0: #first cell
                        dc_no3 = clplus1_no3 - 2*cl_no3 + cpi_no3
                        dc_dz_no3 = (cl_no3 - cpi_no3)/dz
        
                    elif l == (m_usz - 1): #last cell
                        dc_no3 = cl_no3 - 2*cl_no3 + clminus1_no3
                        dc_dz_no3 = (cl_no3 - clminus1_no3)/dz
                        
                    else:
                        dc_no3 = clplus1_no3 - 2*cl_no3 + clminus1_no3
                        dc_dz_no3 = (cl_no3 - clminus1_no3)/dz            
                
                delta_c_no3 = ftransp(teta_sm_i, teta_sm_iplus1, cl_no3, cs_no3, dc_no3, dc_dz_no3, 0, 0, D_no3, UFi_usz, Rxi_2_no3)
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
                Peusz_doc = fPeusz(UFi_usz, D_doc)
                
                if Peusz_doc <= 2:
                    if l == 0: #first cell
                        dc_doc = clplus1_doc - 2*cl_doc + cpi_doc
                        dc_dz_doc = (clplus1_doc - cpi_doc)/(2*dz)
        
                    elif l == (m_usz - 1): #last cell
                        dc_doc = cl_doc - 2*cl_doc + clminus1_doc
                        dc_dz_doc = (cl_doc - clminus1_doc)/dz
                        
                    else:
                        dc_doc = clplus1_doc - 2*cl_doc + clminus1_doc
                        dc_dz_doc = (clplus1_doc - clminus1_doc)/(2*dz)               
                
                else: #Peusz > 2
                    if l == 0: #first cell
                        dc_doc = clplus1_doc - 2*cl_doc + cpi_doc
                        dc_dz_doc = (cl_doc - cpi_doc)/dz
        
                    elif l == (m_usz - 1): #last cell
                        dc_doc = cl_doc - 2*cl_doc + clminus1_doc
                        dc_dz_doc = (cl_doc - clminus1_doc)/dz
                        
                    else:
                        dc_doc = clplus1_doc - 2*cl_doc + clminus1_doc
                        dc_dz_doc = (cl_doc - clminus1_doc)/dz            
                
                delta_c_doc = ftransp(teta_sm_i, teta_sm_iplus1, cl_doc, cs_doc, dc_doc, dc_dz_doc, kads_doc, kdes_doc, D_doc, UFi_usz, Rxi_2_doc)
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

        for j in range(m_sz):
#                 if teta_b_iplus1 > 0:
#                     den_teta = 1/teta_b_iplus1
#                 else:
#                     den_teta = 0

            cj_o2 = cji_o2_list[j]
            cjminus1_o2 = cji_o2_list[j-1]
            if j < (m_sz - 1):
                cjplus1_o2 = cji_o2_list[j+1]
            else:
                cjplus1_o2 = 0
            cmminus1_o2 = cli_o2_list[m_usz-1]

            cj_nh4 = cji_nh4_list[j]
            cjminus1_nh4 = cji_nh4_list[j-1]
            if j < (m_sz - 1):
                cjplus1_nh4 = cji_nh4_list[j+1]
            else:
                cjplus1_nh4 = 0
            cmminus1_nh4 = cli_nh4_list[m_usz-1]

            cj_no3 = cji_no3_list[j]
            cjminus1_no3 = cji_no3_list[j-1]
            if j < (m_sz - 1):
                cjplus1_no3 = cji_no3_list[j+1]
            else:
                cjplus1_no3 = 0
            cmminus1_no3 = cli_no3_list[m_usz-1]

            cj_doc = cji_doc_list[j]
            cjminus1_doc = cji_doc_list[j-1]
            if j < (m_sz - 1):
                cjplus1_doc = cji_doc_list[j+1]
            else:
                cjplus1_doc = 0
            cmminus1_doc = cli_doc_list[m_usz-1]

            cs_o2 = 0

            cs_sz_nh4_a = cs_sz_nh4[t][j]
            cs_nh4 = fcs(cs_sz_nh4_a, teta_b_i, kads2_nh4, cj_nh4, kdes2_nh4, k_nh4_mb)
            csi_sz_nh4.append(cs_nh4)
            if cs_nh4 < 0.00000000000001:
                cs_nh4 = 0

#                 print('====>> t: ', t, ', j: ', j)
#                 print('---- SOIL ----')
#                 print('cs_nh4_a: ', cs_sz_nh4_a, 'teta_b_i: ', teta_b_i, 'kads_nh4: ', kads_nh4, 'cj_nh4: ', cj_nh4, 'kdes_nh4: ', kdes_nh4, 'k_nh4_mb: ', k_nh4_mb)
#                 print('cs_nh4: ', cs_nh4)
#                 print('---- WATER ----')
#                 print('cmminus1_nh4: ', cmminus1_nh4)                
            cs_no3 = 0

            cs_sz_doc_a = cs_sz_doc[t][j]
            cs_doc = fcs(cs_sz_doc_a, teta_b_i, kads2_doc, cj_doc, kdes2_doc, k_doc_mb)
            csi_sz_doc.append(cs_doc)

            UF_sz = []
            alfa2,beta2 = falfa_beta_sz(j)
            #Qet_2 = 0
            UFi_sz = fUF_sz(alfa2, beta2, I2, Qhc, Qet_2, Qorif, Qinf_sz, teta_b_i)
#                 print('I2, Qhc, Qet_2, Qorif, Qinf_sz:', I2, Qhc, Qet_2, Qorif, Qinf_sz)
#                 print('UF: ', UFi_sz, ', alfa2: ', alfa2, ', beta2: ' ,beta2)

            UF_sz.append(UFi_sz)

            Rxi_3_o2 = fR_o2_3(cj_o2, cj_nh4) + fPU_o2_3(cj_o2, c_o2_root, teta_b_i, root_fraction)
            Rxi_3_nh4 = fR_nh4_3() + fPU_nh4_3(cj_nh4, teta_b_i, root_fraction)
            Rxi_3_no3 = fR_no3_3(cj_no3, cj_o2, cj_doc) + fPU_no3_3(cj_no3, teta_b_i, root_fraction)
            Rxi_3_doc = fR_doc_3(cj_doc)

            Rxj_o2.append(Rxi_3_o2*(1/teta_b_iplus1)*dt)
            Rxj_nh4.append(Rxi_3_nh4*(1/teta_b_iplus1)*dt)
            Rxj_no3.append(Rxi_3_no3*(1/teta_b_iplus1)*dt)
            Rxj_doc.append(Rxi_3_doc*(1/teta_b_iplus1)*dt)
            #if t < 20:
                #print('t: ', t, 'j: ', j, 'Rx_sz_o2: ', Rxi_3_o2, 'Rx_sz_nh4: ', Rxi_3_nh4, 'Rx_sz_no3: ', Rxi_3_no3)

### Oxygen
            Pesz_o2 = fPesz(UFi_sz, D_o2)

            if m_usz < (n-1):
                if Pesz_o2 <= 2:
                    if j == 0: #first cell
                        dc_o2 = cjplus1_o2 - 2*cj_o2 + cmminus1_o2
                        dc_dz_o2 = (cjplus1_o2 - cmminus1_o2)/(2*dz)

                    elif j == (m_sz - 1): #last cell
                        dc_o2 = cj_o2 - 2*cj_o2 + cjminus1_o2
                        dc_dz_o2 = (cj_o2 - cjminus1_o2)/dz

                    else:
                        dc_o2 = cjplus1_o2 - 2*cj_o2 + cjminus1_o2
                        dc_dz_o2 = (cjplus1_o2 - cjminus1_o2)/(2*dz)

                else: #Pesz > 2
                    if j == 0: #first cell
                        dc_o2 = cjplus1_o2 - 2*cj_o2 + cmminus1_o2
                        dc_dz_o2 = (cj_o2 - cmminus1_o2)/dz

                    elif j == (m_sz - 1): #last cell
                        dc_o2 = cj_o2 - 2*cj_o2 + cjminus1_o2
                        dc_dz_o2 = (cj_o2 - cjminus1_o2)/dz

                    else:
                        dc_o2 = cjplus1_o2 - 2*cj_o2 + cjminus1_o2
                        dc_dz_o2 = (cj_o2 - cjminus1_o2)/dz


            if m_usz == (n - 1):
                dc_o2 = cj_o2 - 2*cj_o2 + cjminus1_o2
                dc_dz_o2 = (cj_o2 - cjminus1_o2)/(2*dz)

            if m_usz == 0:
                if Pesz_o2 <= 2:
                    if j == 0: #first cell
                        dc_o2 = cjplus1_o2 - 2*cj_o2 + cpi_o2
                        dc_dz_o2 = (cjplus1_o2 - cpi_o2)/(2*dz)

                    elif j == (m_sz - 1): #last cell
                        dc_o2 = cj_o2 - 2*cj_o2 + cjminus1_o2
                        dc_dz_o2 = (cj_o2 - cjminus1_o2)/dz

                    else:
                        dc_o2 = cjplus1_o2 - 2*cj_o2 + cjminus1_o2
                        dc_dz_o2 = (cjplus1_o2 - cjminus1_o2)/(2*dz)

                else: #Pusz > 2
                    if j == 0: #first cell
                        dc_o2 = cjplus1_o2 - 2*cj_o2 + cpi_o2
                        dc_dz_o2 = (cj_o2 - cpi_o2)/dz

                    elif j == (m_sz - 1): #last cell
                        dc_o2 = cj_o2 - 2*cj_o2 + cjminus1_o2
                        dc_dz_o2 = (cj_o2 - cjminus1_o2)/dz

                    else:
                        dc_o2 = cjplus1_o2 - 2*cj_o2 + cjminus1_o2
                        dc_dz_o2 = (cj_o2 - cjminus1_o2)/dz

            delta_c_o2 = ftransp(teta_b_i, teta_b_iplus1, cj_o2, cs_o2, dc_o2, dc_dz_o2, 0, 0, D_o2, UFi_sz, Rxi_3_o2)
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
            Pesz_nh4 = fPesz(UFi_sz, D_nh4)

            if m_usz < (n-1):
                if Pesz_nh4 <= 2:
                    if j == 0: #first cell
                        dc_nh4 = cjplus1_nh4 - 2*cj_nh4 + cmminus1_nh4
                        dc_dz_nh4 = (cjplus1_nh4 - cmminus1_nh4)/(2*dz)

                    elif j == (m_sz - 1): #last cell
                        dc_nh4 = cj_nh4 - 2*cj_nh4 + cjminus1_nh4
                        dc_dz_nh4 = (cj_nh4 - cjminus1_nh4)/dz

                    else:
                        dc_nh4 = cjplus1_nh4 - 2*cj_nh4 + cjminus1_nh4
                        dc_dz_nh4 = (cjplus1_nh4 - cjminus1_nh4)/(2*dz)

                else: #Pesz > 2
                    if j == 0: #first cell
                        dc_nh4 = cjplus1_nh4 - 2*cj_nh4 + cmminus1_nh4
                        dc_dz_nh4 = (cj_nh4 - cmminus1_nh4)/dz

                    elif j == (m_sz - 1): #last cell
                        dc_nh4 = cj_nh4 - 2*cj_nh4 + cjminus1_nh4
                        dc_dz_nh4 = (cj_nh4 - cjminus1_nh4)/dz

                    else:
                        dc_nh4 = cjplus1_nh4 - 2*cj_nh4 + cjminus1_nh4
                        dc_dz_nh4 = (cj_nh4 - cjminus1_nh4)/dz


            if m_usz == (n - 1):
                dc_nh4 = cj_nh4 - 2*cj_nh4 + cjminus1_nh4
                dc_dz_nh4 = (cj_nh4 - cjminus1_nh4)/(2*dz)

            if m_usz == 0:
                if Pesz_nh4 <= 2:
                    if j == 0: #first cell
                        dc_nh4 = cjplus1_nh4 - 2*cj_nh4 + cpi_nh4
                        dc_dz_nh4 = (cjplus1_nh4 - cpi_nh4)/(2*dz)

                    elif j == (m_sz -1): #last cell
                        dc_nh4 = cj_nh4 - 2*cj_nh4 + cjminus1_nh4
                        dc_dz_nh4 = (cj_nh4 - cjminus1_nh4)/dz

                    else:
                        dc_nh4 = cjplus1_nh4 - 2*cj_nh4 + cjminus1_nh4
                        dc_dz_nh4 = (cjplus1_nh4 - cjminus1_nh4)/(2*dz)

                else: #Pusz > 2
                    if j == 0: #first cell
                        dc_nh4 = cjplus1_nh4 - 2*cj_nh4 + cpi_nh4
                        dc_dz_nh4 = (cj_nh4 - cpi_nh4)/dz

                    elif j == (m_sz - 1): #last cell
                        dc_nh4 = cj_nh4 - 2*cj_nh4 + cjminus1_nh4
                        dc_dz_nh4 = (cj_nh4 - cjminus1_nh4)/dz

                    else:
                        dc_nh4 = cjplus1_nh4 - 2*cj_nh4 + cjminus1_nh4
                        dc_dz_nh4 = (cj_nh4 - cjminus1_nh4)/dz

            delta_c_nh4 = ftransp(teta_b_i, teta_b_iplus1, cj_nh4, cs_nh4, dc_nh4, dc_dz_nh4, kads2_nh4, kdes2_nh4, D_nh4, UFi_sz, Rxi_3_nh4)
#                 print('teta_b_i: ', teta_b_i, ', teta_b_iplus1: ', teta_b_iplus1, ', cj_nh4: ', cj_nh4, ', cs_nh4: ', cs_nh4, ', dc_nh4: ', dc_nh4, ', dc_dz_nh4: ', dc_dz_nh4, ', kads2_nh4: ', kads2_nh4, ', kdes2_nh4: ', kdes2_nh4, ', D_nh4: ', D_nh4, ', UFi_sz: ', UFi_sz, ', Rxi_3_nh4: ', Rxi_3_nh4)
#                 print('delta_c_nh4: ', delta_c_nh4)

            if teta_b_iplus1 > 0:
                ci1_nh4 = cj_nh4 + delta_c_nh4
                #print('calculo final - cj_nh4, teta_b_i, teta_b_iplus1, delta_c_nh4: ', cj_nh4, ', ', teta_b_i, ', ', teta_b_iplus1, ', ', delta_c_nh4)
            else:
                ci1_nh4 = 0
            #print('ci1_nh4: ', ci1_nh4)

            if ci1_nh4 <= 0.0000000000000001:
                ci1_nh4 = 0
            else:
                ci1_nh4 = ci1_nh4

            cj_i1_nh4.append(ci1_nh4)

### Nitrate
            Pesz_no3 = fPesz(UFi_sz, D_no3)

            if m_usz < (n-1):
                if Pesz_no3 <= 2:
                    if j == 0: #first cell
                        dc_no3 = cjplus1_no3 - 2*cj_no3 + cmminus1_no3
                        dc_dz_no3 = (cjplus1_no3 - cmminus1_no3)/(2*dz)

                    elif j == (m_sz - 1): #last cell
                        dc_no3 = cj_no3 - 2*cj_no3 + cjminus1_no3
                        dc_dz_no3 = (cj_no3 - cjminus1_no3)/dz

                    else:
                        dc_no3 = cjplus1_no3 - 2*cj_no3 + cjminus1_no3
                        dc_dz_no3 = (cjplus1_no3 - cjminus1_no3)/(2*dz)

                else: #Pesz > 2
                    if j == 0: #first cell
                        dc_no3 = cjplus1_no3 - 2*cj_no3 + cmminus1_no3
                        dc_dz_no3 = (cj_no3 - cmminus1_no3)/dz

                    elif j == (m_sz - 1): #last cell
                        dc_no3 = cj_no3 - 2*cj_no3 + cjminus1_no3
                        dc_dz_no3 = (cj_no3 - cjminus1_no3)/dz

                    else:
                        dc_no3 = cjplus1_no3 - 2*cj_no3 + cjminus1_no3
                        dc_dz_no3 = (cj_no3 - cjminus1_no3)/dz


            if m_usz == (n - 1):
                dc_no3 = cj_no3 - 2*cj_no3 + cjminus1_no3
                dc_dz_no3 = (cj_no3 - cjminus1_no3)/(2*dz)

            if m_usz == 0:
                if Pesz_no3 <= 2:
                    if j == 0: #first cell
                        dc_no3 = cjplus1_no3 - 2*cj_no3 + cpi_no3
                        dc_dz_no3 = (cjplus1_no3 - cpi_no3)/(2*dz)

                    elif j == (m_sz - 1): #last cell
                        dc_no3 = cj_no3 - 2*cj_no3 + cjminus1_no3
                        dc_dz_no3 = (cj_no3 - cjminus1_no3)/dz

                    else:
                        dc_no3 = cjplus1_no3 - 2*cj_no3 + cjminus1_no3
                        dc_dz_no3 = (cjplus1_no3 - cjminus1_no3)/(2*dz)

                else: #Pusz > 2
                    if j == 0: #first cell
                        dc_no3 = cjplus1_no3 - 2*cj_no3 + cpi_no3
                        dc_dz_no3 = (cj_no3 - cpi_no3)/dz

                    elif j == (m_sz - 1): #last cell
                        dc_no3 = cj_no3 - 2*cj_no3 + cjminus1_no3
                        dc_dz_no3 = (cj_no3 - cjminus1_no3)/dz

                    else:
                        dc_no3 = cjplus1_no3 - 2*cj_no3 + cjminus1_no3
                        dc_dz_no3 = (cj_no3 - cjminus1_no3)/dz

            delta_c_no3 = ftransp(teta_b_i, teta_b_iplus1, cj_no3, cs_no3, dc_no3, dc_dz_no3, 0, 0, D_no3, UFi_sz, Rxi_3_no3)
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
            Pesz_doc = fPesz(UFi_sz, D_doc)

            if m_usz < (n-1):
                if Pesz_doc <= 2:
                    if j == 0: #first cell
                        dc_doc = cjplus1_doc - 2*cj_doc + cmminus1_doc
                        dc_dz_doc = (cjplus1_doc - cmminus1_doc)/(2*dz)

                    elif j == (m_sz - 1): #last cell
                        dc_doc = cj_doc - 2*cj_doc + cjminus1_doc
                        dc_dz_doc = (cj_doc - cjminus1_doc)/dz

                    else:
                        dc_doc = cjplus1_doc - 2*cj_doc + cjminus1_doc
                        dc_dz_doc = (cjplus1_doc - cjminus1_doc)/(2*dz)

                else: #Pesz > 2
                    if j == 0: #first cell
                        dc_doc = cjplus1_doc - 2*cj_doc + cmminus1_doc
                        dc_dz_doc = (cj_doc - cmminus1_doc)/dz

                    elif j == (m_sz - 1): #last cell
                        dc_doc = cj_doc - 2*cj_doc + cjminus1_doc
                        dc_dz_doc = (cj_doc - cjminus1_doc)/dz

                    else:
                        dc_doc = cjplus1_doc - 2*cj_doc + cjminus1_doc
                        dc_dz_doc = (cj_doc - cjminus1_doc)/dz


            if m_usz == (n - 1):
                dc_doc = cj_doc - 2*cj_doc + cjminus1_doc
                dc_dz_doc = (cj_doc - cjminus1_doc)/(2*dz)

            if m_usz == 0:
                if Pesz_doc <= 2:
                    if j == 0: #first cell
                        dc_doc = cjplus1_doc - 2*cj_doc + cpi_doc
                        dc_dz_doc = (cjplus1_doc - cpi_doc)/(2*dz)

                    elif j == (m_sz - 1): #last cell
                        dc_doc = cj_doc - 2*cj_doc + cjminus1_doc
                        dc_dz_doc = (cj_doc - cjminus1_doc)/dz

                    else:
                        dc_doc = cjplus1_doc - 2*cj_doc + cjminus1_doc
                        dc_dz_doc = (cjplus1_doc - cjminus1_doc)/(2*dz)

                else: #Pusz > 2
                    if j == 0: #first cell
                        dc_doc = cjplus1_doc - 2*cj_doc + cpi_doc
                        dc_dz_doc = (cj_doc - cpi_doc)/dz

                    elif j == (m_sz - 1): #last cell
                        dc_doc = cj_doc - 2*cj_doc + cjminus1_doc
                        dc_dz_doc = (cj_doc - cjminus1_doc)/dz

                    else:
                        dc_doc = cjplus1_doc - 2*cj_doc + cjminus1_doc
                        dc_dz_doc = (cj_doc - cjminus1_doc)/dz

            delta_c_doc = ftransp(teta_b_i, teta_b_iplus1, cj_doc, cs_doc, dc_doc, dc_dz_doc, kads2_doc, kdes2_doc, D_doc, UFi_sz, Rxi_3_doc)
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
            Mstor_o2_ast = sum(c_usz_o2[t])*Ab*tteta_usz[t]*thusz[t]*1000/m_usz + sum(c_sz_o2[t])*Ab*tteta_sz[t]*thsz[t]*1000/m_sz
            Mstor_o2_ast_list.append(Mstor_o2_ast)
            
            Msoil_o2_a = 0
            Msoil_o2 = 0
            Msoil_o2_acum.append(0)
            
            MRx_o2 = - (sum(Rx_usz_o2[t])*Ab*tteta_usz[t]*thusz[t]*1000/m_usz + sum(Rx_sz_o2[t])*Ab*tteta_sz[t]*thsz[t]*1000/m_sz)
            if t == 0:
                MRx_o2_acum.append(MRx_o2)
            else:
                MRx_o2_acum.append(MRx_o2 + MRx_o2_acum[-1])
    
    
            Min_o2 = tQin[t]*cin_o2_list[t]*dt*1000
            if t == 0:
                Min_o2_acum.append(Min_o2)
            else:
                Min_o2_acum.append(Min_o2 + Min_o2_acum[-1])
    
            Mover_o2 = tQover[t]*cp_o2[t]*dt*1000
            if t == 0:
                Mover_o2_acum.append(Mover_o2)
            else:
                Mover_o2_acum.append(Mover_o2 + Mover_o2_acum[-1])
             
            Mpipe_o2 = tQpipe[t]*c_sz_o2[t][m_sz-1]*dt*1000 
            if t == 0:
                Mpipe_o2_acum.append(Mpipe_o2)
            else:
                Mpipe_o2_acum.append(Mpipe_o2 + Mpipe_o2_acum[-1])
             
            Minfsz_o2 = tQinfsz[t]*c_sz_o2[t][m_sz-1]*dt*1000 
            if t == 0:
                Minfsz_o2_acum.append(Minfsz_o2)
            else:
                Minfsz_o2_acum.append(Minfsz_o2 + Minfsz_o2_acum[-1])
             
            Met_o2 = tQet[t]*cl_i1_o2[0]*dt*1000
            if t == 0:
                Met_o2_acum.append(Met_o2)
            else:
                Met_o2_acum.append(Met_o2 + Met_o2_acum[-1])
             
            Mpz_o2 = thpEND[t]*Ab*1000*cp_o2[t]
            Mpz_o2_list.append(Mpz_o2)
             
        
            Mstor_o2_mb = Min_o2_acum[-1] - Mpz_o2_list[-1] - Mover_o2_acum[-1] - Mpipe_o2_acum[-1] - Minfsz_o2_acum[-1] - Met_o2_acum[-1] - Msoil_o2_acum[-1] - MRx_o2_acum[-1]
            Mstor_o2_mb_list.append(Mstor_o2_mb)        
            
            delta_ast_o2_usz = (Mstor_o2_mb - Mstor_o2_ast)/(Ab*tteta_usz[t]*thusz[t]*1000/m_usz)
            cl_i1_o2[1] = cl_i1_o2[1] + delta_ast_o2_usz
            if cl_i1_o2[1] > 0:
                cl_i1_o2[1] = cl_i1_o2[1]
            else:
                cl_i1_o2[1] = 0
                   
            ### Amonia
            
            Mstor_nh4_ast = sum(c_usz_nh4[t])*Ab*tteta_usz[t]*thusz[t]*1000/m_usz + sum(c_sz_nh4[t])*Ab*tteta_sz[t]*thsz[t]*1000/m_sz
            Mstor_nh4_ast_list.append(Mstor_nh4_ast)
            
            Msoil_nh4_a = cs_usz_nh4_a*ro*Ab*thusz[t]*1000 + cs_sz_nh4_a*ro*Ab*thsz[t]*1000
            Msoil_nh4 = sum(cs_usz_nh4[t])*ro*Ab*thusz[t]*1000/m_usz + sum(cs_sz_nh4[t])*ro*Ab*thsz[t]*1000/m_sz
            Msoil_nh4_acum.append(Msoil_nh4 - Msoil_nh4_a)
            
       
            MRx_nh4 = - (sum(Rx_usz_nh4[t])*Ab*tteta_usz[t]*thusz[t]*1000/m_usz + sum(Rx_sz_nh4[t])*Ab*tteta_sz[t]*thsz[t]*1000/m_sz)
            if t == 0:
                MRx_nh4_acum.append(MRx_nh4)
            else:
                MRx_nh4_acum.append(MRx_nh4 + MRx_nh4_acum[-1])
    
    
            Min_nh4 = tQin[t]*cin_nh4_list[t]*dt*1000
            if t == 0:
                Min_nh4_acum.append(Min_nh4)
            else:
                Min_nh4_acum.append(Min_nh4 + Min_nh4_acum[-1])
    
            Mover_nh4 = tQover[t]*cp_nh4[t]*dt*1000
            if t == 0:
                Mover_nh4_acum.append(Mover_nh4)
            else:
                Mover_nh4_acum.append(Mover_nh4 + Mover_nh4_acum[-1])
             
            Mpipe_nh4 = tQpipe[t]*c_sz_nh4[t][m_sz-1]*dt*1000 
            if t == 0:
                Mpipe_nh4_acum.append(Mpipe_nh4)
            else:
                Mpipe_nh4_acum.append(Mpipe_nh4 + Mpipe_nh4_acum[-1])
             
            Minfsz_nh4 = tQinfsz[t]*c_sz_nh4[t][m_sz-1]*dt*1000 
            if t == 0:
                Minfsz_nh4_acum.append(Minfsz_nh4)
            else:
                Minfsz_nh4_acum.append(Minfsz_nh4 + Minfsz_nh4_acum[-1])
             
            Met_nh4 = tQet[t]*cl_i1_nh4[0]*dt*1000
            if t == 0:
                Met_nh4_acum.append(Met_nh4)
            else:
                Met_nh4_acum.append(Met_nh4 + Met_nh4_acum[-1])
             
            Mpz_nh4 = thpEND[t]*Ab*1000*cp_nh4[t]
            Mpz_nh4_list.append(Mpz_nh4)
             
        
            Mstor_nh4_mb = Min_nh4_acum[-1] - Mpz_nh4_list[-1] - Mover_nh4_acum[-1] - Mpipe_nh4_acum[-1] - Minfsz_nh4_acum[-1] - Met_nh4_acum[-1] - Msoil_nh4_acum[-1] - MRx_nh4_acum[-1]
            Mstor_nh4_mb_list.append(Mstor_nh4_mb)        
             
    
    #         delta_ast_nh4_usz = 0
    #         delta_ast_nh4_sz = 0
            
    #         delta_ast_nh4_usz = (Mstor_nh4_mb - Mstor_nh4_ast)/(n*(Ab*tteta_usz[t]*thusz[t]*1000/m_usz))
    #         delta_ast_nh4_sz = (Mstor_nh4_mb - Mstor_nh4_ast)/(n*(Ab*tteta_sz[t]*thsz[t]*1000/m_sz))
            
            delta_ast_nh4_usz = (Mstor_nh4_mb - Mstor_nh4_ast)/(Ab*tteta_usz[t]*thusz[t]*1000/m_usz)
            cl_i1_nh4[1] = cl_i1_nh4[1] + delta_ast_nh4_usz
            if cl_i1_nh4[1] > 0:
                cl_i1_nh4[1] = cl_i1_nh4[1]
            else:
                cl_i1_nh4[1] = 0
    #         sum_cusz_nh4 = sum(c_usz_nh4[t])
    #         sum_csz_nh4 = sum(c_sz_nh4[t])        
    #         print('t: ', t, 'Mstor_trace: ', Mstor_trace, ', Mstor_mb: ', Mstor_mb, ', dif: ', (Mstor_mb - Mstor_trace))  
    #         print('t: ', t, ', delta_usz: ', delta_trace_nh4_usz, ', delta_sz: ', delta_trace_nh4_sz)
    #         print('teta_sm: ', tteta_usz[t], ', husz: ', thusz[t], ', teta_b: ', tteta_sz[t], ', hsz: ', thsz[t])
    #         print('sum_cusz_nh4: ', sum_cusz_nh4, ', sum_csz_nh4: ', sum_csz_nh4)             
    #         print('Qin: ', tQin[t], ',hpEND_i: ', thpEND[t], ', hpEND_i-1: ', thpEND[t-1], ',Qpipe: ', tQpipe[t], ', Qover: ', tQover[t], ', Qinfsz: ', tQinfsz[t], ', Qet: ', tQet[t])
    #         print('Cin: ', cin_nh4_list[t], ', Cp_i: ', cp_nh4[t], ', Cp_i-1: ', cp_nh4[t-1], ', Csz_ultimo: ', c_pipe)
    #         input()
    #                      
    #         for l in range(m_usz):
    #             cl_i1_nh4[l] = cl_i1_nh4[l] + delta_ast_nh4_usz
    #         for j in range(m_sz):
    #             cj_i1_nh4[j] = cj_i1_nh4[j] + delta_ast_nh4_sz 
        
            ### Nitrate
            
            Mstor_no3_ast = sum(c_usz_no3[t])*Ab*tteta_usz[t]*thusz[t]*1000/m_usz + sum(c_sz_no3[t])*Ab*tteta_sz[t]*thsz[t]*1000/m_sz
            Mstor_no3_ast_list.append(Mstor_no3_ast)
            
            Msoil_no3_a = 0
            Msoil_no3 = 0
            Msoil_no3_acum.append(0)
            
       
            MRx_no3 = - (sum(Rx_usz_no3[t])*Ab*tteta_usz[t]*thusz[t]*1000/m_usz + sum(Rx_sz_no3[t])*Ab*tteta_sz[t]*thsz[t]*1000/m_sz)
            if t == 0:
                MRx_no3_acum.append(MRx_no3)
            else:
                MRx_no3_acum.append(MRx_no3 + MRx_no3_acum[-1])
    
    
            Min_no3 = tQin[t]*cin_no3_list[t]*dt*1000
            if t == 0:
                Min_no3_acum.append(Min_no3)
            else:
                Min_no3_acum.append(Min_no3 + Min_no3_acum[-1])
    
            Mover_no3 = tQover[t]*cp_no3[t]*dt*1000
            if t == 0:
                Mover_no3_acum.append(Mover_no3)
            else:
                Mover_no3_acum.append(Mover_no3 + Mover_no3_acum[-1])
             
            Mpipe_no3 = tQpipe[t]*c_sz_no3[t][m_sz-1]*dt*1000 
            if t == 0:
                Mpipe_no3_acum.append(Mpipe_no3)
            else:
                Mpipe_no3_acum.append(Mpipe_no3 + Mpipe_no3_acum[-1])
             
            Minfsz_no3 = tQinfsz[t]*c_sz_no3[t][m_sz-1]*dt*1000 
            if t == 0:
                Minfsz_no3_acum.append(Minfsz_no3)
            else:
                Minfsz_no3_acum.append(Minfsz_no3 + Minfsz_no3_acum[-1])
             
            Met_no3 = tQet[t]*cl_i1_no3[0]*dt*1000
            if t == 0:
                Met_no3_acum.append(Met_no3)
            else:
                Met_no3_acum.append(Met_no3 + Met_no3_acum[-1])
             
            Mpz_no3 = thpEND[t]*Ab*1000*cp_no3[t]
            Mpz_no3_list.append(Mpz_no3)
             
        
            Mstor_no3_mb = Min_no3_acum[-1] - Mpz_no3_list[-1] - Mover_no3_acum[-1] - Mpipe_no3_acum[-1] - Minfsz_no3_acum[-1] - Met_no3_acum[-1] - Msoil_no3_acum[-1] - MRx_no3_acum[-1]
            Mstor_no3_mb_list.append(Mstor_no3_mb)        
            
            delta_ast_no3_usz = (Mstor_no3_mb - Mstor_no3_ast)/(Ab*tteta_usz[t]*thusz[t]*1000/m_usz)
            cl_i1_no3[1] = cl_i1_no3[1] + delta_ast_no3_usz
            if cl_i1_no3[1] > 0:
                cl_i1_no3[1] = cl_i1_no3[1]
            else:
                cl_i1_no3[1] = 0  
            
            ### DOC
            
            Mstor_doc_ast = sum(c_usz_doc[t])*Ab*tteta_usz[t]*thusz[t]*1000/m_usz + sum(c_sz_doc[t])*Ab*tteta_sz[t]*thsz[t]*1000/m_sz
            Mstor_doc_ast_list.append(Mstor_doc_ast)
            
            Msoil_doc_a = cs_usz_doc_a*ro*Ab*thusz[t]*1000 + cs_sz_doc_a*ro*Ab*thsz[t]*1000
            Msoil_doc = sum(cs_usz_doc[t])*ro*Ab*thusz[t]*1000/m_usz + sum(cs_sz_doc[t])*ro*Ab*thsz[t]*1000/m_sz
            Msoil_doc_acum.append(Msoil_doc - Msoil_doc_a)
            
       
            MRx_doc = - (sum(Rx_usz_doc[t])*Ab*tteta_usz[t]*thusz[t]*1000/m_usz + sum(Rx_sz_doc[t])*Ab*tteta_sz[t]*thsz[t]*1000/m_sz)
            if t == 0:
                MRx_doc_acum.append(MRx_doc)
            else:
                MRx_doc_acum.append(MRx_doc + MRx_doc_acum[-1])
    
    
            Min_doc = tQin[t]*cin_doc_list[t]*dt*1000
            if t == 0:
                Min_doc_acum.append(Min_doc)
            else:
                Min_doc_acum.append(Min_doc + Min_doc_acum[-1])
    
            Mover_doc = tQover[t]*cp_doc[t]*dt*1000
            if t == 0:
                Mover_doc_acum.append(Mover_doc)
            else:
                Mover_doc_acum.append(Mover_doc + Mover_doc_acum[-1])
             
            Mpipe_doc = tQpipe[t]*c_sz_doc[t][m_sz-1]*dt*1000 
            if t == 0:
                Mpipe_doc_acum.append(Mpipe_doc)
            else:
                Mpipe_doc_acum.append(Mpipe_doc + Mpipe_doc_acum[-1])
             
            Minfsz_doc = tQinfsz[t]*c_sz_doc[t][m_sz-1]*dt*1000 
            if t == 0:
                Minfsz_doc_acum.append(Minfsz_doc)
            else:
                Minfsz_doc_acum.append(Minfsz_doc + Minfsz_doc_acum[-1])
             
            Met_doc = tQet[t]*cl_i1_doc[0]*dt*1000
            if t == 0:
                Met_doc_acum.append(Met_doc)
            else:
                Met_doc_acum.append(Met_doc + Met_doc_acum[-1])
             
            Mpz_doc = thpEND[t]*Ab*1000*cp_doc[t]
            Mpz_doc_list.append(Mpz_doc)
             
        
            Mstor_doc_mb = Min_doc_acum[-1] - Mpz_doc_list[-1] - Mover_doc_acum[-1] - Mpipe_doc_acum[-1] - Minfsz_doc_acum[-1] - Met_doc_acum[-1] - Msoil_doc_acum[-1] - MRx_doc_acum[-1]
            Mstor_doc_mb_list.append(Mstor_doc_mb)        
            
            delta_ast_doc_usz = (Mstor_doc_mb - Mstor_doc_ast)/(Ab*tteta_usz[t]*thusz[t]*1000/m_usz)
            cl_i1_doc[1] = cl_i1_doc[1] + delta_ast_doc_usz
            if cl_i1_doc[1] > 0:
                cl_i1_doc[1] = cl_i1_doc[1]
            else:
                cl_i1_doc[1] = 0
        
        
        else: #if hpipe == 0
            
            ### Oxygen
            
            Mstor_o2_ast = sum(c_usz_o2[t])*Ab*tteta_usz[t]*thusz[t]*1000/m_usz
            Mstor_o2_ast_list.append(Mstor_o2_ast)
            
            Msoil_o2_a = 0
            Msoil_o2 = 0
            Msoil_o2_acum.append(0)
            
            MRx_o2 = - (sum(Rx_usz_o2[t])*Ab*tteta_usz[t]*thusz[t]*1000/m_usz)
            if t == 0:
                MRx_o2_acum.append(MRx_o2)
            else:
                MRx_o2_acum.append(MRx_o2 + MRx_o2_acum[-1])
    
    
            Min_o2 = tQin[t]*cin_o2_list[t]*dt*1000
            if t == 0:
                Min_o2_acum.append(Min_o2)
            else:
                Min_o2_acum.append(Min_o2 + Min_o2_acum[-1])
    
            Mover_o2 = tQover[t]*cp_o2[t]*dt*1000
            if t == 0:
                Mover_o2_acum.append(Mover_o2)
            else:
                Mover_o2_acum.append(Mover_o2 + Mover_o2_acum[-1])
             
            Mpipe_o2 = tQpipe[t]*c_usz_o2[t][m_usz-1]*dt*1000 
            if t == 0:
                Mpipe_o2_acum.append(Mpipe_o2)
            else:
                Mpipe_o2_acum.append(Mpipe_o2 + Mpipe_o2_acum[-1])
             
            Minfsz_o2 = tQinfsz[t]*c_usz_o2[t][m_usz-1]*dt*1000 
            if t == 0:
                Minfsz_o2_acum.append(Minfsz_o2)
            else:
                Minfsz_o2_acum.append(Minfsz_o2 + Minfsz_o2_acum[-1])
             
            Met_o2 = tQet[t]*cl_i1_o2[0]*dt*1000
            if t == 0:
                Met_o2_acum.append(Met_o2)
            else:
                Met_o2_acum.append(Met_o2 + Met_o2_acum[-1])
             
            Mpz_o2 = thpEND[t]*Ab*1000*cp_o2[t]
            Mpz_o2_list.append(Mpz_o2)
             
        
            Mstor_o2_mb = Min_o2_acum[-1] - Mpz_o2_list[-1] - Mover_o2_acum[-1] - Mpipe_o2_acum[-1] - Minfsz_o2_acum[-1] - Met_o2_acum[-1] - Msoil_o2_acum[-1] - MRx_o2_acum[-1]
            Mstor_o2_mb_list.append(Mstor_o2_mb)        
            
            delta_ast_o2_usz = (Mstor_o2_mb - Mstor_o2_ast)/(Ab*tteta_usz[t]*thusz[t]*1000/m_usz)
            cl_i1_o2[1] = cl_i1_o2[1] + delta_ast_o2_usz
            if cl_i1_o2[1] > 0:
                cl_i1_o2[1] = cl_i1_o2[1]
            else:
                cl_i1_o2[1] = 0
            
            ### Amonia
        
            Mstor_nh4_ast = sum(c_usz_nh4[t])*Ab*tteta_usz[t]*thusz[t]*1000/m_usz
            Mstor_nh4_ast_list.append(Mstor_nh4_ast)
            
            Msoil_nh4_a = cs_usz_nh4_a*ro*Ab*thusz[t]*1000 
            Msoil_nh4 = sum(cs_usz_nh4[t])*ro*Ab*thusz[t]*1000/m_usz
            Msoil_nh4_acum.append(Msoil_nh4 - Msoil_nh4_a)
            
       
            MRx_nh4 = - (sum(Rx_usz_nh4[t])*Ab*tteta_usz[t]*thusz[t]*1000/m_usz)
            if t == 0:
                MRx_nh4_acum.append(MRx_nh4)
            else:
                MRx_nh4_acum.append(MRx_nh4 + MRx_nh4_acum[-1])
    
    
            Min_nh4 = tQin[t]*cin_nh4_list[t]*dt*1000
            if t == 0:
                Min_nh4_acum.append(Min_nh4)
            else:
                Min_nh4_acum.append(Min_nh4 + Min_nh4_acum[-1])
    
            Mover_nh4 = tQover[t]*cp_nh4[t]*dt*1000
            if t == 0:
                Mover_nh4_acum.append(Mover_nh4)
            else:
                Mover_nh4_acum.append(Mover_nh4 + Mover_nh4_acum[-1])
             
            Mpipe_nh4 = tQpipe[t]*c_usz_nh4[t][m_usz-1]*dt*1000 
            if t == 0:
                Mpipe_nh4_acum.append(Mpipe_nh4)
            else:
                Mpipe_nh4_acum.append(Mpipe_nh4 + Mpipe_nh4_acum[-1])
             
            Minfsz_nh4 = tQinfsz[t]*c_usz_nh4[t][m_usz-1]*dt*1000 
            if t == 0:
                Minfsz_nh4_acum.append(Minfsz_nh4)
            else:
                Minfsz_nh4_acum.append(Minfsz_nh4 + Minfsz_nh4_acum[-1])
             
            Met_nh4 = tQet[t]*cl_i1_nh4[0]*dt*1000
            if t == 0:
                Met_nh4_acum.append(Met_nh4)
            else:
                Met_nh4_acum.append(Met_nh4 + Met_nh4_acum[-1])
             
            Mpz_nh4 = thpEND[t]*Ab*1000*cp_nh4[t]
            Mpz_nh4_list.append(Mpz_nh4)
             
        
            Mstor_nh4_mb = Min_nh4_acum[-1] - Mpz_nh4_list[-1] - Mover_nh4_acum[-1] - Mpipe_nh4_acum[-1] - Minfsz_nh4_acum[-1] - Met_nh4_acum[-1] - Msoil_nh4_acum[-1] - MRx_nh4_acum[-1]
            Mstor_nh4_mb_list.append(Mstor_nh4_mb)        
             
            
            delta_ast_nh4_usz = (Mstor_nh4_mb - Mstor_nh4_ast)/(Ab*tteta_usz[t]*thusz[t]*1000/m_usz)
            cl_i1_nh4[1] = cl_i1_nh4[1] + delta_ast_nh4_usz
            if cl_i1_nh4[1] > 0:
                cl_i1_nh4[1] = cl_i1_nh4[1]
            else:
                cl_i1_nh4[1] = 0

            ### Nitrate
            
            Mstor_no3_ast = sum(c_usz_no3[t])*Ab*tteta_usz[t]*thusz[t]*1000/m_usz
            Mstor_no3_ast_list.append(Mstor_no3_ast)
            
            Msoil_no3_a = 0
            Msoil_no3 = 0
            Msoil_no3_acum.append(0)
            
       
            MRx_no3 = - (sum(Rx_usz_no3[t])*Ab*tteta_usz[t]*thusz[t]*1000/m_usz)
            if t == 0:
                MRx_no3_acum.append(MRx_no3)
            else:
                MRx_no3_acum.append(MRx_no3 + MRx_no3_acum[-1])
    
    
            Min_no3 = tQin[t]*cin_no3_list[t]*dt*1000
            if t == 0:
                Min_no3_acum.append(Min_no3)
            else:
                Min_no3_acum.append(Min_no3 + Min_no3_acum[-1])
    
            Mover_no3 = tQover[t]*cp_no3[t]*dt*1000
            if t == 0:
                Mover_no3_acum.append(Mover_no3)
            else:
                Mover_no3_acum.append(Mover_no3 + Mover_no3_acum[-1])
             
            Mpipe_no3 = tQpipe[t]*c_usz_no3[t][m_usz-1]*dt*1000 
            if t == 0:
                Mpipe_no3_acum.append(Mpipe_no3)
            else:
                Mpipe_no3_acum.append(Mpipe_no3 + Mpipe_no3_acum[-1])
             
            Minfsz_no3 = tQinfsz[t]*c_usz_no3[t][m_usz-1]*dt*1000 
            if t == 0:
                Minfsz_no3_acum.append(Minfsz_no3)
            else:
                Minfsz_no3_acum.append(Minfsz_no3 + Minfsz_no3_acum[-1])
             
            Met_no3 = tQet[t]*cl_i1_no3[0]*dt*1000
            if t == 0:
                Met_no3_acum.append(Met_no3)
            else:
                Met_no3_acum.append(Met_no3 + Met_no3_acum[-1])
             
            Mpz_no3 = thpEND[t]*Ab*1000*cp_no3[t]
            Mpz_no3_list.append(Mpz_no3)
             
        
            Mstor_no3_mb = Min_no3_acum[-1] - Mpz_no3_list[-1] - Mover_no3_acum[-1] - Mpipe_no3_acum[-1] - Minfsz_no3_acum[-1] - Met_no3_acum[-1] - Msoil_no3_acum[-1] - MRx_no3_acum[-1]
            Mstor_no3_mb_list.append(Mstor_no3_mb)        
            
            delta_ast_no3_usz = (Mstor_no3_mb - Mstor_no3_ast)/(Ab*tteta_usz[t]*thusz[t]*1000/m_usz)
            cl_i1_no3[1] = cl_i1_no3[1] + delta_ast_no3_usz
            if cl_i1_no3[1] > 0:
                cl_i1_no3[1] = cl_i1_no3[1]
            else:
                cl_i1_no3[1] = 0  
            
            ### DOC
            
            Mstor_doc_ast = sum(c_usz_doc[t])*Ab*tteta_usz[t]*thusz[t]*1000/m_usz
            Mstor_doc_ast_list.append(Mstor_doc_ast)
            
            Msoil_doc_a = cs_usz_doc_a*ro*Ab*thusz[t]*1000 
            Msoil_doc = sum(cs_usz_doc[t])*ro*Ab*thusz[t]*1000/m_usz
            Msoil_doc_acum.append(Msoil_doc - Msoil_doc_a)
            
       
            MRx_doc = - (sum(Rx_usz_doc[t])*Ab*tteta_usz[t]*thusz[t]*1000/m_usz)
            if t == 0:
                MRx_doc_acum.append(MRx_doc)
            else:
                MRx_doc_acum.append(MRx_doc + MRx_doc_acum[-1])
    
    
            Min_doc = tQin[t]*cin_doc_list[t]*dt*1000
            if t == 0:
                Min_doc_acum.append(Min_doc)
            else:
                Min_doc_acum.append(Min_doc + Min_doc_acum[-1])
    
            Mover_doc = tQover[t]*cp_doc[t]*dt*1000
            if t == 0:
                Mover_doc_acum.append(Mover_doc)
            else:
                Mover_doc_acum.append(Mover_doc + Mover_doc_acum[-1])
             
            Mpipe_doc = tQpipe[t]*c_usz_doc[t][m_usz-1]*dt*1000 
            if t == 0:
                Mpipe_doc_acum.append(Mpipe_doc)
            else:
                Mpipe_doc_acum.append(Mpipe_doc + Mpipe_doc_acum[-1])
             
            Minfsz_doc = tQinfsz[t]*c_usz_doc[t][m_usz-1]*dt*1000 
            if t == 0:
                Minfsz_doc_acum.append(Minfsz_doc)
            else:
                Minfsz_doc_acum.append(Minfsz_doc + Minfsz_doc_acum[-1])
             
            Met_doc = tQet[t]*cl_i1_doc[0]*dt*1000
            if t == 0:
                Met_doc_acum.append(Met_doc)
            else:
                Met_doc_acum.append(Met_doc + Met_doc_acum[-1])
             
            Mpz_doc = thpEND[t]*Ab*1000*cp_doc[t]
            Mpz_doc_list.append(Mpz_doc)
             
        
            Mstor_doc_mb = Min_doc_acum[-1] - Mpz_doc_list[-1] - Mover_doc_acum[-1] - Mpipe_doc_acum[-1] - Minfsz_doc_acum[-1] - Met_doc_acum[-1] - Msoil_doc_acum[-1] - MRx_doc_acum[-1]
            Mstor_doc_mb_list.append(Mstor_doc_mb)        
            
            delta_ast_doc_usz = (Mstor_doc_mb - Mstor_doc_ast)/(Ab*tteta_usz[t]*thusz[t]*1000/m_usz)
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
        if hpipe > 0:
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
    for i in range(m_usz):
        name = 'usz' + str(a)
        column_name.append(name)
        a = a + 1
    data_usz_o2.set_axis(column_name, axis = 'columns', inplace = True)
    
    #print(data_usz_o2)
    
    data_sz_o2 = pd.DataFrame(c_sz_o2)
    a = 0
    column_name = []
    for i in range(m_sz):
        name = 'sz' + str(a)
        column_name.append(name)
        a = a + 1
    data_sz_o2.set_axis(column_name, axis = 'columns', inplace = True)

    data_rx_usz_o2 = pd.DataFrame(Rx_usz_o2)
    a = 0
    column_name = []
    for i in range(m_usz):
        name = 'usz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_usz_o2.set_axis(column_name, axis = 'columns', inplace = True)
    
    data_rx_sz_o2 = pd.DataFrame(Rx_sz_o2)
    a = 0
    column_name = []
    for i in range(m_sz):
        name = 'sz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_sz_o2.set_axis(column_name, axis = 'columns', inplace = True)
    
    frames = [data_usz_o2, data_sz_o2, data_rx_usz_o2, data_rx_sz_o2]
    data_o2 = pd.concat((frames), axis = 1)

    cp_o2.append(0)
    data_o2['pz'] = cp_o2
    
    indice_n = list(range(len(cp_o2)))
    cin_o2_list.append(0)
    c_in = cin_o2_list[:len(indice_n)]
    #print('len c_in:', len(c_in), ', len_indice:', len(indice))
    data_o2['c_in'] = c_in 
    
    data_o2['t'] = indice_n

   
    ## Amonia
    data_usz_nh4 = pd.DataFrame(c_usz_nh4)
    a = 0
    column_name = []
    for i in range(m_usz):
        name = 'usz' + str(a)
        column_name.append(name)
        a = a + 1
    data_usz_nh4.set_axis(column_name, axis = 'columns', inplace = True)
    
    data_sz_nh4 = pd.DataFrame(c_sz_nh4)
    a = 0
    column_name = []
    for i in range(m_sz):
        name = 'sz' + str(a)
        column_name.append(name)
        a = a + 1
    data_sz_nh4.set_axis(column_name, axis = 'columns', inplace = True)
    
    data_s_usz_nh4 = pd.DataFrame(cs_usz_nh4)
    a = 0
    column_name = []
    for i in range(m_usz):
        name = 'usz_soil' + str(a)
        column_name.append(name)
        a = a + 1
    data_s_usz_nh4.set_axis(column_name, axis = 'columns', inplace = True)
    
    data_s_sz_nh4 = pd.DataFrame(cs_sz_nh4)
    a = 0
    column_name = []
    for i in range(m_sz):
        name = 'sz_soil' + str(a)
        column_name.append(name)
        a = a + 1
    data_s_sz_nh4.set_axis(column_name, axis = 'columns', inplace = True)
    
    data_rx_usz_nh4 = pd.DataFrame(Rx_usz_nh4)
    a = 0
    column_name = []
    for i in range(m_usz):
        name = 'usz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_usz_nh4.set_axis(column_name, axis = 'columns', inplace = True)
    
    data_rx_sz_nh4 = pd.DataFrame(Rx_sz_nh4)
    a = 0
    column_name = []
    for i in range(m_sz):
        name = 'sz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_sz_nh4.set_axis(column_name, axis = 'columns', inplace = True)
    
    
    frames = [data_usz_nh4, data_sz_nh4, data_s_usz_nh4, data_s_sz_nh4, data_rx_usz_nh4, data_rx_sz_nh4]
    data_nh4 = pd.concat((frames), axis = 1)
    cp_nh4.append(0)
    cin_nh4_list.append(0)
    c_in = cin_nh4_list[:len(indice_n)]
    #print('len c_in:', len(c_in), ', len_indice:', len(indice))
    data_nh4['c_in'] = c_in
    
    data_nh4['pz'] = cp_nh4
    #indice.append(len(indice)+1)
    data_nh4['t'] = indice_n
    
    
    ## Nitrate
    data_usz_no3 = pd.DataFrame(c_usz_no3)
    a = 0
    column_name = []
    for i in range(m_usz):
        name = 'usz' + str(a)
        column_name.append(name)
        a = a + 1
    data_usz_no3.set_axis(column_name, axis = 'columns', inplace = True)
    
    data_sz_no3 = pd.DataFrame(c_sz_no3)
    a = 0
    column_name = []
    for i in range(m_sz):
        name = 'sz' + str(a)
        column_name.append(name)
        a = a + 1
    data_sz_no3.set_axis(column_name, axis = 'columns', inplace = True)
    
    data_rx_usz_no3 = pd.DataFrame(Rx_usz_no3)
    a = 0
    column_name = []
    for i in range(m_usz):
        name = 'usz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_usz_no3.set_axis(column_name, axis = 'columns', inplace = True)
    
    data_rx_sz_no3 = pd.DataFrame(Rx_sz_no3)
    a = 0
    column_name = []
    for i in range(m_sz):
        name = 'sz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_sz_no3.set_axis(column_name, axis = 'columns', inplace = True)    
    
    frames = [data_usz_no3, data_sz_no3, data_rx_usz_no3, data_rx_sz_no3]
    data_no3 = pd.concat((frames), axis = 1)
    cp_no3.append(0)
    data_no3['pz'] = cp_no3
    
    cin_no3_list.append(0)
    c_in = cin_no3_list[:len(indice_n)]
    data_no3['c_in'] = c_in

    data_no3['t'] = indice_n
    
    ## DOC
    data_usz_doc = pd.DataFrame(c_usz_doc)
    a = 0
    column_name = []
    for i in range(m_usz):
        name = 'usz' + str(a)
        column_name.append(name)
        a = a + 1
    data_usz_doc.set_axis(column_name, axis = 'columns', inplace = True)
    
    data_sz_doc = pd.DataFrame(c_sz_doc)
    a = 0
    column_name = []
    for i in range(m_sz):
        name = 'sz' + str(a)
        column_name.append(name)
        a = a + 1
    data_sz_doc.set_axis(column_name, axis = 'columns', inplace = True)

    data_rx_usz_doc = pd.DataFrame(Rx_usz_doc)
    a = 0
    column_name = []
    for i in range(m_usz):
        name = 'usz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_usz_doc.set_axis(column_name, axis = 'columns', inplace = True)
    
    data_rx_sz_doc = pd.DataFrame(Rx_sz_doc)
    a = 0
    column_name = []
    for i in range(m_sz):
        name = 'sz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_sz_doc.set_axis(column_name, axis = 'columns', inplace = True)
    
    frames = [data_usz_doc, data_sz_doc, data_rx_usz_doc, data_rx_sz_doc]
    data_doc = pd.concat((frames), axis = 1)
    cp_doc.append(0)
    data_doc['pz'] = cp_doc
    
    cin_doc_list.append(0)
    c_in = cin_doc_list[:len(indice_n)]
    data_doc['c_in'] = c_in
    
    data_doc['t'] = indice_n
    return data_o2, data_nh4, data_no3, data_doc
    
if __name__ == '__main__':
    inicio = datetime.datetime.now()
    
    data_o2, data_nh4, data_no3, data_doc = run_Kin()

    data_nh4.to_csv('results/results_40/results_Kin_pf_nh4_2.csv', index = False)
    data_o2.to_csv('results/results_40/results_Kin_pf_o2.csv', index = False)
    data_no3.to_csv('results/results_40/results_Kin_pf_no3.csv', index = False)
    data_doc.to_csv('results/results_40/results_Kin_pf_doc.csv', index = False)
    
    fim = datetime.datetime.now()
    print ('Elapsed time: ', fim - inicio)
    print('Done!')
