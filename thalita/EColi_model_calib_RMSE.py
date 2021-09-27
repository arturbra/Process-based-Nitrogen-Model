__author__ = 'thalita'

# !/usr/bin/env python
# coding: utf-8

import logging
import datetime
import numpy as np
import pandas as pd  # data processing, CSV file I/O (e.g. pd.read_csv)
import configparser
import math
import matplotlib.pyplot as plt

plt.style.use('ggplot')

print('Starting..')
inicio = datetime.datetime.now()

setup_file = "parametros_quali_RMSE.ini"
setup = configparser.ConfigParser()
setup.read(setup_file)

# **1. Parameters and variables**

# Soil characteristics
ro_pd = float(setup['SOIL']['ro_pd'])
lamta1 = float(setup['SOIL']['lamta1'])
d50 = float(setup['SOIL']['d50'])

# lamta2 = float(setup['SOIL']['lamta2'])

# E. Coli constants
kads1 = float(setup['ECOLI']['kads1'])
kdes1 = float(setup['ECOLI']['kdes1'])
# kads2 = float(setup['ECOLI']['kads2'])
# kdes2 = float(setup['ECOLI']['kdes2'])

sita = float(setup['ECOLI']['sita'])
mue1 = float(setup['ECOLI']['mue1'])
# mue2 = float(setup['ECOLI']['mue2'])

# kstr = float(setup['ECOLI']['kstr'])
kstr = 0
etta = float(setup['ECOLI']['etta'])
b1 = float(setup['ECOLI']['b1'])

kads2 = kads1
mue2 = mue1
kdes2 = kdes1
lamta2 = lamta1


####  Parameters to calibrate  ####
def read_from_ini():
    param = {}

   # param['lamta1'] = lamta1
    param['kads1'] = kads1
    param['kdes1'] = kdes1
    param['mue1'] = mue1
    param['sita'] = sita

    # param['lamta1'] = lamta1_ini
    # param['kads1'] = kads1_ini
    # param['kdes1'] = kdes1_ini
    # param['etta'] = etta_ini

    # param['lamta2'] = lamta2_ini
    # param['kads2'] = kads2_ini
    # param['kdes2'] = kdes2_ini

    return param


def get_from_individual(individual):
    param = {}
    #param['lamta1'] = individual[0]
    param['kads1'] = individual[0]
    param['kdes1'] = individual[1]
    param['mue1'] = individual[2]
    param['sita'] = individual[3]

    # #param['lamta1'] = individual[0]
    # param['kads1'] = individual[1]
    # param['kdes1'] = individual[2]
    # param['mue1'] = individual[3]
    # param['sita'] = individual[4]

    # param['lamta2'] = individual[1]
    # param['kads2'] = individual[4]
    # param['kdes2'] = individual[5]
    return param


# **2. Input variables**
from water_flow_model_corrigido import *

wf_run()

cin_file = pd.read_csv('01_Cin_file2.csv')
cin_list = cin_file['Cin e coli [mg/L]'].tolist()

Temp_file = pd.read_csv('01_Temperature_file2.csv')
tTemp = Temp_file['Temperature [celsius]'].tolist()

Qin_file = pd.read_csv('01_Qin_file2.csv')
Qin_list = Qin_file['Qin'].tolist()

# cin_file = pd.read_csv('input_file_Cin.csv')
# cin_list = cin_file['Cin e coli [mg/L]'].tolist()
#
# Temp_file = pd.read_csv('input_file_Temperature.csv')
# tTemp = Temp_file['Temperature [celsius]'].tolist()
#
# Qin_file = pd.read_csv('input_file_Qin.csv')
# Qin_list = Qin_file['Qin'].tolist()

# **3. Definitions of the equations**

n = 11  # number of cells
dz = L / n
m_usz = round((L - hpipe) / dz)
m_sz = n - m_usz
ro = (1 - nusz_ini) * ro_pd


##### Unitary flow to transport equations #####
def falfa_beta_usz(l):
    alfa1 = (m_usz - 1 - l) / (m_usz - 1)
    beta1 = l / (m_usz - 1)

    return alfa1, beta1


def falfa_beta_sz(j):
    if hpipe <= 0.03:
        alfa2 = 0.5
        beta2 = 0.5
    else:
        alfa2 = (m_sz - 1 - j) / (m_sz - 1)
        beta2 = j / (m_sz - 1)

    return alfa2, beta2


def fV_usz(alfa1, beta1, Qpf, Qet1, Qfs, Qhc, Qpipe, Qinf_sz, teta_i_usz):
    if hpipe > 0.03:
        V_usz = (alfa1 * (Qpf - Qet1) + beta1 * (Qfs - Qhc)) / (Ab * teta_i_usz)

    else:
        V_usz = (alfa1 * (Qpf - Qet1) + beta1 * (Qpipe + Qinf_sz - Qhc)) / (Ab * teta_i_usz)

    return V_usz


def fV_sz(alfa2, beta2, Qfs, Qhc, Qet2, Qpipe, Qinf_sz, teta_i_sz):
    V_sz = (alfa2 * (Qfs - Qhc - Qet2) + beta2 * (Qpipe + Qinf_sz)) / (Ab * teta_i_sz)

    return V_sz


def fD_x(V_x):
    D_x = lamta1 * V_x
    return D_x


def fPe_usz(V_usz, D_usz):
    if D_usz > 0:
        Pe_usz = V_usz * dz / D_usz
    else:
        Pe_usz = 100

    return Pe_usz


def fPe_sz(V_sz, D_sz):
    D_sz = lamta2 * V_sz
    if D_sz > 0:
        Pe_sz = V_sz * dz / D_sz
    else:
        Pe_sz = 100

    return Pe_sz


##### Water mass balance to ponding zone #####
def fcpz(cin, Qin, cpz_a, Qpf, Qv, Rxi, hp, hp_a):
    # delta_cpz = ((cin*Qin_p - cpz_a*(Qpf + Qv))*dt)/(hp*Ab) + Rxi*dt
    # cpz = cpz_a + delta_cpz

    # cpz = Rxi*dt + (cin*Qin_p*dt)/(hp*Ab + dt*(Qpf + Qv))

    cpz = (cpz_a * hp_a * Ab + (cin * Qin - cpz_a * (Qpf + Qv) + Rxi * hp * Ab) * dt) / (hp * Ab)

    return cpz


##### Transport equations #####
def ftransport(teta_i, teta_iplus1, ci, csoil_i, dc, dc_dz, kads, kdes, D_x, V_x, Rx):
    if teta_i == 0:
        delta_c_i1 = 0
    elif teta_iplus1 == 0:
        delta_c_i1 = 0
    else:
        # delta_c_i1 = ((1/teta_iplus1)*dt*(-teta_i*kads*ci + ro*kdes*cs_i + teta_i*(D*f*(dc/dz**2) - UF*dc_dz) + Rx))
        delta_c_i1 = ((1 / teta_iplus1) * dt * (-teta_i * kads * ci
                                                + ro * kdes * csoil_i
                                                + teta_i * (D_x * (dc / dz ** 2)
                                                            - V_x * dc_dz)
                                                + Rx))
    return delta_c_i1


def fcsoil(csoil_a, teta, kads, ci, kdes, Rxs):
    # Rxs = mues*csoil_a + kstr*etta*teta*csoil_a # die-off no solo?

    csoil_abs = csoil_a + ((teta / ro) * kads * ci - kdes * csoil_a + Rxs) * dt

    if csoil_abs <= 0:
        csoil = 0
    else:
        csoil = csoil_abs

    return csoil


# def fcsoil(csoil_a, teta, kads, ci, kdes, mues):  #ORIGINAL
#
#     Rxs = mues*csoil_a + kstr*etta*teta*csoil_a # die-off no solo?
#
#     csoil_abs = csoil_a + ((teta / ro) * kads * ci - kdes * csoil_a - Rxs) * dt
#
#     if csoil_abs <= 0:
#         csoil = 0
#     else:
#         csoil = csoil_abs
#
#     return csoil


def fR_dieoff_l(teta_i, muel, ci):
    dieoff_l = -teta_i * muel * ci

    return dieoff_l


def fR_str(teta_i, ci, kstr, etta):
    # R_str = -teta_i*ci*kstr*etta
    R_str = 0
    return R_str


def fR_dieoff_s(mues, csoil):
    dieoff_s = -mues * csoil

    # dieoff_s = 0

    return dieoff_s


# def fR_dieoff_s(teta_i, mues, csoil):
#
#     dieoff_s= -teta_i*mues*csoil
#
#     #dieoff_s = 0
#
#     return dieoff_s

def fRx_pz(muel, cpz_a):
    Rx_pz = -muel * cpz_a
    return Rx_pz


####### **4. Model routine**
def EColi_run(param):
#    lamta1 = param['lamta1']
    kads1 = param['kads1']
    kdes1 = param['kdes1']
    mue1 = param['mue1']
    sita = param['sita']

    # lamta2 = param['lamta2']
    # kads2 = param['kads2']
    # kdes2 = param['kdes2']

    ###Listas de E. coli

    c_usz_list = []
    csoil_usz_list = []
    c_sz_list = []
    csoil_sz_list = []

    Rx_usz_list = []
    Rx_sz_list = []

    # E_coli #PZ
    cpz_a = 0
    cpz = []

    # # E_coli #USZ                  LISTAS ANTES
    # c0_usz = [0]
    # c0_usz = c0_usz * m_usz
    # c_usz_list.append(c0_usz)
    # csoil_usz_list.append(c0_usz)
    # Rx_usz_list.append(c0_usz)
    #
    # # E_coli #SZ
    # c0_sz = [0]
    # c0_sz = c0_sz*m_sz
    # c_sz_list.append(c0_sz)
    # csoil_sz_list.append(c0_sz)
    # Rx_sz_list.append(c0_sz)

    csoil_usz_a = 0
    csoil_sz_a = 0

    # LISTAS DEPOIS
    # E_coli #USZ
    c0_usz = [0]
    c0_usz = c0_usz * m_usz
    csoil0_usz = [csoil_usz_a]
    csoil0_usz = csoil0_usz * m_usz
    c_usz_list.append(c0_usz)
    csoil_usz_list.append(csoil0_usz)
    Rx_usz_list.append(c0_usz)
    # E_coli #SZ
    c0_sz = [0]
    c0_sz = c0_sz * m_sz
    c_sz_list.append(c0_sz)
    csoil_sz_list.append(c0_sz)
    Rx_sz_list.append(c0_sz)

    # corrective step
    Mstor_ast_list = []
    Msoil_acum = []
    MRx_acum = []
    M_in_acum = []
    M_over_acum = []
    Mpipe_acum = []
    Minf_sz_acum = []
    M_et_acum = []
    M_pz_list = []
    Mstor_mb_list = []



    indice = list(range(0, len(tQrain)))

    ##### Ponding Zone

    for t in range(len(indice) - 1):
        cin = cin_list[t]
        hp = thpEND[t]

        if hp < 0.001:
            hp = 0
        else:
            hp = hp

        if t == 0:
            hp_a = 0
        else:
            hp_a = thpEND[t - 1]

        Qin = tQin[t]
        Qpf = tQpf[t]
        Qv = tQover[t]
        Qet1 = tQet1[t]
        Qfs = tQfs[t]
        Qhc = tQhc[t]
        Qet2 = tQet2[t]
        if t == 0:
            Qpipe = 0
        else:
            Qpipe = tQpipe[t]

        Qinf_sz = tQinf_sz[t]
        teta_i_usz = tteta_usz[t]
        teta_i_sz = tteta_sz[t]

        ### Reacao na ponding zone ###
        if sita == 0.0:
            sita = 0.0000000000001

        muel = mue1 * sita ** (tTemp[t] - 20)
        Rx_pz = fRx_pz(muel, cpz_a)
        # Rx_pz=0

        if t < (len(indice) - 1):
            teta_usz_iplus1 = tteta_usz[t + 1]
            teta_sz_iplus1 = tteta_sz[t + 1]
        else:
            teta_usz_iplus1 = tteta_usz[t]
            teta_sz_iplus1 = tteta_sz[t]

        if hp == 0:
            cpz_i = 0
        else:
            cpz_i = fcpz(cin, Qin, cpz_a, Qpf, Qv, Rx_pz, hp, hp_a)

        if cpz_i < 0.00002:
            cpz_i = 0
        else:
            cpz_i = cpz_i

        cpz.append(cpz_i)
        cpz_a = cpz[-1]

        # print('cin', cin)
        # print('Qin', Qin)
        # print('Qpf', Qpf)
        # print('Qv', Qv)
        # print('hp_a', hp_a)
        # print('hp', hp)
        # print('rx_pz', Rx_pz)
        # print('.')

        # print('cpz_i', cpz_i)
        # input()

        #### SOIL MIX

        ## die-off in USZ, according to s[t]
        mues1 = mue1 * sita ** (tTemp[t] - 20)
        muel1 = mue1 * sita ** (tTemp[t] - 20)
        muel2 = mue2 * sita ** (tTemp[t] - 20)
        mues2 = mue2 * sita ** (tTemp[t] - 20)

        # USZ
        cli_list = c_usz_list[t].copy()

        # SZ
        if hpipe >= 0.03:
            cji_list = c_sz_list[t].copy()  # copiar lista

        if m_usz != 0:

            cl_i1 = []
            csoil_i_usz = []
            Rxl = []

            # Predictive step

            ###   USZ   ####
            for l in range(m_usz):

                cl_usz = cli_list[l]
                cl_minus1 = cli_list[l - 1]
                if l < (m_usz - 1):
                    cl_plus1 = cli_list[l + 1]
                else:
                    cl_plus1 = 0

                ## Concentracao no solo - usz ##

                csoil_usz_a = csoil_usz_list[t][l]

                Rxi_soil_usz = fR_str(teta_i_usz, csoil_usz_a, kstr, etta)

                csoil_usz = fcsoil(csoil_usz_a, teta_i_usz, kads1, cl_usz, kdes1, Rxi_soil_usz)

                if csoil_usz < 0.00000000000001:
                    csoil_usz = 0

                csoil_i_usz.append(csoil_usz)

                V_usz = []
                # Qet1 = 0
                alfa1, beta1 = falfa_beta_usz(l)

                Vi_usz = fV_usz(alfa1, beta1, Qpf, Qet1, Qfs, Qhc, Qpipe, Qinf_sz, teta_i_usz)

                V_usz.append(Vi_usz)

                ## Reacoes ## dieoff_l (parte liquida) + straining + dieoff (parte solida)
                # Rxi_usz = 0
                Rxi_usz = fR_dieoff_l(teta_i_usz, muel1, cl_usz) + fR_str(teta_i_usz, cl_usz, kstr, etta) + fR_dieoff_s(
                    mues1, csoil_usz)

                Rxl.append(Rxi_usz * (1 / teta_usz_iplus1) * dt)

                ## Coeficiente de difusao ##
                Di_usz = fD_x(Vi_usz)
                Pe_usz = fPe_usz(Vi_usz, Di_usz)

                if Pe_usz <= 2:

                    if l == 0:  # first cell
                        dc = cl_plus1 - 2 * cl_usz + cpz_i
                        dc_dz = (cl_plus1 - cpz_i) / (2 * dz)

                    elif l == m_usz - 1:  # last cell
                        dc = cl_usz - 2 * cl_usz + cl_minus1
                        dc_dz = (cl_usz - cl_minus1) / dz

                    else:
                        dc = cl_plus1 - 2 * cl_usz + cl_minus1
                        dc_dz = (cl_plus1 - cl_minus1) / (2 * dz)

                else:  # Peusz > 2
                    if l == 0:  # first cell
                        dc = cl_plus1 - 2 * cl_usz + cpz_i
                        dc_dz = (cl_usz - cpz_i) / dz

                    elif l == (m_usz - 1):  # last cell
                        dc = cl_usz - 2 * cl_usz + cl_minus1
                        dc_dz = (cl_usz - cl_minus1) / dz

                    else:
                        dc = cl_plus1 - 2 * cl_usz + cl_minus1
                        dc_dz = (cl_usz - cl_minus1) / dz

                ### Transport usz ###
                # delta_c_usz = ftransport(teta_i_usz, teta_usz_iplus1, cl_usz, csoil_usz, dc, dc_dz, kads, kdes1,
                # Di_usz, Vi_usz, Rxi_usz)
                delta_c_usz = ftransport(teta_i_usz, teta_usz_iplus1, cl_usz, csoil_usz, dc, dc_dz, kads1, kdes1,
                                         Di_usz, Vi_usz, Rxi_usz)
                #
                # delta_c_usz1 = ((1 / teta_usz_iplus1) * dt * (-teta_i_usz * kads1 * cl_usz
                #                                         + ro * kdes1 * csoil_usz
                #                                         + teta_i_usz * (Di_usz * (dc / dz ** 2) - Vi_usz * dc_dz)
                #                                         + Rxi_usz))

                if teta_usz_iplus1 > 0:
                    ci1 = cl_usz + delta_c_usz
                else:
                    ci1 = 0

                if ci1 < 0.01:
                    ci1 = 0
                else:
                    ci1 = ci1

                cl_i1.append(ci1)

                # if l == 3 or l == 4 or l == 5:
                #     print('t = ', t, 'min')
                #     print('layer = ', l)
                #     print('ci1 =', ci1)
                #
                #
                #     input()

            ### SZ ###

            cj_i1 = []
            csoil_i_sz = []
            Rxj = []

            for j in range(m_sz):
                cj_sz = cji_list[j]
                cjminus1 = cji_list[j - 1]
                if j < (m_sz - 1):
                    cjplus1 = cji_list[j + 1]
                else:
                    cjplus1 = 0
                cmminus1 = cli_list[m_usz - 1]

                ## Concentracao no solo - sz ##
                csoil_sz_a = csoil_sz_list[t][j]  # chamar a lista mais ampla

                Rxi_soil_sz = 0

                csoil_sz = fcsoil(csoil_sz_a, teta_i_sz, kads2, cj_sz, kdes2, Rxi_soil_sz)

                if csoil_sz < 0.00000000000001:
                    csoil_sz = 0

                csoil_i_sz.append(csoil_sz)

                V_sz = []
                # Qet1 = 0
                alfa2, beta2 = falfa_beta_sz(j)

                Vi_sz = fV_sz(alfa2, beta2, Qfs, Qhc, Qet2, Qpipe, Qinf_sz, teta_i_sz)
                V_sz.append(Vi_sz)

                ## Reacoes ## dieoff_l (parte liquida) + dieoff_s (parte solida)
                Rxi_sz = fR_dieoff_l(teta_i_sz, muel2, cj_sz) + fR_dieoff_s(mues2, csoil_sz)

                # print('Rxi_sz', Rxi_sz)
                # input()
                # Rxi_sz = 0

                Rxj.append(Rxi_sz * (1 / teta_sz_iplus1) * dt)

                ## Coeficiente de difusao ##
                Di_sz = fD_x(Vi_sz)
                Pe_sz = fPe_sz(Vi_sz, Di_sz)

                if m_usz < (n - 1):
                    if Pe_sz <= 2:
                        if j == 0:  # first cell
                            dc = cjplus1 - 2 * cj_sz + cmminus1
                            dc_dz = (cjplus1 - cmminus1) / (2 * dz)

                        elif j == (m_sz - 1):  # last cell
                            dc = cj_sz - 2 * cj_sz + cjminus1
                            dc_dz = (cj_sz - cjminus1) / dz

                        else:
                            dc = cjplus1 - 2 * cj_sz + cjminus1
                            dc_dz = (cjplus1 - cjminus1) / (2 * dz)

                    else:  # Pesz > 2
                        if j == 0:  # first cell
                            dc = cjplus1 - 2 * cj_sz + cmminus1
                            dc_dz = (cj_sz - cmminus1) / dz

                        elif j == (m_sz - 1):  # last cell
                            dc = cj_sz - 2 * cj_sz + cjminus1
                            dc_dz = (cj_sz - cjminus1) / dz

                        else:
                            dc = cjplus1 - 2 * cj_sz + cjminus1
                            dc_dz = (cj_sz - cjminus1) / dz

                if m_usz == (n - 1):
                    dc = cj_sz - 2 * cj_sz + cjminus1
                    dc_dz = (cj_sz - cjminus1) / (2 * dz)

                if m_usz == 0:
                    if Pe_sz <= 2:
                        if j == 0:  # first cell
                            dc = cjplus1 - 2 * cj_sz + cpz_i
                            dc_dz = (cjplus1 - cpz_i) / (2 * dz)

                        elif j == (m_sz - 1):  # last cell
                            dc = cj_sz - 2 * cj_sz + cjminus1
                            dc_dz = (cj_sz - cjminus1) / dz

                        else:
                            dc = cjplus1 - 2 * cj_sz + cjminus1
                            dc_dz = (cjplus1 - cjminus1) / (2 * dz)

                    else:  # Pe_usz > 2
                        if j == 0:  # first cell
                            dc = cjplus1 - 2 * cj_sz + cpz_i
                            dc_dz = (cj_sz - cpz_i) / dz

                        elif j == (m_sz - 1):  # last cell
                            dc = cj_sz - 2 * cj_sz + cjminus1
                            dc_dz = (cj_sz - cjminus1) / dz

                        else:
                            dc = cjplus1 - 2 * cj_sz + cjminus1
                            dc_dz = (cj_sz - cjminus1) / dz

                ### Transport sz ###
                # delta_c_sz = ftransport(teta_i_sz, teta_sz_iplus1, cj_sz, csoil_sz, dc, dc_dz, kads2, kdes2,
                #                         Di_sz, Vi_sz,Rxi_sz)

                delta_c_sz = ftransport(teta_i_sz, teta_sz_iplus1, cj_sz, csoil_sz, dc, dc_dz, kads2, kdes2,
                                        Di_sz, Vi_sz, Rxi_sz)

                if teta_sz_iplus1 > 0:
                    ci1 = cj_sz + delta_c_sz
                else:
                    ci1 = 0

                if ci1 < 0.01:
                    ci1 = 0
                else:
                    ci1 = ci1

                cj_i1.append(ci1)

        # Corrective step

        if hpipe > 0.03:

            Mstor_ast = sum(c_usz_list[t]) * Ab * tteta_usz[t] * thusz[t] * 1 / m_usz + sum(c_sz_list[t]) * Ab * \
                        tteta_sz[t] * thsz[t] * 1 / m_sz
            Mstor_ast_list.append(Mstor_ast)

            Msoil_a = csoil_usz_a * ro * Ab * thusz[t] * 1 + csoil_sz_a * ro * Ab * thsz[t] * 1
            Msoil = sum(csoil_usz_list[t]) * ro * Ab * thusz[t] * 1 / m_usz + sum(csoil_usz_list[t]) * ro * Ab * thsz[
                t] * 1 / m_sz
            Msoil_acum.append(Msoil - Msoil_a)

            MRx = - (sum(Rx_usz_list[t]) * Ab * tteta_usz[t] * thusz[t] * 1 / m_usz + sum(Rx_sz_list[t]) * Ab *
                     tteta_sz[t] * thsz[t] * 1 / m_sz)

            if t == 0:
                MRx_acum.append(MRx)
            else:
                MRx_acum.append(MRx + MRx_acum[-1])

            M_in = tQin[t] * cin_list[t] * dt * 1
            if t == 0:
                M_in_acum.append(M_in)
            else:
                M_in_acum.append(M_in + M_in_acum[-1])

            M_over = tQover[t] * cpz[t] * dt * 1
            if t == 0:
                M_over_acum.append(M_over)
            else:
                M_over_acum.append(M_over + M_over_acum[-1])

            Mpipe = tQpipe[t] * c_sz_list[t][m_sz - 1] * dt * 1
            if t == 0:
                Mpipe_acum.append(Mpipe)
            else:
                Mpipe_acum.append(Mpipe + Mpipe_acum[-1])

            Minf_sz = tQinf_sz[t] * c_sz_list[t][m_sz - 1] * dt * 1
            if t == 0:
                Minf_sz_acum.append(Minf_sz)
            else:
                Minf_sz_acum.append(Minf_sz + Minf_sz_acum[-1])

            M_et = tQet[t] * cl_i1[0] * dt * 1
            if t == 0:
                M_et_acum.append(M_et)
            else:
                M_et_acum.append(M_et + M_et_acum[-1])

            M_pz = thpEND[t] * Ab * 1 * cpz[t]
            M_pz_list.append(M_pz)

            Mstor_mb = M_in_acum[-1] - M_pz_list[-1] - M_over_acum[-1] - Mpipe_acum[-1] - \
                       Minf_sz_acum[-1] - M_et_acum[-1] - Msoil_acum[-1] - MRx_acum[-1]
            Mstor_mb_list.append(Mstor_mb)

            #         delta_ast_usz = 0
            #         delta_ast_sz = 0

            #         delta_ast_usz = (Mstor_mb - Mstor_ast)/(n*(Ab*tteta_usz[t]*thusz[t]*1000/m_usz))
            #         delta_ast_sz = (Mstor_mb - Mstor_ast)/(n*(Ab*tteta_sz[t]*thsz[t]*1000/m_sz))

            delta_ast_usz = (Mstor_mb - Mstor_ast) / (Ab * tteta_usz[t] * thusz[t] * 1 / m_usz)
            cl_i1[1] = cl_i1[1] + delta_ast_usz
            if cl_i1[1] > 0:
                cl_i1[1] = cl_i1[1]
            else:
                cl_i1[1] = 0

        if hpipe == 0:
            # if hpipe == 0.03:

            Mstor_ast = sum(c_usz_list[t]) * Ab * tteta_usz[t] * thusz[t] * 1 / m_usz
            Mstor_ast_list.append(Mstor_ast)

            Msoil_a = csoil_usz_a * ro * Ab * thusz[t] * 1
            Msoil = sum(csoil_usz_list[t]) * ro * Ab * thusz[t] * 1 / m_usz
            Msoil_acum.append(Msoil - Msoil_a)

            MRx = - (sum(Rx_usz_list[t]) * Ab * tteta_usz[t] * thusz[t] * 1 / m_usz)
            if t == 0:
                MRx_acum.append(MRx)
            else:
                MRx_acum.append(MRx + MRx_acum[-1])

            M_in = tQin[t] * cin_list[t] * dt * 1
            if t == 0:
                M_in_acum.append(M_in)
            else:
                M_in_acum.append(M_in + M_in_acum[-1])

            M_over = tQover[t] * cpz[t] * dt * 1
            if t == 0:
                M_over_acum.append(M_over)
            else:
                M_over_acum.append(M_over + M_over_acum[-1])

            Mpipe = tQpipe[t] * c_usz_list[t][m_usz - 1] * dt * 1
            if t == 0:
                Mpipe_acum.append(Mpipe)
            else:
                Mpipe_acum.append(Mpipe + Mpipe_acum[-1])

            Minf_sz = tQinf_sz[t] * c_usz_list[t][m_usz - 1] * dt * 1
            if t == 0:
                Minf_sz_acum.append(Minf_sz)
            else:
                Minf_sz_acum.append(Minf_sz + Minf_sz_acum[-1])

            M_et = tQet[t] * cl_i1[0] * dt * 1
            if t == 0:
                M_et_acum.append(M_et)
            else:
                M_et_acum.append(M_et + M_et_acum[-1])

            M_pz = thpEND[t] * Ab * 1 * cpz[t]
            M_pz_list.append(M_pz)

            Mstor_mb = M_in_acum[-1] - M_pz_list[-1] - M_over_acum[-1] - Mpipe_acum[-1] - Minf_sz_acum[-1] - M_et_acum[
                -1] - Msoil_acum[-1] - MRx_acum[-1]
            Mstor_mb_list.append(Mstor_mb)

            delta_ast_usz = (Mstor_mb - Mstor_ast) / (Ab * tteta_usz[t] * thusz[t] * 1 / m_usz)
            cl_i1[1] = cl_i1[1] + delta_ast_usz
            if cl_i1[1] > 0:
                cl_i1[1] = cl_i1[1]
            else:
                cl_i1[1] = 0

        csoil_usz_list.append(csoil_i_usz)
        c_usz_list.append(cl_i1)
        Rx_usz_list.append(Rxl)

        if hpipe >= 0.03:
            csoil_sz_list.append(csoil_i_sz)
            c_sz_list.append(cj_i1)
            Rx_sz_list.append(Rxj)

    # **5. Transforming in dataframe **
    ## E Coli
    data_usz = pd.DataFrame(c_usz_list)
    a = 0
    column_name = []
    for i in range(m_usz):
        name = 'usz' + str(a)
        column_name.append(name)
        a = a + 1
    data_usz.set_axis(column_name, axis='columns', inplace=True)

    data_sz = pd.DataFrame(c_sz_list)
    a = 0
    column_name = []
    for i in range(m_sz):
        name = 'sz' + str(a)
        column_name.append(name)
        a = a + 1
    data_sz.set_axis(column_name, axis='columns', inplace=True)

    # print('data_usz', data_usz)

    data_soil_usz = pd.DataFrame(csoil_usz_list)
    a = 0
    column_name = []
    for i in range(m_usz):
        name = 'usz_soil' + str(a)
        column_name.append(name)
        a = a + 1
    data_soil_usz.set_axis(column_name, axis='columns', inplace=True)

    data_soil_sz = pd.DataFrame(csoil_sz_list)
    a = 0
    column_name = []
    for i in range(m_sz):
        name = 'sz_soil' + str(a)
        column_name.append(name)
        a = a + 1
    data_soil_sz.set_axis(column_name, axis='columns', inplace=True)

    data_rx_usz = pd.DataFrame(Rx_usz_list)
    a = 0
    column_name = []
    for i in range(m_usz):
        name = 'usz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_usz.set_axis(column_name, axis='columns', inplace=True)

    data_rx_sz = pd.DataFrame(Rx_sz_list)
    a = 0
    column_name = []
    for i in range(m_sz):
        name = 'sz_rx' + str(a)
        column_name.append(name)
        a = a + 1
    data_rx_sz.set_axis(column_name, axis='columns', inplace=True)

    frames = [data_usz, data_sz, data_soil_usz, data_soil_sz, data_rx_usz, data_rx_sz]
    data_EColi = pd.concat((frames), axis=1)

    # ----- teste
    cpz_list_calib = cpz[:]
    cpz_list_calib.append(0)
    data_EColi['cpz'] = cpz_list_calib
    indice_list = indice[:]
    # indice_list.append(len(indice_list)+1)
    cin = cin_list[:]
    cin.append(0)
    c_in = cin[:len(indice_list)]
    # print('len_cp_o2_list', len(cp_o2_list), 'len c_in:', len(c_in), ', len_indice:', len(indice_list))
    data_EColi['c_in'] = c_in

    data_EColi['t'] = indice_list

    Qorif_list = tQpipe[:len(indice_list)]
    # Qorif_list.append(0)
    # print(len(indice), len(Qorif_list))

    data_EColi['Qorif'] = Qorif_list

    if hpipe > 0:
        data_EColi['Morif'] = data_EColi['Qorif'] * data_EColi['sz3'] * 1000  # m3/s
    else:
        data_EColi['Morif'] = data_EColi['Qorif'] * data_EColi['usz10'] * 1000  # m3/s

    # ------
    # cpz.append(0)
    # data_EColi['cpz'] = cpz
    #
    # indice_n = list(range(len(cpz)))
    # cin_list.append(0)
    # c_in = cin_list[:len(indice_n)]
    #
    # data_EColi['c_in'] = c_in
    #
    # data_EColi['t'] = indice_n

    return data_EColi


if __name__ == '__main__':
    print('Starting water quality module...')
    inicio = datetime.datetime.now()

    param = read_from_ini()

    data_EColi = EColi_run(param)

    data_EColi.to_csv('EColi_results_calib.csv', index=False, sep=';', decimal=',')

    fim = datetime.datetime.now()
    print('Elapsed time: ', fim - inicio)
    print('Done!')
