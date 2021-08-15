__author__ = 'pengfeishen' and 'marinab'

import datetime
import parameters
import configparser
import math
import pandas as pd
import results_tests

WF_INFLOW = parameters.WaterInflow("water_inflow.csv")
WFR = parameters.WaterFlowResults(WF_INFLOW.tQrain, WF_INFLOW.tQin, WF_INFLOW.tEmax)
GENERAL_PARAMETERS = parameters.GeneralParameters("parameters.ini")
PZ = parameters.PondingZone("parameters.ini")
USZ = parameters.UnsaturatedZone('parameters.ini', GENERAL_PARAMETERS.L, GENERAL_PARAMETERS.hpipe, GENERAL_PARAMETERS.dz)
SZ = parameters.SaturatedZone('parameters.ini', GENERAL_PARAMETERS.n, USZ.m_usz)

def run_W():

    hpEND = 0
    sEST = 0
    hszEST = 0
    hp = 0.0
    husz = USZ.husz_ini

    if GENERAL_PARAMETERS.hpipe > 0:
        hsz = GENERAL_PARAMETERS.hpipe
    else:
        hsz = GENERAL_PARAMETERS.dpipe/1000 + 0.03

    nusz = USZ.nusz_ini
    nsz = GENERAL_PARAMETERS.ng
    s = USZ.sw
    for t in range(len(WFR.tQrain)):
        Qin = WFR.tQin[t]
        Qrain = WFR.tQrain[t]
        Emax = WFR.tEmax[t]

        #PZ#
        Qover = PZ.f_weir_overflow(hpEND, GENERAL_PARAMETERS.dt, Qin, Qrain)
        hp = max(hpEND+GENERAL_PARAMETERS.dt/PZ.Ab*(WFR.tQrain[t]+Qin-Qover),0)   #beginning
        Qinfp = PZ.f_infiltration_to_surrounding(USZ.Kf, USZ.A, hp)
        Qpf = PZ.f_infiltration_to_filter_material(USZ.Ks, hp, husz, USZ.A, GENERAL_PARAMETERS.dt, s, nusz, Qinfp)
        hpEND = max(hp-GENERAL_PARAMETERS.dt/PZ.Ab*(Qpf+Qinfp), 0) #end

        #USZ#
        sEST = max(min(s + Qpf*GENERAL_PARAMETERS.dt/(nusz*USZ.A*husz), 1) , 0)
        Qhc = USZ.f_capillary_rise(Emax, sEST, GENERAL_PARAMETERS.Kc)
        Qfs = USZ.f_infiltration_to_sz(hpEND, husz, nusz, GENERAL_PARAMETERS.dt, sEST)
        sEST2 = (sEST*nusz*husz + nsz*hsz)/(nusz*husz + nsz*hsz)
        Qet = PZ.f_evapotranspiration(USZ.sw, USZ.sh, USZ.ss, GENERAL_PARAMETERS.Kc, Emax, USZ.A, sEST2, GENERAL_PARAMETERS.dt)
        Qet1 = Qet * (sEST*nusz*husz)/(sEST*nusz*husz + nsz*hsz)
        Qet2 = Qet - Qet1    
        
        #SZ#
        hszEST = hsz+GENERAL_PARAMETERS.dt*(Qfs - Qhc - Qet2)/USZ.A/nsz
        Qinfsz = SZ.f_infiltration_to_surround(USZ.Kf, USZ.A, PZ.Cs, hszEST)
        Qpipe = SZ.f_underdrain_flow(GENERAL_PARAMETERS.hpipe, USZ.A, nsz, GENERAL_PARAMETERS.dt, Qinfsz, GENERAL_PARAMETERS.Apipe, hszEST, GENERAL_PARAMETERS.Cd)
        hsz = hsz+GENERAL_PARAMETERS.dt*(Qfs - Qhc - Qinfsz- Qpipe - Qet2)/USZ.A/nsz
        husz = GENERAL_PARAMETERS.L - hsz

        #Porosity#
        nsz = SZ.f_porosity(hsz, GENERAL_PARAMETERS.L, GENERAL_PARAMETERS.Dt, GENERAL_PARAMETERS.Dg, GENERAL_PARAMETERS.nf, GENERAL_PARAMETERS.nt, GENERAL_PARAMETERS.ng)
        nusz = USZ.f_porosity(husz, hsz, GENERAL_PARAMETERS.ng, GENERAL_PARAMETERS.Dg, GENERAL_PARAMETERS.Df)
        
        if t == 0:
            husz_a = USZ.husz_ini
            nusz_a = USZ.nusz_ini
            s_a = USZ.sw
        else:
            husz_a = WFR.thusz[t-1]
            nusz_a = WFR.tnusz[t-1]
            s_a = WFR.ts[t-1]
            
        s = max(min(1.0, (s_a*husz_a*nusz_a*USZ.A + GENERAL_PARAMETERS.dt*(Qpf + Qhc - Qfs - Qet1))/(USZ.A*husz*nusz)), USZ.sh)

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


if __name__ == '__main__':
    inicio = datetime.datetime.now()
    run_W()
    #Uncomment this line to generate the .csv with the results
    #data.to_csv('water_flow_results.csv', index = False)
    WFR.water_balance(GENERAL_PARAMETERS.dt)
    fim = datetime.datetime.now()
    print('Elapsed time: ',fim - inicio)
    # wf_test = results_tests.water_flow_comparison_test("results/water_flow_results.csv", WFR)
    # print(len(wf_test))
