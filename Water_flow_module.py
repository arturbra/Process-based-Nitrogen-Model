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

###  3. Definition of functions

## Flows in m3/s

#Total evapotranspiration
def cQet(sw,sh,ss, Kc, Emax, A, sEST):
    if sEST<=sh:
        Qet = 0.0
    elif sEST<=sw:
        Qet = A*Emax*Kc*(sEST-sh)/(sw-sh)
    elif sEST<=ss:
        Qet = A*Emax*Kc*(sEST-sw)/(ss-sw)
    else:
        Qet = A*Emax*Kc

    Qet = Qet/(dt*1000)

    return Qet

#Ponding zone
# Overflow from weir
def cQover(kWeir,wWeir,hp,Hover,expWeir, Ab, dt, Qin, Qrain):
    Vcheck = hp*Ab + dt*(Qin + Qrain)
    if Vcheck > Hover*Ab:
            Hcheck = Vcheck/Ab
            #Qover = min ((Hcheck-Hover)*Ab/dt, kWeir*wWeir*(2*9.81)**0.5*(Hcheck-Hover)**expWeir)
            Qover = min((Hcheck-Hover)*Ab/dt, kWeir*(Hcheck-Hover)**expWeir)
            #Qover = kWeir*(Hcheck-Hover)**expWeir

    else:
            Qover = 0

    return Qover

# Qinfp: Infiltration from the pond to the surrounding soil
def cQinfp(Kf,Ab,A,Cs,Pp,flagp,hpEST):
    if flagp == 1:
        Qinfp = 0
    else:
        Qinfp = Kf*((Ab-A)+Cs*Pp*hpEST)
    return Qinfp

#Qpf: Infiltration from the pond to the filter material
def cQpf(Ks,hp,husz,A,Ab,dt,s,nusz,Qinfp):
    Qpf = min (Ks*A*(hp+husz)/husz, hp*Ab/dt-Qinfp, (1.0-s)*nusz*husz*A/dt)

    return Qpf

#Unsaturated zone
#Capilary rise
def cQhc(A,ss,sfc,Emax,sEST,Kc):
    s2= sEST #min(s + (Qpf-Qet-Qfs)*dt/(nusz*A*husz), 1)
    den = sfc - ss
    if den == 0:
        den = 0.000001

    Cr = 4*Emax*Kc/(2.5*(den)**2)
    if s2>=ss and s2<=sfc:
        Qhc=A*Cr*(s2-ss)*(sfc-s2)
    else:
        Qhc=0

    b = s2 - USZ.ss
    c = sfc - s2

    #logger.debug(f' ss: {ss} \t sfc: {sfc} \t Emax: {Emax} \t s2: {s2} \t sfc - ss: {den} \t s2 - ss: {b} \t sfc - s2: {c} \t Cr: {Cr} \t Qhc: {Qhc}')


    return Qhc

#Infiltration from USZ to SZ
def cQfs(A, Ks, hp, husz, gama, nusz, dt,sfc, sEST):
        if sEST >= sfc:
            Qfs = min((A*Ks*(hp+husz)/husz)*sEST**gama, (sEST-sfc)*nusz*A*husz/dt)
        else:
            Qfs = 0
        return Qfs



#Saturated zone
#Infiltration to surrounding soil
def cQinfsz(Kf,A,Cs,Psz,flagsz,hszEST):
    if flagsz == 1: #lined
        Qinfsz = 0.0
    else:
        Qinfsz = Kf*(A+Cs*Psz*hszEST)
    return Qinfsz

#Underdrain flow
def cQpipe(hpipe,A,nsz,dt,Qinfsz,Apipe,hszEST,Cd):
    if hszEST<=hpipe:
        Qpipe = 0

    else:
        Qpipemax = (hszEST-hpipe)*A*nsz/dt - Qinfsz
        Qpipemax = max(0,Qpipemax)
        Qpipepossible = Cd*Apipe*((hszEST- hpipe)*2*9.81)**0.5
        Qpipepossible = max(0,Qpipepossible)
        Qpipe = min (Qpipemax, Qpipepossible)

    return Qpipe

#Porosity
def cnsz(hsz,L,Dt,Dg,nf,nt,ng):
    if hsz > Dt+Dg and hsz <= L:
        nsz = ((ng*Dg+nt*Dt+nf*(hsz-Dg-Dt)))/hsz
    elif hsz > Dg and hsz <= Dg+Dt:
        nsz = (ng*Dg+nt*(hsz-Dg))/hsz
    else:
        nsz = ng
    return nsz


def cnusz(husz,hsz,nusz_ini, ng, Dg, Df):
    if hsz < Dg:
        nusz = (nusz_ini*Df + ng*(Dg - hsz))/husz
    else:
        nusz = nusz_ini
    return nusz

#if hsz < Dg, nusz = (nf*Df + ng*(Dg-hsz))/husz

###   4. Model routine

WFR.indice = list(range(0, len(WFR.tQrain)))

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
        Qover = cQover(PZ.Kweir, PZ.wWeir, hpEND, PZ.Hover, PZ.expWeir, PZ.Ab, GENERAL_PARAMETERS.dt, Qin, Qrain)
#         print('Qover: ', Qover)
#         input()

        hp = max(hpEND+GENERAL_PARAMETERS.dt/PZ.Ab*(WFR.tQrain[t]+Qin-Qover),0)   #beginning
#         print('hp:', hp)
#         input()

        Qinfp = cQinfp(USZ.Kf, PZ.Ab, USZ.A, PZ.Cs, PZ.Pp, PZ.flagp, hp)
        Qpf = cQpf(USZ.Ks, hp, husz, USZ.A, PZ.Ab, GENERAL_PARAMETERS.dt, s, nusz, Qinfp)
#         print('Qpf: ', Qpf, ', Ks:', Ks, ', hp:', hp, ', husz: ', husz, ', A: ',A, ', Ab: ', Ab, ', GENERAL_PARAMETERS.dt:', GENERAL_PARAMETERS.dt, ', s: ', s, ', nusz: ', nusz, ', Qinfp: ', Qinfp)
#         input()

        hpEND = max(hp-GENERAL_PARAMETERS.dt/PZ.Ab*(Qpf+Qinfp), 0) #end
#         print('hpEND:', hpEND)
#         input()

        #USZ#
        sEST = max(min(s + Qpf*GENERAL_PARAMETERS.dt/(nusz*USZ.A*husz), 1) , 0)
        Qhc = cQhc(USZ.A, USZ.ss, USZ.sfc, Emax, sEST, GENERAL_PARAMETERS.Kc)
        #print('Qhc: ', Qhc ,', A:', A, ', USZ.ss: ', USZ.ss, ', sfc:', sfc, ', Emax:', Emax, ', sEST:', sEST, ', GENERAL_PARAMETERS.Kc:', GENERAL_PARAMETERS.Kc)
        #input()

        Qfs = cQfs(USZ.A, USZ.Ks, hpEND, husz, USZ.gama, nusz, GENERAL_PARAMETERS.dt, USZ.sfc, sEST)
        #print('Qfs: ', Qfs, ', A:', A, ', Ks:', Ks, ', hpEND: ', hpEND, ', husz: ', husz, ', gama:', gama, ', nusz:', nusz, ', GENERAL_PARAMETERS.dt: ', GENERAL_PARAMETERS.dt, ', sfc: ', sfc, ', sEST: ', sEST)
        #input()

        sEST2 = (sEST*nusz*husz + nsz*hsz)/(nusz*husz + nsz*hsz)

        Qet = cQet(USZ.sw, USZ.sh, USZ.ss, GENERAL_PARAMETERS.Kc, Emax, USZ.A, sEST2)
        
        Qet1 = Qet * (sEST*nusz*husz)/(sEST*nusz*husz + nsz*hsz)
        Qet2 = Qet - Qet1    
        
        #SZ#
        hszEST = hsz+GENERAL_PARAMETERS.dt*(Qfs - Qhc - Qet2)/USZ.A/nsz
        #print('hsz: ', hsz, ', hszEST: ', hszEST)
        Qinfsz = cQinfsz(USZ.Kf, USZ.A, PZ.Cs, SZ.Psz, SZ.flagsz, hszEST)
        Qpipe = cQpipe(GENERAL_PARAMETERS.hpipe, USZ.A, nsz, GENERAL_PARAMETERS.dt, Qinfsz, GENERAL_PARAMETERS.Apipe, hszEST, GENERAL_PARAMETERS.Cd)
        
        hsz = hsz+GENERAL_PARAMETERS.dt*(Qfs - Qhc - Qinfsz- Qpipe - Qet2)/USZ.A/nsz
        #print('hsz: ', hsz, ', Qfs: ', Qfs, ', Qhc: ', Qhc, ', Qinfsz: ', Qinfsz, ', Qpipe: ', Qpipe, ', Qet: ', Qet)
        husz = GENERAL_PARAMETERS.L - hsz
#         print('husz: ', husz)
#         input()

        #porosity#
        nsz = cnsz(hsz, GENERAL_PARAMETERS.L, GENERAL_PARAMETERS.Dt, GENERAL_PARAMETERS.Dg, GENERAL_PARAMETERS.nf, GENERAL_PARAMETERS.nt, GENERAL_PARAMETERS.ng)
        nusz = cnusz(husz, hsz, USZ.nusz_ini, GENERAL_PARAMETERS.ng, GENERAL_PARAMETERS.Dg, GENERAL_PARAMETERS.Df)
        
        if t == 0:
            husz_a = USZ.husz_ini
            nusz_a = USZ.nusz_ini
            s_a = USZ.sw
        else:
            husz_a = WFR.thusz[t-1]
            nusz_a = WFR.tnusz[t-1]
            s_a = WFR.ts[t-1]
            
        s = max(min(1.0, (s_a*husz_a*nusz_a*USZ.A + GENERAL_PARAMETERS.dt*(Qpf + Qhc - Qfs - Qet1))/(USZ.A*husz*nusz)), USZ.sh)
        #print('s:', s, ', s_a: ', s_a, ', husz_a: ', husz_a, ', nusz_a: ', nusz_a, ', A: ', A, ', Qpf: ', Qpf, ', Qhc: ', Qhc, ', Qfs: ', Qfs, ', Qet1: ', Qet1, ', husz: ', husz, ', nusz: ' , nusz)
        #input()

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
    
    ###   6. Saving the results in CSV
    
    dict_data = {
        't': WFR.tt,
        'Qin': WFR.tQin[:len(WFR.tt)],
        'Qet': WFR.tQet[:len(WFR.tt)],
        'hpEND': WFR.thpEND[:len(WFR.tt)],
        'Qpf': WFR.tQpf[:len(WFR.tt)],
        'Qover': WFR.tQover[:len(WFR.tt)],
        'Qfs': WFR.tQfs[:len(WFR.tt)],
        'Qet_1': WFR.tQet1[:len(WFR.tt)],
        'Qhc': WFR.tQhc[:len(WFR.tt)],
        'Qpipe': WFR.tQpipe[:len(WFR.tt)],
        'Qet_2': WFR.tQet2[:len(WFR.tt)],
        'teta_usz': WFR.tteta_usz[:len(WFR.tt)],
        'teta_sz': WFR.tteta_sz[:len(WFR.tt)],
        'Qrain': WFR.tQrain[:len(WFR.tt)],
        'Qinfp': WFR.tQinfp[:len(WFR.tt)],
        'Qinfsz': WFR.tQinfsz[:len(WFR.tt)],
        'hp': WFR.thp[:len(WFR.tt)],
        's': WFR.ts[:len(WFR.tt)],
        'husz': WFR.thusz[:len(WFR.tt)],
        'hsz': WFR.thsz[:len(WFR.tt)],
        'nsz': WFR.tnsz[:len(WFR.tt)],
        'nusz': WFR.tnusz[:len(WFR.tt)],
        'hszEST': WFR.thszEST[:len(WFR.tt)]
             }
    
    
    data = pd.DataFrame(dict_data)
    
#     data[['Qin','Qover','Qpipe', 'Qinfsz']].plot(figsize=(15,8), linewidth=1)
#     plt.show()
      
    #Uncomment this line to generate the .csv with the results
    #data.to_csv('water_flow_results.csv', index = False)
    
    ###  5. Water balance
    
    Qin_total = data['Qin'].sum()
    Vtotal_in = Qin_total * GENERAL_PARAMETERS.dt
    
    Qover_total = data['Qover'].sum()
    Vtotal_over = Qover_total * GENERAL_PARAMETERS.dt
    
#     Qinf_sz_total = data['Qinfsz'].sum()
#     Vtotal_inf_sz = Qinf_sz_total * GENERAL_PARAMETERS.dt
    
    Qpipe_total = data['Qpipe'].sum()
    Vtotal_pipe = Qpipe_total * GENERAL_PARAMETERS.dt
    
    #Vtotal_prec = (P * Ac)/1000
    
    #Vtotal_in_2 = Vtotal_prec*C
    
    Qpeak_over = data['Qover'].max()
    
    Qpf_total = data['Qpf'].sum()
    Vtotal_pf = Qpf_total * GENERAL_PARAMETERS.dt
    
    Qfs_total = data['Qfs'].sum()
    Vtotal_fs = Qfs_total * GENERAL_PARAMETERS.dt
    
    Smax = data['s'].max()
    
    hmax = data['hp'].max()
    
    tpeak = data.loc[data['Qover'] == Qpeak_over, 't'].iloc[0]
    #tpeak = dados[dados['Qv'] == Qpeak_over]['t'].values
    
    print('Vtotal_in :', Vtotal_in) #m3
    print('Vtotal_over :', Vtotal_over) #m3
    print('Vtotal_pipe (m3):', Vtotal_pipe) #m3
    #print('Vtotal_inf_sz (m3):', Vtotal_inf_sz) #m3
    print('Vtotal_pf :', Vtotal_pf) #m3
    print('Vtotal_fs :', Vtotal_fs) #m3
    #print('Vtotal_prec :', Vtotal_prec) #m3
    #print('Vtotal_in_2 :', Vtotal_in_2) #m3
    print('Qpeak_over :', Qpeak_over*1000) #L/s
    print('Smax :', Smax) #m3
    print('hmax :', hmax) #m
    #print('Qmax :', Qmax*1000) #L/s 
    print('tpeak :', tpeak) #min


    fim = datetime.datetime.now()
    print ('Elapsed time: ',fim - inicio)
    wf_test = results_tests.water_flow_comparison_test("results/water_flow_results.csv", WFR)
    print(len(wf_test))
