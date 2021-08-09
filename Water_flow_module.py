__author__ = 'pengfeishen' and 'marinab'

import datetime
import numpy as np 
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
import configparser
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import math
import logging

rodada = 'teste_1'

# logging.basicConfig(filename = 'teste_Qhc_pf.txt', level=logging.DEBUG, format=f'{rodada}: - %(asctime)s - %(levelname)s - %(message)s')
# logger = logging.getLogger(__name__)


###  1. Import input files
## All files should have the same length!!

#Rainfall file
Qrain_file = pd.read_csv('Qrain_file.csv')
tQrain = Qrain_file['Qrain'].tolist()

#Evapotranspiration file  ## verificar unidade
Emax_file = pd.read_csv('ET.csv')
tEmax = Emax_file['ET'].tolist()

#Inflow file
Qin_file = pd.read_csv('inflow.csv')
tQin = Qin_file['Qin'].tolist()


###  2. Definition of parameters

setup_file = "Water_flow_module.ini"
setup = configparser.ConfigParser()
setup.read(setup_file)

#General
#Ew = float(setup['GENERAL']['Ew'])
Kc = float(setup['GENERAL']['Kc']) #verificar esse valor antes de rodar
Df = float(setup['GENERAL']['Df'])
Dw = float(setup['GENERAL']['Dw'])
Dt = float(setup['GENERAL']['Dt'])
Dg = float(setup['GENERAL']['Dg'])
L = float(setup['GENERAL']['L'])
nf = float(setup['GENERAL']['nf'])
nw = float(setup['GENERAL']['nw'])
nt = float(setup['GENERAL']['nt'])
ng = float(setup['GENERAL']['ng'])
nn = float(setup['GENERAL']['nn'])
rwv = float(setup['GENERAL']['rwv'])

#Ponding zone
Ab = float(setup['PONDING_ZONE']['Ap'])
Hover = float(setup['PONDING_ZONE']['Hover'])
Kweir = float(setup['PONDING_ZONE']['Kweir'])
wWeir = float(setup['PONDING_ZONE']['wWeir'])
expWeir = float(setup['PONDING_ZONE']['expWeir'])
Cs = float(setup['PONDING_ZONE']['Cs'])
Pp = float(setup['PONDING_ZONE']['Pp'])
flagp = float(setup['PONDING_ZONE']['flagp'])

#Unsaturated zone
A = float(setup['UNSATURATED_ZONE']['A'])
husz_ini = float(setup['UNSATURATED_ZONE']['husz'])
nusz_ini = float(setup['UNSATURATED_ZONE']['nusz'])
Ks = float(setup['UNSATURATED_ZONE']['Ks'])
sh = float(setup['UNSATURATED_ZONE']['sh'])
sw = float(setup['UNSATURATED_ZONE']['sw'])
sfc = float(setup['UNSATURATED_ZONE']['sfc'])
ss = float(setup['UNSATURATED_ZONE']['ss'])
gama = float(setup['UNSATURATED_ZONE']['gama'])
Kf = float(setup['UNSATURATED_ZONE']['Kf'])

#Saturated zone
Psz = float(setup['SATURATED_ZONE']['Psz'])
hpipe = float(setup['SATURATED_ZONE']['hpipe'])
flagsz = float(setup['SATURATED_ZONE']['flagsz'])
dpipe = float(setup['SATURATED_ZONE']['dpipe'])
Cd = float(setup['SATURATED_ZONE']['Cd'])
Apipe = math.pi*(dpipe/(1000*2))**2
 

#Timestep
dt = float(setup['TIMESTEP']['dt'])


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
    
    b = s2 - ss
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

tt = []
tQover = []
tQpf = []
tQinfp = []
tQfs = []
tQhc = []
tQet = []
tQinfsz = []
tQpipe = []
tQet1 = []
tQet2 = []

thp = []
ts = []
thsz = []
thszEST = []
thusz =[]
tnsz = []
tnusz = []
thpEND = []
tteta_usz = []
tteta_sz = []

indice = list(range(0, len(tQrain)))

def run_W():
    
    hpEND = 0
    sEST = 0
    hszEST = 0

    hp = 0.0
    husz = husz_ini
    if hpipe > 0:
        hsz = hpipe
    else:
        hsz = dpipe/1000 + 0.03
    nusz = nusz_ini
    nsz = ng
    s = sw

    for t in range(len(tQrain)):
        Qin = tQin[t]
        Qrain = tQrain[t]
        Emax = tEmax[t]
        
        #PZ#
        Qover = cQover(Kweir, wWeir, hpEND, Hover, expWeir, Ab, dt, Qin, Qrain)
#         print('Qover: ', Qover)
#         input()
        
        hp = max(hpEND+dt/Ab*(tQrain[t]+Qin-Qover),0)   #beginning
#         print('hp:', hp)
#         input()
        
        Qinfp = cQinfp(Kf, Ab, A, Cs, Pp, flagp, hp)
        Qpf = cQpf(Ks, hp, husz, A, Ab, dt, s, nusz, Qinfp)
#         print('Qpf: ', Qpf, ', Ks:', Ks, ', hp:', hp, ', husz: ', husz, ', A: ',A, ', Ab: ', Ab, ', dt:', dt, ', s: ', s, ', nusz: ', nusz, ', Qinfp: ', Qinfp)
#         input()
        
        hpEND = max(hp-dt/Ab*(Qpf+Qinfp), 0) #end
#         print('hpEND:', hpEND)
#         input()
            
        #USZ#
        sEST = max(min(s + Qpf*dt/(nusz*A*husz), 1) , 0)
        Qhc = cQhc(A, ss, sfc, Emax, sEST, Kc)
        #print('Qhc: ', Qhc ,', A:', A, ', ss: ', ss, ', sfc:', sfc, ', Emax:', Emax, ', sEST:', sEST, ', Kc:', Kc)
        #input()
        
        Qfs = cQfs(A, Ks, hpEND, husz, gama, nusz, dt, sfc, sEST)
        #print('Qfs: ', Qfs, ', A:', A, ', Ks:', Ks, ', hpEND: ', hpEND, ', husz: ', husz, ', gama:', gama, ', nusz:', nusz, ', dt: ', dt, ', sfc: ', sfc, ', sEST: ', sEST)
        #input()
        
        sEST2 = (sEST*nusz*husz + nsz*hsz)/(nusz*husz + nsz*hsz)
        
        Qet = cQet(sw, sh, ss, Kc, Emax, A, sEST2)
        
        Qet1 = Qet * (sEST*nusz*husz)/(sEST*nusz*husz + nsz*hsz)
        Qet2 = Qet - Qet1    
        
        #SZ#
        hszEST = hsz+dt*(Qfs - Qhc - Qet2)/A/nsz
        #print('hsz: ', hsz, ', hszEST: ', hszEST)
        Qinfsz = cQinfsz(Kf, A, Cs, Psz, flagsz, hszEST)
        Qpipe = cQpipe(hpipe, A, nsz, dt, Qinfsz, Apipe, hszEST, Cd)
        
        hsz = hsz+dt*(Qfs - Qhc - Qinfsz- Qpipe - Qet2)/A/nsz  
        #print('hsz: ', hsz, ', Qfs: ', Qfs, ', Qhc: ', Qhc, ', Qinfsz: ', Qinfsz, ', Qpipe: ', Qpipe, ', Qet: ', Qet)
        husz = L - hsz
#         print('husz: ', husz)
#         input()

        #porosity#
        nsz = cnsz(hsz, L, Dt, Dg, nf, nt, ng)
        nusz = cnusz(husz, hsz, nusz_ini, ng, Dg, Df)
        
        if t == 0:
            husz_a = husz_ini
            nusz_a = nusz_ini
            s_a = sw
        else:
            husz_a = thusz[t-1]
            nusz_a = tnusz[t-1]
            s_a = ts[t-1]
            
        s = max(min(1.0, (s_a*husz_a*nusz_a*A + dt*(Qpf + Qhc - Qfs - Qet1))/(A*husz*nusz)), sh)
        #print('s:', s, ', s_a: ', s_a, ', husz_a: ', husz_a, ', nusz_a: ', nusz_a, ', A: ', A, ', Qpf: ', Qpf, ', Qhc: ', Qhc, ', Qfs: ', Qfs, ', Qet1: ', Qet1, ', husz: ', husz, ', nusz: ' , nusz)
        #input()
        
        #save all results
        tt.append(t)
        tQover.append(Qover)
        tQpf.append(Qpf)
        tQinfp.append(Qinfp)
        tQfs.append(Qfs)
        tQhc.append(Qhc)
        tQet.append(Qet)
        tQinfsz.append(Qinfsz)
        tQpipe.append(Qpipe)
        tQet1.append(Qet1)
        tQet2.append(Qet2)
    
        thp.append(hp)
        ts.append(s)
        thusz.append(husz)
        thsz.append(hsz)
        thszEST.append(hszEST)
        tnsz.append(nsz)
        tnusz.append(nusz)
        thpEND.append(hpEND)
        tteta_usz.append(s*nusz)
        tteta_sz.append(nsz)




if __name__ == '__main__':
    inicio = datetime.datetime.now()
    
    run_W()         
    
    ###   6. Saving the results in CSV
    
    dict_data = {
        't': tt,
        'Qin': tQin[:len(tt)],
        'Qet': tQet[:len(tt)],
        'hpEND': thpEND[:len(tt)],
        'Qpf': tQpf[:len(tt)],
        'Qover': tQover[:len(tt)],
        'Qfs': tQfs[:len(tt)],
        'Qet_1': tQet1[:len(tt)],
        'Qhc': tQhc[:len(tt)],
        'Qpipe': tQpipe[:len(tt)],
        'Qet_2': tQet2[:len(tt)],
        'teta_usz': tteta_usz[:len(tt)],
        'teta_sz': tteta_sz[:len(tt)],
        'Qrain': tQrain[:len(tt)],      
        'Qinfp': tQinfp[:len(tt)],        
        'Qinfsz': tQinfsz[:len(tt)],
        'hp': thp[:len(tt)],
        's': ts[:len(tt)],
        'husz': thusz[:len(tt)],
        'hsz': thsz[:len(tt)],
        'nsz': tnsz[:len(tt)],
        'nusz': tnusz[:len(tt)],
        'hszEST' : thszEST[:len(tt)]
             }
    
    
    data = pd.DataFrame(dict_data)
    
#     data[['Qin','Qover','Qpipe', 'Qinfsz']].plot(figsize=(15,8), linewidth=1)
#     plt.show()
      
    
    data.to_csv('water_flow_results.csv', index = False)
    
    ###  5. Water balance
    
    Qin_total = data['Qin'].sum()
    Vtotal_in = Qin_total * dt 
    
    Qover_total = data['Qover'].sum()
    Vtotal_over = Qover_total * dt 
    
#     Qinf_sz_total = data['Qinfsz'].sum()
#     Vtotal_inf_sz = Qinf_sz_total * dt 
    
    Qpipe_total = data['Qpipe'].sum()
    Vtotal_pipe = Qpipe_total * dt
    
    #Vtotal_prec = (P * Ac)/1000
    
    #Vtotal_in_2 = Vtotal_prec*C
    
    Qpeak_over = data['Qover'].max()
    
    Qpf_total = data['Qpf'].sum()
    Vtotal_pf = Qpf_total * dt
    
    Qfs_total = data['Qfs'].sum()
    Vtotal_fs = Qfs_total * dt  
    
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
    
