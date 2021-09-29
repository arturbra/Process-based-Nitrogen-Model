__author__ = 'pengfeishen' and 'marinab'


import numpy as np 
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
import datetime
import configparser
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import math

###  1. Import input files
## All files should have the same length!!

#Rainfall file
Qrain_file = pd.read_csv('thalita/01_Qrain_file2.csv')
tQrain = Qrain_file['Qrain'].tolist()

#Evapotranspiration file
Emax_file = pd.read_csv('thalita/01_Emax_file2.csv')
tEmax = Emax_file['Emax'].tolist()

#Inflow file
Qin_file = pd.read_csv('thalita/01_Qin_file2.csv')
tQin = Qin_file['Qin'].tolist()


#Rainfall file
# Qrain_file = pd.read_csv('input_file_Qrain.csv')
# tQrain = Qrain_file['Qrain'].tolist()
#
# #Evapotranspiration file
# Emax_file = pd.read_csv('input_file_Emax.csv')
# tEmax = Emax_file['Emax'].tolist()
#
# #Inflow file
# Qin_file = pd.read_csv('input_file_Qin.csv')
# tQin = Qin_file['Qin'].tolist()

###  2. Definition of parameters

setup_file = "thalita/parametros_quanti_RSUP.ini"
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
Ab = float(setup['PONDING_ZONE']['Ab'])
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

# Total evapotranspiration
def cQet(sw, sh, ss, Kc, Emax, A, sEST):
    if sEST <= sh:
        Qet = 0.0
    elif sEST <= sw:
        Qet = A * Emax * Kc * (sEST - sh) / (sw - sh)
    elif sEST <= ss:
        Qet = A * Emax * Kc * (sEST - sw) / (ss - sw)
    else:
        Qet = A * Emax * Kc

    Qet = Qet / (dt * 1000)

    return Qet

# Ponding zone
# Overflow from weir
def cQover(kWeir, wWeir, hp, Hover, expWeir, Ab, dt, Qin, Qrain):
    Vcheck = hp * Ab + dt * (Qin + Qrain)
    if Vcheck > Hover * Ab:
        Hcheck = Vcheck / Ab
        Qover = min((Hcheck - Hover) * Ab / dt, kWeir * wWeir * (2 * 9.81) ** 0.5 * (Hcheck - Hover) ** expWeir)
        # Qover = kWeir*(Hcheck-Hover)**expWeir

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

# Qpf: Infiltration from the pond to the filter material
def cQpf(Ks, hp, husz, A, Ab, dt, s, nusz, Qinfp):
    Qpf = min(Ks * A * (hp + husz) / husz, hp * Ab / dt - Qinfp, (1.0 - s) * nusz * husz * A / dt)

    return Qpf



#Unsaturated zone
# Capilary rise
def cQhc(A, ss, sfc, Emax, sEST, Kc):
    s2 = sEST  # min(s + (Qpf-Qet-Qfs)*dt/(nusz*A*husz), 1)
    den = sfc - ss
    if den == 0:
        den = 0.000001

    Cr = 4 * Emax * Kc / (2.5 * (den) ** 2)
    if s2 >= ss and s2 <= sfc:
        Qhc = A * Cr * (s2 - ss) * (sfc - s2)
    else:
        Qhc = 0

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
def cQinf_sz(Kf,A,Cs,Psz,flagsz,hszEST):
    if flagsz == 1: #lined
        Qinf_sz = 0.0
    else:
        Qinf_sz = Kf*(A+Cs*Psz*hszEST)

    return Qinf_sz

#Underdrain flow
def cQpipe(hpipe,A,nsz,dt,Qinf_sz,Apipe,hszEST, Cd):
    if hszEST <= hpipe:
        Qpipe = 0
    
    else:
        Qpipemax = (hszEST-hpipe)*A*nsz/dt - Qinf_sz
        Qpipemax = max(0, Qpipemax)
        Qpipepossible = Cd*Apipe*((hszEST - hpipe)*2*9.81)**0.5
        Qpipepossible = max(0, Qpipepossible)
        Qpipe = min(Qpipemax, Qpipepossible)

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


###  4. Initial values

tt = []
tQover = []
tQpf = []
tQinfp = []
tQfs = []
tQhc = []
tQet = []
tQinf_sz = []
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

###   5. Model routine

def wf_run():
    hpEND = 0
    sEST = 0
    hszEST = 0

    hp = 0.0
    husz = husz_ini


    if hpipe > 0:
        hsz = hpipe
    else:
        #hsz = dpipe / 1000 + 0.03
        hsz = 0

    nusz = nusz_ini
    nsz = ng
    s = sw

    for t in range(len(tQrain)):
        Qin = tQin[t]
        Qrain = tQrain[t]
        Emax = tEmax[t]

        #PONDING ZONE#
        Qover = cQover(Kweir, wWeir, hpEND, Hover, expWeir, Ab, dt, Qin, Qrain)

        hp = max(hpEND + (dt / Ab) * (tQrain[t] + Qin - Qover), 0) #beginning

        Qinfp = cQinfp(Kf, Ab, A, Cs, Pp, flagp, hp)
        Qpf = cQpf(Ks, hp, husz, A, Ab, dt, s, nusz, Qinfp)

        hpEND = max(hp - dt / Ab * (Qpf + Qinfp), 0)  # end


        #UNSATURATED ZONE#
        sEST = max(min(s + Qpf * dt / (nusz * A * husz), 1), 0)

        if hpipe <= 0.03:
            Qhc = 0
        else:
            Qhc = cQhc(A, ss, sfc, Emax, sEST, Kc)

        Qfs = cQfs(A, Ks, hpEND, husz, gama, nusz, dt, sfc, sEST)

        sEST2 = (sEST * nusz * husz + nsz * hsz) / (nusz * husz + nsz * hsz)

        Qet = cQet(sw, sh, ss, Kc, Emax, A, sEST2)

        Qet1 = Qet * (sEST * nusz * husz) / (sEST * nusz * husz + nsz * hsz)
        Qet2 = Qet - Qet1
        if hpipe > 0.03:                                                ###### acrescentei isso
            Qet2 = Qet - Qet1
        else:
            Qet2 = 0


        #SATURATED ZONE#
        hszEST = hsz + dt * (Qfs - Qhc - Qet2) / A / nsz

        Qinf_sz = cQinf_sz(Kf, A, Cs, Psz, flagsz, hszEST)
        Qpipe = cQpipe(hpipe, A, nsz, dt, Qinf_sz, Apipe, hszEST, Cd)

        hsz = hsz + dt * (Qfs - Qhc - Qinf_sz - Qpipe - Qet2) / A / nsz
        if hsz < 0.0001:                                                   ####acrescentei isso
            hsz = 0

        husz = L - hsz
    
    
        #porosity#
        nsz = cnsz(hsz, L, Dt, Dg, nf, nt, ng)
        nusz = cnusz(husz, hsz, nusz_ini, ng, Dg, Df)

        if t == 0:
            husz_a = husz_ini
            nusz_a = nusz_ini
            s_a = sw
        else:
            husz_a = thusz[t - 1]
            nusz_a = tnusz[t - 1]
            s_a = ts[t - 1]

        s = max(min(1.0, (s_a * husz_a * nusz_a * A + dt * (Qpf + Qhc - Qfs - Qet1)) / (A * husz * nusz)), sh)

        #save all results
        tt.append(t)
        tQover.append(Qover)
        tQpf.append(Qpf)
        tQinfp.append(Qinfp)
        tQfs.append(Qfs)
        tQhc.append(Qhc)
        tQet.append(Qet)
        tQinf_sz.append(Qinf_sz)
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
        tteta_usz.append(s * nusz)
        tteta_sz.append(nsz)



if __name__ == '__main__':

    print('Starting..')
    inicio = datetime.datetime.now()

    wf_run()

    ###   6. Saving the results in CSV

    dict_data = {
        't': tt,
        'Qrain': tQrain[:len(tt)],
        'Qin': tQin[:len(tt)],
        'Qet': tQet[:len(tt)],
        'Qet_1': tQet1[:len(tt)],
        'Qet_2': tQet2[:len(tt)],
        'hpEND': thpEND[:len(tt)],
        'Qpf': tQpf[:len(tt)],
        'Qover': tQover[:len(tt)],
        'Qfs': tQfs[:len(tt)],
        'Qhc': tQhc[:len(tt)],
        'Qpipe': tQpipe[:len(tt)],
        'teta_usz': tteta_usz[:len(tt)],
        'teta_sz': tteta_sz[:len(tt)],
        'Qinfp': tQinfp[:len(tt)],
        'Qinf_sz': tQinf_sz[:len(tt)],
        'hp': thp[:len(tt)],
        's': ts[:len(tt)],
        'husz': thusz[:len(tt)],
        'hsz': thsz[:len(tt)],
        'nsz': tnsz[:len(tt)],
        'nusz': tnusz[:len(tt)],
        'hszEST': thszEST[:len(tt)]
        }
    

    data = pd.DataFrame(dict_data)

    data[['Qin', 'Qover', 'Qpipe', 'Qinf_sz']].plot(figsize=(9, 5), linewidth=1)
    plt.show()
      
    #data.to_csv('water_flow_model_results_valid.csv', index=False, sep=';', decimal=',')

    ###  7. Water balance
    
    Qin_total = data['Qin'].sum()
    Vtotal_in = Qin_total * dt * 1000
    
    Qover_total = data['Qover'].sum()
    Vtotal_over = Qover_total * dt * 1000

    Qet_total = data['Qet'].sum()
    Vtotal_et = Qet_total*dt*1000
    
    Qpipe_total = data['Qpipe'].sum()
    Vtotal_pipe = Qpipe_total * dt*1000
    
    # Qinf_sz_total = data['Qinf_sz'].sum()
    # Vtotal_inf_sz = Qinf_sz_total * dt

    # Vtotal_prec = (P * Ac)/1000
    # Vtotal_in_2 = Vtotal_prec*C

    Qpeak_over = data['Qover'].max()

    Qpf_total = data['Qpf'].sum()
    Vtotal_pf = Qpf_total * dt

    Qfs_total = data['Qfs'].sum()
    Vtotal_fs = Qfs_total * dt

    Smax = data['s'].max()
    
    hmax = data['hp'].max()
    
    tpeak = data.loc[data['Qover'] == Qpeak_over, 't'].iloc[0]


    # calculo do balanço de massa teorico

    V_storage_teorico = (Vtotal_in - Vtotal_pipe - Vtotal_et - Vtotal_over) #litros


    #calculo do balanço de massa real

    n = 11  # number of cells
    dz = L / n
    m_usz = round((L - hpipe) / dz)
    m_sz = n - m_usz

    #usz
    Vlayer_usz = (Ab * (L - hpipe)) / m_usz              # mudei isso, antes era: Vlayer_usz = (Ab*(Df))/m_usz

    usz_tot_vol_inicial = Vlayer_usz*tteta_usz[0]*m_usz
    total_linhas_usz = len(tteta_usz) - 1
    usz_tot_vol_final = Vlayer_usz*tteta_usz[total_linhas_usz]*m_usz

    # sz
    if hpipe > 0:

        Vlayer_sz = (Ab * hpipe) / m_sz                  # mudei isso, antes era: Vlayer_sz = (Ab*(Dg))/m_sz

        sz_tot_vol_inicial = Vlayer_sz * tteta_sz[0] * m_sz
        total_linhas_sz = len(tteta_sz) - 1
        sz_tot_vol_final = Vlayer_sz * tteta_sz[total_linhas_sz] * m_sz

        V_storage_real = (usz_tot_vol_final + sz_tot_vol_final - usz_tot_vol_inicial - sz_tot_vol_inicial) * 1000

    else:
        V_storage_real = (usz_tot_vol_final - usz_tot_vol_inicial) * 1000


    verificacao_balanco = (V_storage_teorico - V_storage_real)


    print('.')
    print("\033[1:30m........ Informações ........ \033[m")
    print('Vtotal_in: {:.2f} L'.format(Vtotal_in)) #m3
    print('Vtotal_pipe: {:.2f} L'.format(Vtotal_pipe)) #m3
    print('Vtotal_et: {:.2f} L'.format(Vtotal_et))  # m3
    print('.'*30)
    print('.')
    print("\033[1:30m...... Balanço de massa ...... \033[m")
    print('V_storage_teorico: {:.2f} L'.format(V_storage_teorico))
    print('V_storage_real: {:.2f} L'.format(V_storage_real))
    print('Verificação do balanco: {:.2f} L'.format(verificacao_balanco))


    porcentagem = ((verificacao_balanco*100)/V_storage_teorico)

    if porcentagem <= 1 and porcentagem >= - 1:
        print("\033[1:34mBalanço de massa aceitável! {:.2f} % de erro \033[m".format(porcentagem))
    else:
        print("\033[1:31mBalanço de massa não aceitável: {:.2f} % de erro\033[m".format(porcentagem))





    #print('Vtotal_over :', Vtotal_over) #m3
    #print('Vtotal_inf_sz (m3):', Vtotal_inf_sz) #m3
    #print('Vtotal_pf :', Vtotal_pf) #m3
    #print('Vtotal_fs :', Vtotal_fs) #m3
    #print('Vtotal_prec :', Vtotal_prec) #m3
    #print('Vtotal_in_2 :', Vtotal_in_2) #m3
    #print('Qpeak_over :', Qpeak_over*1000) #L/s
    #print('Smax :', Smax) #m3
    #print('hmax :', hmax) #m
    #print('Qmax :', Qmax*1000) #L/s
    #print('tpeak :', tpeak) #min

    print('.' * 30)
    print('.')
    fim = datetime.datetime.now()
    print('Elapsed time: ', fim - inicio)
    print('Done!')