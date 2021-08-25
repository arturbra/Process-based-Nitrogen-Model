import configparser
from math import pi
import pandas as pd




class WaterInflow:
    def __init__(self, file):
        csv_file = pd.read_csv(file, sep=";")
        self.tQin = csv_file["Qin"]
        self.tEmax = csv_file["ET"]
        self.tQrain = csv_file["Qrain"]


class WaterFlowResults:
    def __init__(self, tQrain, tQin, tEmax):
        self.indice = list(range(0, len(tQrain)))
        self.tt = []
        self.tQover = []
        self.tQpf = []
        self.tQinfp = []
        self.tQfs = []
        self.tQhc = []
        self.tQet = []
        self.tQinfsz = []
        self.tQpipe = []
        self.tQet1 = []
        self.tQet2 = []
        self.thp = []
        self.ts = []
        self.thusz = []
        self.thsz = []
        self.thszEST = []
        self.tnsz = []
        self.tnusz = []
        self.thpEND = []
        self.tteta_usz = []
        self.tteta_sz = []
        self.tQin = tQin
        self.tQrain = tQrain
        self.tEmax = tEmax
        self.indice = list(range(0, len(self.tQrain)))

    def water_balance(self, dt):
        Qin_total = self.tQin.sum()
        Qover_total = pd.array(self.tQover).sum()
        Qpipe_total = pd.array(self.tQpipe).sum()
        Qpeak_over = pd.array(self.tQover).max()
        Qpf_total = pd.array(self.tQpf).sum()
        Qfs_total = pd.array(self.tQfs).sum()

        Vtotal_in = Qin_total * dt
        Vtotal_over = Qover_total * dt
        Vtotal_pipe = Qpipe_total * dt
        Vtotal_pf = Qpf_total * dt
        Vtotal_fs = Qfs_total * dt

        Smax = pd.array(self.ts).max()
        hmax = pd.array(self.thp).max()

        print('Vtotal_in :', Vtotal_in)  # m3
        print('Vtotal_over :', Vtotal_over)  # m3
        print('Vtotal_pipe (m3):', Vtotal_pipe)  # m3
        print('Vtotal_pf :', Vtotal_pf)  # m3
        print('Vtotal_fs :', Vtotal_fs)  # m3
        print('Qpeak_over :', Qpeak_over * 1000)  # L/s
        print('Smax :', Smax)  # m3
        print('hmax :', hmax)  # m


class GeneralParameters:
    def __init__(self, setup_file):
        setup = configparser.ConfigParser()
        setup.read(setup_file)
        self.Ew = float(setup['GENERAL']['Ew'])
        self.Kc = float(setup['GENERAL']['Kc'])  # verificar esse valor antes de rodar
        self.Df = float(setup['GENERAL']['Df'])
        self.Dw = float(setup['GENERAL']['Dw'])
        self.Dt = float(setup['GENERAL']['Dt'])
        self.Dg = float(setup['GENERAL']['Dg'])
        self.L = float(setup['GENERAL']['L'])
        self.nf = float(setup['GENERAL']['nf'])
        self.nw = float(setup['GENERAL']['nw'])
        self.nt = float(setup['GENERAL']['nt'])
        self.ng = float(setup['GENERAL']['ng'])
        self.nn = float(setup['GENERAL']['nn'])
        self.rwv = float(setup['GENERAL']['rwv'])
        self.dt = float(setup['TIMESTEP']['dt'])
        self.n = int(setup['MODEL']['n'])
        self.dz = self.L / self.n
        self.d50 = float(setup['MODEL']['d50'])
        self.k_nit = float(setup['NITRIFICATION']['k_nit'])
        self.k_denit = float(setup['DENITRIFICATION']['k_denit'])
        self.hpipe = float(setup['SATURATED_ZONE']['hpipe'])
        self.dpipe = float(setup['SATURATED_ZONE']['dpipe'])
        self.Cd = float(setup['SATURATED_ZONE']['Cd'])
        self.Apipe = pi * (self.dpipe / (1000 * 2)) ** 2


class PondingZone:
    def __init__(self, setup_file):
        setup = configparser.ConfigParser()
        setup.read(setup_file)
        self.Ab = float(setup['PONDING_ZONE']['Ap'])
        self.Hover = float(setup['PONDING_ZONE']['Hover'])
        self.kWeir = float(setup['PONDING_ZONE']['Kweir'])
        self.wWeir = float(setup['PONDING_ZONE']['wWeir'])
        self.expWeir = float(setup['PONDING_ZONE']['expWeir'])
        self.Cs = float(setup['PONDING_ZONE']['Cs'])
        self.Pp = float(setup['PONDING_ZONE']['Pp'])
        self.flagp = float(setup['PONDING_ZONE']['flagp'])
        self.k_denit = float(setup['DENITRIFICATION']['k_denit_pz'])

    def f_concentration(self, cin, Qin_p, cp_a, Qpf, Qv, Rxi, height_pz, height_pz_before, Ab, dt):
        # delta_cp = ((cin*Qin_p - cp_a*(Qpf + Qv))*GENERAL_PARAMETERS.dt)/(height_pz*PZ.Ab) + Rxi*GENERAL_PARAMETERS.dt
        # cp = cp_a + delta_cp

        # cp = Rxi*GENERAL_PARAMETERS.dt + (cin*Qin_p*GENERAL_PARAMETERS.dt)/(height_pz*PZ.Ab + GENERAL_PARAMETERS.dt*(Qpf + Qv))

        cp = (cp_a * height_pz_before * Ab + (cin * Qin_p - cp_a * (Qpf + Qv) + Rxi * height_pz * Ab) * dt) / (
                    height_pz * Ab)

        return cp

    def f_evapotranspiration(self, sw, sh, ss, Kc, Emax, A, sEST, dt):
        if sEST <= sh:
            evapotranspiration = 0.0
        elif sEST <= sw:
            evapotranspiration = A * Emax * Kc * (sEST - sh) / (sw - sh)
        elif sEST <= ss:
            evapotranspiration = A * Emax * Kc * (sEST - sw) / (ss - sw)
        else:
            evapotranspiration = A * Emax * Kc

        Qet = evapotranspiration / (dt * 1000)
        return evapotranspiration

    def f_weir_overflow(self, height_pz, dt, Qin, Qrain):
        volume = height_pz * self.Ab + dt * (Qin + Qrain)
        if volume > self.Hover * self.Ab:
            height = volume / self.Ab
            weir_overflow = min((height - self.Hover) * self.Ab / dt, self.kWeir * (height - self.Hover) ** self.expWeir)
        else:
            weir_overflow = 0
        return weir_overflow

    def f_infiltration_to_surrounding(self, Kf, A, hpEST):
        if self.flagp == 1:
            infiltration = 0
        else:
            infiltration = Kf * ((self.Ab - A) + self.Cs * self.Pp * hpEST)
        return infiltration

    def f_infiltration_to_filter_material(self, Ks, height_pz, husz, A, dt, s, nusz, Qinfp):
        infiltration = min(Ks * A * (height_pz + husz) / husz, height_pz * self.Ab / dt - Qinfp, (1.0 - s) * nusz * husz * A / dt)
        return infiltration


class UnsaturatedZone:
    def __init__(self, setup_file, L, hpipe, dz):
        setup = configparser.ConfigParser()
        setup.read(setup_file)
        self.A = float(setup['UNSATURATED_ZONE']['A'])
        self.husz_ini = float(setup['UNSATURATED_ZONE']['husz'])
        self.nusz_ini = float(setup['UNSATURATED_ZONE']['nusz'])
        self.Ks = float(setup['UNSATURATED_ZONE']['Ks'])
        self.sh = float(setup['UNSATURATED_ZONE']['sh'])
        self.sw = float(setup['UNSATURATED_ZONE']['sw'])
        self.sfc = float(setup['UNSATURATED_ZONE']['sfc'])
        self.ss = float(setup['UNSATURATED_ZONE']['ss'])
        self.gama = float(setup['UNSATURATED_ZONE']['gama'])
        self.Kf = float(setup['UNSATURATED_ZONE']['Kf'])
        self.m_usz = round((L - hpipe) / dz)

    def f_capillary_rise(self, Emax, sEST, Kc):
        s2 = sEST
        den = self.sfc - self.ss
        if den == 0:
            den = 0.000001
        Cr = 4 * Emax * Kc / (2.5 * (den) ** 2)
        if s2 >= self.ss and s2 <= self.sfc:
            capillary_rise = self.A * Cr * (s2 - self.ss) * (self.sfc - s2)
        else:
            capillary_rise = 0
        return capillary_rise

    def f_infiltration_to_sz(self, height_pz, husz, nusz, dt, sEST):
        if sEST >= self.sfc:
            Qfs = min((self.A * self.Ks * (height_pz + husz) / husz) * sEST ** self.gama, (sEST - self.sfc) * nusz * self.A * husz / dt)
        else:
            Qfs = 0
        return Qfs

    def f_porosity(self, husz, hsz, ng, Dg, Df):
        if hsz < Dg:
            nusz = (self.nusz_ini * Df + ng * (Dg - hsz)) / husz
        else:
            nusz = self.nusz_ini
        return nusz

    def f_alfa_beta(self, l):
        alfa = (self.m_usz - 1 - l) / (self.m_usz - 1)
        beta = l / (self.m_usz - 1)
        return alfa, beta

    def f_unit_flux(self, alfa, beta, Qpf, Qet_1, Qfs, Qhc, Qorif, Qinf_sz, teta_sm_i, hpipe, Ab):
        if hpipe > 0:
            UF_usz = (alfa * (Qpf - Qet_1) + beta * (Qfs - Qhc)) / (Ab * teta_sm_i)

        else:
            UF_usz = (alfa * (Qpf - Qet_1) + beta * (Qorif + Qinf_sz - Qhc)) / (Ab * teta_sm_i)

        return UF_usz

    def f_peclet(self, UF_usz, D, dz):
        if D > 0:
            Peusz = UF_usz * dz / D
        else:
            Peusz = 100

        return Peusz


class SaturatedZone:
    def __init__(self, setup_file, n, m_usz):
        setup = configparser.ConfigParser()
        setup.read(setup_file)
        self.Psz = float(setup['SATURATED_ZONE']['Psz'])
        self.flagsz = float(setup['SATURATED_ZONE']['flagsz'])
        self.m_sz = n - m_usz

    def f_infiltration_to_surround(self, Kf, A, Cs, hszEST):
        if self.flagsz == 1:  # lined
            infiltration_to_surround = 0.0
        else:
            infiltration_to_surround = Kf * (A + Cs * self.Psz * hszEST)
        return infiltration_to_surround

    def f_underdrain_flow(self, hpipe, A, nsz, dt, Qinfsz, Apipe, hszEST, Cd):
        if hszEST <= hpipe:
            underdrain_flow = 0
        else:
            underdrain_flow_max = (hszEST - hpipe) * A * nsz / dt - Qinfsz
            underdrain_flow_max = max(0, underdrain_flow_max)
            underdrain_flow_possible = Cd * Apipe * ((hszEST - hpipe) * 2 * 9.81) ** 0.5
            underdrain_flow_possible = max(0, underdrain_flow_possible)
            underdrain_flow = min(underdrain_flow_max, underdrain_flow_possible)

        return underdrain_flow

    def f_porosity(self, hsz, L, Dt, Dg, nf, nt, ng):
        if hsz > Dt + Dg and hsz <= L:
            nsz = ((ng * Dg + nt * Dt + nf * (hsz - Dg - Dt))) / hsz
        elif hsz > Dg and hsz <= Dg + Dt:
            nsz = (ng * Dg + nt * (hsz - Dg)) / hsz
        else:
            nsz = ng
        return nsz

    def f_alfa_beta(self, j):
        alfa2 = (self.m_sz - 1 - j) / (self.m_sz - 1)
        beta2 = j / (self.m_sz - 1)
        return alfa2, beta2

    def f_unit_flux(self, alfa2, beta2, Qfs, Qhc, Qet_2, Qorif, Qinf_sz, teta_b_i, Ab):
        UF_sz = (alfa2 * (Qfs - Qhc - Qet_2) + beta2 * (Qorif + Qinf_sz)) / (Ab * teta_b_i)

        return UF_sz

    def f_peclet(self, UF_sz, D, dz):
        if D > 0:
            Pesz = UF_sz * dz / D
        else:
            Pesz = 100

        return Pesz


class Calibration:
    def __init__(self, setup_file):
        setup = configparser.ConfigParser()
        setup.read(setup_file)
        self.Kc_min = float(setup['FLOW_CALIBRATION']['Kc_min'])
        self.Kc_max = float(setup['FLOW_CALIBRATION']['Kc_max'])
        self.Ks_min = float(setup['FLOW_CALIBRATION']['Ks_min'])
        self.Ks_max = float(setup['FLOW_CALIBRATION']['Ks_max'])
        self.sh_min = float(setup['FLOW_CALIBRATION']['sh_min'])
        self.sh_max = float(setup['FLOW_CALIBRATION']['sh_max'])
        self.sw_min = float(setup['FLOW_CALIBRATION']['sw_min'])
        self.sw_max = float(setup['FLOW_CALIBRATION']['sw_max'])
        self.sfc_min = float(setup['FLOW_CALIBRATION']['sfc_min'])
        self.sfc_max = float(setup['FLOW_CALIBRATION']['sfc_max'])
        self.ss_min = float(setup['FLOW_CALIBRATION']['ss_min'])
        self.ss_max = float(setup['FLOW_CALIBRATION']['ss_max'])
        self.Kf_min = float(setup['FLOW_CALIBRATION']['Kf_min'])
        self.Kf_max = float(setup['FLOW_CALIBRATION']['Kf_max'])
        self.Cd_min = float(setup['FLOW_CALIBRATION']['Cd_min'])
        self.Cd_max = float(setup['FLOW_CALIBRATION']['Cd_max'])
        self.pop_init = int(setup['FLOW_CALIBRATION']['pop'])
        self.gen_init = int(setup['FLOW_CALIBRATION']['gen'])
        self.obs_file = setup['FLOW_CALIBRATION']['obs_file']
        self.lamta_min = float(setup['QUALITY_CALIBRATION']['lamta_min'])
        self.lamta_max = float(setup['QUALITY_CALIBRATION']['lamta_max'])
        self.kads_nh4_min = float(setup['QUALITY_CALIBRATION']['kads_nh4_min'])
        self.kads_nh4_max = float(setup['QUALITY_CALIBRATION']['kads_nh4_max'])
        self.kdes_nh4_min = float(setup['QUALITY_CALIBRATION']['kdes_nh4_min'])
        self.kdes_nh4_max = float(setup['QUALITY_CALIBRATION']['kdes_nh4_max'])
        self.kads2_nh4_min = float(setup['QUALITY_CALIBRATION']['kads2_nh4_min'])
        self.kads2_nh4_max = float(setup['QUALITY_CALIBRATION']['kads2_nh4_max'])
        self.kdes2_nh4_min = float(setup['QUALITY_CALIBRATION']['kdes2_nh4_min'])
        self.kdes2_nh4_max = float(setup['QUALITY_CALIBRATION']['kdes2_nh4_max'])
        self.k_nit_min = float(setup['QUALITY_CALIBRATION']['k_nit_min'])
        self.k_nit_max = float(setup['QUALITY_CALIBRATION']['k_nit_max'])
        self.k_denit_min = float(setup['QUALITY_CALIBRATION']['k_denit_min'])
        self.k_denit_max = float(setup['QUALITY_CALIBRATION']['k_denit_max'])
        self.D_nh4_min = float(setup['QUALITY_CALIBRATION']['D_nh4_min'])
        self.D_nh4_max = float(setup['QUALITY_CALIBRATION']['D_nh4_max'])
        self.D_no3_min = float(setup['QUALITY_CALIBRATION']['D_no3_min'])
        self.D_no3_max = float(setup['QUALITY_CALIBRATION']['D_no3_max'])
        self.pop_init = int(setup['QUALITY_CALIBRATION']['pop'])
        self.gen_init = int(setup['QUALITY_CALIBRATION']['gen'])
        self.obs_file_nh4 = setup['QUALITY_CALIBRATION']['obs_file_nh4']
        self.obs_file_no3 = setup['QUALITY_CALIBRATION']['obs_file_no3']
        self.show_summary = bool(setup['QUALITY_CALIBRATION']['show_summary'])


class SoilPlant:
    def __init__(self, setup_file, nusz_ini):
        setup = configparser.ConfigParser()
        setup.read(setup_file)
        self.ro_pd = float(setup['SOIL_PLANT']['ro'])
        self.f = float(setup['SOIL_PLANT']['f'])
        self.lamda = float(setup['SOIL_PLANT']['lamda'])
        self.lamta = float(setup['SOIL_PLANT']['lamta'])
        self.root_fraction = float(setup['SOIL_PLANT']['root_fraction'])
        self.c_o2_root = float(setup['SOIL_PLANT']['c_o2_root'])
        self.ro = (1 - nusz_ini) * self.ro_pd


class Pollutant:
    def __init__(self, m_usz, m_sz):
        self.initial_concentration_usz = [0] * m_usz
        self.initial_concentration_sz = [0] * m_sz
        self.concentration_PZ = []
        self.concentration_USZ = [self.initial_concentration_usz]
        self.concentration_soil_USZ = [self.initial_concentration_usz]
        self.concentration_SZ = [self.initial_concentration_sz]
        self.concentration_soil_SZ = [self.initial_concentration_sz]
        self.reaction_rate_USZ = [self.initial_concentration_usz]
        self.reaction_rate_SZ = [self.initial_concentration_sz]
        self.concentration_PZ_before = 0
        self.concentration_soil_USZ_before = 0
        self.concentration_soil_SZ_before = 0
        self.mass_stormwater = []
        self.mass_soil = []
        self.mass_accumulated_reaction = []
        self.mass_accumulated_inflow = []
        self.mass_accumulated_overflow = []
        self.mass_accumulated_pipe_outflow = []
        self.mass_accumulated_infiltrated_to_SZ = []
        self.mass_accumulated_evapotranspiration = []
        self.mass_PZ = []
        self.mass_balance = []

    def f_transport(self, teta_i, teta_iplus1, ci, cs_i, dc, dc_dz, kads, kdes, D, UF, Rx, ro, f, dt, dz):
        if teta_i == 0:
            delta_c_i1 = 0
        elif teta_iplus1 == 0:
            delta_c_i1 = 0
        else:
            delta_c_i1 = ((1 / teta_iplus1) * dt * (
                    -teta_i * kads * ci + ro * kdes * cs_i + teta_i * (
                        D * f * (dc / dz ** 2) - UF * dc_dz) + Rx))
            # delta_c_i1 = ((1/teta_i)*GENERAL_PARAMETERS.dt*(-teta_i*kads*ci + ro*kdes*cs_i + teta_i*(D*SOIL_PLANT.f*(dc/GENERAL_PARAMETERS.dz**2) - UF*dc_dz) + Rx))
        return delta_c_i1

    def f_concentration_soil(self, cs_a, teta, kads, ci, kdes, kmicro, ro, dt):
        # Rxs = kmicro*cs_a
        Rxs = 0

        # Rxs = Um*(teta*cs_a/(Km + teta*cs_a))
        cs_abs = cs_a + ((teta / ro) * kads * ci - kdes * cs_a - Rxs) * dt

        if cs_abs <= 0:
            cs = 0
        else:
            cs = cs_abs

        return cs


class Ammonia(Pollutant):
    def __init__(self, m_usz, m_sz, setup_file, inflow_file):
        super(Ammonia, self).__init__(m_usz, m_sz)
        setup = configparser.ConfigParser()
        setup.read(setup_file)
        self.D = float(setup['NH4']['D_nh4'])
        self.kads = float(setup['NH4']['kads_nh4'])
        self.kdes = float(setup['NH4']['kdes_nh4'])
        self.kads2 = float(setup['NH4']['kads2_nh4'])
        self.kdes2 = float(setup['NH4']['kdes2_nh4'])
        self.k_mb = float(setup['NH4']['k_nh4_mb'])
        self.Fm = float(setup['NH4']['Fm_nh4'])
        self.Km = float(setup['NH4']['Km_nh4'])
        csv_file = pd.read_csv(inflow_file, sep=';')
        self.concentration_inflow = csv_file['nh4'].tolist()

    def f_reaction_pz(self):
        return 0

    def f_reaction_usz(self, C_nh4_iminus1, k_nit):
        R_nit = -k_nit * C_nh4_iminus1
        return R_nit

    def f_reaction_sz(self):
        return 0

    def f_plant_uptake_usz(self, C_nh4_2, teta_sm, root_fraction):
        PU_nh4_2 = -root_fraction * (self.Fm * teta_sm * C_nh4_2 / (self.Km + teta_sm * C_nh4_2))
        return PU_nh4_2

    def f_plant_uptake_sz(self, C_nh4_3, teta_b, root_fraction):
        PU_nh4_3 = -root_fraction * (self.Fm * teta_b * C_nh4_3 / (self.Km + teta_b * C_nh4_3))
        return PU_nh4_3

class Nitrate(Pollutant):
    def __init__(self, m_usz, m_sz, setup_file, inflow_file):
        super(Nitrate, self).__init__(m_usz, m_sz)
        setup = configparser.ConfigParser()
        setup.read(setup_file)
        self.D = float(setup['NO3']['D_no3'])
        self.Fm = float(setup['NO3']['Fm_no3'])
        self.Km = float(setup['NO3']['Km_no3'])
        self.kads = 0
        self.kdes = 0
        self.kads2 = 0
        self.kdes2 = 0
        csv_file = pd.read_csv(inflow_file, sep=';')
        self.concentration_inflow = csv_file['no3'].tolist()

    def f_reaction_pz(self):
        return 0

    def f_reaction_usz(self, C_nh4_iminus1, k_nit):
        R_nit = k_nit * C_nh4_iminus1
        return R_nit

    def f_reaction_sz(self, C_no3_iminus1, C_o2_i, C_doc_iminus1, k_denit):
        #     Of = O2.K/(O2.K+C_o2_i)
        #     bDOCf = (C_doc_iminus1 + DOC.bDOCd*GENERAL_PARAMETERS.dt)/(C_doc_iminus1 + DOC.bDOCd*GENERAL_PARAMETERS.dt + DOC.KbDOC)
        #
        #     k2 = GENERAL_PARAMETERS.k_denit*Of*bDOCf
        k2 = k_denit  ###testando sem influencia de DOC e O2
        R_denit = - k2 * C_no3_iminus1
        return R_denit


    def f_plant_uptake_usz(self, C_no3_2, teta_sm, root_fraction):
        PU_no3_2 = -root_fraction * (self.Fm * teta_sm * C_no3_2 / (self.Km + teta_sm * C_no3_2))

        return PU_no3_2

    def f_plant_uptake_sz(self, C_no3_3, teta_b, root_fraction):
        PU_no3_3 = -root_fraction * (self.Fm * teta_b * C_no3_3 / (self.Km + teta_b * C_no3_3))

        return PU_no3_3


class Oxygen(Pollutant):
    def __init__(self, m_usz, m_sz, setup_file, inflow_file):
        super(Oxygen, self).__init__(m_usz, m_sz)
        setup = configparser.ConfigParser()
        setup.read(setup_file)
        self.D = float(setup['O2']['D_o2'])
        self.K = float(setup['O2']['k_inib_o2'])  # perguntar pq K no lugar de k_inib_02
        self.k = float(setup['O2']['k_o2'])
        self.kads = 0
        self.kdes = 0
        self.kads2 = 0
        self.kdes2 = 0
        csv_file = pd.read_csv(inflow_file, sep=';')
        self.concentration_inflow = csv_file['o2'].tolist()
        
    def f_reaction_pz(self):
        return 0

    def f_reaction_usz(self, C_o2_iminus1, C_nh4_iminus1, k_nit):
        R_o2 = -self.k * C_o2_iminus1 - k_nit * C_nh4_iminus1 / 2
        return R_o2

    def f_reaction_sz(self, C_o2_iminus1, C_nh4_iminus1, k_nit):
        R_o2 = -self.k * C_o2_iminus1 - k_nit * C_nh4_iminus1 / 2

        return R_o2


    def f_plant_uptake_usz(self, C_o2_2, C_o2_root, teta_sm, root_fraction, lamda):
        PU_o2_2 = -root_fraction * (lamda * (teta_sm * C_o2_root - teta_sm * C_o2_2))
        return PU_o2_2

    def f_plant_uptake_sz(self, C_o2_3, C_o2_root, teta_b, root_fraction, lamda):
        PU_o2_3 = -root_fraction * (lamda * (teta_b * C_o2_root - teta_b * C_o2_3))
        return PU_o2_3


class DissolvedOrganicCarbon(Pollutant):
    def __init__(self, m_usz, m_sz, setup_file, inflow_file):
        super(DissolvedOrganicCarbon, self).__init__(m_usz, m_sz)
        setup = configparser.ConfigParser()
        setup.read(setup_file)
        self.D = float(setup['DOC']['D_doc'])
        self.fb = float(setup['DOC']['fb_doc'])
        self.bDOCd = float(setup['DOC']['bDOCd'])
        self.KbDOC = float(setup['DOC']['KbDOC'])
        self.k = float(setup['DOC']['k_doc'])
        self.kads = float(setup['DOC']['kads_doc'])
        self.kdes = float(setup['DOC']['kdes_doc'])
        self.kads2 = float(setup['DOC']['kads2_doc'])
        self.kdes2 = float(setup['DOC']['kdes2_doc'])
        self.k_mb = float(setup['DOC']['k_doc_mb'])
        csv_file = pd.read_csv(inflow_file, sep=';')
        self.concentration_inflow = csv_file['doc'].tolist()
        
    def f_reaction_pz(self):
        return 0

    def f_reaction_usz(self, C_doc_iminus1):
        R_doc = -self.k * C_doc_iminus1 + self.bDOCd
        return R_doc

    def f_reaction_sz(self, C_doc_iminus1):
        R_doc = -self.k * C_doc_iminus1 + self.bDOCd
        return R_doc