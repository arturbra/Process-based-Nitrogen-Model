import configparser
from math import pi
import pandas as pd


class ConcentrationInflow:
    def __init__(self, file):
        csv_file = pd.read_csv(file, sep=';')
        self.cin_nh4 = csv_file['nh4'].tolist()
        self.cin_no3 = csv_file['no3'].tolist()
        self.cin_o2 = csv_file['o2'].tolist()
        self.cin_doc = csv_file['doc'].tolist()


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
        self.Kweir = float(setup['PONDING_ZONE']['Kweir'])
        self.wWeir = float(setup['PONDING_ZONE']['wWeir'])
        self.expWeir = float(setup['PONDING_ZONE']['expWeir'])
        self.Cs = float(setup['PONDING_ZONE']['Cs'])
        self.Pp = float(setup['PONDING_ZONE']['Pp'])
        self.flagp = float(setup['PONDING_ZONE']['flagp'])
        self.k_denit_pz = float(setup['DENITRIFICATION']['k_denit_pz'])

    def f_concentration(self, cin, Qin_p, cp_a, I1, Qv, Rxi, hp, hp_a, Ab, dt, threshold):
        concentration = (cp_a * hp_a * Ab + (cin * Qin_p - cp_a * (I1 + Qv) + Rxi * hp * Ab) * dt) / (hp * Ab)
        if concentration < threshold:
            concentration = 0
        return concentration

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

    def f_weir_overflow(self, hp, dt, Qin, Qrain):
        volume = hp * self.Ab + dt * (Qin + Qrain)
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

    def f_infiltration_to_filter_material(self, Ks, hp, husz, A, dt, s, nusz, Qinfp):
        infiltration = min(Ks * A * (hp + husz) / husz, hp * self.Ab / dt - Qinfp, (1.0 - s) * nusz * husz * A / dt)
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

    def f_infiltration_to_sz(self, hp, husz, nusz, dt, sEST):
        if sEST >= self.sfc:
            Qfs = min((self.A * self.Ks * (hp + husz) / husz) * sEST ** self.gama, (sEST - self.sfc) * nusz * self.A * husz / dt)
        else:
            Qfs = 0
        return Qfs

    def f_porosity(self, husz, hsz, ng, Dg, Df):
        if hsz < Dg:
            nusz = (self.nusz_ini * Df + ng * (Dg - hsz)) / husz
        else:
            nusz = self.nusz_ini
        return nusz

    def f_alfa_beta(self, layer):
        alfa = (self.m_usz - 1 - layer) / (self.m_usz - 1)
        beta = layer / (self.m_usz - 1)
        return alfa, beta

    def f_unit_flux(self, l, Qpf, Qet_1, Qfs, Qhc, Qorif, Qinf_sz, teta_sm_i, Ab, hpipe):
        alfa, beta = self.f_alfa_beta(l)
        if hpipe > 0:
            unit_flux = (alfa * (Qpf - Qet_1) + beta * (Qfs - Qhc)) / (Ab * teta_sm_i)
        else:
            unit_flux = (alfa * (Qpf - Qet_1) + beta * (Qorif + Qinf_sz - Qhc)) / \
                        (Ab * teta_sm_i)
        return unit_flux

    def f_peclet(self, unit_flux, D, dz):
        if D > 0:
            peclet = unit_flux * dz / D
        else:
            peclet = 100
        return peclet


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

    def f_alfa_beta(self, layer):
        alfa = (self.m_sz - 1 - layer) / (self.m_sz - 1)
        beta = layer / (self.m_sz - 1)
        return alfa, beta

    def f_unit_flux(self, sz_layer, Qfs, Qhc, Qet_2, Qorif, Qinf_sz, teta_b_i, Ab):
        alfa, beta = self.f_alfa_beta(sz_layer)
        unit_flux = (alfa * (Qfs - Qhc - Qet_2) + beta * (Qorif + Qinf_sz)) / (Ab * teta_b_i)
        return unit_flux

    def f_peclet(self, unit_flux, D, dz):
        if D > 0:
            peclet = unit_flux * dz / D
        else:
            peclet = 100
        return peclet


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

    def f_plant_uptake_usz(self, concentration, teta_sm, parameter, root_influx=0, michaelis_uptake_constant=0):
        if parameter == "O2":
            plant_uptake = -self.root_fraction * (self.lamda * (teta_sm * self.c_o2_root - teta_sm * concentration))
        else:
            plant_uptake = -self.root_fraction * (root_influx * teta_sm * concentration / (michaelis_uptake_constant + teta_sm * concentration))
        return plant_uptake

    def f_plant_uptake_sz(self, concentration, teta_b, parameter, root_influx=0, michaleis_uptake_constant=0):
        if parameter == "O2":
            plant_uptake = -self.root_fraction * (self.lamda * (teta_b * self.c_o2_root - teta_b * concentration))
        else:
            plant_uptake = -self.root_fraction * (root_influx * teta_b * concentration / (michaleis_uptake_constant + teta_b * concentration))
        return plant_uptake


class Nutrient:
    def __init__(self, m_usz, m_sz):
        self.cp_a = 0
        self.cs_usz_a = 0
        self.cs_sz_a = 0
        self.cp = []
        self.cpi = 0
        self.c0_usz = [0] * m_usz
        self.c0_sz = [0] * m_sz
        self.cs0_usz = [self.cs_usz_a] * m_usz
        self.cs0_sz = [self.cs_sz_a] * m_sz
        self.c_usz = [self.c0_usz]
        self.cs_usz = [self.c0_usz]
        self.c_sz = [self.c0_sz]
        self.cl_i1 = []
        self.csi_usz = []
        self.Rxl = []
        self.cs_sz = [self.c0_sz]
        self.Rx_pz = 0
        self.Rx_usz = [self.c0_usz]
        self.Rxi_usz = 0
        self.Rx_sz = [self.c0_sz]
        self.mass_storm_asterisk_accumulated = []  # O passo corretivo é feito com a massa, M indica a massa que vai ser usada no passo corretivo. Ast é asterisco.
        self.mass_reaction_rate_accumulated = []
        self.mass_inflow_accumulated = []
        self.mass_overflow_accumulated = []
        self.mass_pipe_outflow_accumulated = []
        self.mass_infiltration_sz_accumulated = []
        self.mass_evapotranspiration_accumulated = []
        self.mass_pz_accumulated = []
        self.mass_soil_accumulated = []
        self.mass_balance = []

    def f_transport(self, teta_i, teta_iplus1, ci, cs_i, dc, dc_dz, kads, kdes, D, UF, Rx, dt, ro, f, dz):
        if teta_i == 0:
            delta_c_i1 = 0
        elif teta_iplus1 == 0:
            delta_c_i1 = 0
        else:
            delta_c_i1 = ((1 / teta_iplus1) * dt * (
                        -teta_i * kads * ci + ro * kdes * cs_i + teta_i * (D * f * (dc / dz ** 2) - UF * dc_dz) + Rx))
        return delta_c_i1

    def f_concentration_soil(self, cs_a, teta, kads, ci, kdes, ro, dt, method="FO", kmicro=0, Um=0, Km=0):
        if method == "FO":
            Rxs = kmicro * cs_a
        if method == "MM":
            Rxs = Um * (teta * cs_a / (Km + teta * cs_a))

        cs_abs = cs_a + ((teta / ro) * kads * ci - kdes * cs_a - Rxs) * dt

        if cs_abs <= 0:
            cs = 0
        else:
            cs = cs_abs
        return cs

    def concentration_water_phase_usz(self, m_usz, layer):
        if layer < (m_usz - 1):
            clplus1 = self.cli[layer + 1]
        else:
            clplus1 = 0
        return clplus1

    def concentration_water_phase_sz(self, m_sz, layer):
        if layer < (m_sz - 1):
            cjplus1 = self.cji[layer + 1]
        else:
            cjplus1 = 0
        return cjplus1

    def concentration_soil_phase_usz(self, usz_layer, time, teta, ci, ro, dt, threshold, method="FO", kmicro=0, Um=0, Km=0):
        self.cs = self.cs_usz[time][usz_layer]
        cs_next_iteration = self.f_concentration_soil(self.cs, teta, self.kads, ci, self.kdes, ro, dt, method, kmicro, Um, Km)
        if cs_next_iteration < threshold:
            cs_next_iteration = 0
        return cs_next_iteration

    def concentration_soil_phase_sz(self, sz_layer, time, teta, ci, ro, dt, threshold, method= "FO", kmicro = 0, Um=0, Km=0):
        self.cs =self.cs_sz[time][sz_layer]
        cs_next_iteration = self.f_concentration_soil(self.cs, teta, self.kads2, ci, self.kdes2, ro, dt, method, kmicro, Um, Km)
        if cs_next_iteration < threshold:
            cs_next_iteration = 0
        return(cs_next_iteration)

    def concentration_delta_usz(self, peclet, usz_layer, m_usz, dz, teta_usz, teta_sm_iplus1, uf, dt, ro, f, threshold):
        if peclet <= 2:
            if usz_layer == 0:
                dc = self.cliplus1 - 2 * self.cli[usz_layer] + self.cpi
                dc_dz = (self.cliplus1 - self.cpi) / (2 * dz)
            elif usz_layer == (m_usz - 1):
                dc = self.cli[usz_layer] - 2 * self.cli[usz_layer] + self.cli[usz_layer - 1]
                dc_dz = (self.cli[usz_layer] - self.cli[usz_layer - 1]) / dz
            else:
                dc = self.cliplus1 - 2 * self.cli[usz_layer] + self.cli[usz_layer - 1]
                dc_dz = (self.cliplus1 - self.cli[usz_layer - 1]) / (2 * dz)
        else:
            if usz_layer == 0:
                dc = self.cliplus1 - 2 * self.cli[usz_layer] + self.cpi
                dc_dz = (self.cli[usz_layer] - self.cpi) / dz
            elif usz_layer == (m_usz - 1):
                dc = self.cli[usz_layer] - 2 * self.cli[usz_layer] + self.cli[usz_layer - 1]
                dc_dz = (self.cli[usz_layer] - self.cli[usz_layer - 1]) / dz
            else:
                dc = self.cliplus1 - 2 * self.cli[usz_layer] + self.cli[usz_layer - 1]
                dc_dz = (self.cli[usz_layer] - self.cli[usz_layer - 1]) / dz

        delta_concentration = self.f_transport(teta_usz, teta_sm_iplus1, self.cli[usz_layer], self.cs, dc, dc_dz,
                                               self.kads, self.kdes, self.D, uf, self.Rxi_usz, dt, ro, f, dz)
        if teta_sm_iplus1 > 0:
            concentration = self.cli[usz_layer] + delta_concentration
        else:
            concentration = 0

        if concentration <= threshold:
            concentration = 0

        return concentration

    def concentration_delta_sz(self, m_usz, m_sz, n, peclet, sz_layer, dz, teta_sz, teta_b_iplus1, uf, dt, ro, f, threshold):
        if m_usz < (n - 1):
            if peclet <= 2:
                if sz_layer == 0:
                    dc = self.cjplus1 - 2 * self.cji[sz_layer] + self. cli[m_usz - 1]
                    dc_dz = (self.cjplus1 - self.cli[m_usz - 1]) / (2 * dz)
                elif sz_layer == (m_sz - 1):
                    dc = self.cji[sz_layer] - 2 * self.cji[sz_layer] + self.cji[sz_layer - 1]
                    dc_dz = (self.cji[sz_layer] - self.cji[sz_layer - 1]) / dz
                else:
                    dc = self.cjplus1 - 2 * self.cji[sz_layer] + self.cji[sz_layer - 1]
                    dc_dz = (self.cjplus1 - self.cji[sz_layer - 1]) / (2 * dz)
            else:
                if sz_layer == 0:
                    dc = self.cjplus1 - 2 * self.cji[sz_layer] + self.cli[m_usz - 1]
                    dc_dz = (self.cji[sz_layer] - self.cli[m_usz - 1]) / dz
                elif sz_layer == (m_sz - 1):
                    dc = self.cji[sz_layer] - 2 * self.cji[sz_layer] + self.cji[sz_layer - 1]
                    dc_dz = (self.cji[sz_layer] - self.cji[sz_layer - 1]) / dz
                else:
                    dc = self.cjplus1 - 2 * self.cji[sz_layer] + self.cji[sz_layer - 1]
                    dc_dz = (self.cji[sz_layer] - self.cji[sz_layer - 1]) / dz
                    
        if m_usz == (n - 1):
            dc = self.cji[sz_layer] - 2 * self.cji[sz_layer] + self.cji[sz_layer - 1]
            dc_dz = (self.cji[sz_layer] - self.cji[sz_layer - 1]) / (2 * dz)
        
        if m_usz == 0:
            if peclet <= 2:
                if sz_layer == 0:
                    dc = self.cjplus1 - 2 * self.cji[sz_layer] + self.cpi
                    dc_dz = (self.cjplus1 - self.cpi) / (2 * dz)
                elif sz_layer == (m_sz - 1):
                    dc = self.cji[sz_layer] - 2 * self.cji[sz_layer] + self.cji[sz_layer - 1]
                    dc_dz = (self.cji[sz_layer] - self.cji[sz_layer - 1]) / dz
                else:
                    dc = self.cjplus1 - 2 * self.cji[sz_layer] + self.cji[sz_layer - 1]
                    dc_dz = (self.cjplus1 - self.cji[sz_layer - 1]) / (2 * dz)
            else: 
                if sz_layer == 0:  # first cell
                    dc = self.cjplus1 - 2 * self.cji[sz_layer] + self.cpi
                    dc_dz = (self.cji[sz_layer] - self.cpi) / dz

                elif sz_layer == (SZ.m_sz - 1):  # last cell
                    dc = self.cji[sz_layer] - 2 * self.cji[sz_layer] + self.cji[sz_layer - 1]
                    dc_dz = (self.cji[sz_layer] - self.cji[sz_layer - 1]) / dz

                else:
                    dc = self.cjplus1 - 2 * self.cji[sz_layer] + self.cji[sz_layer - 1]
                    dc_dz = (self.cji[sz_layer] - self.cji[sz_layer - 1]) / dz
        delta_concentration = self.f_transport(teta_sz, teta_b_iplus1, self.cji[sz_layer], self.cs, dc, dc_dz,
                                               self.kads, self.kdes, self.D, uf, self.Rxi_sz, dt, ro, f, dz)
        if teta_b_iplus1 > 0:
            concentration = self.cji[sz_layer] + delta_concentration
        else:
            concentration = 0

        if concentration <= threshold:
            concentration = 0

        return concentration

    def f_mass_storm_asterisk(self, time, Ab, teta_usz, teta_sz, m_usz, m_sz, husz, hsz, hpipe):
        usz_therm = sum(self.c_usz[time]) * Ab * teta_usz[time] * husz[time] * 1000 / m_usz
        if hpipe > 0:
            sz_therm = sum(self.c_sz[time]) * Ab * teta_sz[time] * hsz[time] * 1000 / m_sz
        else:
            sz_therm = 0
        mass_storm_ast = usz_therm + sz_therm
        return mass_storm_ast

    def f_mass_soil(self, time, ro, Ab, husz, hsz, m_usz, m_sz, hpipe):
        mass_soil_before_usz = self.cs_usz_a * ro * Ab * husz[time] * 1000
        mass_soil_now_usz = sum(self.cs_usz[time]) * ro * Ab * husz[time] * 1000 / m_usz
        if hpipe > 0:
            mass_soil_before_sz = self.cs_sz_a * ro * Ab * hsz * 1000
            mass_soil_now_sz = sum(self.cs_sz[time]) * ro * Ab * hsz * 1000 / m_sz
        else:
            mass_soil_before_sz = 0
            mass_soil_now_sz = 0
        mass_soil_before = mass_soil_before_usz + mass_soil_before_sz
        mass_soil_now = mass_soil_now_usz + mass_soil_now_sz
        mass_soil = mass_soil_now - mass_soil_before
        return mass_soil

    def f_mass_reaction_rate(self, time, Ab, teta_usz, teta_sz, m_usz, m_sz, husz, hsz, hpipe):
        usz_therm = sum(self.Rx_usz[time]) * Ab * teta_usz[time] * husz[time] * 1000 / m_usz
        if hpipe > 0:
            sz_therm = sum(self.Rx_sz[time]) * Ab * teta_sz[time] * hsz[time] * 1000 / m_sz
        else:
            sz_therm = 0
        mass_reaction_rate = - (usz_therm + sz_therm)
        return mass_reaction_rate

    def f_mass_inflow(self, time, water_inflow, inflow_concentration, dt):
        mass_inflow = water_inflow[time] * inflow_concentration[time] * dt * 1000
        return mass_inflow

    def f_mass_overflow(self, time, overflow, dt):
        mass_overflow = overflow[time] * self.cp[time] * dt * 1000
        return mass_overflow

    def f_mass_pipe_outflow(self, time, pipe_outflow, m_usz, m_sz, dt, hpipe):
        if hpipe > 0:
            mass_pipe_outflow = pipe_outflow[time] * self.c_sz[time][m_sz - 1] * dt * 1000
        else:
            mass_pipe_outflow = pipe_outflow[time] * self.c_usz[time][m_usz - 1] * dt * 1000

        return mass_pipe_outflow

    def f_mass_infiltration_sz(self, time, water_infiltration, m_usz, m_sz, dt, hpipe):
        if hpipe > 0:
            mass_infiltration_sz = water_infiltration[time] * self.c_sz[time][m_sz - 1] * dt * 1000
        else:
            mass_infiltration_sz = water_infiltration[time] * self.c_usz[time][m_usz - 1] * dt * 1000
        return mass_infiltration_sz

    def f_mass_evapotranspiration(self, time, water_evapotranspiration, dt):
        mass_evapotranspiration = water_evapotranspiration[time] * self.cl_i1[0] * dt * 1000
        return mass_evapotranspiration

    def f_mass_pz(self, time, hpz, Ab):
        mass_pz = hpz[time] * Ab * self.cp[time] * 1000
        return mass_pz

    def f_mass_balance(self, time, mass_storm_asterisk, mass_soil, mass_reaction_rate, mass_inflow, mass_overflow,
                       mass_pipe_outflow, mass_infiltration_sz, mass_evapotranspiration, mass_pz):
        self.mass_pz_accumulated.append(mass_pz)
        self.mass_storm_asterisk_accumulated.append(mass_storm_asterisk)
        self.mass_soil_accumulated.append(mass_soil)

        if time == 0:
            self.mass_reaction_rate_accumulated.append(mass_reaction_rate)
            self.mass_inflow_accumulated.append(mass_inflow)
            self.mass_overflow_accumulated.append(mass_overflow)
            self.mass_pipe_outflow_accumulated.append(mass_pipe_outflow)
            self.mass_infiltration_sz_accumulated.append(mass_infiltration_sz)
            self.mass_evapotranspiration_accumulated.append(mass_evapotranspiration)
        else:
            self.mass_reaction_rate_accumulated.append(mass_reaction_rate + self.mass_reaction_rate_accumulated[-1])
            self.mass_inflow_accumulated.append(mass_inflow + self.mass_inflow_accumulated[-1])
            self.mass_overflow_accumulated.append(mass_overflow + self.mass_overflow_accumulated[-1])
            self.mass_pipe_outflow_accumulated.append(mass_pipe_outflow + self.mass_pipe_outflow_accumulated[-1])
            self.mass_infiltration_sz_accumulated.append(mass_infiltration_sz + self.mass_infiltration_sz_accumulated[-1])
            self.mass_evapotranspiration_accumulated.append(mass_evapotranspiration + self.mass_evapotranspiration_accumulated[-1])

        mass_balance = self.mass_inflow_accumulated[-1] - self.mass_pz_accumulated[-1] - self.mass_overflow_accumulated[-1] - \
                       self.mass_pipe_outflow_accumulated[-1] - self.mass_infiltration_sz_accumulated[-1] - \
                       self.mass_evapotranspiration_accumulated[-1] - self.mass_soil_accumulated[-1] - \
                       self.mass_reaction_rate_accumulated[-1]

        return mass_balance

    def f_delta_mass_balance_usz(self, time, Ab, teta_usz, husz, m_usz):
        delta_mass_balance_usz = (self.mass_balance[-1] - self.mass_storm_asterisk_accumulated[-1]) / (Ab * teta_usz[time] * husz[time] * 1000 / m_usz)
        concentration = self.cl_i1[1] + delta_mass_balance_usz
        if concentration < 0:
            concentration = 0
        return concentration

    def water_quality_results(self, m_usz, m_sz, inflow_concentration):
        columns_name = ["usz" + str(num) for num in range(m_usz)]
        columns_name += ["sz" + str(num) for num in range(m_sz)]
        if len(self.cs_usz) != 0:
            columns_name += ["usz_soil" + str(num) for num in range(m_usz)]
        if len(self.cs_sz) != 0:
            columns_name += ["sz_soil" + str(num) for num in range(m_sz)]
        columns_name += ["usz_rx" + str(num) for num in range(m_usz)]
        columns_name += ["sz_rx" + str(num) for num in range(m_sz)]


        df_usz = pd.DataFrame(self.c_usz)
        df_sz = pd.DataFrame(self.c_sz)
        df_soil_usz = pd.DataFrame(self.cs_usz)
        df_soil_sz = pd.DataFrame(self.cs_sz)
        df_rx_usz = pd.DataFrame(self.Rx_usz)
        df_rx_sz = pd.DataFrame(self.Rx_sz)
        frames = [df_usz, df_sz, df_soil_usz, df_soil_sz, df_rx_usz, df_rx_sz]
        df = pd.concat(frames, axis=1)
        df.columns = columns_name

        self.cp.append(0)
        df['pz'] = self.cp
        df['c_in'] = inflow_concentration
        df['t'] = list(range(len(self.c_usz)))
        return df


class Ammonia(Nutrient):
    def __init__(self, m_usz, m_sz, setup_file):
        super(Ammonia, self).__init__(m_usz, m_sz)
        setup = configparser.ConfigParser()
        setup.read(setup_file)
        self.D = float(setup['NH4']['D_nh4'])
        self.kads = float(setup['NH4']['kads_nh4'])
        self.kdes = float(setup['NH4']['kdes_nh4'])
        self.kads2 = float(setup['NH4']['kads2_nh4'])
        self.kdes2 = float(setup['NH4']['kdes2_nh4'])
        self.k_nh4_mb = float(setup['NH4']['k_nh4_mb'])
        self.Fm_nh4 = float(setup['NH4']['Fm_nh4'])
        self.Km_nh4 = float(setup['NH4']['Km_nh4'])


    def f_reaction_pz(self):
        return 0

    def f_reaction_usz(self, C_nh4_iminus1, k_nit):
        reaction = -k_nit * C_nh4_iminus1
        return reaction

    def f_reaction_sz(self):
        return 0

class Nitrate(Nutrient):
    def __init__(self, m_usz, m_sz, setup_file):
        super(Nitrate, self).__init__(m_usz, m_sz)
        setup = configparser.ConfigParser()
        setup.read(setup_file)
        self.cs_usz = []
        self.cs_sz = []
        self.D = float(setup['NO3']['D_no3'])
        self.Fm_no3 = float(setup['NO3']['Fm_no3'])
        self.Km_no3 = float(setup['NO3']['Km_no3'])
        self.kads = 0
        self.kdes = 0

    def f_reaction_pz(self):
        return 0

    def f_reaction_usz(self, C_nh4_iminus1, k_nit):
        reaction = k_nit * C_nh4_iminus1
        return reaction

    def f_reaction_sz(self, C_no3_iminus1, C_o2_i, C_doc_iminus1, k_denit):
        #     Of = K_o2/(K_o2+C_o2_i)
        #     bDOCf = (C_doc_iminus1 + bDOCd*dt)/(C_doc_iminus1 + bDOCd*dt + KbDOC)
        #
        #     k2 = k_denit*Of*bDOCf
        ###testando sem influencia de DOC e O2
        reaction = -k_denit * C_no3_iminus1
        return reaction

    def concentration_soil_phase_usz(self):
        return 0

    def f_mass_soil(self):
        return 0


class Oxygen(Nutrient):
    def __init__(self, m_usz, m_sz, setup_file):
        super(Oxygen, self).__init__(m_usz, m_sz)
        setup = configparser.ConfigParser()
        setup.read(setup_file)
        self.cs_usz = []
        self.cs_sz = []
        self.D = float(setup['O2']['D_o2'])
        self.K_o2 = float(setup['O2']['k_inib_o2'])  # perguntar pq K no lugar de k_inib_02
        self.k_o2 = float(setup['O2']['k_o2'])
        self.kads = 0
        self.kdes = 0

        
    def f_reaction_pz(self):
        return 0
    
    def f_reaction_usz(self, C_o2_iminus1, C_nh4_iminus1, k_nit):
        reaction = -self.k_o2 * C_o2_iminus1 - k_nit * C_nh4_iminus1 / 2
        return reaction
        
    def f_reaction_sz(self, C_o2_iminus1, C_nh4_iminus1, k_nit):
        reaction = -self.k_o2 * C_o2_iminus1 - k_nit * C_nh4_iminus1 / 2
        return reaction

    def concentration_soil_phase_usz(self):
        return 0

    def concentration_soil_phase_sz(self):
        return 0

    def f_mass_soil(self):
        return 0


class DissolvedOrganicCarbon(Nutrient):
    def __init__(self, m_usz, m_sz, setup_file):
        super(DissolvedOrganicCarbon, self).__init__(m_usz, m_sz)
        setup = configparser.ConfigParser()
        setup.read(setup_file)
        self.D = float(setup['DOC']['D_doc'])
        self.fb_doc = float(setup['DOC']['fb_doc'])
        self.bDOCd = float(setup['DOC']['bDOCd'])
        self.KbDOC = float(setup['DOC']['KbDOC'])
        self.k_doc = float(setup['DOC']['k_doc'])
        self.kads = float(setup['DOC']['kads_doc'])
        self.kdes = float(setup['DOC']['kdes_doc'])
        self.kads2 = float(setup['DOC']['kads2_doc'])
        self.kdes2 = float(setup['DOC']['kdes2_doc'])
        self.k_doc_mb = float(setup['DOC']['k_doc_mb'])
        
    def f_reaction_pz(self):
        return 0
    
    def f_reaction_usz(self, C_doc_iminus1):
        reaction = -self.k_doc * C_doc_iminus1 + self.bDOCd
        return reaction

    def f_reaction_sz(self, C_doc_iminus1):
        reaction = -self.k_doc * C_doc_iminus1 + self.bDOCd
        return reaction