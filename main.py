import parameters
from modified_mpire import water_flow_module, nitrogen_module, ecoli_module
import water_flow_calibration

if __name__ == "__main__":
    while True:
        action = input("Type f to run the water_flow_module\n"
                       "Type n for the nitrogen_module\n"
                       "Type e for the ecoli_module\n"
                       "Type cf to run the calibration of water_flow\n"
                       "Type cn to run the calibration of nitrogen (NOT YET IMPLEMENTED)\n"
                       "Type ce to run the calibration of ecoli (NOT YET IMPLEMENTED\n"
                       "Option: "
                       )

        if action == "f":
            SETUP_FILE = "parameters.ini"
            WATER_FLOW_INPUT_FILE = "C:\\Users\\artur\\OneDrive\\Mestrado\\TC\\Dissertação\\evento_qualidade\\eventos_nitrogenio_ecoli\\evento_03\\concentration_inflow.csv"
            GP = parameters.GeneralParameters(SETUP_FILE, WATER_FLOW_INPUT_FILE)
            USZ = parameters.UnsaturatedZone(SETUP_FILE, GP.L, GP.hpipe, GP.dz)
            PZ = parameters.PondingZone(SETUP_FILE)
            SZ = parameters.SaturatedZone(SETUP_FILE, GP, USZ)
            water_flow_module(GP, USZ, PZ, SZ)
            print("Done")
            MB = parameters.MassBalance(GP, PZ, USZ, SZ)
            break

        if action == 'rtc_f':
            SETUP_FILE = "parameters.ini"
            WATER_FLOW_INPUT_FILE = "C:\\Users\\artur\\OneDrive\\Mestrado\\TC\\Dissertação\\evento_qualidade\\eventos_nitrogenio_ecoli\\evento_03\\concentration_inflow.csv"
            WATER_QUALITY_INPUT_FILE = WATER_FLOW_INPUT_FILE
            GP = parameters.GeneralParameters(SETUP_FILE, WATER_FLOW_INPUT_FILE)
            USZ = parameters.UnsaturatedZone(SETUP_FILE, GP.L, GP.hpipe, GP.dz)
            PZ = parameters.PondingZone(SETUP_FILE)
            SZ = parameters.SaturatedZone(SETUP_FILE, GP, USZ)
            RTC = parameters.RealTimeControl(SETUP_FILE)

            SOIL_PLANT = parameters.SoilPlant(SETUP_FILE, USZ.nusz_ini)
            NH4 = parameters.Ammonia(USZ.m_usz, SZ.m_sz, SETUP_FILE, WATER_QUALITY_INPUT_FILE)
            NO3 = parameters.Nitrate(USZ.m_usz, SZ.m_sz, SETUP_FILE, WATER_QUALITY_INPUT_FILE)
            O2 = parameters.Oxygen(USZ.m_usz, SZ.m_sz, SETUP_FILE, WATER_QUALITY_INPUT_FILE)
            DOC = parameters.DissolvedOrganicCarbon(USZ.m_usz, SZ.m_sz, SETUP_FILE, WATER_QUALITY_INPUT_FILE)

            water_flow_module(GP, USZ, PZ, SZ, RTC)
            data_o2, data_nh4, data_no3, data_doc = nitrogen_module(GP, PZ, USZ, SZ, SOIL_PLANT, NH4, NO3, O2, DOC)
            print("Done")
            MB = parameters.MassBalance(GP, PZ, USZ, SZ)

        elif action == "n":
            SETUP_FILE = "parameters.ini"
            WATER_FLOW_INPUT_FILE = "C:\\Users\\artur\\OneDrive\\Mestrado\\TC\\Dissertação\\evento_qualidade\\eventos_nitrogenio_ecoli\\evento_03\\concentration_inflow.csv"
            WATER_QUALITY_INPUT_FILE = WATER_FLOW_INPUT_FILE
            GP = parameters.GeneralParameters(SETUP_FILE, WATER_FLOW_INPUT_FILE)
            USZ = parameters.UnsaturatedZone(SETUP_FILE, GP.L, GP.hpipe, GP.dz)
            PZ = parameters.PondingZone(SETUP_FILE)
            SZ = parameters.SaturatedZone(SETUP_FILE, GP, USZ)
            SOIL_PLANT = parameters.SoilPlant(SETUP_FILE, USZ.nusz_ini)
            NH4 = parameters.Ammonia(USZ.m_usz, SZ.m_sz, SETUP_FILE, WATER_QUALITY_INPUT_FILE)
            NO3 = parameters.Nitrate(USZ.m_usz, SZ.m_sz, SETUP_FILE, WATER_QUALITY_INPUT_FILE)
            O2 = parameters.Oxygen(USZ.m_usz, SZ.m_sz, SETUP_FILE, WATER_QUALITY_INPUT_FILE)
            DOC = parameters.DissolvedOrganicCarbon(USZ.m_usz, SZ.m_sz, SETUP_FILE, WATER_QUALITY_INPUT_FILE)
            water_flow_module(GP, USZ, PZ, SZ)
            data_o2, data_nh4, data_no3, data_doc = nitrogen_module(GP, PZ, USZ, SZ, SOIL_PLANT, NH4, NO3, O2, DOC)
            print("Done")
            MB = parameters.MassBalance(GP, PZ, USZ, SZ)
            break

        elif action == "e":
            SETUP_FILE = "parameters.ini"
            WATER_FLOW_INPUT_FILE = "C:\\Users\\artur\\OneDrive\\Mestrado\\TC\\Dissertação\\evento_qualidade\\eventos_nitrogenio_ecoli\\evento_03\\concentration_inflow.csv"
            WATER_QUALITY_INPUT_FILE = WATER_FLOW_INPUT_FILE
            GP = parameters.GeneralParameters(SETUP_FILE, WATER_FLOW_INPUT_FILE)
            USZ = parameters.UnsaturatedZone(SETUP_FILE, GP.L, GP.hpipe, GP.dz)
            PZ = parameters.PondingZone(SETUP_FILE)
            SZ = parameters.SaturatedZone(SETUP_FILE, GP, USZ)
            SOIL_PLANT = parameters.SoilPlant(SETUP_FILE, USZ.nusz_ini)
            EC = parameters.Ecoli(USZ.nusz_ini, USZ.m_usz, SZ.m_sz, SETUP_FILE, WATER_QUALITY_INPUT_FILE)
            water_flow_module(GP, USZ, PZ, SZ)
            data_ecoli = ecoli_module(GP, PZ, USZ, SZ, SOIL_PLANT, EC)
            print("Done")
            MB = parameters.MassBalance(GP, PZ, USZ, SZ)
            break

        elif action == "cf":
            setup_file = "parameters.ini"
            CALIBRATION = parameters.Calibration("parameters.ini")
            water_flow_calibration.water_flow_calibration(CALIBRATION)
            print("Done")
            break

        elif action == "cq":
            pass

        else:
            print("Wrong option\n")
