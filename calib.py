#!/usr/bin/env python
# coding: utf-8

import random
from deap import base, creator
from deap import tools
from model_pf import *

# !!!! para calibrar para mais de um evento, fazer um arquivo continuo com os eventos !!!!

setup_file = "tc_pf.ini"
setup = configparser.ConfigParser()
setup.read(setup_file)

# **1. Parameters to calibrate**
Kc_min = float(setup['CALIBRATION']['Kc_min'])
Kc_max = float(setup['CALIBRATION']['Kc_max'])
Ks_min = float(setup['CALIBRATION']['Ks_min'])
Ks_max = float(setup['CALIBRATION']['Ks_max'])
sh_min = float(setup['CALIBRATION']['sh_min'])
sh_max = float(setup['CALIBRATION']['sh_max'])
sw_min = float(setup['CALIBRATION']['sw_min'])
sw_max = float(setup['CALIBRATION']['sw_max'])
sfc_min = float(setup['CALIBRATION']['sfc_min'])
sfc_max = float(setup['CALIBRATION']['sfc_max'])
ss_min = float(setup['CALIBRATION']['ss_min'])
ss_max = float(setup['CALIBRATION']['ss_max'])
Kf_min = float(setup['CALIBRATION']['Kf_min'])
Kf_max = float(setup['CALIBRATION']['Kf_max'])
Cd_min = float(setup['CALIBRATION']['Cd_min'])
Cd_max = float(setup['CALIBRATION']['Cd_max'])

pop_init = int(setup['CALIBRATION']['pop'])
gen_init = int(setup['CALIBRATION']['gen'])
obs_file = setup['CALIBRATION']['obs_file']

show_summary = setup['CALIBRATION']['show_summary']

obs_df = pd.read_csv(obs_file) #csv file must have the same size of df
obs_df.set_index('t', inplace = True)

# **2. Definition of basic functions**
#creating evaluation function
def evaluate(individual):
    param = get_from_individual(individual)
    
    #param = read_from_ini()
    pen_total = penalty(individual)
    
    dados = run_W(param)
       

    #print('Vin:', dados['Qin'].sum(), 'Vv:', dados['Qv'].sum(), 'Vorif:', dados['Qorif'].sum())
    #print('---------------')
    
    dados.set_index('t', inplace = True)
     
    merge_df = dados.merge(obs_df, left_index=True, right_index=True)
    #merge_df.to_csv('verif_calib.csv', index = False)
      
    #NSE for Orifice
    avg_obs = merge_df['Qorif_obs'].mean() #media dos valores observados
      
    merge_df['dif_sim_obs'] = merge_df['Qorif_obs'] - merge_df['Qpipe'] #dif obs e simulado
    merge_df['dif_sim_obs_square'] = merge_df['dif_sim_obs']*merge_df['dif_sim_obs'] #quadrado da diferenca
    num = merge_df['dif_sim_obs_square'].sum() #soma dos quadrados da diferencas
      
    merge_df['dif_obs_avg'] = merge_df['Qorif_obs'] - avg_obs #dif obs e media
    merge_df['dif_obs_avg_square'] = merge_df['dif_obs_avg']*merge_df['dif_obs_avg'] #quadrado da diferenca
    den = merge_df['dif_obs_avg_square'].sum()
      
    NSE_1 = 1 - num/den

    #NSE for PD
    avg_obs_2 = merge_df['h_obs'].mean() #media dos valores observados
      
    merge_df['dif_sim_obs_2'] = merge_df['h_obs'] - merge_df['hpEND'] #dif obs e simulado
    merge_df['dif_sim_obs_square_2'] = merge_df['dif_sim_obs_2']*merge_df['dif_sim_obs_2'] #quadrado da diferenca
    num_2 = merge_df['dif_sim_obs_square_2'].sum() #soma dos quadrados da diferencas
      
    merge_df['dif_obs_avg_2'] = merge_df['h_obs'] - avg_obs_2 #dif obs e media
    merge_df['dif_obs_avg_square_2'] = merge_df['dif_obs_avg_2']*merge_df['dif_obs_avg_2'] #quadrado da diferenca
    den_2 = merge_df['dif_obs_avg_square_2'].sum()
      
    NSE_2 = 1 - num_2/den_2
    
    NSE_final = (NSE_1 + NSE_2)/2 + pen_total
    

    return NSE_final,

#creating ranges for penality
def penalty(individual):
    if Kc_min <= individual[0] <= Kc_max:
        pen0 = 0
    else:
        pen0 = -10
    
    if Ks_min <= individual[1] <= Ks_max:
        pen1 = 0
    else:
        pen1 = -10
    
    if sh_min <= individual[2] <= sh_max:
        pen2 = 0
    else:
        pen2 = -10    
    
    if sw_min <= individual[3] <= sw_max:
        pen3 = 0
    else:
        pen3 = -10
    
    if sfc_min <= individual[4] <= sfc_max:
        pen4 = 0
    else:
        pen4 = -10
    
    if ss_min <= individual[5] <= ss_max:
        pen5 = 0
    else:
        pen5 = -10

    if Kf_min <= individual[6] <= Kf_max:
        pen6 = 0
    else:
        pen6 = -10
    
    if Cd_min <= individual[7] <= Cd_max:
        pen7 = 0
    else:
        pen7 = -10
    
        
    pen_total = pen0 + pen1 + pen2 + pen3 + pen4 + pen5 + pen6 + pen7
    
    return pen_total

#creating functions to save results
def get_individual_values(individual):
    global Kc
    global Ks
    global sh
    global sw
    global sfc
    global ss
    global Kf
    global Cd
                
    Kc = individual[0]
    Ks = individual[1]
    sh = individual[2]
    sw = individual[3]
    sfc = individual[4]
    ss = individual[5]
    Kf = individual[6]
    Cd = individual[7]
    
# def save_config(individual):
#     get_individual_values(individual)
#     setup.set('GREEN_AMPT_INFILTRATION', 'psi', str(psi))
#     setup.set('GREEN_AMPT_INFILTRATION', 'ksat', str(ksat))
#     setup.set('GREEN_AMPT_INFILTRATION', 'te', str(te))
#     setup.set('GREEN_AMPT_INFILTRATION', 'teta_w', str(teta_w))
#     setup.set('GREEN_AMPT_INFILTRATION', 'teta_s', str(teta_s))
#     setup.set('UNDERDRAIN', 'Cd', str(Cd))
#     setup.set('UNDERDRAIN', 'dH', str(dH))
#     setup.set('PLANTS', 'kc', str(kc))
#     setup.set('BOTTOM_INFILTRATION', 'Ksat_ss', str(Ksat_ss))
#     setup.set('BOTTOM_INFILTRATION', 'psi_ss', str(psi_ss))
# #     setup.set('BIORETENTION_WEIR', 'kWeir', str(kWeir))
# #     setup.set('BIORETENTION_WEIR', 'expWeir', str(expWeir))     
#         
#     with open('calibrated_quanti.ini', 'wt') as configfile:
#         setup.write(configfile)
        
# **3. GA setup**
creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", list, fitness=creator.FitnessMax)

toolbox = base.Toolbox()
toolbox.register("Kc", random.uniform, Kc_min, Kc_max) 
toolbox.register("Ks", random.uniform, Ks_min, Ks_max)
toolbox.register("sh", random.uniform, sh_min, sh_max)
toolbox.register("sw", random.uniform, sw_min, sw_max)
toolbox.register("sfc", random.uniform, sfc_min, sfc_max)
toolbox.register("ss", random.uniform, ss_min, ss_max)
toolbox.register("Kf", random.uniform, Kf_min, Kf_max)
toolbox.register("Cd", random.uniform, Cd_min, Cd_max)

toolbox.register("individual", tools.initCycle, creator.Individual, (toolbox.Kc, toolbox.Ks, toolbox.sh, toolbox.sw, toolbox.sfc, toolbox.ss, toolbox.Kf, toolbox.Cd))
toolbox.register("population", tools.initRepeat, list, toolbox.individual)
toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)
toolbox.register("select", tools.selTournament, tournsize=3)
toolbox.register("evaluate", evaluate)

def calibrate():
    random.seed()
    pop = toolbox.population(n=pop_init)
    #print(pop)
    cx_prob, mut_prob, ngen = 0.5, 0.2, gen_init
    
    fitnesses = list(map(toolbox.evaluate, pop))
    for ind, fit in zip(pop, fitnesses):
        ind.fitness.values = fit

    print("Evaluated {} individuals".format(len(pop)))
    
    contador = 0
    for g in range(ngen):
        print("=========================")
        print("Generation {}".format(g))
        
        # choose the next generation
        offspring = toolbox.select(pop, len(pop))

        # clone selected individuals
        offspring = list(map(toolbox.clone, offspring))

        # compute the crossover
        for child1, child2 in zip(offspring[::2], offspring[1::2]):
            if random.random() < cx_prob:
                toolbox.mate(child1, child2)
                del child1.fitness.values
                del child2.fitness.values

        # compute the mutation
        for mutant in offspring:
            if random.random() < mut_prob:
                toolbox.mutate(mutant)
                del mutant.fitness.values

        # evaluate invalid individuals
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        print("Evaluated {} individuals".format (len(invalid_ind)))

        # offsprint complete replacement strategy
        pop[:] = offspring
        fits = [ind.fitness.values[0] for ind in pop]
        length = len(pop)
        mean = sum(fits) / length
        sum2 = sum(x*x for x in fits)
        std = abs(sum2 / length - mean**2)**0.5
        #nse = max(fits)
        best_ind = tools.selBest(pop, 1, fit_attr='fitness')[0]
        if show_summary == 'True':
            print("  Min %s" % min(fits))
            print("  Max %s" % max(fits))
            print("  Avg %s" % mean)
            print("  Std %s" % std)
            print("  => Best_ind:", best_ind)
            print("  => Best_fit:", best_ind.fitness.values)
        
        contador += 1
        if contador == 10:
            partial_results_file = open("partial_results.txt", "a")
            A = ['Generation: ', str(g), "\n"]
            B = ["best_ind: ", str(best_ind), "\n"]
            C = ["best_ind_NSE: ", str(best_ind.fitness.values), "\n"]
            partial_results_file.writelines(A)
            partial_results_file.writelines(B)
            partial_results_file.writelines(C)
            partial_results_file.close()
            contador = 0
            
        #stop criteria:
        # when |f(x_i) - f(x_{i-1})|/f(x_i) stops decreasing
        #during n generations.
        #  n = 3.
        # if abs((nse - self.previousError)/nse) < 0.0000001:
        #     self.repeatCounter += 1
        # else:
        #     self.repeatCounter = 0
        # self.previousError = nse
        # if self.repeatCounter == 10:
        #     break
        #saving the best result at each generation
        
    best_ind = tools.selBest(pop, 1)[0]
    #save_config(best_ind)
    
    print("Done!")    
    print("Best_ind:", best_ind, "Best_fit:", best_ind.fitness.values)
    
if __name__ == "__main__":
    inicio = datetime.datetime.now()
    print(inicio)
    calibrate()
    fim = datetime.datetime.now()
    print ('Elapsed time: ',fim - inicio)  
    