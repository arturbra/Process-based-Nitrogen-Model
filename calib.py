#!/usr/bin/env python
# coding: utf-8

import random
from deap import base, creator
from deap import tools

print('Start importing model_Kin..')
from model_Kin_pf_new import *

print('Finished importing.')


setup_file = "tc_quali_Kin.ini"
setup = configparser.ConfigParser()
setup.read(setup_file)

# **1. Parameters to calibrate**
lamta_min = float(setup['CALIBRATION']['lamta_min'])
lamta_max = float(setup['CALIBRATION']['lamta_max'])
kads_nh4_min = float(setup['CALIBRATION']['kads_nh4_min'])
kads_nh4_max = float(setup['CALIBRATION']['kads_nh4_max'])
kdes_nh4_min = float(setup['CALIBRATION']['kdes_nh4_min'])
kdes_nh4_max = float(setup['CALIBRATION']['kdes_nh4_max'])
kads2_nh4_min = float(setup['CALIBRATION']['kads2_nh4_min'])
kads2_nh4_max = float(setup['CALIBRATION']['kads2_nh4_max'])
kdes2_nh4_min = float(setup['CALIBRATION']['kdes2_nh4_min'])
kdes2_nh4_max = float(setup['CALIBRATION']['kdes2_nh4_max'])
k_nit_min = float(setup['CALIBRATION']['k_nit_min'])
k_nit_max = float(setup['CALIBRATION']['k_nit_max'])
k_denit_min = float(setup['CALIBRATION']['k_denit_min'])
k_denit_max = float(setup['CALIBRATION']['k_denit_max'])
D_nh4_min = float(setup['CALIBRATION']['D_nh4_min'])
D_nh4_max = float(setup['CALIBRATION']['D_nh4_max'])
D_no3_min = float(setup['CALIBRATION']['D_no3_min'])
D_no3_max = float(setup['CALIBRATION']['D_no3_max'])

pop_init = int(setup['CALIBRATION']['pop'])
gen_init = int(setup['CALIBRATION']['gen'])
obs_file_nh4 = setup['CALIBRATION']['obs_file_nh4']
obs_file_no3 = setup['CALIBRATION']['obs_file_no3']
show_summary = setup['CALIBRATION']['show_summary']

obs_df_nh4 = pd.read_csv(obs_file_nh4) #csv file must have the same size of df
obs_df_nh4.set_index('t', inplace = True)
obs_df_no3 = pd.read_csv(obs_file_no3) #csv file must have the same size of df
obs_df_no3.set_index('t', inplace = True)

# **2. Definition of basic functions**
#creating evaluation function
def evaluate(individual):
    param = get_from_individual(individual)
    #param = read_from_ini()
    
    pen_total = penalty(individual)
    
    #print('individuo:', individual)
    
    data_o2, data_nh4, data_no3, data_doc = run_Kin(param)
    
    data_nh4.set_index('t', inplace = True)
    data_no3.set_index('t', inplace = True)
    
#     data_nh4.to_csv('nh4_evaluate.csv', index = False)
#     data_no3.to_csv('no3_evaluate.csv', index = False)
     
    merge_df_nh4 = data_nh4.merge(obs_df_nh4, left_index=True, right_index=True)
    merge_df_no3 = data_no3.merge(obs_df_nh4, left_index=True, right_index=True)
    
#     merge_df_nh4.to_csv('nh4_merge.csv', index = False)
#     merge_df_no3.to_csv('no3_merge.csv', index = False)
      
    #NSE for NH4   
    avg_obs_nh4 = merge_df_nh4['Morif_obs'].mean() #media dos valores observados
      
    merge_df_nh4['dif_sim_obs'] = merge_df_nh4['Morif_obs'] - merge_df_nh4['Morif'] #dif obs e simulado
    merge_df_nh4['dif_sim_obs_square'] = merge_df_nh4['dif_sim_obs']*merge_df_nh4['dif_sim_obs'] #quadrado da diferenca
    num = merge_df_nh4['dif_sim_obs_square'].sum() #soma dos quadrados da diferencas
      
    merge_df_nh4['dif_obs_avg'] = merge_df_nh4['Morif_obs'] - avg_obs_nh4 #dif obs e media
    merge_df_nh4['dif_obs_avg_square'] = merge_df_nh4['dif_obs_avg']*merge_df_nh4['dif_obs_avg'] #quadrado da diferenca
    den = merge_df_nh4['dif_obs_avg_square'].sum()
      
    NSE_nh4 = 1 - num/den

    #NSE for NO3
    avg_obs_no3 = merge_df_no3['Morif_obs'].mean() #media dos valores observados
      
    merge_df_no3['dif_sim_obs_2'] = merge_df_no3['Morif_obs'] - merge_df_no3['Morif'] #dif obs e simulado
    merge_df_no3['dif_sim_obs_square_2'] = merge_df_no3['dif_sim_obs_2']*merge_df_no3['dif_sim_obs_2'] #quadrado da diferenca
    num_2 = merge_df_no3['dif_sim_obs_square_2'].sum() #soma dos quadrados da diferencas
      
    merge_df_no3['dif_obs_avg_2'] = merge_df_no3['Morif_obs'] - avg_obs_no3 #dif obs e media
    merge_df_no3['dif_obs_avg_square_2'] = merge_df_no3['dif_obs_avg_2']*merge_df_no3['dif_obs_avg_2'] #quadrado da diferenca
    den_2 = merge_df_no3['dif_obs_avg_square_2'].sum()
    
    NSE_no3 = 1 - num_2/den_2
    
    NSE_final = (NSE_nh4 + NSE_no3)/2 + pen_total
    
    #print('ver arquivos csv e depois continuar para proxima rodada')
    #input()
                      
    return NSE_final,

#creating ranges for penality
def penalty(individual):
    if lamta_min < individual[0] < lamta_max:
        pen0 = 0
    else:
        pen0 = -10 
    
    if kads_nh4_min < individual[1] < kads_nh4_max:
        pen1 = 0
    else:
        pen1 = -10
    
    if kdes_nh4_min < individual[2] < kdes_nh4_max:
        pen2 = 0
    else:
        pen2 = -10    
    
    if kads2_nh4_min < individual[3] < kads2_nh4_max:
        pen3 = 0
    else:
        pen3 = -10
    
    if kdes2_nh4_min < individual[4] < kdes2_nh4_max:
        pen4 = 0
    else:
        pen4 = -10
    
    if k_nit_min < individual[5] < k_nit_max:
        pen5 = 0
    else:
        pen5 = -10

    if k_denit_min < individual[6] < k_denit_max:
        pen6 = 0
    else:
        pen6 = -10

    if D_nh4_min < individual[7] < D_nh4_max:
        pen7 = 0
    else:
        pen7 = -10
        
    if D_no3_min < individual[8] < D_no3_max:
        pen8 = 0
    else:
        pen8 = -10    

    pen_total = pen0 + pen1 + pen2 + pen3 + pen4 + pen5 + pen6 + pen7 + pen8
    
    return pen_total

#creating functions to save results
def get_individual_values(individual):
    global lamta
    global kads_nh4
    global kdes_nh4
    global kads2_nh4
    global kdes2_nh4
    global k_nit
    global k_denit
    global D_nh4
    global D_no3
               
    lamta = individual[0]
    kads_nh4 = individual[1]
    kdes_nh4 = individual[2]
    kads2_nh4 = individual[3]
    kdes2_nh4 = individual[4]
    k_nit = individual[5]
    k_denit = individual[6]
    D_nh4 = individual[7]
    D_no3 = individual[8]
    
def save_config(individual):
    get_individual_values(individual)
    setup.set('SOIL_PLANT', 'lamta', str(lamta))
    setup.set('NH4', 'kads_nh4', str(kads_nh4))
    setup.set('NH4', 'kdes_nh4', str(kdes_nh4))
    setup.set('NH4', 'kads2_nh4', str(kads2_nh4))
    setup.set('NH4', 'kdes2_nh4', str(kdes2_nh4))
    setup.set('NITRIFICATION', 'k_nit', str(k_nit))
    setup.set('DENITRIFICATION', 'k_denit', str(k_denit))
    setup.set('NH4', 'D_nh4', str(D_nh4))
    setup.set('NO3', 'D_no3', str(D_no3))
       
    with open('calibrated_quanli_Kin.ini', 'wt') as configfile:
        setup.write(configfile)
        
# **3. GA setup**
creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", list, fitness=creator.FitnessMax)

toolbox = base.Toolbox()
toolbox.register("lamta", random.uniform, lamta_min, lamta_max) 
toolbox.register("kads_nh4", random.uniform, kads_nh4_min, kads_nh4_max)
toolbox.register("kdes_nh4", random.uniform, kdes_nh4_min, kdes_nh4_max)
toolbox.register("kads2_nh4", random.uniform, kads2_nh4_min, kads2_nh4_max)
toolbox.register("kdes2_nh4", random.uniform, kdes2_nh4_min, kdes2_nh4_max)
toolbox.register("k_nit", random.uniform, k_nit_min, k_nit_max)
toolbox.register("k_denit", random.uniform, k_denit_min, k_denit_max)
toolbox.register("D_nh4", random.uniform, D_nh4_min, D_nh4_max)
toolbox.register("D_no3", random.uniform, D_no3_min, D_no3_max)
toolbox.register("individual", tools.initCycle, creator.Individual, (toolbox.lamta, toolbox.kads_nh4, toolbox.kdes_nh4, toolbox.kads2_nh4, toolbox.kdes2_nh4, toolbox.k_nit, toolbox.k_denit, toolbox.D_nh4, toolbox.D_no3))
toolbox.register("population", tools.initRepeat, list, toolbox.individual)
toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)
toolbox.register("select", tools.selTournament, tournsize=3)
toolbox.register("evaluate", evaluate)

def calibrate():
    random.seed()
    pop = toolbox.population(n=pop_init)
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
        for x in fits:
            #print('x:', x)
            a = x*x
            #print('a:', a)
            
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
    save_config(best_ind)
    
    print("Done!")    
    print("Best_ind:", best_ind, "Best_fit:", best_ind.fitness.values)
    
if __name__ == "__main__":
    print('Starting calibrator...')
    inicio = datetime.datetime.now()
    print(inicio)
    calibrate()
    fim = datetime.datetime.now()
    print ('Elapsed time: ',fim - inicio)  
    print('Finish calibration.')