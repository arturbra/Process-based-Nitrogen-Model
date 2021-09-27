#!/usr/bin/env python
# coding: utf-8

import random
from deap import base, creator
from deap import tools

print('Start importing Ecoli_model..')
from EColi_model_calib_NSE_log import * #ajeitar aqui

print('Finished importing.')


setup_file = "parametros_quali_NSE.ini"
setup = configparser.ConfigParser()
setup.read(setup_file)

# **1. Parameters to calibrate**
#lamta1_min = float(setup['CALIBRATION']['lamta1_min'])
#lamta1_max = float(setup['CALIBRATION']['lamta1_max'])
# lamta2_min = float(setup['CALIBRATION']['lamta2_min'])
# lamta2_max = float(setup['CALIBRATION']['lamta2_max'])
kads1_min = float(setup['CALIBRATION']['kads1_min'])
kads1_max = float(setup['CALIBRATION']['kads1_max'])
kdes1_min = float(setup['CALIBRATION']['kdes1_min'])
kdes1_max = float(setup['CALIBRATION']['kdes1_max'])
# kads2_min = float(setup['CALIBRATION']['kads2_min'])
# kads2_max = float(setup['CALIBRATION']['kads2_max'])
# kdes2_min = float(setup['CALIBRATION']['kdes2_min'])
# kdes2_max = float(setup['CALIBRATION']['kdes2_max'])
#etta_min = float(setup['CALIBRATION']['etta_min'])
#etta_max = float(setup['CALIBRATION']['etta_max'])
mue1_min = float(setup['CALIBRATION']['mue1_min'])
mue1_max = float(setup['CALIBRATION']['mue1_max'])
sita_min = float(setup['CALIBRATION']['sita_min'])
sita_max = float(setup['CALIBRATION']['sita_max'])


pop_init = int(setup['CALIBRATION']['pop'])
gen_init = int(setup['CALIBRATION']['gen'])
obs_file_EColi = setup['CALIBRATION']['obs_file_EColi']
show_summary = setup['CALIBRATION']['show_summary']

obs_df_EColi = pd.read_csv(obs_file_EColi) #csv file must have the same size of df
obs_df_EColi.set_index('t', inplace = True)
# obs_df_no3 = pd.read_csv(obs_file_no3) #csv file must have the same size of df
# obs_df_no3.set_index('t', inplace = True)

# **2. Definition of basic functions**
#creating evaluation function
def evaluate(individual):
    param = get_from_individual(individual)
    #param = read_from_ini()
    
    pen_total = penalty(individual)
    
    #print('individuo:', individual)
    
    data_EColi_log = EColi_run(param) #colocar o dataframe de ecoli
    
    # data_EColi_log.set_index('t', inplace = True)
    # data_EColi_log.set_index('t', inplace = True)
    
    data_EColi_log.to_csv('ecoli_evaluate.csv', index = False, sep=';', decimal=',')
    # data_no3.to_csv('no3_evaluate.csv', index = False)


    merge_df_EColi = data_EColi_log.merge(obs_df_EColi, left_index=True, right_index=True)

    # print(data_EColi_log)
    # #print(merge_df_EColi['Morif'])
    # input()
    # merge_df_EColi.to_csv('nh4_merge.csv', index = False)
    # merge_df_no3.to_csv('no3_merge.csv', index = False)
      
    #NSE for EColi
    avg_obs_EColi = merge_df_EColi['Morif_obs'].mean() #media dos valores observados
    # print(avg_obs_EColi)
    # input()
      
    merge_df_EColi['dif_sim_obs'] = merge_df_EColi['Morif_obs'] - merge_df_EColi['Morif'] #dif obs e simulado
    merge_df_EColi['dif_sim_obs_square'] = merge_df_EColi['dif_sim_obs']*merge_df_EColi['dif_sim_obs'] #quadrado da diferenca
    num = merge_df_EColi['dif_sim_obs_square'].sum() #soma dos quadrados da diferencas
      
    merge_df_EColi['dif_obs_avg'] = merge_df_EColi['Morif_obs'] - avg_obs_EColi #dif obs e media
    merge_df_EColi['dif_obs_avg_square'] = merge_df_EColi['dif_obs_avg']*merge_df_EColi['dif_obs_avg'] #quadrado da diferenca
    den = merge_df_EColi['dif_obs_avg_square'].sum()
      
    NSE_EColi = 1 - num/den
    NSE_final = NSE_EColi + pen_total

    return NSE_final,


    # #NSE for NO3
    # avg_obs_no3 = merge_df_no3['Morif_obs'].mean() #media dos valores observados
    #   
    # merge_df_no3['dif_sim_obs_2'] = merge_df_no3['Morif_obs'] - merge_df_no3['Morif'] #dif obs e simulado
    # merge_df_no3['dif_sim_obs_square_2'] = merge_df_no3['dif_sim_obs_2']*merge_df_no3['dif_sim_obs_2'] #quadrado da diferenca
    # num_2 = merge_df_no3['dif_sim_obs_square_2'].sum() #soma dos quadrados da diferencas
    #   
    # merge_df_no3['dif_obs_avg_2'] = merge_df_no3['Morif_obs'] - avg_obs_no3 #dif obs e media
    # merge_df_no3['dif_obs_avg_square_2'] = merge_df_no3['dif_obs_avg_2']*merge_df_no3['dif_obs_avg_2'] #quadrado da diferenca
    # den_2 = merge_df_no3['dif_obs_avg_square_2'].sum()
    # 
    # NSE_no3 = 1 - num_2/den_2
    #NSE_final = (NSE_nh4 + NSE_no3) / 2 + pen_total

    #print('ver arquivos csv e depois continuar para proxima rodada')
    #input()
                      

#creating ranges for penality
def penalty(individual):            #mudar essa parte, corrigir as variaveis
    # if lamta1_min < individual[0] < lamta1_max:
    #     pen0 = 0
    # else:
    #     pen0 = -10

    if kads1_min < individual[0] < kads1_max:
        pen1 = 0
    else:
        pen1 = -10
    
    if kdes1_min < individual[1] < kdes1_max:
        pen2 = 0
    else:
        pen2 = -10
    

    # if etta_min < individual[3] < etta_max:
    #     pen3 = 0
    # else:
    #     pen3 = -10

    if mue1_min < individual[2] < mue1_max:
        pen3 = 0
    else:
        pen3 = -10

    if sita_min < individual[3] < sita_max:
        pen4 = 0
    else:
        pen4 = -10

    # if k_denit_min < individual[6] < k_denit_max:
    #     pen6 = 0
    # else:
    #     pen6 = -10
    #
    # if D_nh4_min < individual[7] < D_nh4_max:
    #     pen7 = 0
    # else:
    #     pen7 = -10
    #
    # if D_no3_min < individual[8] < D_no3_max:
    #     pen8 = 0
    # else:
    #     pen8 = -10

    pen_total = pen1 + pen2 + pen3 + pen4 #+ pen5 + pen6 + pen7 + pen8 +pen0
    
    return pen_total

#creating functions to save results
def get_individual_values(individual):   #ajustar os parametros
    #global lamta1
    #global lamta2
    global kads1
    global kdes1
    # global kads2
    # global kdes2
    #global etta
    global mue1
    global sita

    #lamta1 = individual[0]
    kads1 = individual[0]
    kdes1 = individual[1]
    #etta = individual[3]
    mue1 = individual[2]
    sita = individual[3]

    # lamta2 = individual[1]
    # kads2 = individual[4]
    # kdes2 = individual[5]
#PAREI AQUI
def save_config(individual):    #ajustar os parametros
    get_individual_values(individual)
    #setup.set('SOIL', 'lamta1', str(lamta1))
    #setup.set('SOIL', 'lamta2', str(lamta2))
    setup.set('ECOLI', 'kads1', str(kads1))
    setup.set('ECOLI', 'kdes1', str(kdes1))
    # setup.set('ECOLI', 'kads2', str(kads2))
    # setup.set('ECOLI', 'kdes2', str(kdes2))
    #setup.set('ECOLI', 'etta', str(etta))
    setup.set('ECOLI', 'mue1', str(mue1))
    setup.set('ECOLI', 'sita', str(sita))


    with open('calibrated_quanli_Kin.ini', 'wt') as configfile:  #???
        setup.write(configfile)
        
# **3. GA setup**
creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", list, fitness=creator.FitnessMax)

toolbox = base.Toolbox()  #mudar os nomes dos parametros!
#toolbox.register("lamta1", random.uniform, lamta1_min, lamta1_max)
# toolbox.register("lamta2", random.uniform, lamta2_min, lamta2_max)
toolbox.register("kads1", random.uniform, kads1_min, kads1_max)
toolbox.register("kdes1", random.uniform, kdes1_min, kdes1_max)
# toolbox.register("kads2", random.uniform, kads2_min, kads2_max)
# toolbox.register("kdes2", random.uniform, kdes2_min, kdes2_max)
#toolbox.register("etta", random.uniform, etta_min, etta_max)
toolbox.register("mue1", random.uniform, mue1_min, mue1_max)
toolbox.register("sita", random.uniform, sita_min, sita_max)
toolbox.register("individual", tools.initCycle, creator.Individual, (toolbox.kads1, toolbox.kdes1, toolbox.mue1, toolbox.sita))
# toolbox.register("individual", tools.initCycle, creator.Individual, (toolbox.lamta1, toolbox.lamta2, toolbox.kads1, toolbox.kdes1, toolbox.kads2, toolbox.kdes2, toolbox.etta))
toolbox.register("population", tools.initRepeat, list, toolbox.individual)
toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)
toolbox.register("select", tools.selTournament, tournsize=3)
toolbox.register("evaluate", evaluate)

def calibrate():
    print('Inicio da função calibrate()')
    random.seed()
    pop = toolbox.population(n=pop_init)
    cx_prob, mut_prob, ngen = 0.5, 0.2, gen_init

    fitnesses = list(map(toolbox.evaluate, pop))

    print('fit')

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
            
        #stop criteria:       #parada manual quando para de evoluir! por enquanto +/- 1h
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