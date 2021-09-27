#!/usr/bin/env python
# coding: utf-8

import random
from deap import base, creator
from deap import tools
from math import log10

print('Start importing EColi_model_calib_emc..')
from EColi_model_calib_RMSE import *

print('Finished importing.')


setup_file = "parametros_quali_RMSE.ini"
setup = configparser.ConfigParser()
setup.read(setup_file)

# **1. Parameters to calibrate**
#lamta1_min = float(setup['CALIBRATION']['lamta1_min'])
#lamta1_max = float(setup['CALIBRATION']['lamta1_max'])
kads1_min = float(setup['CALIBRATION']['kads1_min'])
kads1_max = float(setup['CALIBRATION']['kads1_max'])
kdes1_min = float(setup['CALIBRATION']['kdes1_min'])
kdes1_max = float(setup['CALIBRATION']['kdes1_max'])
mue1_min = float(setup['CALIBRATION']['mue1_min'])
mue1_max = float(setup['CALIBRATION']['mue1_max'])
sita_min = float(setup['CALIBRATION']['sita_min'])
sita_max = float(setup['CALIBRATION']['sita_max'])

pop_init = int(setup['CALIBRATION']['pop'])
gen_init = int(setup['CALIBRATION']['gen'])
obs_file_EColi = setup['CALIBRATION']['obs_file_EColi']
show_summary = setup['CALIBRATION']['show_summary']

obs_df_EColi = pd.read_csv(obs_file_EColi) #csv file must have the same size of df
#obs_df_EColi.set_index('t', inplace = True)


# **2. Definition of basic functions**
#creating evaluation function

def event_separator(Qin_mod):
    eventos = []
    
    contador_zero = 0
    evento = 0
    
    for q in Qin_mod: ## nao tem o mesmo t
        # print('q=', q)
        # print('evento=', evento)
        # print('contador_zero', contador_zero)


        inflow = q
        if inflow == 0:
            contador_zero += 1

        if contador_zero > 200 and inflow > 0: #novo evento
            evento +=1
            contador_zero = 0
        
        eventos.append(evento)
        # print('evento de novo', evento)
        # print('contador_zero de novo', contador_zero)
        #
        # print(".")

    return(eventos)

def calculate_EMC(df_poluente):    
    Mtotal = df_poluente['Morif'].sum()
    Vtotal = df_poluente['Qorif'].sum()*1000
    EMC = Mtotal / Vtotal
    
    return EMC
    
def calculate_EMCs(df_poluentes): 
    neventos = df_poluentes['eventos'].unique()

    EMCs = []
    
    for evento in neventos:
        if evento == 0:
            df_poluente = df_poluentes[df_poluentes['eventos'] == evento]
            EMC = calculate_EMC(df_poluente)
            EMC = log10(EMC)
            EMCs.append(EMC)

            #pass
        else:
            df_poluente = df_poluentes[df_poluentes['eventos'] == evento]
            EMC = calculate_EMC(df_poluente)
            EMC = log10(EMC)
            EMCs.append(EMC)

        #print(EMC)
    
    return EMCs


def evaluate(individual):
    
    param = get_from_individual(individual)

    pen_total = penalty(individual)

    # print('individuo:', individual)

    data_EColi_log = EColi_run(param)  # colocar o dataframe de ecoli

    M_ecoli_list = data_EColi_log['Morif'].tolist()
    Qin_mod = tQin[:len(M_ecoli_list)]   #troquei Qin por tQin


    eventos = event_separator(Qin_mod)
    
    data_EColi_log['eventos'] = eventos
    
    EMC_EColi_model_list = calculate_EMCs(data_EColi_log)

    #RMSE for Ecoli
    EMC_EColi_list = obs_df_EColi['EMC'].tolist()
    #print(EMC_EColi_list)
    #input()

    error_EColi_list = []


    for n in range(len(EMC_EColi_list)):

        #print('for n...', n)
        EMC_EColi_obs = EMC_EColi_list[n]
        #print('1- obs', EMC_EColi_obs)

        #print('n:', n)
        #print(' EMC_EColi_model_list:', EMC_EColi_model_list)




        EMC_EColi_model = EMC_EColi_model_list[n]
        #print('2- model', EMC_EColi_model)

        # print('EMC_EColi_model',EMC_EColi_model)
        # print('EMC_EColi_obs',EMC_EColi_obs)
        # input()
        error = (EMC_EColi_model - EMC_EColi_obs)**2
        error_EColi_list.append(error)

    
    N = len(EMC_EColi_list)
    RMSE_EColi = (sum(error_EColi_list)/N)**(1/2)
    # print(RMSE_EColi)
    # print(10**RMSE_EColi)
    # input()

    # #RMSE for NO3
    # EMC_no3_list = obs_df_no3['EMC'].tolist()
    # error_no3_list = []
    #
    # for n in range(len(EMC_no3_list)):
    #     EMC_no3_obs = EMC_no3_list[n]
    #     EMC_no3_model = EMC_no3_model_list[n]
    #     error = (EMC_no3_model - EMC_no3_obs)**2
    #     error_no3_list.append(error)
    #
    # N = len(EMC_no3_list)
    # RMSE_no3 = (sum(error_no3_list)/N)**(1/2)
    #
    #
    RMSE_final = RMSE_EColi - pen_total
    
    #input()
                      
    return RMSE_final,

#creating ranges for penality
def penalty(individual):
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

    if mue1_min < individual[2] < mue1_max:
        pen3 = 0
    else:
        pen3 = -10

    if sita_min < individual[3] < sita_max:
        pen4 = 0
    else:
        pen4 = -10
    

    pen_total = pen1 + pen2 + pen3 + pen4
    #pen0
    
    return pen_total

#creating functions to save results
def get_individual_values(individual):
    #global lamta1
    global kads1
    global kdes1
    global mue1
    global sita

    # lamta1 = individual[0]
    kads1 = individual[0]
    kdes1 = individual[1]
    mue1 = individual[2]
    sita = individual[3]

   # # lamta1 = individual[0]
   #  kads1 = individual[1]
   #  kdes1 = individual[2]
   #  mue1 = individual[3]
   #  sita = individual[4]
    
def save_config(individual):
    get_individual_values(individual)
    #setup.set('SOIL', 'lamta1', str(lamta1))
    setup.set('ECOLI', 'kads1', str(kads1))
    setup.set('ECOLI', 'kdes1', str(kdes1))
    setup.set('ECOLI', 'mue1', str(mue1))
    setup.set('ECOLI', 'sita', str(sita))


    with open('calibrated_quali_EColi.ini', 'wf_run') as configfile: #???
        setup.write(configfile)
        
# **3. GA setup**
creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
creator.create("Individual", list, fitness=creator.FitnessMin)

toolbox = base.Toolbox()
#toolbox.register("lamta1", random.uniform, lamta1_min, lamta1_max)
toolbox.register("kads1", random.uniform, kads1_min, kads1_max)
toolbox.register("kdes1", random.uniform, kdes1_min, kdes1_max)
toolbox.register("mue1", random.uniform, mue1_min, mue1_max)
toolbox.register("sita", random.uniform, sita_min, sita_max) #toolbox.lamta1
toolbox.register("individual", tools.initCycle, creator.Individual, (toolbox.kads1, toolbox.kdes1, toolbox.mue1, toolbox.sita))
toolbox.register("population", tools.initRepeat, list, toolbox.individual)
toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)
toolbox.register("select", tools.selTournament, tournsize=3)
toolbox.register("evaluate", evaluate)

def calibrate():
    print("Inicio calibrate()")
    random.seed()
    pop = toolbox.population(n=pop_init)
    cx_prob, mut_prob, ngen = 0.5, 0.2, gen_init    #MUDEI AQUI antes: "cx_prob, mut_prob, ngen = 0.5, 0.2, gen_init"

    fitnesses = list(map(toolbox.evaluate, pop))
    for ind, fit in zip(pop, fitnesses):
        ind.fitness.values = fit
    print("fit")

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