import parameters
from deap import base, creator, tools
import random
from scoop import futures


SETUP_FILE = 'parameters.ini'
CALIBRATION = parameters.Calibration(SETUP_FILE)
coefficient_limits = {0: (CALIBRATION.Kc_min, CALIBRATION.Kc_max), 1: (CALIBRATION.Ks_min, CALIBRATION.Ks_max),
                      2: (CALIBRATION.sh_min, CALIBRATION.sh_max), 3: (CALIBRATION.sw_min, CALIBRATION.sw_max),
                      4: (CALIBRATION.sfc_min, CALIBRATION.sfc_max), 5: (CALIBRATION.ss_min, CALIBRATION.ss_max),
                      6: (CALIBRATION.Kf_min, CALIBRATION.Kf_max), 7: (CALIBRATION.Cd_min, CALIBRATION.Cd_max)}


creator.create('FitnessMax', base.Fitness, weights=(1.0,))
creator.create('Individual', list, fitness=creator.FitnessMax)
toolbox = base.Toolbox()
toolbox.register("Kc", random.uniform, CALIBRATION.Kc_min, CALIBRATION.Kc_max)
toolbox.register("Ks", random.uniform, CALIBRATION.Ks_min, CALIBRATION.Ks_max)
toolbox.register("sh", random.uniform, CALIBRATION.sh_min, CALIBRATION.sh_max)
toolbox.register("sw", random.uniform, CALIBRATION.sw_min, CALIBRATION.sw_max)
toolbox.register("sfc", random.uniform, CALIBRATION.sfc_min, CALIBRATION.sfc_max)
toolbox.register("ss", random.uniform, CALIBRATION.ss_min, CALIBRATION.ss_max)
toolbox.register("Kf", random.uniform, CALIBRATION.Kf_min, CALIBRATION.Kf_max)
toolbox.register("Cd", random.uniform, CALIBRATION.Cd_min, CALIBRATION.Cd_max)
toolbox.register("individual", tools.initCycle, creator.Individual,
                 (toolbox.Kc, toolbox.Ks, toolbox.sh, toolbox.sw, toolbox.sfc, toolbox.ss, toolbox.Kf, toolbox.Cd))
toolbox.register('population', tools.initRepeat, list, toolbox.individual)
individual = toolbox.individual()
toolbox.register('evaluate', CALIBRATION.evaluate_by_nash)
toolbox.register('mate', tools.cxTwoPoint)
toolbox.register('mutate', tools.mutGaussian, mu=0, sigma=0.1, indpb=0.05)
toolbox.register('select', tools.selTournament, tournsize=3)


def water_flow_calibration(calib):
    print("Calibration started:")
    random.seed()
    pop = toolbox.population(n=calib.pop_init)
    cx_prob, mut_prob, ngen = 0.5, 0.05, calib.gen_init
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
                
        for mutant in offspring:
            if random.random() < mut_prob:
                toolbox.mutate(mutant)
                del mutant.fitness.values
        
        # evaluate invalid individuals
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = []
        for i in range(len(invalid_ind)):
            fitnesses.append(toolbox.evaluate(invalid_ind[i]))
        fitnesses = map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit
            
        print("Evaluated {} individuals".format(len(invalid_ind)))
        
        # offsprint complete replacement strategy
        pop[:] = offspring
        fits = [ind.fitness.values[0] for ind in pop]
        length = len(pop)
        mean = sum(fits) / length
        sum2 = sum(x * x for x in fits)
        std = abs(sum2 / length - mean ** 2) ** 0.5
        # nse = max(fits)
        best_ind = tools.selBest(pop, 1, fit_attr='fitness')[0]
        if calib.show_summary == 'True':
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
            

    best_ind = tools.selBest(pop, 1)[0]
    print("Done!")
    print("Best_ind:", best_ind, "Best_fit:", best_ind.fitness.values)

if __name__ == "__main__":
    toolbox.register("map", futures.map)
    water_flow_calibration(CALIBRATION)
