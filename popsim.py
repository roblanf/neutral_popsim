#Copyright (C) 2012 Robert Lanfear
#
#This program is free software: you can redistribute it and/or modify it
#under the terms of the GNU General Public License as published by the
#Free Software Foundation, either version 3 of the License, or (at your
#option) any later version.
#
#This program is distributed in the hope that it will be useful, but
#WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#General Public License for more details. You should have received a copy
#of the GNU General Public License along with this program.  If not, see
#<http://www.gnu.org/licenses/>. 

'''This little script is a population genetics simulation of a single locus
    populations start as a list of zeros, and a single iteration ends when the population 
    is fixed for 1's. The idea is to flexibly simulate overlapping generations and 
    fluctuating population size.
    
    Parameters:
        Ss   -   a list of probabilities that an individual survives a given timestep
                S=0.0 is non-overlapping generations
        N   -   a list of lists of population sizes to switch between 
        pNs  - a list of  probabilities (per timestep) of switching to a new population size from N
        us  -   a list of mutation rates from 0->1. NB this should be <<1/Ne    
'''

#TODO allow parameters to be set as lists, so that one can do grid searches easily.
#or just use range()


import numpy as np
uniform = np.random.uniform
poisson = np.random.poisson
import random 
from itertools import product


Ss = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]     #survival probability per generation. 0.0 is non-overlapping generations
Ns = [[10, 20]]                                           #list of population sizes
pNs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,  0.8, 0.9, 1.0]                        #probability of changing popn size each generation (just choose a new popsize from the list)
us = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5]                   #mutation rate from 0->1. NB this should be <<1/Ne
R = 10000                                                     #number of reps
outfilename = "results.txt"

def new_popsize(n, N):
    '''Choose a new population size
        N is just a list of possible population sizes, it could be of any length
        This function works by picking a new population size value from N
        and making sure that it's different from the current population size (n)
    '''
    #get the position of the current population size in the list of possible sizes
    n_index = N.index(n)
    #make a copy of the list of population sizes
    N_copy = list(N)    
    #and remove the current population size from it
    N_copy.pop(n_index)
    #now pick the new population size by randomly picking from the remaining sizes.
    n_new = random.choice(N_copy)
    return n_new

def mutate(population, u):
    '''Add mutations to a population. NB, we only mutate from 0 to 1
    '''
    #first we figure out how many mutants we're going to generate
    num_mutants = poisson(u*len(population))
    
    #now we mutate some of the population
    #we allow for multiple hits here - we might choose a 1 and mutate to 1, but that's OK
    #but we don't allow for the same individual to get hit twice in one generation
    indices = range(len(population))
    random.shuffle(indices)
    mutate = indices[0:num_mutants]

    for i in mutate:
        population[i] = 1

    return population
    
def update_one_gen(pop, n, S, u):
    '''Take a list of 0s and 1s that represents a population
        - kill stuff off according to S
        - then grow or shrink the population to make it equal n
    '''

    #1. if the max. population size changed, reduce it by randomly killing individuals
    diff = len(pop)-n
    if diff>0: #our current pop is bigger than we want, so kill of diff individuals
        random.shuffle(pop) #randomise the order of the individuals
        pop = pop[diff:] #and delete the ones at the start    

    #2. copy the parental population, so we can add offspring later
    original_pop = list(pop)

    #3. Now randomly kill some of the parents
    number_to_kill = poisson(S*len(pop))
    random.shuffle(pop)
    pop = pop[number_to_kill:]

    #3. now add new offspring from the original parents 
    diff = len(pop)-n
    if diff<0: #our current pop is smaller than we want, so grow it
        newborns = []
        for i in range(-1*diff):
            newborns.append(random.choice(original_pop)) #we sample the newborns from the original list of parents
        newborns = mutate(newborns, u) #newborns might mutate...
        pop = pop+newborns
    return pop

def run_pop(pop, S, N, pN, u):
    ''' run a population until we hit a termination condition'''
    gen = 0
    n = len(pop)
    while True:
        pop = update_one_gen(pop, n, S, u)
        
        #test termination conditions
        zeros = pop.count(0)
        if pop.count(0)==0:     #fixed for 1
            break
            
        #if we didn't terminate, update pop size, and generations
        if uniform()<pN:
            n = new_popsize(n, N)
        gen += 1
    
    return pop, gen

outfile = open(outfilename, "w")
header =  "S\tN\tpN\tR\tu\tfix_times\tmean_fix\tstd_fix\tk\n"
outfile.write(header)
outfile.close()
print header

for params in product(Ss, Ns, pNs, us):
    S = params[0]
    N = params[1]
    pN = params[2]
    u = params[3]
    
    fixation_rates = []
    for i in range(R):
        if i%100 == 0: print i
    
        start_n = N[0] #start with the first population size in the list
            
        #start with all zeros
        pop = [0]*start_n
        
        #keep updating until we reach T, or until we fix for one allele
        pop, gen = run_pop(pop, S, N, pN, u)    
        
        if pop.count(1)==len(pop): #fixed for 1's
            fixation_rates.append(gen)
            
    result = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(str(S), str(N), str(pN), str(R), str(u), str(fixation_rates), str(np.mean(fixation_rates)), str(np.std(fixation_rates)), str(1.0/np.mean(fixation_rates)))
    print result
    outfile = open(outfilename, "a")
    outfile.write(result)
    outfile.close()