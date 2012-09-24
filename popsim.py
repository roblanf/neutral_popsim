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
#<http://www.gnu.org/licenses/>. PartitionFinder also includes the PhyML
#program and the PyParsing library both of which are protected by their
#own licenses and conditions, using PartitionFinder implies that you
#agree with those licences and conditions as well.

'''This little script is a population genetics simulation of a single locus
    populations start as a list of zeros, and a single iteration ends when the population 
    is fixed for 1's. The idea is to flexibly simulate overlapping generations and 
    fluctuating population size.
    
    Parameters:
        S   -   probability that an individual survives a given timestep
                S=0.0 is non-overlapping generations
        N   -   a list of population sizes to switch between 
        pN  -   probability (per timestep) of switching to a new population size from N
        u   -   mutation rate from 0->1. NB this should be <<1/Ne
    
'''


import numpy as np
uniform = np.random.uniform
import random 


S = 0.0         #survival probability per generation. 0.0 is non-overlapping generations
N = [1]        #list of population sizes
pN = 0.0        #probability of changing popn size each generation (just choose a new popsize from the list)
R = 1000       #number of reps
T = -1          #time to run simulation, if -1, it goes until a termination step is reached
u = 0.00001   #mutation rate from 0->1. NB this should be <<1/Ne

def get_start_pop(n):
    pop = [0]*n         #population of zeros
    pop[0] = 1          #a new mutant at the first pos.
    return pop

def new_popsize(n):
    '''Choose a new population size from N, that isn't the same as n
    '''
    n_index = N.index(n)
    N_copy = list(N)    
    N_copy.pop(n_index)
    n_new = random.choice(N_copy)
    return n_new

def mutate(population):
    '''Add mutations to a population
        NB, we only mutate from 0 to 1
        And, there's probably a much smarter way to do this, e.g. draw the number of muts
        from a poisson distribution then apply them. This one is computationally VERY 
        expensive
    '''

    for i in range(len(population)):
        if population[i] == 0 and uniform()<u: 
            population[i] = 1
            #print "MUTANT"
    
    return population
    
def update_one_gen(pop, n):
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
    ordering = range(len(pop))
    ordering.reverse()
    for i in ordering:    
        if uniform()>S:
            pop.pop(i) #pop.pop(), haha


    #3. now add new offspring from the original parents 
    diff = len(pop)-n
    if diff<0: #our current pop is smaller than we want, so grow it
        newborns = []
        for i in range(-1*diff):
            newborns.append(random.choice(original_pop)) #we sample the newborns from the original list of parents

        newborns = mutate(newborns) #newborns might mutate...
        pop = pop+newborns

    return pop

def run_pop(pop):
    ''' run a population until we hit a termination condition'''
    gen = 0
    n = len(pop)
    while True:
        #print pop, n
        pop = update_one_gen(pop, n)
        #print "1s = ", pop.count(1)
        
        #test termination conditions
        if gen==T:
            #print "TIMED OUT"
            break
        zeros = pop.count(0)
        if pop.count(0)==0:     #fixed for 1
            #print "FIXED FOR ONES"
            break
        #if pop.count(0)==n:      #fixed for 0
            #print "FIXED FOR ZEROS"
            #break
            
            
        #if we didn't terminate, update pop size, and generations
        if uniform()<pN:
            #print "CHANGING POPSIZE"
            n = new_popsize(n)       
        gen += 1
    
    return pop, gen


fixation_rates = []
extinction_rates = []
for i in range(R):
    if i%10 == 0: print i

    start_n = N[0] #start with the first population size in the list
    
    #start with a population with 1 mutant in it
    #pop = get_start_pop(start_n)
    print "********* NEW RUN %d *************" %i
    
    #start with all zeros
    pop = [0]*start_n
    
    #keep updating until we reach T, or until we fix for one allele
    pop, gen = run_pop(pop)    
    
    if pop.count(1)==len(pop): #fixed for 1's
        fixation_rates.append(gen)
    else:
        extinction_rates.append(gen)
    
    

print fixation_rates

print "S\tN\tpN\tR\tu\tfix_times\tmean_fix\tstd_fix\n"
print S, "\t", N, "\t", pN, "\t", R, "\t", u, "\t", fixation_rates, "\t", np.mean(fixation_rates), "\t", np.std(fixation_rates)
