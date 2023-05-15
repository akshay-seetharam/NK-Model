import numpy as np
import random
from matplotlib import pyplot as plt


class Genome:
    def __init__(self, n, a, k, genotype=None, mutation_proportion=0.1):
        """ Constructor for state
        Arguments:
            n: length of genome
            a: options for choice at each locus in genome
            k: number of alleles that affect fitness from one allele
            genotype: predefined genotype, None by default
            mutation_proportion: proportion of possible mutants that form
        Returns:
            New instance of State """
        self.n = n
        self.a = a
        self.k = k
        self.mutation_proportion = mutation_proportion
        self.genotype = [random.randint(0, a - 1) for _ in range(n)] if genotype is None else genotype
        self.hash = g2h(self.genotype, a)
        # self.k_others = k_others
        self.ws = hash_ws[self.hash]

    def __str__(self):
        return f'Genome with n = {self.n}, a = {self.a}, k = {self.k}, fitness = {self.fitness()}'

    def fitness(self):
        """ Fitness of Genome object
        Returns:
            fitness: fitness of Genome instance """
        fitnesses = [0] * self.n
        for i in range(self.n):
            if self.k != 0:
                fitnesses[i] = sum([*[self.ws[j] for j in k_others[i]], self.ws[i]]) / self.k
            else:
                fitnesses[i] = self.ws[i]

        return sum(fitnesses) / self.n

    def get_mutants(self):
        """ Gets mutants
        Returns:
                mutants: all genotypes one mutation away, taking mutation_porportion into account """
        mutants = []
        for i in range(self.n):
            for j in range(self.a):
                if self.genotype[i] == j or np.random.random() > self.mutation_proportion:
                    continue
                new_genotype = self.genotype.copy()
                new_genotype[i] = j
                mutants.append(new_genotype)

        return mutants


def next_gen(genome_dict):
    """ Iterates to generate next generation
    Arguments:
        genome_dict: dictionary whose elements are (genome object, portion of population)
    Returns:
        new_genome_dict: dictionary following format of genome_dict after evolving through one cycle """
    new_genome_dict = {}
    for old_genome, portion in genome_dict.items():
        # Calculate percentage that survived and didn't mutate
        new_portion_raw = portion * old_genome.fitness()
        new_portion = new_portion_raw * (1 - old_genome.mutation_proportion)
        new_genome_dict[old_genome] = new_portion

        # Mutations
        new_genotypes = old_genome.get_mutants()
        for new_genotype in new_genotypes:
            new_genome_dict[Genome(old_genome.n,
                                   old_genome.a,
                                   old_genome.k,
                                   genotype=new_genotype,
                                   mutation_proportion=old_genome.mutation_proportion
                                   )] = new_portion_raw * old_genome.mutation_proportion / len(new_genotypes)

    deletable = []
    for new_genome, new_portion in new_genome_dict.items():
        if new_portion < 0.01: # pruning threshold
            deletable.append(new_genome)
    for to_be_deleted in deletable:
        del new_genome_dict[to_be_deleted]
    # Rescale
    total_pop = sum([portion for genome, portion in new_genome_dict.items()])
    for genome in new_genome_dict.keys():
        new_genome_dict[genome] /= total_pop

    return new_genome_dict

def avg_fitness(gd):
    """Returns average fitness for genome dictionary"""
    return sum([k.fitness() * v for k, v in gd.items()])

def evolve_over_time(genome_dict, gens):
    """Evolves genome_dict over generations
        Arguments:
            genome_dict: dictionary of genomes
            gens: number of generations over which to evolve
        Returns:
            historical_genome_dict: all genome dictionaries
            fitnesses: average fitness over time"""
    historical_genome_dict = [genome_dict]
    for i in range(gens):
        historical_genome_dict.append(next_gen(historical_genome_dict[-1]))
    fitnesses = [avg_fitness(gd) for gd in historical_genome_dict]

    return historical_genome_dict, fitnesses


def generate_hash_ws(N, A):
    ''' Assigns a list of ws to each genotype hash '''
    d = {}
    for i in range(A ** N):
        d[i] = [random.random() for _ in range(N)]
    return d


def g2h(genotype, A):
    ''' Takes each genotype and outputs its hash '''
    h = 0
    for j, v in enumerate(genotype):
        h += (v * A ** j)
    return h


def generate_k_others(N, K):
    ''' Assigns a (static global) set of K other genes affecting fitness of each gene. '''
    l = []
    for i in range(N):
        band = [j for j in range(N) if i != j]
        choices = random.choices(band, k=K)
        l.append(choices)
    return l


if __name__ == '__main__':
    """ BEGIN PARAMETERS """
    N = 100
    A = 2
    K = 3
    mutation_proportion = 0.05
    generations = 20
    populations = 20
    """ END PARAMETERS"""

    k_others = generate_k_others(N, K)
    hash_ws = generate_hash_ws(N, A)
    from pprint import pprint

    print('K-Others:')
    # pprint(k_others)
    print('Hash -> Ws:')
    # pprint(hash_ws)

    all_fitnesses = []
    for i in range(populations):
        genome = Genome(N, A, K, None, mutation_proportion)
        print(f'run {i} starting with', genome, genome.genotype)
        my_genome_dict = {genome: 1.0}
        historical_genomes, fitnesses = evolve_over_time(my_genome_dict, generations)
        all_fitnesses.append(fitnesses)

    for i in all_fitnesses:
        plt.plot(range(len(fitnesses)), i)
    plt.xlabel('Generation')
    plt.ylabel('Average Fitness')
    plt.show()
