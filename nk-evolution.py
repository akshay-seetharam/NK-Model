import numpy as np
import random
from matplotlib import pyplot as plt


class Genome:
    def __init__(self, n, a, k, genotype=None):
        """ Constructor for state
        Arguments:
            n: length of genome
            a: options for choice at each locus in genome
            k: number of alleles that affect fitness from one allele
            genotype: predefined genotype, None by default
        Returns:
            New instance of State """
        self.n = n
        self.a = a
        self.k = k
        self.genotype = [random.randint(0, a - 1) for _ in range(n)] if genotype is None else genotype
        self.k_alleles = []
        for i in range(n):
            self.k_alleles.append([])
            for j in range(k):
                other_allele = random.randint(0, n - 1)
                while other_allele == i or other_allele in self.k_alleles[-1]:
                    other_allele = random.randint(0, n - 1)
                self.k_alleles[-1].append(other_allele)
        self.ws = [np.random.random() for _ in range(n)]

    def __str__(self):
        return f'Genome with n = {self.n}, a = {self.a}, k = {self.k}, fitness = {self.fitness()}'

    def fitness(self):
        """ Fitness of Genome object
        Returns:
            fitness: fitness of Genome instance """
        fitnesses = [0] * self.n
        for i in range(self.n):
            if self.k != 0:
                fitnesses[i] = sum([*[self.ws[j] for j in self.k_alleles[i]], self.ws[i]]) / self.k
            else:
                fitnesses[i] = self.ws[i]

        return sum(fitnesses) / self.n

    def get_mutants(self):
        """ Gets mutants
        Returns:
                mutants: all genotypes one mutation away """
        mutants = []
        for i in range(self.n):
            for j in range(self.a):
                if self.genotype[i] == j:
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
        new_genotypes = old_genome.get_mutants()
        for new_genotype in new_genotypes:
            new_genome_dict[Genome(old_genome.n, old_genome.a, old_genome.k, genotype=new_genotype)] = 1.0 * portion / len(new_genotypes)

    return new_genome_dict

if __name__ == '__main__':
    """ BEGIN PARAMETERS """
    N = 2
    A = 2
    K = 1
    generations = 1000
    """ END PARAMETERS"""

    genome = Genome(N, A, K)
    print(genome, genome.genotype)
    my_genome_dict = {genome: 1.0}
    print([(i.genotype, j) for i, j in next_gen(my_genome_dict).items()])
