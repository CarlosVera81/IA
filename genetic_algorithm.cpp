#include "genetic_algorithm.h"
#include "decoder.h"
#include "evaluator.h"
#include <algorithm>
#include <iostream>

const Solution& tournamentSelect(const std::vector<Solution> &pop,
                                 int tournamentK,
                                 std::mt19937 &rng) {
    std::uniform_int_distribution<int> dist(0, (int)pop.size() - 1);

    int bestIdx = dist(rng);
    double bestFit = pop[bestIdx].fitness;

    for (int i = 1; i < tournamentK; ++i) {
        int idx = dist(rng);
        if (pop[idx].fitness < bestFit) {
            bestFit = pop[idx].fitness;
            bestIdx = idx;
        }
    }
    return pop[bestIdx];
}

Solution crossover(const Instance &inst,
                   const Solution &p1,
                   const Solution &p2,
                   std::mt19937 &rng) {
    int L = (int)p1.genome.size();

    if (L < 2 || (int)p2.genome.size() != L) {
        return (p1.fitness <= p2.fitness) ? p1 : p2;
    }
    Solution child;
    child.genome.resize(L);
    std::uniform_int_distribution<int> cutDist(1, L - 1);
    int cut = cutDist(rng);
    for (int i = 0; i < cut; ++i) {
        child.genome[i] = p1.genome[i];
    }
    for (int i = cut; i < L; ++i) {
        child.genome[i] = p2.genome[i];
    }
    decodeGenomeToRects(inst, child.genome, child.zonas, rng);
    evaluateSolution(inst, child);
    return child;
}

Solution mutate(const Instance &inst,const Solution &original,std::mt19937 &rng, double strongMutationRate) {
    Solution mutated = original;
    std::uniform_real_distribution<double> prob01(0.0, 1.0);
    std::normal_distribution<double> gauss(0.0, 0.1); 
    std::bernoulli_distribution flipOrient(0.1);
    std::bernoulli_distribution newOrient(0.5);
    std::uniform_real_distribution<double> newGamma(0.2, 0.8);
    std::uniform_real_distribution<double> newBeta(0.0, 1.0);
    bool strongMutation = (prob01(rng) < strongMutationRate); 
    if (mutated.genome.empty()) {
        return generateRandomStripSolution(inst, rng);
    }
    int L = (int)mutated.genome.size();
    if (strongMutation) {
        mutated.genome.resize(L);
        for (int i = 0; i < L; ++i) {
            mutated.genome[i].horizontal = newOrient(rng);
            mutated.genome[i].gamma      = newGamma(rng);
            mutated.genome[i].beta       = newBeta(rng);
        }
    } else {
        int g = std::uniform_int_distribution<int>(0, L-1)(rng);
        mutated.genome[g].gamma += gauss(rng);
        mutated.genome[g].beta  += gauss(rng);
        if (flipOrient(rng)) mutated.genome[g].horizontal = !mutated.genome[g].horizontal;
    }
    decodeGenomeToRects(inst, mutated.genome, mutated.zonas, rng);
    evaluateSolution(inst, mutated);
    return mutated;
}

Solution generateRandomStripSolution(const Instance &inst, std::mt19937 &rng) {
    Solution sol;
    sol.genome.clear();
    sol.zonas.clear();

    int p = inst.p;
    int totalCells = inst.N * inst.M;

    if (p <= 0 || p > totalCells) {
        sol.fitness = 1e18;
        return sol;
    }

    int numGenes = p - 1;
    sol.genome.resize(numGenes);

    std::bernoulli_distribution orientation(0.5);
    std::uniform_real_distribution<double> gammaDist(0.2, 0.8);
    std::uniform_real_distribution<double> betaDist(0.0, 1.0);

    for (int g = 0; g < numGenes; ++g) {
        sol.genome[g].horizontal = orientation(rng);
        sol.genome[g].gamma      = gammaDist(rng);
        sol.genome[g].beta       = betaDist(rng);
    }
    decodeGenomeToRects(inst, sol.genome, sol.zonas, rng);
    evaluateSolution(inst, sol);
    return sol;
}

Solution runGA(const Instance &inst,
               const GAParams &params,
               std::mt19937 &rng) {
    std::vector<Solution> population;
    population.reserve(params.populationSize);

   
    for (int i = 0; i < params.populationSize; ++i) {
        Solution s = generateRandomStripSolution(inst, rng);
        population.push_back(std::move(s));
    }

    // Mejor de todos (para tracking global)
    Solution bestOverall = population[0];
    for (const auto &s : population) {
        if (s.fitness < bestOverall.fitness) {
            bestOverall = s;
        }
    }

    std::uniform_real_distribution<double> prob01(0.0, 1.0);

    // ---- Bucle principal GA ----
    for (int gen = 0; gen < params.generations; ++gen) {
        // Ordenar población por fitness (mejor primero)
        std::sort(population.begin(), population.end(),
                  [](const Solution &a, const Solution &b) {
                      return a.fitness < b.fitness;
                  });

        // Elitismo: mantenemos el mejor
        Solution elite = population[0];
        if (elite.fitness < bestOverall.fitness) {
            bestOverall = elite;
        }

        std::vector<Solution> newPop;
        newPop.reserve(params.populationSize);
        newPop.push_back(elite); // copiamos el elitista

        // Rellenamos el resto de la población
        while ((int)newPop.size() < params.populationSize) {
            const Solution &parent1 = tournamentSelect(population, params.tournamentK, rng);
            const Solution &parent2 = tournamentSelect(population, params.tournamentK, rng);

            Solution child;

            // Crossover
            if (prob01(rng) < params.crossoverRate) {
                child = crossover(inst, parent1, parent2, rng);
            } else {
                // Sin crossover: clon del padre1
                child = parent1;
            }

            // Mutación
            if (prob01(rng) < params.mutationRate) {
                child = mutate(inst, child, rng, params.strongMutationRate);
            }

            newPop.push_back(std::move(child));
        }

        population = std::move(newPop);

        // Log opcional
        if (gen % 100 == 0) {
            std::cout << "[GA] Gen " << gen
                      << " mejor fitness = " << bestOverall.fitness << "\n";
        }
    }

    return bestOverall;
}
