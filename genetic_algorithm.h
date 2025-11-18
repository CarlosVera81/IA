#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

#include "instance.h"
#include "solution.h"
#include <random>

// Parámetros del algoritmo genético
struct GAParams {
    int populationSize = 50;
    int generations    = 200;
    int tournamentK    = 3;
    double crossoverRate = 0.8;
    double mutationRate  = 0.3;
};

// Selección por torneo
const Solution& tournamentSelect(const std::vector<Solution> &pop,
                                 int tournamentK,
                                 std::mt19937 &rng);

// Operador de cruce
Solution crossover(const Instance &inst,
                   const Solution &p1,
                   const Solution &p2,
                   std::mt19937 &rng);

// Operador de mutación
Solution mutate(const Instance &inst,
                const Solution &original,
                std::mt19937 &rng);

// Genera una solución aleatoria inicial
Solution generateRandomStripSolution(const Instance &inst, std::mt19937 &rng);

// Ejecuta el algoritmo genético completo
Solution runGA(const Instance &inst,
               const GAParams &params,
               std::mt19937 &rng);

#endif // GENETIC_ALGORITHM_H
