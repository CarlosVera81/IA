#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

#include "instance.h"
#include "solution.h"
#include <random>

// Parámetros del algoritmo genético
struct GAParams {
    int populationSize = 50;   ///< Tamaño de la población
    int generations    = 200;  ///< Número de generaciones
    int tournamentK    = 3;    ///< Tamaño del torneo para selección
    double crossoverRate = 0.8;///< Probabilidad de cruce
    double mutationRate  = 0.3;///< Probabilidad de mutación
    std::string instancePath = "Instancias/Pequeñas"; ///< Ruta de las instancias a procesar
};

// Selección por torneo
/**
 * @brief Selecciona un individuo de la población mediante torneo.
 * 
 * @param pop Población actual.
 * @param tournamentK Tamaño del torneo.
 * @param rng Generador de números aleatorios.
 * @return Referencia constante al individuo seleccionado.
 */
const Solution& tournamentSelect(const std::vector<Solution> &pop,
                                 int tournamentK,
                                 std::mt19937 &rng);

// Operador de cruce
/**
 * @brief Realiza el cruce entre dos padres para generar un hijo.
 * 
 * @param inst Instancia del problema.
 * @param p1 Primer padre.
 * @param p2 Segundo padre.
 * @param rng Generador de números aleatorios.
 * @return Nueva solución hija.
 */
Solution crossover(const Instance &inst,
                   const Solution &p1,
                   const Solution &p2,
                   std::mt19937 &rng);

// Operador de mutación
/**
 * @brief Aplica mutación a una solución.
 * 
 * @param inst Instancia del problema.
 * @param original Solución original.
 * @param rng Generador de números aleatorios.
 * @return Nueva solución mutada.
 */
Solution mutate(const Instance &inst,
                const Solution &original,
                std::mt19937 &rng);

/**
 * @brief Genera una solución inicial aleatoria basada en franjas (strips).
 * 
 * @param inst Instancia del problema.
 * @param rng Generador de números aleatorios.
 * @return Solución generada.
 */
Solution generateRandomStripSolution(const Instance &inst, std::mt19937 &rng);

// Ejecuta el algoritmo genético completo
/**
 * @brief Ejecuta el Algoritmo Genético completo.
 * 
 * @param inst Instancia del problema.
 * @param params Parámetros del GA.
 * @param rng Generador de números aleatorios.
 * @return La mejor solución encontrada.
 */
Solution runGA(const Instance &inst,
               const GAParams &params,
               std::mt19937 &rng);

#endif // GENETIC_ALGORITHM_H
