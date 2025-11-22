#ifndef DECODER_H
#define DECODER_H

#include "instance.h"
#include "solution.h"
#include <random>
#include <vector>

// Decodifica un genoma (secuencia de genes) a rectángulos concretos
/**
 * @brief Decodifica un genoma (secuencia de genes) a rectángulos concretos.
 * 
 * @param inst Instancia del problema.
 * @param genome Genoma a decodificar.
 * @param rects Vector donde se almacenarán los rectángulos resultantes.
 * @param rng Generador de números aleatorios.
 */
void decodeGenomeToRects(const Instance &inst,
                         const std::vector<Gene> &genome,
                         std::vector<Rect> &rects,
                         std::mt19937 &rng);

// Par de índices para rectángulos adyacentes
using RectPair = std::pair<int, int>;

// Busca pares de rectángulos adyacentes verticalmente
/**
 * @brief Busca pares de rectángulos adyacentes verticalmente.
 * 
 * @param sol Solución a analizar.
 * @return Vector de pares de índices de rectángulos adyacentes.
 */
std::vector<RectPair> findVerticalAdjacencies(const Solution &sol);

// Busca pares de rectángulos adyacentes horizontalmente
/**
 * @brief Busca pares de rectángulos adyacentes horizontalmente.
 * 
 * @param sol Solución a analizar.
 * @return Vector de pares de índices de rectángulos adyacentes.
 */
std::vector<RectPair> findHorizontalAdjacencies(const Solution &sol);

#endif // DECODER_H
