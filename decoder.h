#ifndef DECODER_H
#define DECODER_H

#include "instance.h"
#include "solution.h"
#include <random>
#include <vector>

// Decodifica un genoma (secuencia de genes) a rectángulos concretos
void decodeGenomeToRects(const Instance &inst,
                         const std::vector<Gene> &genome,
                         std::vector<Rect> &rects,
                         std::mt19937 &rng);

// Par de índices para rectángulos adyacentes
using RectPair = std::pair<int, int>;

// Busca pares de rectángulos adyacentes verticalmente
std::vector<RectPair> findVerticalAdjacencies(const Solution &sol);

// Busca pares de rectángulos adyacentes horizontalmente
std::vector<RectPair> findHorizontalAdjacencies(const Solution &sol);

#endif // DECODER_H
