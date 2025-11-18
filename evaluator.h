#ifndef EVALUATOR_H
#define EVALUATOR_H

#include "instance.h"
#include "solution.h"
#include <vector>

// Calcula la media de los valores en un rectángulo
double computeZoneMean(const Instance &inst, const Rect &rect);

// Calcula estadísticas de una zona: media, varianza y SSE
void computeZoneStats(const Instance &inst,
                      const Rect &rect,
                      double &zoneMean,
                      double &zoneVar,
                      double &zoneSSE);

// Verifica que la partición sea válida (sin solapamiento, cobertura total)
bool checkPartition(const Instance &inst,
                    const Solution &sol,
                    std::vector<std::vector<int>> *labelsOut = nullptr);

// Evalúa una solución calculando SSE total y aplicando penalización si viola homogeneidad
double evaluateSolution(const Instance &inst, Solution &sol);

// Imprime la matriz de etiquetas de zonas
void printLabelMatrix(const Instance &inst, const Solution &sol);

// Imprime las varianzas de cada zona
void printZoneVariances(const Instance &inst, const Solution &sol);

#endif // EVALUATOR_H
