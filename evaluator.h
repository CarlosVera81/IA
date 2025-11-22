#ifndef EVALUATOR_H
#define EVALUATOR_H

#include "instance.h"
#include "solution.h"
#include <vector>

// Calcula la media de los valores en un rectángulo
/**
 * @brief Calcula la media de los valores en un rectángulo.
 * 
 * @param inst Instancia del problema.
 * @param rect Rectángulo a evaluar.
 * @return Media de los valores.
 */
double computeZoneMean(const Instance &inst, const Rect &rect);

// Calcula estadísticas de una zona: media, varianza y SSE
/**
 * @brief Calcula estadísticas de una zona: media, varianza y SSE.
 * 
 * @param inst Instancia del problema.
 * @param rect Rectángulo a evaluar.
 * @param zoneMean Referencia para devolver la media.
 * @param zoneVar Referencia para devolver la varianza.
 * @param zoneSSE Referencia para devolver la Suma de Errores al Cuadrado.
 */
void computeZoneStats(const Instance &inst,
                      const Rect &rect,
                      double &zoneMean,
                      double &zoneVar,
                      double &zoneSSE);

// Verifica que la partición sea válida (sin solapamiento, cobertura total)
/**
 * @brief Verifica que la partición sea válida (sin solapamiento, cobertura total).
 * 
 * @param inst Instancia del problema.
 * @param sol Solución a verificar.
 * @param labelsOut Puntero opcional para devolver la matriz de etiquetas.
 * @return true si es válida, false en caso contrario.
 */
bool checkPartition(const Instance &inst,
                    const Solution &sol,
                    std::vector<std::vector<int>> *labelsOut = nullptr);

// Evalúa una solución calculando SSE total y aplicando penalización si viola homogeneidad
/**
 * @brief Evalúa una solución calculando SSE total y aplicando penalización si viola homogeneidad.
 * 
 * @param inst Instancia del problema.
 * @param sol Solución a evaluar (se actualiza su fitness).
 * @return El fitness calculado.
 */
double evaluateSolution(const Instance &inst, Solution &sol);

// Imprime la matriz de etiquetas de zonas
/**
 * @brief Imprime la matriz de etiquetas de zonas.
 * 
 * @param inst Instancia del problema.
 * @param sol Solución a imprimir.
 */
void printLabelMatrix(const Instance &inst, const Solution &sol);

// Imprime las varianzas de cada zona
/**
 * @brief Imprime las varianzas de cada zona.
 * 
 * @param inst Instancia del problema.
 * @param sol Solución a imprimir.
 */
void printZoneVariances(const Instance &inst, const Solution &sol);

#endif // EVALUATOR_H
