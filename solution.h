#ifndef SOLUTION_H
#define SOLUTION_H

#include <vector>


/**
 * @brief Representa un gen en el algoritmo genético.
 */
struct Gene {
    bool horizontal; ///< Orientación del corte (true: horizontal, false: vertical)
    double gamma;    ///< Posición relativa del corte
    double beta;     ///< Selección de sub-rectángulo
};

// Rectángulo en la grilla
/**
 * @brief Representa un rectángulo en la grilla.
 */
struct Rect {
    int r1, c1; ///< Fila y columna de la esquina superior izquierda
    int r2, c2; ///< Fila y columna de la esquina inferior derecha
};


/**
 * @brief Representa una solución completa (individuo).
 */
struct Solution {
    std::vector<Gene> genome; ///< Genoma del individuo
    std::vector<Rect> zonas;  ///< Zonas decodificadas
    double fitness = 1e18;    ///< Valor de fitness (Error Total)
};

#endif // SOLUTION_H
