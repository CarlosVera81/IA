#ifndef INSTANCE_H
#define INSTANCE_H

#include <vector>
#include <string>

struct Instance {
    int N = 0;
    int M = 0;
    int p = 0;        // cantidad de zonas/sensores
    double alpha = 0; // nivel de homogeneidad
    std::vector<std::vector<double>> S; ///< Matriz de datos de la instancia

    double globalMean = 0.0;
    double globalVar  = 0.0;

    /**
     * @brief Carga una instancia desde un archivo de texto.
     * 
     * @param path Ruta al archivo de la instancia.
     * @return true si la carga fue exitosa, false en caso contrario.
     */
    bool loadFromFile(const std::string &path);

    /**
     * @brief Calcula estadísticas globales de la instancia (media y varianza).
     */
    void computeGlobalStats();
};

#endif // INSTANCE_H
