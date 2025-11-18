#ifndef UTILS_H
#define UTILS_H

#include "genetic_algorithm.h"
#include "instance.h"
#include <string>
#include <random>

// Carga parámetros del GA desde archivo de configuración
bool loadGAConfig(const std::string &configPath, GAParams &params);

// Solicita al usuario el valor de p (cantidad de zonas)
int pedirP();

// Solicita al usuario el valor de alpha (nivel de homogeneidad)
double pedirAlpha();

// Procesa todas las instancias .spp en una carpeta dada
void processInstancesInFolder(const std::string &folderPath,
                              int p,
                              double alpha,
                              const GAParams &params,
                              std::mt19937 &rng);

#endif // UTILS_H
