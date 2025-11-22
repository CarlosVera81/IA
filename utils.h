#ifndef UTILS_H
#define UTILS_H

#include "genetic_algorithm.h"
#include "instance.h"
#include <string>
#include <random>

/**
 * @brief Carga los parámetros del Algoritmo Genético desde un archivo de configuración.
 * 
 * @param configPath Ruta al archivo de configuración.
 * @param params Estructura donde se almacenarán los parámetros cargados.
 * @return true si la carga fue exitosa, false en caso contrario.
 */
bool loadGAConfig(const std::string &configPath, GAParams &params);

/**
 * @brief Solicita al usuario el valor de p (cantidad de zonas).
 * 
 * @return El valor entero de p ingresado por el usuario.
 */
int pedirP();

/**
 * @brief Solicita al usuario el valor de alpha (nivel de homogeneidad).
 * 
 * @return El valor double de alpha ingresado por el usuario.
 */
double pedirAlpha();

/**
 * @brief Procesa todas las instancias .spp encontradas en una carpeta.
 * 
 * @param folderPath Ruta de la carpeta que contiene las instancias.
 * @param p Cantidad de zonas deseadas.
 * @param alpha Parámetro de homogeneidad.
 * @param params Parámetros del algoritmo genético.
 * @param rng Generador de números aleatorios.
 */
void processInstancesInFolder(const std::string &folderPath,
                              int p,
                              double alpha,
                              const GAParams &params,
                              std::mt19937 &rng);

#endif // UTILS_H
