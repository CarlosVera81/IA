#include "utils.h"
#include "genetic_algorithm.h"
#include <iostream>
#include <random>

/**
 * @brief Punto de entrada principal del programa.
 * 
 * Carga la configuración, solicita parámetros al usuario y procesa las instancias
 * ubicadas en la carpeta especificada.
 * 
 * @return 0 si la ejecución fue exitosa.
 */
int main() {
    // Cargar configuración del GA desde archivo
    GAParams gaParams;
    loadGAConfig("ga_config.txt", gaParams);

    // Solicitar parámetros del problema
    int p = pedirP();
    double alpha = pedirAlpha();

    // Ruta de la carpeta de instancias (desde config)
    std::string folderPath = gaParams.instancePath;
    
    // Opcional: descomentar para pedir la ruta al usuario
    /*
    std::cout << "Ingrese ruta de la carpeta de instancias "
              << "(por ejemplo, Instancias/Grandes): ";
    if (!std::getline(std::cin, folderPath)) {
        std::cerr << "Error leyendo la ruta de la carpeta.\n";
        return 1;
    }
    
    if (folderPath.empty()) {
        std::cerr << "Error: la ruta de la carpeta no puede ser vacía.\n";
        return 1;
    }
    */

    // Inicializar generador de números aleatorios
    std::random_device rd;
    std::mt19937 rng(rd());

    // Procesar todas las instancias en la carpeta
    processInstancesInFolder(folderPath, p, alpha, gaParams, rng);

    return 0;
}
