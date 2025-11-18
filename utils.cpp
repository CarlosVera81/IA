#include "utils.h"
#include "evaluator.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <filesystem>

namespace fs = std::filesystem;

bool loadGAConfig(const std::string &configPath, GAParams &params) {
    std::ifstream in(configPath);
    if (!in.is_open()) {
        std::cerr << "No se pudo abrir el archivo de configuración: " << configPath << "\n";
        std::cerr << "Usando parámetros por defecto.\n";
        return false;
    }

    std::string line;
    while (std::getline(in, line)) {
        // Ignorar líneas vacías y comentarios
        if (line.empty() || line[0] == '#') continue;

        // Buscar '='
        size_t pos = line.find('=');
        if (pos == std::string::npos) continue;

        std::string key = line.substr(0, pos);
        std::string value = line.substr(pos + 1);

        // Quitar espacios en blanco
        auto trim = [](std::string &s) {
            s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
                return !std::isspace(ch);
            }));
            s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
                return !std::isspace(ch);
            }).base(), s.end());
        };

        trim(key);
        trim(value);

        // Asignar valores
        if (key == "populationSize") {
            params.populationSize = std::stoi(value);
        } else if (key == "generations") {
            params.generations = std::stoi(value);
        } else if (key == "tournamentK") {
            params.tournamentK = std::stoi(value);
        } else if (key == "crossoverRate") {
            params.crossoverRate = std::stod(value);
        } else if (key == "mutationRate") {
            params.mutationRate = std::stod(value);
        }
    }

    std::cout << "Configuración del GA cargada desde " << configPath << ":\n";
    std::cout << "  populationSize = " << params.populationSize << "\n";
    std::cout << "  generations = " << params.generations << "\n";
    std::cout << "  tournamentK = " << params.tournamentK << "\n";
    std::cout << "  crossoverRate = " << params.crossoverRate << "\n";
    std::cout << "  mutationRate = " << params.mutationRate << "\n\n";

    return true;
}

int pedirP() {
    while (true) {
        std::cout << "Ingrese p (cantidad de zonas/sensores, entero positivo): ";
        std::string line;
        if (!std::getline(std::cin, line)) {
            std::cerr << "Error leyendo entrada. Saliendo.\n";
            exit(1);
        }

        // Quitar espacios al inicio/fin
        auto notSpace = [](int ch){ return !std::isspace(ch); };
        line.erase(line.begin(), std::find_if(line.begin(), line.end(), notSpace));
        line.erase(std::find_if(line.rbegin(), line.rend(), notSpace).base(), line.end());

        if (line.empty()) {
            std::cout << "p no puede ser vacío. Intente nuevamente.\n";
            continue;
        }

        std::stringstream ss(line);
        int p;
        char extra;
        if (!(ss >> p) || (ss >> extra)) {
            std::cout << "Entrada inválida. p debe ser un entero positivo.\n";
            continue;
        }

        if (p <= 0) {
            std::cout << "p debe ser un entero positivo (> 0).\n";
            continue;
        }

        return p;
    }
}

double pedirAlpha() {
    while (true) {
        std::cout << "Ingrese alpha (nivel de homogeneidad, 0 < alpha <= 1): ";
        std::string line;
        if (!std::getline(std::cin, line)) {
            std::cerr << "Error leyendo entrada. Saliendo.\n";
            exit(1);
        }

        // Quitar espacios al inicio/fin
        auto notSpace = [](int ch){ return !std::isspace(ch); };
        line.erase(line.begin(), std::find_if(line.begin(), line.end(), notSpace));
        line.erase(std::find_if(line.rbegin(), line.rend(), notSpace).base(), line.end());

        if (line.empty()) {
            std::cout << "alpha no puede ser vacío. Intente nuevamente.\n";
            continue;
        }

        std::stringstream ss(line);
        double alpha;
        char extra;
        if (!(ss >> alpha) || (ss >> extra)) {
            std::cout << "Entrada inválida. alpha debe ser un número real.\n";
            continue;
        }

        if (alpha <= 0.0 || alpha > 1.0) {
            std::cout << "alpha debe cumplir 0 < alpha <= 1.\n";
            continue;
        }

        return alpha;
    }
}

void processInstancesInFolder(const std::string &folderPath,
                              int p,
                              double alpha,
                              const GAParams &params,
                              std::mt19937 &rng) {
    if (!fs::exists(folderPath) || !fs::is_directory(folderPath)) {
        std::cerr << "La ruta '" << folderPath << "' no existe o no es un directorio.\n";
        return;
    }

    std::vector<std::string> files;
    for (const auto &entry : fs::directory_iterator(folderPath)) {
        if (entry.is_regular_file() && entry.path().extension() == ".spp") {
            files.push_back(entry.path().string());
        }
    }

    if (files.empty()) {
        std::cerr << "No se encontraron archivos .spp en la carpeta '" << folderPath << "'.\n";
        return;
    }

    std::sort(files.begin(), files.end());

    std::cout << "\nSe encontraron " << files.size()
              << " instancias .spp en '" << folderPath << "'.\n";

    for (const std::string &filePath : files) {
        std::cout << "\n========================================\n";
        std::cout << "Procesando instancia: " << filePath << "\n";

        Instance inst;
        if (!inst.loadFromFile(filePath)) {
            std::cerr << "  -> Error cargando la instancia. Se omite.\n";
            continue;
        }

        inst.p = p;
        inst.alpha = alpha;

        std::cout << "  N = " << inst.N << ", M = " << inst.M
                  << ", p = " << inst.p << ", alpha = " << inst.alpha << "\n";
        std::cout << "  Varianza global = " << inst.globalVar << "\n";

        // Ejecutar el GA
        Solution best = runGA(inst, params, rng);

        std::cout << "  Mejor fitness encontrado = " << best.fitness << "\n";
        std::cout << "  Matriz Z (etiquetas de zonas) del mejor individuo:\n";
        printLabelMatrix(inst, best);
        printZoneVariances(inst, best);
    }
}
