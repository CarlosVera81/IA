#include "instance.h"
#include <fstream>
#include <iostream>

bool Instance::loadFromFile(const std::string &path) {
    std::ifstream in(path);
    if (!in.is_open()) {
        std::cerr << "Error al abrir archivo de instancia: " << path << "\n";
        return false;
    }

    in >> N >> M;
    if (!in) {
        std::cerr << "Error leyendo N y M en " << path << "\n";
        return false;
    }

    S.assign(N, std::vector<double>(M, 0.0));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            in >> S[i][j];
            if (!in) {
                std::cerr << "Error leyendo S[" << i << "][" << j << "] en " << path << "\n";
                return false;
            }
        }
    }

    computeGlobalStats();
    return true;
}

void Instance::computeGlobalStats() {
    int totalCells = N * M;
    double sum = 0.0;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < M; ++j)
            sum += S[i][j];

    globalMean = sum / totalCells;

    double sumSq = 0.0;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < M; ++j) {
            double diff = S[i][j] - globalMean;
            sumSq += diff * diff;
        }

    globalVar = sumSq / totalCells;
}
