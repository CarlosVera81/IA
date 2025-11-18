#include "evaluator.h"
#include <iostream>
#include <algorithm>

const double INF = 1e18;
const double PENALIZACION_GRANDE = 1e12;

double computeZoneMean(const Instance &inst, const Rect &rect) {
    double sum = 0.0;
    int count = 0;
    for (int i = rect.r1; i <= rect.r2; ++i) {
        for (int j = rect.c1; j <= rect.c2; ++j) {
            sum += inst.S[i][j];
            ++count;
        }
    }
    return sum / std::max(1, count);
}

void computeZoneStats(const Instance &inst,const Rect &rect,double &zoneMean,double &zoneVar,double &zoneSSE) {
    zoneMean = computeZoneMean(inst, rect);
    double sse = 0.0;
    int count = 0;
    for (int i = rect.r1; i <= rect.r2; ++i) {
        for (int j = rect.c1; j <= rect.c2; ++j) {
            double diff = inst.S[i][j] - zoneMean;
            sse += diff * diff;
            ++count;
        }
    }
    zoneSSE = sse;
    zoneVar = (count > 0) ? (sse / count) : 0.0;
}

bool checkPartition(const Instance &inst,
                    const Solution &sol,
                    std::vector<std::vector<int>> *labelsOut) {
    int N = inst.N, M = inst.M;
    std::vector<std::vector<int>> labels(N, std::vector<int>(M, -1));

    if ((int)sol.zonas.size() != inst.p) {
        return false;
    }

    for (int k = 0; k < (int)sol.zonas.size(); ++k) {
        const Rect &r = sol.zonas[k];

        // chequeo de límites
        if (r.r1 < 0 || r.c1 < 0 || r.r2 >= N || r.c2 >= M ||
            r.r1 > r.r2 || r.c1 > r.c2) {
            return false;
        }

        // chequeo de solapamiento
        for (int i = r.r1; i <= r.r2; ++i) {
            for (int j = r.c1; j <= r.c2; ++j) {
                if (labels[i][j] != -1) {
                    // ya estaba cubierto → solape
                    return false;
                }
                labels[i][j] = k;
            }
        }
    }

    // chequeo de cobertura exacta: todas las celdas deben estar cubiertas
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            if (labels[i][j] == -1) {
                return false;
            }
        }
    }

    if (labelsOut) {
        *labelsOut = std::move(labels);
    }
    return true;
}

double evaluateSolution(const Instance &inst, Solution &sol) {

    std::vector<std::vector<int>> labels;
    if (!checkPartition(inst, sol, &labels)) {
        sol.fitness = INF;
        return sol.fitness;
    }

    double totalSSE = 0.0;
    bool homogeneityOK = true;

    for (const Rect &r : sol.zonas) {
        double zoneMean = 0.0;
        double zoneVar  = 0.0;
        double zoneSSE  = 0.0;
        computeZoneStats(inst, r, zoneMean, zoneVar, zoneSSE);
        totalSSE += zoneSSE;
        if (zoneVar > inst.alpha * inst.globalVar + 1e-9) {
            homogeneityOK = false;
        }
    }

    if (!homogeneityOK) {
        totalSSE += PENALIZACION_GRANDE;
    }
    sol.fitness = totalSSE;
    return sol.fitness;
}

void printLabelMatrix(const Instance &inst, const Solution &sol) {
    std::vector<std::vector<int>> labels;
    if (!checkPartition(inst, sol, &labels)) {
        std::cout << "Solución no es una partición válida.\n";
        return;
    }

    std::cout << "Matriz Z (etiquetas de zonas):\n";
    for (int i = 0; i < inst.N; ++i) {
        for (int j = 0; j < inst.M; ++j) {
            std::cout << labels[i][j] + 1 << " "; 
        }
        std::cout << "\n";
    }
}

void printZoneVariances(const Instance &inst, const Solution &sol) {
    std::cout << "Varianza intra-zona de la mejor solución:\n";

    for (size_t k = 0; k < sol.zonas.size(); ++k) {
        const Rect &r = sol.zonas[k];

        double zoneMean = 0.0;
        double zoneVar  = 0.0;
        double zoneSSE  = 0.0;

        computeZoneStats(inst, r, zoneMean, zoneVar, zoneSSE);

        std::cout << "  Zona " << (k + 1)
                  << " -> media = " << zoneMean
                  << ", varianza = " << zoneVar
                  << ", SSE = " << zoneSSE
                  << "\n";
    }
}
