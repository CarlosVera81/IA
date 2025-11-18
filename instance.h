#ifndef INSTANCE_H
#define INSTANCE_H

#include <vector>
#include <string>

struct Instance {
    int N = 0;
    int M = 0;
    int p = 0;        // cantidad de zonas/sensores
    double alpha = 0; // nivel de homogeneidad
    std::vector<std::vector<double>> S;

    double globalMean = 0.0;
    double globalVar  = 0.0;

    bool loadFromFile(const std::string &path);
    void computeGlobalStats();
};

#endif // INSTANCE_H
