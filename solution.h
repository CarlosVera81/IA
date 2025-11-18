#ifndef SOLUTION_H
#define SOLUTION_H

#include <vector>


struct Gene {
    bool horizontal; // orientación
    double gamma;    // posición
    double beta;     // selección
};

// Rectángulo en la grilla
struct Rect {
    int r1, c1; // esquina superior izquierda
    int r2, c2; // esquina inferior derecha
};


struct Solution {
    std::vector<Gene> genome; 
    std::vector<Rect> zonas;  
    double fitness = 1e18;    
};

#endif // SOLUTION_H
