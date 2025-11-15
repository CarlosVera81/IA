#include <bits/stdc++.h>
using namespace std;
namespace fs = std::filesystem;

const double INF = 1e18;
const double PENALIZACION_GRANDE = 1e12;


struct GAParams {
    int populationSize = 50;    // tamaño de población
    int generations    = 200;   // número de generaciones
    int tournamentK    = 3;     // tamaño de torneo para selección
    double crossoverRate = 0.8; // prob. de aplicar crossover
    double mutationRate  = 0.3; // prob. de mutar un hijo
};




// -------------------------
// Datos de la instancia
// -------------------------
struct Instance {
    int N = 0;
    int M = 0;
    int p = 0;        // se fijan DESPUÉS por consola
    double alpha = 0; // idem
    vector<vector<double>> S;

    double globalMean = 0.0;
    double globalVar  = 0.0;

    bool loadFromFile(const string &path) {
        ifstream in(path);
        if (!in.is_open()) {
            cerr << "Error al abrir archivo de instancia: " << path << "\n";
            return false;
        }

        in >> N >> M;
        if (!in) {
            cerr << "Error leyendo N y M en " << path << "\n";
            return false;
        }

        S.assign(N, vector<double>(M, 0.0));
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                in >> S[i][j];
                if (!in) {
                    cerr << "Error leyendo S[" << i << "][" << j << "] en " << path << "\n";
                    return false;
                }
            }
        }

        computeGlobalStats();
        return true;
    }

    void computeGlobalStats() {
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
};



// -------------------------
// Representación de solución
// -------------------------
// Gen que representa un corte guillotinable
struct Gene {
    bool horizontal; // true = corte horizontal, false = vertical
    double alpha;    // posición relativa del corte en [0,1]
    double beta; // nuevo

};

struct Rect {
    int r1, c1; // esquina superior izquierda
    int r2, c2; // esquina inferior derecha
};

// Ahora la solución tiene "genoma" (secuencia de cortes) + rectángulos decodificados
struct Solution {
    vector<Gene> genome; // tamaño = p-1
    vector<Rect> zonas;  // rectángulos resultantes de aplicar el genoma
    double fitness = INF;
};

const Solution& tournamentSelect(const vector<Solution> &pop,
                                 int tournamentK,
                                 mt19937 &rng) {
    std::uniform_int_distribution<int> dist(0, (int)pop.size() - 1);

    int bestIdx = dist(rng);
    double bestFit = pop[bestIdx].fitness;

    for (int i = 1; i < tournamentK; ++i) {
        int idx = dist(rng);
        if (pop[idx].fitness < bestFit) {
            bestFit = pop[idx].fitness;
            bestIdx = idx;
        }
    }
    return pop[bestIdx];
}


// -------------------------
// Utilidades para evaluación
// -------------------------

// Calcula media de los valores en un rectángulo
double computeZoneMean(const Instance &inst, const Rect &rect) {
    double sum = 0.0;
    int count = 0;
    for (int i = rect.r1; i <= rect.r2; ++i) {
        for (int j = rect.c1; j <= rect.c2; ++j) {
            sum += inst.S[i][j];
            ++count;
        }
    }
    return sum / max(1, count);
}

// Calcula SSE y varianza de la zona
void computeZoneStats(const Instance &inst,
                      const Rect &rect,
                      double &zoneMean,
                      double &zoneVar,
                      double &zoneSSE) {
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

// Pares de índices de rectángulos

// -------------------------------------------------------------
// Nuevo: decodificar un genoma (secuencia de cortes) a rectángulos
// -------------------------------------------------------------
void decodeGenomeToRects(const Instance &inst,
                         const vector<Gene> &genome,
                         vector<Rect> &rects,
                         mt19937 &rng)
{
    rects.clear();
    // Empezamos con el rectángulo completo
    rects.push_back(Rect{0, 0, inst.N - 1, inst.M - 1});

    // Aplicamos cada gen en orden
    for (const Gene &g : genome) {
        // Buscar rectángulo con mayor área que aún se pueda partir
        int bestIdx = -1;
        int bestArea = -1;
        for (int i = 0; i < (int)rects.size(); ++i) {
            int h = rects[i].r2 - rects[i].r1 + 1;
            int w = rects[i].c2 - rects[i].c1 + 1;
            int area = h * w;
            if (area > bestArea && area >= 2) {
                bestArea = area;
                bestIdx = i;
            }
        }
        if (bestIdx == -1) break; // no se puede seguir partiendo

        Rect r = rects[bestIdx];
        int h = r.r2 - r.r1 + 1;
        int w = r.c2 - r.c1 + 1;

        bool horizontal = g.horizontal;
        // Si orientación pedida no es viable, forzamos la otra
        if (horizontal && h < 2 && w >= 2) horizontal = false;
        if (!horizontal && w < 2 && h >= 2) horizontal = true;

        Rect rA, rB;

        if (horizontal && h >= 2) {
            int minRow = r.r1;
            int maxRow = r.r2 - 1; // al menos 1 fila en cada lado

            double alpha = std::max(0.0, std::min(1.0, g.alpha));
            int offset = (int)std::round(alpha * (h - 1));
            int cut = minRow + offset;
            cut = std::max(minRow, std::min(cut, maxRow));

            rA = Rect{r.r1, r.c1, cut,     r.c2};
            rB = Rect{cut + 1, r.c1, r.r2, r.c2};
        } else if (!horizontal && w >= 2) {
            int minCol = r.c1;
            int maxCol = r.c2 - 1; // al menos 1 col en cada lado

            double alpha = std::max(0.0, std::min(1.0, g.alpha));
            int offset = (int)std::round(alpha * (w - 1));
            int cut = minCol + offset;
            cut = std::max(minCol, std::min(cut, maxCol));

            rA = Rect{r.r1, r.c1, r.r2, cut};
            rB = Rect{r.r1, cut + 1, r.r2, r.c2};
        } else {
            // Caso patológico: no se puede cortar en ninguna dirección
            continue;
        }

        rects[bestIdx] = rA;
        rects.push_back(rB);
    }

    // Si por redondeos no llegamos a p rectángulos, completamos partiendo al azar
    std::uniform_int_distribution<int> coin(0, 1);
    while ((int)rects.size() < inst.p) {
        int idx = -1;
        for (int i = 0; i < (int)rects.size(); ++i) {
            int h = rects[i].r2 - rects[i].r1 + 1;
            int w = rects[i].c2 - rects[i].c1 + 1;
            if (h * w >= 2) {
                idx = i;
                break;
            }
        }
        if (idx == -1) break;

        Rect r = rects[idx];
        int h = r.r2 - r.r1 + 1;
        int w = r.c2 - r.c1 + 1;

        bool splitHor;
        if (h >= 2 && w >= 2) {
            splitHor = (coin(rng) == 0);
        } else if (h >= 2) {
            splitHor = true;
        } else if (w >= 2) {
            splitHor = false;
        } else {
            break;
        }

        Rect rA, rB;
        if (splitHor) {
            std::uniform_int_distribution<int> distRow(r.r1, r.r2 - 1);
            int cut = distRow(rng);
            rA = Rect{r.r1, r.c1, cut,     r.c2};
            rB = Rect{cut + 1, r.c1, r.r2, r.c2};
        } else {
            std::uniform_int_distribution<int> distCol(r.c1, r.c2 - 1);
            int cut = distCol(rng);
            rA = Rect{r.r1, r.c1, r.r2, cut};
            rB = Rect{r.r1, cut + 1, r.r2, r.c2};
        }

        rects[idx] = rA;
        rects.push_back(rB);
    }
}


using RectPair = pair<int,int>;

// Busca pares de rectángulos adyacentes verticalmente:
// mismo rango de filas y columnas contiguas.
vector<RectPair> findVerticalAdjacencies(const Solution &sol) {
    vector<RectPair> res;
    int n = (int)sol.zonas.size();
    for (int i = 0; i < n; ++i) {
        const Rect &a = sol.zonas[i];
        for (int j = i + 1; j < n; ++j) {
            const Rect &b = sol.zonas[j];
            // mismo rango de filas
            if (a.r1 == b.r1 && a.r2 == b.r2) {
                // a a la izquierda de b
                if (a.c2 + 1 == b.c1) {
                    res.emplace_back(i, j);
                }
                // b a la izquierda de a
                if (b.c2 + 1 == a.c1) {
                    res.emplace_back(j, i);
                }
            }
        }
    }
    return res;
}

// Busca pares de rectángulos adyacentes horizontalmente:
// mismo rango de columnas y filas contiguas.
vector<RectPair> findHorizontalAdjacencies(const Solution &sol) {
    vector<RectPair> res;
    int n = (int)sol.zonas.size();
    for (int i = 0; i < n; ++i) {
        const Rect &a = sol.zonas[i];
        for (int j = i + 1; j < n; ++j) {
            const Rect &b = sol.zonas[j];
            // mismo rango de columnas
            if (a.c1 == b.c1 && a.c2 == b.c2) {
                // a arriba de b
                if (a.r2 + 1 == b.r1) {
                    res.emplace_back(i, j);
                }
                // b arriba de a
                if (b.r2 + 1 == a.r1) {
                    res.emplace_back(j, i);
                }
            }
        }
    }
    return res;
}


// Verifica que los rectángulos:
// - estén dentro de la grilla,
// - no se solapen,
// - cubran todas las celdas exactamente una vez.
// Devuelve true si la partición es válida y además construye
// una matriz de etiquetas (opcionalmente útil para depurar/imprimir).
bool checkPartition(const Instance &inst,
                    const Solution &sol,
                    vector<vector<int>> *labelsOut = nullptr) {
    int N = inst.N, M = inst.M;
    vector<vector<int>> labels(N, vector<int>(M, -1));

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

// Evalúa una solución: calcula SSE total y aplica penalización
// si alguna zona viola la homogeneidad.
double evaluateSolution(const Instance &inst, Solution &sol) {
    // Primero, verificamos partición rectangular válida
    vector<vector<int>> labels;
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

        // Restricción de homogeneidad: Var(z_k) <= alpha * Var(S)
        if (zoneVar > inst.alpha * inst.globalVar + 1e-9) {
            homogeneityOK = false;
        }
    }

    //if (!homogeneityOK) {
        // Penalizar fuerte
    //   totalSSE += PENALIZACION_GRANDE;
   // }

    sol.fitness = totalSSE;
    return sol.fitness;
}

// -------------------------
// Generación de solución inicial
// -------------------------


Solution crossover(const Instance &inst,
                   const Solution &p1,
                   const Solution &p2,
                   mt19937 &rng) {
    int L = (int)p1.genome.size();

    // Por seguridad (por si algo quedó raro)
    if (L == 0 || (int)p2.genome.size() != L) {
        // Fallback: comportamiento viejo
        return (p1.fitness <= p2.fitness) ? p1 : p2;
    }

    Solution child;
    child.genome.resize(L);

    // Crossover 1-point sobre el genoma
    std::uniform_int_distribution<int> cutDist(1, L - 1);
    int cut = cutDist(rng);

    for (int i = 0; i < cut; ++i) {
        child.genome[i] = p1.genome[i];
    }
    for (int i = cut; i < L; ++i) {
        child.genome[i] = p2.genome[i];
    }

    // Decodificar genoma → rectángulos y evaluar
    decodeGenomeToRects(inst, child.genome, child.zonas, rng);
    evaluateSolution(inst, child);

    return child;
}


// Genera una partición de la grilla en p "tiras" horizontales o verticales,
// todas rectangulares, sin solaparse y cubriendo todo el mapa.
// Es una forma sencilla pero 100% factible de solución inicial.
Solution generateRandomStripSolution(const Instance &inst, mt19937 &rng) {
    Solution sol;
    sol.genome.clear();
    sol.zonas.clear();

    int p = inst.p;
    int totalCells = inst.N * inst.M;

    if (p <= 0 || p > totalCells) {
        sol.fitness = INF;
        return sol;
    }

    int numGenes = p - 1;
    sol.genome.resize(numGenes);

    std::bernoulli_distribution orientation(0.5);
    std::uniform_real_distribution<double> alphaDist(0.2, 0.8);

    for (int g = 0; g < numGenes; ++g) {
        sol.genome[g].horizontal = orientation(rng);
        sol.genome[g].alpha      = alphaDist(rng);
    }

    // Decodificar genoma → rectángulos y evaluar
    decodeGenomeToRects(inst, sol.genome, sol.zonas, rng);
    evaluateSolution(inst, sol);

    return sol;
}






// Intenta aplicar una mutación local moviendo una frontera entre dos rectángulos.
// Devuelve true si logró mutar algo, false si no encontró pares adyacentes.
bool localBoundaryMutation(const Instance &inst, Solution &sol, mt19937 &rng) {
    (void)inst; // por ahora no lo usamos aquí

    // Buscamos pares adyacentes
    auto vertPairs = findVerticalAdjacencies(sol);
    auto horPairs  = findHorizontalAdjacencies(sol);

    if (vertPairs.empty() && horPairs.empty()) {
        return false; // no hay nada que podamos mover localmente
    }

    std::uniform_int_distribution<int> coin(0, 1);

    bool useVertical;
    if (!vertPairs.empty() && !horPairs.empty()) {
        useVertical = (coin(rng) == 0);
    } else if (!vertPairs.empty()) {
        useVertical = true;
    } else {
        useVertical = false;
    }

    if (useVertical) {
        // Elegir un par vertical al azar
        std::uniform_int_distribution<int> dist(0, (int)vertPairs.size() - 1);
        auto [idxL, idxR] = vertPairs[dist(rng)];
        Rect &left  = sol.zonas[idxL];
        Rect &right = sol.zonas[idxR];

        // left: [r1..r2] x [c1..c2]
        // right: [r1..r2] x [c2+1..c2'] (por construcción)
        int minBoundary = left.c1;          // left no puede quedar con menos de 1 col
        int maxBoundary = right.c2 - 1;     // right no puede quedar con menos de 1 col

        int oldBoundary = left.c2;

        if (minBoundary == maxBoundary) {
            // solo un posible corte → no hay movimiento real
            return false;
        }

        // Permitimos mover a lo más ±2 columnas alrededor del corte actual
        std::uniform_int_distribution<int> deltaDist(-2, 2);
        int newBoundary = oldBoundary + deltaDist(rng);
        newBoundary = max(newBoundary, minBoundary);
        newBoundary = min(newBoundary, maxBoundary);

        if (newBoundary == oldBoundary) {
            // no cambió nada
            return false;
        }

        left.c2  = newBoundary;
        right.c1 = newBoundary + 1;
        return true;
    } else {
        // Horizontal
        std::uniform_int_distribution<int> dist(0, (int)horPairs.size() - 1);
        auto [idxTop, idxBot] = horPairs[dist(rng)];
        Rect &top = sol.zonas[idxTop];
        Rect &bot = sol.zonas[idxBot];

        // top: [r1..r2] x [c1..c2]
        // bot: [r2+1..r2'] x [c1..c2]
        int minBoundary = top.r1;          // top al menos 1 fila
        int maxBoundary = bot.r2 - 1;      // bot al menos 1 fila

        int oldBoundary = top.r2;

        if (minBoundary == maxBoundary) {
            return false;
        }

        std::uniform_int_distribution<int> deltaDist(-2, 2);
        int newBoundary = oldBoundary + deltaDist(rng);
        newBoundary = max(newBoundary, minBoundary);
        newBoundary = min(newBoundary, maxBoundary);

        if (newBoundary == oldBoundary) {
            return false;
        }

        top.r2 = newBoundary;
        bot.r1 = newBoundary + 1;
        return true;
    }
}


Solution mutate(const Instance &inst,
                const Solution &original,
                mt19937 &rng) {
    Solution mutated = original;

    std::uniform_real_distribution<double> prob01(0.0, 1.0);
    std::normal_distribution<double> gauss(0.0, 0.1);           // perturbación suave de alpha
    std::bernoulli_distribution flipOrient(0.1);                // 10% cambiar orientación
    std::bernoulli_distribution newOrient(0.5);
    std::uniform_real_distribution<double> newAlpha(0.2, 0.8);  // para mutación fuerte

    bool strongMutation = (prob01(rng) < 0.1); // 10%: mutación global fuerte

    if (mutated.genome.empty()) {
        // Si algo salió mal, regeneramos
        return generateRandomStripSolution(inst, rng);
    }

    if (strongMutation) {
        // Reseteamos completamente el genoma (diversificación fuerte)
        int L = inst.p - 1;
        mutated.genome.resize(L);
        for (int i = 0; i < L; ++i) {
            mutated.genome[i].horizontal = newOrient(rng);
            mutated.genome[i].alpha      = newAlpha(rng);
        }
    } else {
        // Mutación local: intensificación
        int L = (int)mutated.genome.size();
        for (int i = 0; i < L; ++i) {
            // Con cierta probabilidad, perturbamos alpha
            if (prob01(rng) < 0.3) {
                mutated.genome[i].alpha += gauss(rng);
                // Acotar para evitar cortes extremos
                mutated.genome[i].alpha = std::max(0.05, std::min(0.95, mutated.genome[i].alpha));
            }
            // Con menor probabilidad, cambiamos la orientación
            if (flipOrient(rng)) {
                mutated.genome[i].horizontal = !mutated.genome[i].horizontal;
            }
        }
    }

    // Volvemos a decodificar y evaluar
    decodeGenomeToRects(inst, mutated.genome, mutated.zonas, rng);
    evaluateSolution(inst, mutated);
    return mutated;
}

Solution runGA(const Instance &inst,
               const GAParams &params,
               mt19937 &rng) {
    vector<Solution> population;
    population.reserve(params.populationSize);

    // ---- Inicialización de la población ----
    for (int i = 0; i < params.populationSize; ++i) {
        // generateRandomStripSolution ya:
        // - genera un genoma aleatorio
        // - decodifica a rectángulos
        // - llama a evaluateSolution
        Solution s = generateRandomStripSolution(inst, rng);
        population.push_back(std::move(s));
    }

    // Mejor de todos (para tracking global)
    Solution bestOverall = population[0];
    for (const auto &s : population) {
        if (s.fitness < bestOverall.fitness) {
            bestOverall = s;
        }
    }

    std::uniform_real_distribution<double> prob01(0.0, 1.0);

    // ---- Bucle principal GA ----
    for (int gen = 0; gen < params.generations; ++gen) {
        // Ordenar población por fitness (mejor primero)
        std::sort(population.begin(), population.end(),
                  [](const Solution &a, const Solution &b) {
                      return a.fitness < b.fitness;
                  });

        // Elitismo: mantenemos el mejor
        Solution elite = population[0];
        if (elite.fitness < bestOverall.fitness) {
            bestOverall = elite;
        }

        vector<Solution> newPop;
        newPop.reserve(params.populationSize);
        newPop.push_back(elite); // copiamos el elitista

        // Rellenamos el resto de la población
        while ((int)newPop.size() < params.populationSize) {
            const Solution &parent1 = tournamentSelect(population, params.tournamentK, rng);
            const Solution &parent2 = tournamentSelect(population, params.tournamentK, rng);

            Solution child;

            // Crossover
            if (prob01(rng) < params.crossoverRate) {
                // crossover:
                // - mezcla genome de p1 y p2
                // - decodifica a rectángulos
                // - llama a evaluateSolution
                child = crossover(inst, parent1, parent2, rng);
            } else {
                // Sin crossover: clon del padre1 (fitness ya válido)
                child = parent1;
            }

            // Mutación
            if (prob01(rng) < params.mutationRate) {
                // mutate:
                // - altera el genome
                // - decodifica a rectángulos
                // - llama a evaluateSolution
                child = mutate(inst, child, rng);
            }

            // OJO: ya no hace falta llamar a evaluateSolution aquí,
            // porque crossover/mutate/generación inicial ya lo hacen.

            newPop.push_back(std::move(child));
        }

        population = std::move(newPop);

        // Log opcional
        if (gen % 10 == 0) {
            cout << "[GA] Gen " << gen
                 << " mejor fitness = " << bestOverall.fitness << "\n";
        }
    }

    return bestOverall;
}



void printLabelMatrix(const Instance &inst, const Solution &sol) {
    vector<vector<int>> labels;
    if (!checkPartition(inst, sol, &labels)) {
        cout << "Solución no es una partición válida.\n";
        return;
    }

    cout << "Matriz Z (etiquetas de zonas):\n";
    for (int i = 0; i < inst.N; ++i) {
        for (int j = 0; j < inst.M; ++j) {
            cout << labels[i][j] + 1 << " "; // +1 para usar etiquetas 1..p
        }
        cout << "\n";
    }
}
void printZoneVariances(const Instance &inst, const Solution &sol) {
    cout << "Varianza intra-zona de la mejor solución:\n";

    for (size_t k = 0; k < sol.zonas.size(); ++k) {
        const Rect &r = sol.zonas[k];

        double zoneMean = 0.0;
        double zoneVar  = 0.0;
        double zoneSSE  = 0.0;

        computeZoneStats(inst, r, zoneMean, zoneVar, zoneSSE);

        cout << "  Zona " << (k + 1)
             << " -> media = " << zoneMean
             << ", varianza = " << zoneVar
             << ", SSE = " << zoneSSE
             << "\n";
    }
}




// Procesa todas las instancias .spp en una carpeta dada
void processInstancesInFolder(const string &folderPath,
                              int p,
                              double alpha,
                              mt19937 &rng) {
    if (!fs::exists(folderPath) || !fs::is_directory(folderPath)) {
        cerr << "La ruta '" << folderPath << "' no existe o no es un directorio.\n";
        return;
    }

    vector<string> files;
    for (const auto &entry : fs::directory_iterator(folderPath)) {
        if (entry.is_regular_file() && entry.path().extension() == ".spp") {
            files.push_back(entry.path().string());
        }
    }

    if (files.empty()) {
        cerr << "No se encontraron archivos .spp en la carpeta '" << folderPath << "'.\n";
        return;
    }

    sort(files.begin(), files.end());

    cout << "\nSe encontraron " << files.size()
         << " instancias .spp en '" << folderPath << "'.\n";

    // Parámetros del GA (puedes ajustarlos o incluso leerlos por consola)
    GAParams gaParams;
    gaParams.populationSize = 500;
    gaParams.generations    = 10000;
    gaParams.tournamentK    = 2;
    gaParams.crossoverRate  = 0.5;
    gaParams.mutationRate   = 0.8;

    for (const string &filePath : files) {
        cout << "\n========================================\n";
        cout << "Procesando instancia: " << filePath << "\n";

        Instance inst;
        if (!inst.loadFromFile(filePath)) {
            cerr << "  -> Error cargando la instancia. Se omite.\n";
            continue;
        }

        inst.p = p;
        inst.alpha = alpha;

        cout << "  N = " << inst.N << ", M = " << inst.M
             << ", p = " << inst.p << ", alpha = " << inst.alpha << "\n";
        cout << "  Varianza global = " << inst.globalVar << "\n";

        // Ejecutar el GA
        Solution best = runGA(inst, gaParams, rng);

        cout << "  Mejor fitness encontrado = " << best.fitness << "\n";
        cout << "  Matriz Z (etiquetas de zonas) del mejor individuo:\n";
        printLabelMatrix(inst, best);
        printZoneVariances(inst, best);
    }
}



int pedirP() {
    while (true) {
        cout << "Ingrese p (cantidad de zonas/sensores, entero positivo): ";
        string line;
        if (!std::getline(cin, line)) {
            cerr << "Error leyendo entrada. Saliendo.\n";
            exit(1);
        }

        // Quitar espacios al inicio/fin
        auto notSpace = [](int ch){ return !std::isspace(ch); };
        line.erase(line.begin(), std::find_if(line.begin(), line.end(), notSpace));
        line.erase(std::find_if(line.rbegin(), line.rend(), notSpace).base(), line.end());

        if (line.empty()) {
            cout << "p no puede ser vacío. Intente nuevamente.\n";
            continue;
        }

        std::stringstream ss(line);
        int p;
        char extra;
        if (!(ss >> p) || (ss >> extra)) {
            cout << "Entrada inválida. p debe ser un entero positivo.\n";
            continue;
        }

        if (p <= 0) {
            cout << "p debe ser un entero positivo (> 0).\n";
            continue;
        }

        return p;
    }
}

double pedirAlpha() {
    while (true) {
        cout << "Ingrese alpha (nivel de homogeneidad, 0 < alpha <= 1): ";
        string line;
        if (!std::getline(cin, line)) {
            cerr << "Error leyendo entrada. Saliendo.\n";
            exit(1);
        }

        // Quitar espacios al inicio/fin
        auto notSpace = [](int ch){ return !std::isspace(ch); };
        line.erase(line.begin(), std::find_if(line.begin(), line.end(), notSpace));
        line.erase(std::find_if(line.rbegin(), line.rend(), notSpace).base(), line.end());

        if (line.empty()) {
            cout << "alpha no puede ser vacío. Intente nuevamente.\n";
            continue;
        }

        std::stringstream ss(line);
        double alpha;
        char extra;
        if (!(ss >> alpha) || (ss >> extra)) {
            cout << "Entrada inválida. alpha debe ser un número real.\n";
            continue;
        }

        if (alpha <= 0.0 || alpha > 1.0) {
            cout << "alpha debe cumplir 0 < alpha <= 1.\n";
            continue;
        }

        return alpha;
    }
}




// -------------------------
// main: prueba mínima
// -------------------------
int main() {
    

    int p = pedirP();
    double alpha = pedirAlpha();

    string folderPath = "Instancias/Pequeñas";
    /*
    cout << "Ingrese ruta de la carpeta de instancias "
            "(por ejemplo, Instancias/Grandes): ";
    if (!std::getline(cin, folderPath)) {
        cerr << "Error leyendo la ruta de la carpeta.\n";
        return 1;
    }
    */
    if (folderPath.empty()) {
        cerr << "Error: la ruta de la carpeta no puede ser vacía.\n";
        return 1;
    }

    // RNG global
    std::random_device rd;
    std::mt19937 rng(rd());

    // Procesar todas las instancias .spp en esa carpeta
    processInstancesInFolder(folderPath, p, alpha, rng);

    return 0;
}
