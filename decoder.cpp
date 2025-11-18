#include "decoder.h"
#include <algorithm>
#include <numeric>
#include <cmath>

void decodeGenomeToRects(const Instance &inst,
                         const std::vector<Gene> &genome,
                         std::vector<Rect> &rects,
                         std::mt19937 &rng)
{
    rects.clear();
    // Rectángulo completo
    rects.push_back(Rect{0, 0, inst.N - 1, inst.M - 1});

    // Aplicamos cada gen
    for (const Gene &g : genome) {
        if ((int)rects.size() >= inst.p) break;

        if (rects.empty()) break;

        // ---- Seleccionar rectángulo usando g.beta ----
        std::vector<int> idx(rects.size());
        std::iota(idx.begin(), idx.end(), 0);

        // Ordenar por área descendente
        std::sort(idx.begin(), idx.end(), [&](int a, int b) {
            int ha = rects[a].r2 - rects[a].r1 + 1;
            int wa = rects[a].c2 - rects[a].c1 + 1;
            int hb = rects[b].r2 - rects[b].r1 + 1;
            int wb = rects[b].c2 - rects[b].c1 + 1;
            return ha * wa > hb * wb;
        });

        double betaClamped = std::max(0.0, std::min(1.0, g.beta));
        int pos = (int)std::floor(betaClamped * idx.size());
        if (pos >= (int)idx.size()) pos = (int)idx.size() - 1;
        if (pos < 0) pos = 0;

        int chosen = -1;
        // Intentar desde pos, luego pos+1, etc. hasta encontrar algo que se pueda cortar
        for (int k = 0; k < (int)idx.size(); ++k) {
            int cand = idx[(pos + k) % idx.size()];
            Rect rTest = rects[cand];
            int h = rTest.r2 - rTest.r1 + 1;
            int w = rTest.c2 - rTest.c1 + 1;
            if (h * w >= 2) {
                chosen = cand;
                break;
            }
        }

        if (chosen == -1) {
            // Nada para cortar
            break;
        }

        Rect r = rects[chosen];
        int h = r.r2 - r.r1 + 1;
        int w = r.c2 - r.c1 + 1;

        bool horizontal = g.horizontal;
        // Ajustar orientación si no es viable
        if (horizontal && h < 2 && w >= 2) horizontal = false;
        if (!horizontal && w < 2 && h >= 2) horizontal = true;

        Rect rA, rB;

        if (horizontal && h >= 2) {
            int minRow = r.r1;
            int maxRow = r.r2 - 1; // al menos 1 fila en cada lado

            double gammaClamped = std::max(0.0, std::min(1.0, g.gamma));
            int offset = (int)std::round(gammaClamped * (h - 1));
            int cutRow = minRow + offset;
            if (cutRow < minRow) cutRow = minRow;
            if (cutRow > maxRow) cutRow = maxRow;

            rA = Rect{r.r1, r.c1, cutRow, r.c2};
            rB = Rect{cutRow + 1, r.c1, r.r2, r.c2};
        } else if (!horizontal && w >= 2) {
            int minCol = r.c1;
            int maxCol = r.c2 - 1; // al menos 1 col en cada lado

            double gammaClamped = std::max(0.0, std::min(1.0, g.gamma));
            int offset = (int)std::round(gammaClamped * (w - 1));
            int cutCol = minCol + offset;
            if (cutCol < minCol) cutCol = minCol;
            if (cutCol > maxCol) cutCol = maxCol;

            rA = Rect{r.r1, r.c1, r.r2, cutCol};
            rB = Rect{r.r1, cutCol + 1, r.r2, r.c2};
        } else {
            // No se pudo cortar
            continue;
        }

        rects[chosen] = rA;
        rects.push_back(rB);
    }

    // Si por algún motivo no alcanzamos p rectángulos, rellenamos partiendo al azar
    std::uniform_int_distribution<int> coin(0, 1);
    while ((int)rects.size() < inst.p) {
        int idxRect = -1;
        for (int i = 0; i < (int)rects.size(); ++i) {
            int h = rects[i].r2 - rects[i].r1 + 1;
            int w = rects[i].c2 - rects[i].c1 + 1;
            if (h * w >= 2) {
                idxRect = i;
                break;
            }
        }
        if (idxRect == -1) break;

        Rect r = rects[idxRect];
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
            int cutRow = distRow(rng);
            rA = Rect{r.r1, r.c1, cutRow, r.c2};
            rB = Rect{cutRow + 1, r.c1, r.r2, r.c2};
        } else {
            std::uniform_int_distribution<int> distCol(r.c1, r.c2 - 1);
            int cutCol = distCol(rng);
            rA = Rect{r.r1, r.c1, r.r2, cutCol};
            rB = Rect{r.r1, cutCol + 1, r.r2, r.c2};
        }

        rects[idxRect] = rA;
        rects.push_back(rB);
    }
}

std::vector<RectPair> findVerticalAdjacencies(const Solution &sol) {
    std::vector<RectPair> res;
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

std::vector<RectPair> findHorizontalAdjacencies(const Solution &sol) {
    std::vector<RectPair> res;
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
