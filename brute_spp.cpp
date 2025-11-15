#include <bits/stdc++.h>
using namespace std;

struct Rect {
    int r1, c1, r2, c2; // [r1..r2], [c1..c2]
};

int N, M, P;
double alpha;
vector<vector<double>> val;

// Mejor solución global SIN restricción de alpha
double bestSSE = 1e100;
vector<Rect> bestRects;

// Mejor solución que SÍ cumple alpha (si existe)
double bestSSEAlpha = 1e100;
vector<Rect> bestRectsAlpha;

// Stats globales
double globalVar = 0.0;

// Calcula SSE y varianza (var = SSE / n) de un rectángulo
pair<double,double> rectStats(const Rect &r) {
    double sum = 0.0;
    double sum2 = 0.0;
    int cnt = 0;
    for (int i = r.r1; i <= r.r2; ++i) {
        for (int j = r.c1; j <= r.c2; ++j) {
            double x = val[i][j];
            sum += x;
            sum2 += x * x;
            ++cnt;
        }
    }
    if (cnt == 0) return {0.0, 0.0};
    double mean = sum / cnt;
    double sse = 0.0;
    for (int i = r.r1; i <= r.r2; ++i) {
        for (int j = r.c1; j <= r.c2; ++j) {
            double x = val[i][j];
            double d = x - mean;
            sse += d * d;
        }
    }
    double var = sse / cnt; // var poblacional
    return {sse, var};
}

// Evalúa una partición completa: devuelve (esFactibleAlpha, SSEtotal)
pair<bool,double> evalPartition(const vector<Rect> &rects) {
    double totalSSE = 0.0;
    bool okAlpha = true;
    double limitVar = alpha * globalVar; // alpha * Var(S)

    for (const Rect &r : rects) {
        auto [sse, var] = rectStats(r);
        totalSSE += sse;
        if (var > limitVar) {
            okAlpha = false;
        }
    }
    return {okAlpha, totalSSE};
}

// Backtracking: vamos partiendo rectángulos hasta tener P
void backtrack(vector<Rect> &rects) {
    if ((int)rects.size() == P) {
        auto [okAlpha, totalSSE] = evalPartition(rects);

        // Siempre actualizamos la mejor solución SIN restricción de alpha
        if (totalSSE < bestSSE) {
            bestSSE = totalSSE;
            bestRects = rects;
            cout << "[MEJOR GLOBAL] SSE = " << fixed << setprecision(6)
                 << bestSSE
                 << (okAlpha ? " (cumple alpha)\n" : " (NO cumple alpha)\n");
        }

        // Y si además cumple alpha, actualizamos el mejor con alpha
        if (okAlpha && totalSSE < bestSSEAlpha) {
            bestSSEAlpha = totalSSE;
            bestRectsAlpha = rects;
            cout << "  [MEJOR CON ALPHA] SSE = " << fixed << setprecision(6)
                 << bestSSEAlpha << " (cumple alpha)\n";
        }

        return;
    }

    int currentSize = rects.size();
    for (int idx = 0; idx < currentSize; ++idx) {
        Rect r = rects[idx];
        int h = r.r2 - r.r1 + 1;
        int w = r.c2 - r.c1 + 1;

        // Cortes horizontales posibles
        if (h >= 2) {
            for (int cutRow = r.r1; cutRow < r.r2; ++cutRow) {
                Rect rA{r.r1, r.c1, cutRow,   r.c2};
                Rect rB{cutRow + 1, r.c1, r.r2, r.c2};

                rects[idx] = rA;
                rects.push_back(rB);
                backtrack(rects);
                rects.pop_back();
                rects[idx] = r; // restaurar
            }
        }

        // Cortes verticales posibles
        if (w >= 2) {
            for (int cutCol = r.c1; cutCol < r.c2; ++cutCol) {
                Rect rA{r.r1, r.c1, r.r2, cutCol};
                Rect rB{r.r1, cutCol + 1, r.r2, r.c2};

                rects[idx] = rA;
                rects.push_back(rB);
                backtrack(rects);
                rects.pop_back();
                rects[idx] = r; // restaurar
            }
        }
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // Leer N,M
    if (!(cin >> N >> M)) {
        cerr << "Error leyendo N y M\n";
        return 1;
    }
    val.assign(N, vector<double>(M));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            cin >> val[i][j];
        }
    }

    // Leer p y alpha
    cin >> P >> alpha; // ejemplo: 5 0.2

    // Calcular varianza global
    {
        double sum = 0.0, sum2 = 0.0;
        int cnt = N * M;
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < M; ++j) {
                double x = val[i][j];
                sum += x;
                sum2 += x * x;
            }
        double mean = sum / cnt;
        double sse = 0.0;
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < M; ++j) {
                double x = val[i][j];
                double d = x - mean;
                sse += d * d;
            }
        globalVar = sse / cnt;
        cout << fixed << setprecision(6);
        cout << "Varianza global = " << globalVar
             << ", alpha*Var = " << alpha * globalVar << "\n";
    }

    // Partición inicial: todo el mapa
    vector<Rect> rects;
    rects.push_back(Rect{0, 0, N - 1, M - 1});

    backtrack(rects);

    cout << "\n===== RESULTADOS FINALES =====\n";

    // Siempre mostramos la mejor solución SIN alpha
    if (bestSSE < 1e99) {
        cout << "Mejor SSE (sin restringir alpha) = " << bestSSE << "\n";
        cout << "Rectángulos (r1,c1,r2,c2):\n";
        for (size_t i = 0; i < bestRects.size(); ++i) {
            auto &r = bestRects[i];
            cout << "Zona " << (i+1) << ": "
                 << "(" << r.r1 << "," << r.c1 << ") - "
                 << "(" << r.r2 << "," << r.c2 << ")\n";
        }

        // Matriz de etiquetas Z para esa mejor SSE
        vector<vector<int>> Z(N, vector<int>(M, 0));
        for (int k = 0; k < (int)bestRects.size(); ++k) {
            auto &r = bestRects[k];
            for (int i = r.r1; i <= r.r2; ++i)
                for (int j = r.c1; j <= r.c2; ++j)
                    Z[i][j] = k + 1;
        }
        cout << "Matriz Z (mejor SSE):\n";
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                cout << Z[i][j] << " ";
            }
            cout << "\n";
        }
    } else {
        cout << "No se encontró ninguna partición (esto no debería pasar para tamaños pequeños)\n";
    }

    // Y adicionalmente, si existe alguna que cumpla alpha, la mostramos
    if (bestSSEAlpha < 1e99) {
        cout << "\nMejor SSE que SÍ cumple alpha = " << bestSSEAlpha << "\n";
        cout << "Rectángulos (r1,c1,r2,c2):\n";
        for (size_t i = 0; i < bestRectsAlpha.size(); ++i) {
            auto &r = bestRectsAlpha[i];
            cout << "Zona " << (i+1) << ": "
                 << "(" << r.r1 << "," << r.c1 << ") - "
                 << "(" << r.r2 << "," << r.c2 << ")\n";
        }

        vector<vector<int>> Z2(N, vector<int>(M, 0));
        for (int k = 0; k < (int)bestRectsAlpha.size(); ++k) {
            auto &r = bestRectsAlpha[k];
            for (int i = r.r1; i <= r.r2; ++i)
                for (int j = r.c1; j <= r.c2; ++j)
                    Z2[i][j] = k + 1;
        }
        cout << "Matriz Z (mejor que cumple alpha):\n";
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                cout << Z2[i][j] << " ";
            }
            cout << "\n";
        }
    } else {
        cout << "\nNo se encontró ninguna partición que cumpla la restricción de alpha.\n";
    }

    return 0;
}
