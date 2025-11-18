# Algoritmo Genético para el Problema SPP

Este proyecto implementa un algoritmo genético para resolver el Sensor Placement Problem (SPP) con cortes guillotinables.

## Estructura del Proyecto

### Archivos de Código

- **instance.h/cpp**: Manejo de datos de instancias del problema
  - Carga de instancias desde archivos .spp
  - Cálculo de estadísticas globales (media y varianza)

- **solution.h/cpp**: Estructuras de datos para soluciones
  - `Gene`: Representa un corte guillotinable (orientación, posición gamma, selector beta)
  - `Rect`: Representa un rectángulo en la grilla
  - `Solution`: Genoma + rectángulos decodificados + fitness

- **decoder.h/cpp**: Decodificación de genomas a particiones
  - Convierte secuencia de genes en rectángulos concretos
  - Detección de adyacencias entre rectángulos

- **evaluator.h/cpp**: Evaluación y validación de soluciones
  - Cálculo de SSE (Sum of Squared Errors)
  - Verificación de restricciones de homogeneidad
  - Validación de particiones
  - Impresión de resultados

- **genetic_algorithm.h/cpp**: Operadores y lógica del algoritmo genético
  - Selección por k torneo
  - Operador de cruce (crossover)
  - Operador de mutación (local e intensiva)
  - Bucle principal del GA con elitismo

- **utils.h/cpp**: Utilidades y entrada/salida
  - Carga de configuración desde archivo
  - Solicitud de parámetros por consola
  - Procesamiento por lotes de instancias

- **spp_main.cpp**: Punto de entrada del programa

### Archivos de Configuración

- **ga_config.txt**: Parámetros del algoritmo genético
  - populationSize: Tamaño de la población
  - generations: Número de generaciones
  - tournamentK: Tamaño del torneo para selección
  - crossoverRate: Probabilidad de aplicar crossover
  - mutationRate: Probabilidad de mutar un hijo

### Otros Archivos

- **Makefile**: Compilación automatizada del proyecto

## Compilación

```bash
make
```

Para recompilar desde cero:
```bash
make rebuild
```

Para limpiar archivos generados:
```bash
make clean
```

## Ejecución

```bash
./spp_solver
```

El programa solicitará:
1. **p**: Cantidad de zonas/sensores (entero positivo)
2. **alpha**: Nivel de homogeneidad (0 < alpha ≤ 1)

Los parámetros del GA se cargan automáticamente desde `ga_config.txt`.

