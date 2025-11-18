# Estructura del Proyecto SPP - Algoritmo Genético

```
/home/carlos/IA/
│
├── Header Files (.h)
│   ├── instance.h           # Definición de estructura Instance
│   ├── solution.h           # Estructuras Gene, Rect, Solution
│   ├── decoder.h            # Decodificación de genomas
│   ├── evaluator.h          # Evaluación y validación
│   ├── genetic_algorithm.h  # Operadores del GA
│   └── utils.h              # Utilidades de I/O y configuración
│
├── Implementation Files (.cpp)
│   ├── instance.cpp         # Carga de instancias y estadísticas
│   ├── solution.cpp         # (Vacío - estructura simple)
│   ├── decoder.cpp          # Lógica de decodificación de genomas
│   ├── evaluator.cpp        # Cálculo de fitness y validación
│   ├── genetic_algorithm.cpp# Selección, cruce, mutación, GA principal
│   └── utils.cpp            # Carga config, entrada usuario, procesamiento
│
├── Main Program
│   └── spp_main.cpp         # Punto de entrada del programa
│
├── Configuration
│   └── ga_config.txt        # Parámetros del algoritmo genético
│
├── Build System
│   └── Makefile             # Compilación automatizada
│
├── Documentation
│   └── README.md            # Documentación del proyecto
│
├── Version Control
│   └── .gitignore           # Archivos a ignorar en git
│
└── Legacy Files (mantener como backup)
    ├── main.cpp
    ├── main2.cpp
    └── main3.cpp            # Versión original (modificada)
```

## Diagrama de Dependencias

```
spp_main.cpp
    ├── utils.h
    │   ├── genetic_algorithm.h
    │   │   ├── instance.h
    │   │   ├── solution.h
    │   │   ├── decoder.h
    │   │   │   ├── instance.h
    │   │   │   └── solution.h
    │   │   └── evaluator.h
    │   │       ├── instance.h
    │   │       └── solution.h
    │   ├── instance.h
    │   └── evaluator.h
    └── genetic_algorithm.h
```

## Flujo de Ejecución

1. **Inicio** (spp_main.cpp)
   - Cargar configuración del GA desde ga_config.txt
   - Solicitar parámetros p y alpha al usuario

2. **Procesamiento** (utils.cpp)
   - Iterar sobre archivos .spp en la carpeta especificada
   - Para cada instancia:

3. **Carga de Instancia** (instance.cpp)
   - Leer dimensiones N x M
   - Leer matriz de valores S
   - Calcular estadísticas globales

4. **Algoritmo Genético** (genetic_algorithm.cpp)
   - Generar población inicial aleatoria
   - Para cada generación:
     - Evaluar fitness de cada individuo
     - Seleccionar padres por torneo
     - Aplicar crossover y mutación
     - Mantener elitismo

5. **Decodificación** (decoder.cpp)
   - Convertir genoma (secuencia de genes) en rectángulos

6. **Evaluación** (evaluator.cpp)
   - Validar partición rectangular
   - Calcular SSE total
   - Verificar restricciones de homogeneidad
   - Aplicar penalización si es necesario

7. **Salida** (evaluator.cpp)
   - Imprimir matriz de etiquetas
   - Mostrar varianzas por zona
   - Reportar fitness final
```
