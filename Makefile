# Makefile para el proyecto SPP con Algoritmo Genético

CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O3 -march=native
LDFLAGS = -lstdc++fs

# Archivos fuente
SOURCES = spp_main.cpp instance.cpp evaluator.cpp decoder.cpp genetic_algorithm.cpp utils.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = spp_solver

# Target principal
all: $(EXECUTABLE)

# Compilar el ejecutable
$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Compilar archivos objeto
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Limpiar archivos generados
clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

# Ejecutar el programa
run: $(EXECUTABLE)
	./$(EXECUTABLE)

# Regla para recompilar todo
rebuild: clean all

.PHONY: all clean run rebuild
