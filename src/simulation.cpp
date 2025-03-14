#include <string>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <thread>
#include <chrono>
#include <mpi.h>
#include <vector>
#include <algorithm>

#include "model.hpp"
#include "display.hpp"

using namespace std::string_literals;
using namespace std::chrono_literals;

struct ParamsType {
    double length{1.};
    unsigned discretization{20u};
    std::array<double,2> wind{0.,0.};
    Model::LexicoIndices start{10u,10u};
};

// Funciones de análisis de argumentos (se mantienen iguales)
ParamsType parse_arguments(int nargs, char* args[]);
bool check_params(ParamsType& params);
void display_params(ParamsType const& params);

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Parse de parámetros (todos ven el mismo dominio global)
    auto params = parse_arguments(argc-1, &argv[1]);
    display_params(params);
    if (!check_params(params)) {
        MPI_Finalize();
        return EXIT_FAILURE;
    }
    
    int N = params.discretization; // tamaño global (filas y columnas)
    
    // Suponemos que el proceso 0 solo visualiza y los trabajadores (ranks 1...size-1) realizan los cálculos.
    int workerCount = size - 1;
    if(workerCount < 1) {
        std::cerr << "Se requieren al menos 2 procesos (1 visualizador + 1 trabajador)." << std::endl;
        MPI_Finalize();
        return EXIT_FAILURE;
    }
    
    // Dividir el dominio entre los procesos trabajadores.
    int rows_per_worker = N / workerCount;
    int extra_rows = N % workerCount;
    int local_start = 0, local_end = 0;
    if(rank != 0) {
        int worker_rank = rank - 1; // trabajadores numerados de 0 a workerCount-1
        local_start = worker_rank * rows_per_worker + std::min(worker_rank, extra_rows);
        local_end = local_start + rows_per_worker + (worker_rank < extra_rows ? 1 : 0);
    }
    int local_rows = (rank == 0) ? 0 : (local_end - local_start);
    
    // Creamos una instancia global del modelo (se asume que la misma configuración se utiliza en todos)
    Model simu(params.length, N, params.wind, params.start);
    
    // Proceso 0: reserva arreglos globales para visualización.
    std::vector<std::uint8_t> global_fire, global_vegetal;
    if (rank == 0) {
        global_fire.resize(N * N, 0u);
        global_vegetal.resize(N * N, 255u);
    }
    
    // Procesos trabajadores: cada uno reserva arreglos locales extendidos (con celdas fantasma)
    // Tamaño: (local_rows + 2) x N
    std::vector<std::uint8_t> local_fire, local_vegetal;
    if (rank != 0) {
        int extended_rows = local_rows + 2; // fila 0: halo superior, filas 1..local_rows: datos internos, fila local_rows+1: halo inferior
        local_fire.resize(extended_rows * N, 0u);
        local_vegetal.resize(extended_rows * N, 255u);
        // Inicializar la parte interna con el estado inicial extraído de simu
        for (int r = 0; r < local_rows; r++) {
            std::copy(simu.fire_map().begin() + (local_start + r) * N,
                      simu.fire_map().begin() + (local_start + r + 1) * N,
                      local_fire.begin() + (r + 1) * N);
            std::copy(simu.vegetal_map().begin() + (local_start + r) * N,
                      simu.vegetal_map().begin() + (local_start + r + 1) * N,
                      local_vegetal.begin() + (r + 1) * N);
        }
    }
    
    // Bucle principal de simulación
    bool active = true;
    while (active) {
        if (rank != 0) {
            // 1. Intercambio de halos
            // Halo superior: si existe (rank > 1)
            if (rank > 1) {
                MPI_Sendrecv(local_fire.data() + 1 * N, N, MPI_UNSIGNED_CHAR,
                             rank - 1, 0,
                             local_fire.data(), N, MPI_UNSIGNED_CHAR,
                             rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Sendrecv(local_vegetal.data() + 1 * N, N, MPI_UNSIGNED_CHAR,
                             rank - 1, 0,
                             local_vegetal.data(), N, MPI_UNSIGNED_CHAR,
                             rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            // Halo inferior: si existe (rank < size-1)
            if (rank < size - 1) {
                MPI_Sendrecv(local_fire.data() + local_rows * N, N, MPI_UNSIGNED_CHAR,
                             rank + 1, 1,
                             local_fire.data() + (local_rows + 1) * N, N, MPI_UNSIGNED_CHAR,
                             rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Sendrecv(local_vegetal.data() + local_rows * N, N, MPI_UNSIGNED_CHAR,
                             rank + 1, 1,
                             local_vegetal.data() + (local_rows + 1) * N, N, MPI_UNSIGNED_CHAR,
                             rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            
            // 2. Actualización local: se actualiza la parte interna usando la función update_local
            // Esta función debe estar implementada en Model (como se mostró en el ejemplo anterior)
            active = simu.update_local(local_fire, local_vegetal, local_rows, N);
            
            // 3. Enviar al proceso 0 las filas internas (excluyendo halos)
            int count = local_rows * N;
            MPI_Send(local_fire.data() + N, count, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD);
            MPI_Send(local_vegetal.data() + N, count, MPI_UNSIGNED_CHAR, 0, 1, MPI_COMM_WORLD);
        } else {
            // En el proceso 0: solo se encarga de la visualización
            // (Podrías también actualizar un estado global, si lo deseas)
            // Recibir las porciones de cada trabajador
            for (int i = 1; i < size; i++) {
                int worker_rank = i - 1;
                int start_row_i = worker_rank * rows_per_worker + std::min(worker_rank, extra_rows);
                int end_row_i = start_row_i + rows_per_worker + (worker_rank < extra_rows ? 1 : 0);
                int count = (end_row_i - start_row_i) * N;
                MPI_Recv(global_fire.data() + start_row_i * N, count,
                         MPI_UNSIGNED_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(global_vegetal.data() + start_row_i * N, count,
                         MPI_UNSIGNED_CHAR, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            // Actualizar la visualización
            auto displayer = Displayer::instance();
            displayer->update(global_vegetal, global_fire);
            
            // Opcional: avanzar el estado global (por ejemplo, llamar a simu.update())
            simu.update();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    } // Fin del bucle de simulación

    MPI_Finalize();
    return EXIT_SUCCESS;
}
