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

struct ParamsType
{
    double length{1.};
    unsigned discretization{20u};
    std::array<double,2> wind{0.,0.};
    Model::LexicoIndices start{10u,10u};
};

void analyze_arg( int nargs, char* args[], ParamsType& params )
{
    if (nargs ==0) return;
    std::string key(args[0]);
    if (key == "-l"s)
    {
        if (nargs < 2)
        {
            std::cerr << "Manque une valeur pour la longueur du terrain !" << std::endl;
            exit(EXIT_FAILURE);
        }
        params.length = std::stoul(args[1]);
        analyze_arg(nargs-2, &args[2], params);
        return;
    }
    auto pos = key.find("--longueur=");
    if (pos < key.size())
    {
        auto subkey = std::string(key,pos+11);
        params.length = std::stoul(subkey);
        analyze_arg(nargs-1, &args[1], params);
        return;
    }

    if (key == "-n"s)
    {
        if (nargs < 2)
        {
            std::cerr << "Manque une valeur pour le nombre de cases par direction pour la discrétisation du terrain !" << std::endl;
            exit(EXIT_FAILURE);
        }
        params.discretization = std::stoul(args[1]);
        analyze_arg(nargs-2, &args[2], params);
        return;
    }
    pos = key.find("--number_of_cases=");
    if (pos < key.size())
    {
        auto subkey = std::string(key, pos+18);
        params.discretization = std::stoul(subkey);
        analyze_arg(nargs-1, &args[1], params);
        return;
    }

    if (key == "-w"s)
    {
        if (nargs < 2)
        {
            std::cerr << "Manque une paire de valeurs pour la direction du vent !" << std::endl;
            exit(EXIT_FAILURE);
        }
        std::string values =std::string(args[1]);
        params.wind[0] = std::stod(values);
        auto pos = values.find(",");
        if (pos == values.size())
        {
            std::cerr << "Doit fournir deux valeurs séparées par une virgule pour définir la vitesse" << std::endl;
            exit(EXIT_FAILURE);
        }
        auto second_value = std::string(values, pos+1);
        params.wind[1] = std::stod(second_value);
        analyze_arg(nargs-2, &args[2], params);
        return;
    }
    pos = key.find("--wind=");
    if (pos < key.size())
    {
        auto subkey = std::string(key, pos+7);
        params.wind[0] = std::stoul(subkey);
        auto pos = subkey.find(",");
        if (pos == subkey.size())
        {
            std::cerr << "Doit fournir deux valeurs séparées par une virgule pour définir la vitesse" << std::endl;
            exit(EXIT_FAILURE);
        }
        auto second_value = std::string(subkey, pos+1);
        params.wind[1] = std::stod(second_value);
        analyze_arg(nargs-1, &args[1], params);
        return;
    }

    if (key == "-s"s)
    {
        if (nargs < 2)
        {
            std::cerr << "Manque une paire de valeurs pour la position du foyer initial !" << std::endl;
            exit(EXIT_FAILURE);
        }
        std::string values =std::string(args[1]);
        params.start.column = std::stod(values);
        auto pos = values.find(",");
        if (pos == values.size())
        {
            std::cerr << "Doit fournir deux valeurs séparées par une virgule pour définir la position du foyer initial" << std::endl;
            exit(EXIT_FAILURE);
        }
        auto second_value = std::string(values, pos+1);
        params.start.row = std::stod(second_value);
        analyze_arg(nargs-2, &args[2], params);
        return;
    }
    pos = key.find("--start=");
    if (pos < key.size())
    {
        auto subkey = std::string(key, pos+8);
        params.start.column = std::stoul(subkey);
        auto pos = subkey.find(",");
        if (pos == subkey.size())
        {
            std::cerr << "Doit fournir deux valeurs séparées par une virgule pour définir la vitesse" << std::endl;
            exit(EXIT_FAILURE);
        }
        auto second_value = std::string(subkey, pos+1);
        params.start.row = std::stod(second_value);
        analyze_arg(nargs-1, &args[1], params);
        return;
    }
}

ParamsType parse_arguments( int nargs, char* args[] )
{
    if (nargs == 0) return {};
    if ( (std::string(args[0]) == "--help"s) || (std::string(args[0]) == "-h") )
    {
        std::cout << 
R"RAW(Usage : simulation [option(s)]
  Lance la simulation d'incendie en prenant en compte les [option(s)].
  Les options sont :
    -l, --longueur=LONGUEUR     Définit la taille LONGUEUR (réel en km) du carré représentant la carte de la végétation.
    -n, --number_of_cases=N     Nombre n de cases par direction pour la discrétisation
    -w, --wind=VX,VY            Définit le vecteur vitesse du vent (pas de vent par défaut).
    -s, --start=COL,ROW         Définit les indices I,J de la case où commence l'incendie (milieu de la carte par défaut)
)RAW";
        exit(EXIT_SUCCESS);
    }
    ParamsType params;
    analyze_arg(nargs, args, params);
    return params;
}

bool check_params(ParamsType& params)
{
    bool flag = true;
    if (params.length <= 0)
    {
        std::cerr << "[ERREUR FATALE] La longueur du terrain doit être positive et non nulle !" << std::endl;
        flag = false;
    }

    if (params.discretization <= 0)
    {
        std::cerr << "[ERREUR FATALE] Le nombre de cellules par direction doit être positive et non nulle !" << std::endl;
        flag = false;
    }

    if ( (params.start.row >= params.discretization) || (params.start.column >= params.discretization) )
    {
        std::cerr << "[ERREUR FATALE] Mauvais indices pour la position initiale du foyer" << std::endl;
        flag = false;
    }
    
    return flag;
}

void display_params(ParamsType const& params)
{
    std::cout << "Parametres définis pour la simulation : \n"
              << "\tTaille du terrain : " << params.length << std::endl 
              << "\tNombre de cellules par direction : " << params.discretization << std::endl 
              << "\tVecteur vitesse : [" << params.wind[0] << ", " << params.wind[1] << "]" << std::endl
              << "\tPosition initiale du foyer (col, ligne) : " << params.start.column << ", " << params.start.row << std::endl;
}


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
    
    // Suponemos que el proceso 0 solo visualiza y los trabajadores (ranks 1...size-1) realizan los cálculos.
    int workerCount = size - 1;
    if(workerCount < 1) {
        std::cerr << "Se requieren al menos 2 procesos (1 visualizador + 1 trabajador)." << std::endl;
        MPI_Finalize();
        return EXIT_FAILURE;
    }
    
    // Dividir el dominio entre los procesos trabajadores.
    int rows_per_worker = params.discretization / workerCount;
    int extra_rows = params.discretization% workerCount;
    int local_start = 0, local_end = 0;
    if(rank != 0) {
        int worker_rank = rank - 1; // trabajadores numerados de 0 a workerCount-1
        local_start = worker_rank * rows_per_worker + std::min(worker_rank, extra_rows);
        local_end = local_start + rows_per_worker + (worker_rank < extra_rows ? 1 : 0);
    }
    int local_rows = (rank == 0) ? 0 : (local_end - local_start);
    
    // Creamos una instancia global del modelo (se asume que la misma configuración se utiliza en todos)
    Model simu(params.length, params.discretization, params.wind, params.start);
    
    // Proceso 0: reserva arreglos globales para visualización.
    std::vector<std::uint8_t> global_fire, global_vegetal;
    if (rank == 0) {
        global_fire.resize(params.discretization* params.discretization, 0u);
        global_vegetal.resize(params.discretization * params.discretization, 255u);
    }
    
    // Procesos trabajadores: cada uno reserva arreglos locales extendidos (con celdas fantasma)
    // Tamaño: (local_rows + 2) x N
    std::vector<std::uint8_t> local_fire, local_vegetal;
    if (rank != 0) {
        int extended_rows = local_rows + 2; // fila 0: halo superior, filas 1..local_rows: datos internos, fila local_rows+1: halo inferior
        local_fire.resize(extended_rows * params.discretization, 0u);
        local_vegetal.resize(extended_rows * params.discretization, 255u);
        // Inicializar la parte interna con el estado inicial extraído de simu
        for (int r = 0; r < local_rows; r++) {
            std::copy(simu.fire_map().begin() + (local_start + r) * params.discretization,
                      simu.fire_map().begin() + (local_start + r + 1) * params.discretization,
                      local_fire.begin() + (r + 1) * params.discretization);
            std::copy(simu.vegetal_map().begin() + (local_start + r) * params.discretization,
                      simu.vegetal_map().begin() + (local_start + r + 1) * params.discretization,
                      local_vegetal.begin() + (r + 1) * params.discretization);
        }
    }
    
    // Bucle principal de simulación
    bool active = true;
    while (active) {
        if (rank != 0) {
            // 1. Intercambio de halos
            // Halo superior: si existe (rank > 1)
            if (rank > 1) {
                MPI_Sendrecv(local_fire.data() + 1 * params.discretization, params.discretization, MPI_UNSIGNED_CHAR,
                             rank - 1, 0,
                             local_fire.data(), params.discretization, MPI_UNSIGNED_CHAR,
                             rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Sendrecv(local_vegetal.data() + 1 * params.discretization, params.discretization, MPI_UNSIGNED_CHAR,
                             rank - 1, 0,
                             local_vegetal.data(), params.discretization, MPI_UNSIGNED_CHAR,
                             rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            // Halo inferior: si existe (rank < size-1)
            if (rank < size - 1) {
                MPI_Sendrecv(local_fire.data() + local_rows * params.discretization, params.discretization, MPI_UNSIGNED_CHAR,
                             rank + 1, 1,
                             local_fire.data() + (local_rows + 1) * params.discretization, params.discretization, MPI_UNSIGNED_CHAR,
                             rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Sendrecv(local_vegetal.data() + local_rows * params.discretization, params.discretization, MPI_UNSIGNED_CHAR,
                             rank + 1, 1,
                             local_vegetal.data() + (local_rows + 1) * params.discretization, params.discretization, MPI_UNSIGNED_CHAR,
                             rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            
            // 2. Actualización local: se actualiza la parte interna usando la función update_local
            // Esta función debe estar implementada en Model (como se mostró en el ejemplo anterior)
            active = simu.update_local(local_fire, local_vegetal, local_rows, params.discretization);
            
            // 3. Enviar al proceso 0 las filas internas (excluyendo halos)
            int count = local_rows * params.discretization;
            MPI_Send(local_fire.data() + params.discretization, count, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD);
            MPI_Send(local_vegetal.data() + params.discretization, count, MPI_UNSIGNED_CHAR, 0, 1, MPI_COMM_WORLD);
        } else {
            // En el proceso 0: solo se encarga de la visualización
            // (Podrías también actualizar un estado global, si lo deseas)
            // Recibir las porciones de cada trabajador
            for (int i = 1; i < size; i++) {
                int worker_rank = i - 1;
                int start_row_i = worker_rank * rows_per_worker + std::min(worker_rank, extra_rows);
                int end_row_i = start_row_i + rows_per_worker + (worker_rank < extra_rows ? 1 : 0);
                int count = (end_row_i - start_row_i) * params.discretization;
                MPI_Recv(global_fire.data() + start_row_i * params.discretization, count,
                         MPI_UNSIGNED_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(global_vegetal.data() + start_row_i * params.discretization, count,
                         MPI_UNSIGNED_CHAR, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            // Actualizar la visualización
            auto displayer = Displayer::instance();
            displayer->update(global_vegetal, global_fire);
            
            // Opcional: avanzar el estado global (por ejemplo, llamar a simu.update())
            // simu.update();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    } // Fin del bucle de simulación

    MPI_Finalize();
    return EXIT_SUCCESS;
}
