#include <string>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <thread>
#include <chrono>

#include <mpi.h>

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

int main( int nargs, char* args[] )
{

    MPI_Init(&nargs, &args);
    int rank , size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size < 2)
    {
        std::cerr << "Le nombre de processus doit être au moins 2" << std::endl;
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    // auto start = std::chrono::high_resolution_clock::now();

    auto params = parse_arguments(nargs-1, &args[1]);
    display_params(params);
    if (!check_params(params)) return EXIT_FAILURE;

    if(rank == 0) {
        auto displayer = Displayer::init_instance(params.discretization, params.discretization);
        // Buffers para recibir los mapas (asumimos que fire_map y vegetal_map tienen tamaño discretization²)
        std::vector<std::uint8_t> global_fire_map(params.discretization * params.discretization);
        std::vector<std::uint8_t> global_vegetal_map(params.discretization * params.discretization);
        bool continuer = true;
        SDL_Event event;
        while(continuer) {
            // Recibimos los mapas enviados por el proceso simulador (rank 1)
            MPI_Request reqs[2];
            MPI_Irecv(global_fire_map.data(), global_fire_map.size(), MPI_UNSIGNED_CHAR, 1, 0, MPI_COMM_WORLD, &reqs[0]);
            MPI_Irecv(global_vegetal_map.data(), global_vegetal_map.size(), MPI_UNSIGNED_CHAR, 1, 1, MPI_COMM_WORLD, &reqs[1]);
            MPI_Waitall(2, reqs, MPI_STATUSES_IGNORE);
            // Actualizamos la visualización
            displayer->update(global_vegetal_map, global_fire_map);
            // Chequeamos si el usuario cierra la ventana
            if(SDL_PollEvent(&event) && event.type == SDL_QUIT) {
                continuer = false;
            }
            // Enviamos el flag de continuación al proceso de simulación
            int flag = continuer ? 1 : 0;
            MPI_Request req_flag;
            MPI_Isend(&flag, 1, MPI_INT, 1, 2, MPI_COMM_WORLD, &req_flag);
            MPI_Wait(&req_flag, MPI_STATUS_IGNORE);
        }
    }
    // Proceso de simulación: rank 1 (únicamente este hace update())
    else if(rank == 1) {
        auto simu = Model(params.length, params.discretization, params.wind, params.start);
        bool continuer = true;
        while(continuer && simu.update()) {
            if ((simu.time_step() & 31) == 0) 
                std::cout << "Time step " << simu.time_step() << "\n===============" << std::endl;
            // Enviamos el estado actual de los mapas al proceso de visualización

            std::vector<std::uint8_t> tmp_fire = simu.fire_map();
            std::vector<std::uint8_t> tmp_vegetal = simu.vegetal_map();
            MPI_Request reqs[2];

            MPI_Isend(tmp_fire.data(), tmp_fire.size(), MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD, &reqs[0]);
            MPI_Isend(tmp_vegetal.data(), tmp_vegetal.size(), MPI_UNSIGNED_CHAR, 0, 1, MPI_COMM_WORLD, &reqs[1]);
            MPI_Waitall(2, reqs, MPI_STATUSES_IGNORE);
            // Esperamos el flag de continuación desde el visualizador
            int flag = 0;
            MPI_Request req_flag;
            MPI_Irecv(&flag, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &req_flag);
            MPI_Wait(&req_flag, MPI_STATUS_IGNORE);
            continuer = (flag != 0);
        }
    }
    MPI_Finalize();
    return EXIT_SUCCESS;
}
