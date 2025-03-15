#include <string>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <thread>
#include <chrono>
#include <fstream>
#include <mpi.h>
#include "model.hpp"
#include "display.hpp"

using namespace std::string_literals;
using namespace std::chrono_literals;

struct ParamsType {
    double length{1.};
    unsigned discretization{20u};
    std::array<double,2> wind{0.,0.};
    Model::LexicoIndices start{10u,10u};
    unsigned num_threads{1};
};

void analyze_arg(int nargs, char* args[], ParamsType& params) {
    if(nargs==0)return;
    std::string key(args[0]);
    if(key=="-l"s) {
        if(nargs<2){ std::cerr<<"Manque une valeur pour la longueur du terrain !"<<std::endl; exit(EXIT_FAILURE); }
        params.length = std::stoul(args[1]);
        analyze_arg(nargs-2, &args[2], params); return;
    }
    auto pos = key.find("--longueur=");
    if(pos < key.size()) {
        auto subkey = std::string(key, pos+11);
        params.length = std::stoul(subkey);
        analyze_arg(nargs-1, &args[1], params); return;
    }
    if(key=="-n"s) {
        if(nargs<2){ std::cerr<<"Manque une valeur pour le nombre de cases par direction pour la discrétisation du terrain !"<<std::endl; exit(EXIT_FAILURE); }
        params.discretization = std::stoul(args[1]);
        analyze_arg(nargs-2, &args[2], params); return;
    }
    pos = key.find("--number_of_cases=");
    if(pos < key.size()) {
        auto subkey = std::string(key, pos+18);
        params.discretization = std::stoul(subkey);
        analyze_arg(nargs-1, &args[1], params); return;
    }
    if(key=="-w"s) {
        if(nargs<2){ std::cerr<<"Manque une paire de valeurs pour la direction du vent !"<<std::endl; exit(EXIT_FAILURE); }
        std::string values = std::string(args[1]);
        params.wind[0] = std::stod(values);
        auto pos = values.find(",");
        if(pos == values.size()){ std::cerr<<"Doit fournir deux valeurs séparées par une virgule pour définir la vitesse"<<std::endl; exit(EXIT_FAILURE); }
        auto second_value = std::string(values, pos+1);
        params.wind[1] = std::stod(second_value);
        analyze_arg(nargs-2, &args[2], params); return;
    }
    pos = key.find("--wind=");
    if(pos < key.size()){
        auto subkey = std::string(key, pos+7);
        params.wind[0] = std::stoul(subkey);
        auto pos = subkey.find(",");
        if(pos == subkey.size()){ std::cerr<<"Doit fournir deux valeurs séparées par une virgule pour définir la vitesse"<<std::endl; exit(EXIT_FAILURE); }
        auto second_value = std::string(subkey, pos+1);
        params.wind[1] = std::stod(second_value);
        analyze_arg(nargs-1, &args[1], params); return;
    }
    if(key=="-s"s) {
        if(nargs<2){ std::cerr<<"Manque une paire de valeurs pour la position du foyer initial !"<<std::endl; exit(EXIT_FAILURE); }
        std::string values = std::string(args[1]);
        params.start.column = std::stod(values);
        auto pos = values.find(",");
        if(pos == values.size()){ std::cerr<<"Doit fournir deux valeurs séparées par une virgule pour définir la position du foyer initial"<<std::endl; exit(EXIT_FAILURE); }
        auto second_value = std::string(values, pos+1);
        params.start.row = std::stod(second_value);
        analyze_arg(nargs-2, &args[2], params); return;
    }
    pos = key.find("--start=");
    if(pos < key.size()){
        auto subkey = std::string(key, pos+8);
        params.start.column = std::stoul(subkey);
        auto pos = subkey.find(",");
        if(pos == subkey.size()){ std::cerr<<"Doit fournir deux valeurs séparées par une virgule pour définir la vitesse"<<std::endl; exit(EXIT_FAILURE); }
        auto second_value = std::string(subkey, pos+1);
        params.start.row = std::stod(second_value);
        analyze_arg(nargs-1, &args[1], params); return;
    }
    if(key=="-t"s) {
        if(nargs<2){ std::cerr<<"Manque une valeur pour le nombre de threads !"<<std::endl; exit(EXIT_FAILURE); }
        params.num_threads = std::stoul(args[1]);
        analyze_arg(nargs-2, &args[2], params); return;
    }
    pos = key.find("--threads=");
    if(pos < key.size()){
        auto subkey = std::string(key, pos+10);
        params.num_threads = std::stoul(subkey);
        analyze_arg(nargs-1, &args[1], params); return;
    }
}

ParamsType parse_arguments(int nargs, char* args[]) {
    if(nargs==0)return {};
    if((std::string(args[0])=="--help"s) || (std::string(args[0])=="-h"s)){
        std::cout<<R"RAW(Usage : simulation [option(s)]
  -l, --longueur=LONGUEUR
  -n, --number_of_cases=N
  -w, --wind=VX,VY
  -s, --start=COL,ROW
  -t, --threads=N
)RAW";
        exit(EXIT_SUCCESS);
    }
    ParamsType params;
    analyze_arg(nargs, args, params);
    return params;
}

bool check_params(ParamsType& params) {
    bool flag = true;
    if(params.length<=0){ std::cerr<<"[ERREUR FATALE] La longueur du terrain doit être positive et non nulle !"<<std::endl; flag=false; }
    if(params.discretization<=0){ std::cerr<<"[ERREUR FATALE] Le nombre de cellules par direction doit être positive et non nulle !"<<std::endl; flag=false; }
    if((params.start.row>=params.discretization)||(params.start.column>=params.discretization)){
        std::cerr<<"[ERREUR FATALE] Mauvais indices pour la position initiale du foyer"<<std::endl; flag=false; }
    if(params.num_threads < 1){ std::cerr<<"[ERREUR FATALE] Le nombre de threads doit être au moins 1 !"<<std::endl; flag=false; }
    return flag;
}

void display_params(ParamsType const& params) {
    std::cout<<"Parametres définis pour la simulation : \n"
             <<"\tTaille du terrain : "<<params.length<<"\n"
             <<"\tNombre de cellules par direction : "<<params.discretization<<"\n"
             <<"\tVecteur vitesse : ["<<params.wind[0]<<", "<<params.wind[1]<<"]\n"
             <<"\tPosition initiale du foyer (col, ligne) : "<<params.start.column<<", "<<params.start.row<<"\n"
             <<"\tNombre de threads utilisés : "<<params.num_threads<<"\n";
}

std::string get_timestamp() {
    auto now = std::chrono::system_clock::now();
    std::time_t now_time = std::chrono::system_clock::to_time_t(now);
    std::tm local_time = *std::localtime(&now_time);
    char buffer[50];
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%d_%H-%M-%S", &local_time);
    return std::string(buffer);
}

int main(int nargs, char* args[])
{
    MPI_Init(&nargs, &args);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if(size < 2) {
        std::cerr << "Le nombre de processus doit être au moins 2" << std::endl;
        MPI_Finalize();
        return EXIT_FAILURE;
    }
    
    
    // Lê parâmetros passados pela linha de comando
    auto params = parse_arguments(nargs - 1, &args[1]);

    // (Opcional) Exibir parâmetros no rank 0 para conferência
    if(rank == 0) {
        display_params(params);
    }

    // Verifica validade dos parâmetros
    if(!check_params(params)) {
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    // Abre o arquivo CSV somente no rank 0
    std::string filename =
        "Resultats/mpi_n" + std::to_string(params.discretization) +
        "_t" + std::to_string(params.num_threads) +
        "_process" + std::to_string(size) + ".csv";

    std::ofstream csvFile;
    if(rank == 0) {
        csvFile.open(filename);
        csvFile << "Paramètres de la simulation\n";
        csvFile << "Taille du terrain," << params.length << " km\n";
        csvFile << "Nombre de cellules par direction," << params.discretization << "\n";
        csvFile << "Vecteur vitesse du vent," << params.wind[0] << "," << params.wind[1] << "\n";
        csvFile << "Position initiale du foyer," << params.start.column << "," << params.start.row << "\n";
        csvFile << "Nombre de threads utilisés," << params.num_threads << "\n\n";

        // Cabeçalho do CSV: TimeStep, Time_Avancement, T_Affichage, T_Total
        csvFile << "TimeStep,Time_Avancement,T_Affichage,T_Total\n";
    }

    // -------------------- Rank 0: exibe e grava no CSV --------------------
    if(rank == 0)
    {
        // Cria janela/gráfico
        auto displayer = Displayer::init_instance(
                             params.discretization, params.discretization);

        // Vetores para receber mapas
        std::vector<std::uint8_t> global_fire_map(
            params.discretization * params.discretization);
        std::vector<std::uint8_t> global_vegetal_map(
            params.discretization * params.discretization);

        SDL_Event event;
        bool continuer = true;

        // Marca tempo inicial
        double global_start = MPI_Wtime();
        int iter_count = 0;

        while(continuer)
        {
            // Tempo no início da iteração
            double iter_start = MPI_Wtime();

            // Recebe mapas de outro rank
            MPI_Request reqs[2];
            MPI_Irecv(global_fire_map.data(), global_fire_map.size(),
                      MPI_UNSIGNED_CHAR, 1, 0, MPI_COMM_WORLD, &reqs[0]);
            MPI_Irecv(global_vegetal_map.data(), global_vegetal_map.size(),
                      MPI_UNSIGNED_CHAR, 1, 1, MPI_COMM_WORLD, &reqs[1]);
            MPI_Waitall(2, reqs, MPI_STATUSES_IGNORE);

            // Recebe flag de continuação (1 = continua, 0 = parou)
            int flag;
            MPI_Request req_flag;
            MPI_Irecv(&flag, 1, MPI_INT, 1, 2, MPI_COMM_WORLD, &req_flag);
            MPI_Wait(&req_flag, MPI_STATUS_IGNORE);

            // Mede tempo de exibição
            double display_start = MPI_Wtime();
            displayer->update(global_vegetal_map, global_fire_map);
            double display_end   = MPI_Wtime();
            double display_time  = display_end - display_start;

            // Se fechar a janela ou se flag==0 -> encerrar
            if(SDL_PollEvent(&event) && event.type == SDL_QUIT) {
                continuer = false;
            }
            if(flag == 0) {
                continuer = false;
            }

            // Envia dummy para sinalizar que terminou exibição
            int dummy = 0;
            MPI_Request req_dummy;
            MPI_Isend(&dummy, 1, MPI_INT, 1, 3, MPI_COMM_WORLD, &req_dummy);
            MPI_Wait(&req_dummy, MPI_STATUS_IGNORE);

            // Tempo total da iteração e total global
            double iter_end   = MPI_Wtime();
            double iter_time  = iter_end - iter_start;   // tempo desta iteração
            double total_time = iter_end - global_start; // tempo desde o início

            iter_count++;

            // Grava dados no CSV
            csvFile << iter_count << ","
                    << iter_time << ","
                    << display_time << ","
                    << total_time << "\n";
        }

        csvFile.close();
    }
    // -------------------- Rank 1: realiza a simulação --------------------
    else if(rank == 1)
    {
        Model simu(params.length, params.discretization,
                   params.wind, params.start);

        while(true)
        {
            // Atualiza simulação
            bool cont = simu.update();

            // Captura mapas
            auto tmp_fire    = simu.fire_map();
            auto tmp_vegetal = simu.vegetal_map();

            // Envia mapas para rank 0
            MPI_Request reqs[2];
            MPI_Isend(tmp_fire.data(),    tmp_fire.size(),
                      MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD, &reqs[0]);
            MPI_Isend(tmp_vegetal.data(), tmp_vegetal.size(),
                      MPI_UNSIGNED_CHAR, 0, 1, MPI_COMM_WORLD, &reqs[1]);
            MPI_Waitall(2, reqs, MPI_STATUSES_IGNORE);

            // Se acabou o fogo, manda flag=0 e encerra
            if(!cont) {
                int term_flag = 0;
                MPI_Send(&term_flag, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
                break;
            }

            // Caso contrário, envia flag=1 (continua)
            int flag = 1;
            MPI_Request req_flag;
            MPI_Isend(&flag, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &req_flag);
            MPI_Wait(&req_flag, MPI_STATUSES_IGNORE);

            // Recebe dummy do rank 0 (para saber que ele terminou a exibição)
            int dummy;
            MPI_Request req_dummy;
            MPI_Irecv(&dummy, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, &req_dummy);
            MPI_Wait(&req_dummy, MPI_STATUS_IGNORE);

            // (Opcional) imprimir algo a cada 100 steps
            if((simu.time_step() % 100) == 0) {
                std::cout << "[Rank 1] Time step " 
                          << simu.time_step() << " still burning.\n";
            }
        }
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}
