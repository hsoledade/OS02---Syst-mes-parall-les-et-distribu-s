#include <fstream>  // Para manipular arquivos
#include <string>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <thread>
#include <chrono>
#include <ctime>  // Para obter data e hora

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
    unsigned num_threads{1};  // Novo parâmetro para o número de threads
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

    // =====================
    // ADICIONADO - PARÂMETRO PARA O NÚMERO DE THREADS
    // =====================
    if (key == "-t"s)  // Parâmetro "-t N" para definir o número de threads
    {
        if (nargs < 2)
        {
            std::cerr << "Manque une valeur pour le nombre de threads !" << std::endl;
            exit(EXIT_FAILURE);
        }
        params.num_threads = std::stoul(args[1]);  // Converte string para número inteiro
        analyze_arg(nargs - 2, &args[2], params);
        return;
    }

    pos = key.find("--threads=");
    if (pos < key.size())
    {
        auto subkey = std::string(key, pos + 10);
        params.num_threads = std::stoul(subkey);
        analyze_arg(nargs - 1, &args[1], params);
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
    -t, --threads=N             Définit le nombre de threads utilisés pour la simulation.
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

    if (params.num_threads < 1)
    {
        std::cerr << "[ERREUR FATALE] Le nombre de threads doit être au moins 1 !" << std::endl;
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
              << "\tPosition initiale du foyer (col, ligne) : " << params.start.column << ", " << params.start.row << std::endl
              << "\tNombre de threads utilisés : " << params.num_threads << std::endl;
}

// Função para obter data e hora formatada para o nome do arquivo
std::string get_timestamp()
{
    auto now = std::chrono::system_clock::now();
    std::time_t now_time = std::chrono::system_clock::to_time_t(now);
    std::tm local_time = *std::localtime(&now_time);  // Converter para o horário local

    char buffer[50];
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%d_%H-%M-%S", &local_time);  // Formato: YYYY-MM-DD_HH-MM-SS
    return std::string(buffer);
}

int main( int nargs, char* args[] )
{
    // Parsear argumentos
    auto params = parse_arguments(nargs - 1, &args[1]);
    
    // Gerar nome do arquivo CSV com tamanho da grid, número de threads e data/hora
    std::string filename = "Resultats/Sequentiel/n" + std::to_string(params.discretization) + "_t" + std::to_string(params.num_threads) + "_sequentiel" + ".csv";
    
    // Criar arquivo CSV
    std::ofstream csvFile(filename);
    
    // Escrever os parâmetros no início do arquivo CSV
    csvFile << "Paramètres de la simulation\n";
    csvFile << "Taille du terrain," << params.length << " km\n";
    csvFile << "Nombre de cellules par direction," << params.discretization << "\n";
    csvFile << "Vecteur vitesse du vent," << params.wind[0] << "," << params.wind[1] << "\n";
    csvFile << "Position initiale du foyer," << params.start.column << "," << params.start.row << "\n";
    csvFile << "Nombre de threads utilisés," << params.num_threads << "\n\n";
    
    // Escrever cabeçalho das medições
    csvFile << "TimeStep,Time_Avancement,T_Affichage,T_total\n";

    display_params(params);
    if (!check_params(params)) return EXIT_FAILURE;

    auto displayer = Displayer::init_instance(params.discretization, params.discretization);
    auto simu = Model(params.length, params.discretization, params.wind, params.start);

    SDL_Event event;
    while (true)
    {   
        while (SDL_PollEvent(&event))
            if (event.type == SDL_QUIT)
                return EXIT_SUCCESS;

        if ((simu.time_step() & 31) == 0) 
            std::cout << "Time step " << simu.time_step() << "\n===============" << std::endl;

        auto start_total = std::chrono::high_resolution_clock::now();  // Início da simulação

        // Medir tempo de simu.update()
        auto step_start = std::chrono::high_resolution_clock::now();
        bool continua = simu.update();
        auto step_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> step_time = step_end - step_start;

        if (!continua)  // Se a simulação terminou, saímos do loop
            break;
        
        // Medir tempo de displayer.update()
        auto step_start_display = std::chrono::high_resolution_clock::now();
        displayer->update(simu.vegetal_map(), simu.fire_map());
        auto step_end_display = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> step_time_display = step_end_display - step_start_display; 

        // Medir o tempo total
        auto end_total = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> total_time = end_total - start_total;

        // Salvar os tempos no arquivo CSV
        csvFile << simu.time_step() << "," << step_time.count() << "," << step_time_display.count() << "," << total_time.count() << "\n";
    }

    // Fechar o arquivo CSV
    csvFile.close();

    return EXIT_SUCCESS;
}
















