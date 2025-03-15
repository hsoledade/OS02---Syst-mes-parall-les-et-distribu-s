#include <stdexcept>
#include <cmath>
#include <iostream>
#include "model.hpp"

#include <vector>
#include <unordered_map>
#include <omp.h>  // Biblioteca do OpenMP

namespace
{
    double pseudo_random(std::size_t index, std::size_t time_step)
    {
        std::uint_fast32_t xi = std::uint_fast32_t(index*(time_step+10000));
        std::uint_fast32_t r  = (48271*xi)%2147483647;
        return r/2147483646.;
    }

    double log_factor(std::uint8_t value)
    {
        return std::log(1.+value)/std::log(256);
    }
}

// Struct para cada célula local
struct CellUpdate
{
    bool used;          // se essa célula foi atualizada
    std::uint8_t value; // nova intensidade do fogo

    CellUpdate() : used(false), value(0u) {}
};

Model::Model( double t_length, unsigned t_discretization, std::array<double,2> t_wind,
              LexicoIndices t_start_fire_position, double t_max_wind )
    :   m_length(t_length),
        m_distance(-1),
        m_geometry(t_discretization),
        m_wind(t_wind),
        m_wind_speed(std::sqrt(t_wind[0]*t_wind[0] + t_wind[1]*t_wind[1])),
        m_max_wind(t_max_wind),
        m_vegetation_map(t_discretization*t_discretization, 255u),
        m_fire_map(t_discretization*t_discretization, 0u)
{
    if (t_discretization == 0)
    {
        throw std::range_error("Le nombre de cases par direction doit être plus grand que zéro.");
    }
    m_distance = m_length/double(m_geometry);

    // Inicia o fogo na posição solicitada
    auto index = get_index_from_lexicographic_indices(t_start_fire_position);
    m_fire_map[index]   = 255u;
    m_fire_front[index] = 255u;
    constexpr double alpha0 = 4.52790762e-01;
    constexpr double alpha1 = 9.58264437e-04;
    constexpr double alpha2 = 3.61499382e-05;

    if (m_wind_speed < t_max_wind)
        p1 = alpha0 + alpha1*m_wind_speed + alpha2*(m_wind_speed*m_wind_speed);
    else 
        p1 = alpha0 + alpha1*t_max_wind + alpha2*(t_max_wind*t_max_wind);
    p2 = 0.3;

    if (m_wind[0] > 0)
    {
        alphaEastWest = std::abs(m_wind[0]/t_max_wind)+1;
        alphaWestEast = 1.-std::abs(m_wind[0]/t_max_wind);    
    }
    else
    {
        alphaWestEast = std::abs(m_wind[0]/t_max_wind)+1;
        alphaEastWest = 1. - std::abs(m_wind[0]/t_max_wind);
    }

    if (m_wind[1] > 0)
    {
        alphaSouthNorth = std::abs(m_wind[1]/t_max_wind) + 1;
        alphaNorthSouth = 1. - std::abs(m_wind[1]/t_max_wind);
    }
    else
    {
        alphaNorthSouth = std::abs(m_wind[1]/t_max_wind) + 1;
        alphaSouthNorth = 1. - std::abs(m_wind[1]/t_max_wind);
    }
}

bool Model::update()
{
    // Dicionário atual de fogo e uma cópia
    auto next_front = m_fire_front;

    // Lista das células em fogo
    std::vector<std::size_t> fire_cells;
    fire_cells.reserve(m_fire_front.size());
    for (auto& it : m_fire_front)
        fire_cells.push_back(it.first);

    int num_threads = omp_get_max_threads();

    // Ao invés de std::unordered_map por thread, usamos um vetor 2D:
    // local_front[tid][index_da_celula] => struct com 'used' e 'value'
    std::vector<std::vector<CellUpdate>> local_front(num_threads);
    // Para cada thread, alocamos um vetor do tamanho total de células
    for(int t=0; t<num_threads; t++) {
        local_front[t].resize(m_fire_map.size());
    }

    // --------------------------
    // 1) Propagação e decaimento em paralelo
    // --------------------------
    #pragma omp parallel for num_threads(num_threads)
    for (size_t i = 0; i < fire_cells.size(); ++i)
    {
        int tid = omp_get_thread_num();
        std::size_t cell_index = fire_cells[i];
        LexicoIndices coord = get_lexicographic_from_index(cell_index);

        // Intensidade atual do fogo nessa célula
        std::uint8_t old_intensity = m_fire_front[cell_index];
        double power = log_factor(old_intensity);

        // Função lambda para atualizar uma célula no local_front[tid]
        auto update_cell = [&](std::size_t idx, std::uint8_t new_val)
        {
            CellUpdate &cellRef = local_front[tid][idx];
            // Se for a 1a vez => seta
            if(!cellRef.used) {
                cellRef.used  = true;
                cellRef.value = new_val;
            }
            else {
                // Caso já exista valor, pegamos o máximo
                cellRef.value = std::max(cellRef.value, new_val);
            }
        };

        // Função lambda para tentar propagar
        auto try_spread = [&](std::size_t neighbor_index, double factor, std::size_t seedFactor)
        {
            double tirage      = pseudo_random(cell_index * seedFactor + m_time_step, m_time_step);
            double green_power = m_vegetation_map[neighbor_index];
            double correction  = power * log_factor(green_power);

            if (tirage < factor * p1 * correction)
            {
                update_cell(neighbor_index, 255u);
            }
        };

        // ===========================================
        // Propagação
        // ===========================================
        if (coord.row < m_geometry - 1)
        {
            try_spread(cell_index + m_geometry, alphaSouthNorth, 1ul);
        }
        if (coord.row > 0)
        {
            try_spread(cell_index - m_geometry, alphaNorthSouth, 13427ul);
        }
        if (coord.column < m_geometry-1)
        {
            try_spread(cell_index + 1, alphaEastWest, 13427ul*13427ul);
        }
        if (coord.column > 0)
        {
            try_spread(cell_index - 1, alphaWestEast, 13427ul*13427ul*13427ul);
        }

        // ===========================================
        // Decaimento
        // ===========================================
        if (old_intensity == 255u)
        {
            double tirage = pseudo_random(cell_index*52513 + m_time_step, m_time_step);
            if (tirage < p2)
            {
                std::uint8_t half_val = old_intensity >> 1; 
                update_cell(cell_index, half_val);
            }
        }
        else
        {
            std::uint8_t half_val = old_intensity >> 1;
            update_cell(cell_index, half_val);
        }
    }

    // --------------------------
    // 2) Consolidar os dados de local_front no next_front
    // --------------------------
    for(int t=0; t<num_threads; t++)
    {
        for(size_t idx=0; idx<local_front[t].size(); idx++)
        {
            if(local_front[t][idx].used) 
            {
                std::uint8_t val = local_front[t][idx].value;
                // Faz merge com o que está em next_front
                auto it = next_front.find(idx);
                if(it == next_front.end()) 
                {
                    // Não havia fogo nessa célula
                    if(val == 0)
                    {
                        // Se val=0 => nada
                        m_fire_map[idx] = 0;
                    }
                    else
                    {
                        next_front[idx] = val;
                        m_fire_map[idx] = val;
                    }
                }
                else
                {
                    // Já existia => pegamos o max
                    std::uint8_t old_val = it->second; 
                    bool wasIgnition = (old_val == 255u);
                    bool newIgnition = (val == 255u);
                    std::uint8_t merged;

                    if (wasIgnition && !newIgnition)
                    {
                        // Célula estava em 255, agora 'val' é menor => decaimento
                        merged = std::min(old_val, val);
                    }
                    else if (!wasIgnition && newIgnition)
                    {
                        // Agora surgiu ignição => 255
                        merged = std::max(old_val, val);
                    }
                    else if (wasIgnition && newIgnition)
                    {
                        // Ambos ignição => 255
                        merged = 255u;
                    }
                    else
                    {
                        // Ambos decaimento => pega o menor
                        merged = std::min(old_val, val);
                    }

                    // Se ficou 0, apaga do dicionário
                    if (merged == 0) {
                        next_front.erase(idx);
                        m_fire_map[idx] = 0;
                    }
                    else {
                        next_front[idx] = merged;
                        m_fire_map[idx] = merged;
                    }

                }
            }
        }
    }

    // --------------------------
    // 3) Reduz a vegetação
    // --------------------------
    #pragma omp parallel for num_threads(num_threads)
    for(size_t i=0; i<fire_cells.size(); i++)
    {
        std::size_t idx = fire_cells[i];
        if(m_vegetation_map[idx] > 0)
        {
            m_vegetation_map[idx] -= 1;
        }
    }

    // Atualiza
    m_fire_front = std::move(next_front);
    m_time_step += 1;

    return !m_fire_front.empty();
}

std::size_t Model::get_index_from_lexicographic_indices(LexicoIndices t_lexico_indices) const
{
    return t_lexico_indices.row*this->geometry() + t_lexico_indices.column;
}

auto Model::get_lexicographic_from_index(std::size_t t_global_index) const -> LexicoIndices
{
    LexicoIndices ind_coords;
    ind_coords.row    = t_global_index / this->geometry();
    ind_coords.column = t_global_index % this->geometry();
    return ind_coords;
}
