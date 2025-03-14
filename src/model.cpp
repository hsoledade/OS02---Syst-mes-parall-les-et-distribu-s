#include <stdexcept>
#include <cmath>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include "model.hpp"

namespace
{
    double pseudo_random( std::size_t index, std::size_t time_step )
    {
        std::uint_fast32_t xi = std::uint_fast32_t(index*(time_step+1));
        std::uint_fast32_t r  = (48271*xi)%2147483647;
        return r/2147483646.;
    }

    double log_factor( std::uint8_t value )
    {
        return std::log(1.+value)/std::log(256);
    }
}

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
    auto index = get_index_from_lexicographic_indices(t_start_fire_position);
    m_fire_map[index] = 255u;
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

// -----------------------------------------------------------------------------
// Esta función actualiza el bloque local (incluyendo halos) y retorna true si
// aún hay fuego activo en el dominio interno.
bool Model::update_local(std::vector<std::uint8_t>& local_fire,
                           std::vector<std::uint8_t>& local_vegetal,
                           int local_rows, int cols)
{
    // Creamos una copia para calcular el nuevo estado.
    std::vector<std::uint8_t> new_fire = local_fire;
    // Se iteran únicamente las celdas internas: filas 1 .. local_rows
    for (int i = 1; i <= local_rows; i++) {
        for (int j = 0; j < cols; j++) {
            int idx = i * cols + j;
            // Solo se procesan celdas encendidas
            if (local_fire[idx] == 0)
                continue;
            double power = log_factor(local_fire[idx]);

            // Vecino superior (usa fila 0 si i==1)
            if (i > 0) {
                int top_idx = (i - 1) * cols + j;
                double tirage = pseudo_random(idx + 1, m_time_step);
                double green_power = local_vegetal[top_idx];
                double correction = power * log_factor(green_power);
                if (tirage < alphaSouthNorth * p1 * correction)
                    new_fire[top_idx] = 255;
            }

            // Vecino inferior (fila i+1, puede estar en el halo inferior)
            if (i < local_rows + 1) {
                int bottom_idx = (i + 1) * cols + j;
                double tirage = pseudo_random(idx + 2, m_time_step);
                double green_power = local_vegetal[bottom_idx];
                double correction = power * log_factor(green_power);
                if (tirage < alphaNorthSouth * p1 * correction)
                    new_fire[bottom_idx] = 255;
            }

            // Vecino izquierdo
            if (j > 0) {
                int left_idx = i * cols + (j - 1);
                double tirage = pseudo_random(idx + 3, m_time_step);
                double green_power = local_vegetal[left_idx];
                double correction = power * log_factor(green_power);
                if (tirage < alphaEastWest * p1 * correction)
                    new_fire[left_idx] = 255;
            }

            // Vecino derecho
            if (j < cols - 1) {
                int right_idx = i * cols + (j + 1);
                double tirage = pseudo_random(idx + 4, m_time_step);
                double green_power = local_vegetal[right_idx];
                double correction = power * log_factor(green_power);
                if (tirage < alphaWestEast * p1 * correction)
                    new_fire[right_idx] = 255;
            }

            // Decaimiento del fuego
            if (local_fire[idx] == 255) {
                double tirage = pseudo_random(idx * 2, m_time_step);
                if (tirage < p2) {
                    new_fire[idx] >>= 1;
                }
            } else {
                new_fire[idx] >>= 1;
            }
        }
    }

    // Reducir la vegetación en las celdas que están en fuego
    for (int i = 1; i <= local_rows; i++) {
        for (int j = 0; j < cols; j++) {
            int idx = i * cols + j;
            if (new_fire[idx] > 0 && local_vegetal[idx] > 0)
                local_vegetal[idx]--;
        }
    }

    // Actualizar el estado interno y determinar si hay fuego activo
    bool active = false;
    for (int i = 1; i <= local_rows; i++) {
        for (int j = 0; j < cols; j++) {
            int idx = i * cols + j;
            local_fire[idx] = new_fire[idx];
            if (local_fire[idx] > 0)
                active = true;
        }
    }
    m_time_step++; // incrementar el tiempo global
    return active;
}

// -----------------------------------------------------------------------------
// Métodos auxiliares que mapean índices globales (sin celdas fantasma)
std::size_t Model::get_index_from_lexicographic_indices(LexicoIndices t_lexico_indices) const
{
    return t_lexico_indices.row * this->geometry() + t_lexico_indices.column;
}

auto Model::get_lexicographic_from_index(std::size_t t_global_index) const -> LexicoIndices
{
    LexicoIndices ind_coords;
    ind_coords.row = t_global_index / this->geometry();
    ind_coords.column = t_global_index % this->geometry();
    return ind_coords;
}
