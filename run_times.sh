#!/bin/bash

# Definir os valores de grid (-n), processos (np) e threads (-t)
grids=(100 200 300 400 500)
processes=(1 2 4 8)
threads=(1 2 4 8)

# Loop para percorrer todas as combinações
for n in "${grids[@]}"; do
    for p in "${processes[@]}"; do
        for t in "${threads[@]}"; do
            echo "🔹 Executando simulação: Grid=$n | np=$p | threads=$t"

            # Executar a simulação com mpirun -np <p>, passando também -n e -t
            mpirun -np "$p" ./simulation.exe -n "$n" -t "$t"

            echo "✅ Simulação concluída: Grid=$n | np=$p | threads=$t"
            echo "-----------------------------------------"
        done
    done
done

echo "🎉 Todas as simulações foram executadas!"
