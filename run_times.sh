#!/bin/bash

# Definir os valores de grid (-n) e threads (-t)
grids=(100 200 300 400 500)
threads=(1 2 4 8)

# Loop para percorrer todas as combinações
for n in "${grids[@]}"; do
    for t in "${threads[@]}"; do
        echo "🔹 Executando simulação: Grid=$n | Threads=$t"
        
        # Executar a simulação
        ./simulation.exe -n "$n" -t "$t"
        
        echo "✅ Simulação concluída: Grid=$n | Threads=$t"
        echo "-----------------------------------------"
    done
done

echo "🎉 Todas as simulações foram executadas!"
