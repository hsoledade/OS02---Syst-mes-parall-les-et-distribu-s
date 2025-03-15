import pandas as pd
import matplotlib.pyplot as plt
import glob
import os

# Criar a pasta de saída se não existir
output_folder = "graphs/Amdahl/"
os.makedirs(output_folder, exist_ok=True)

# Dicionário onde as chaves serão "n100", "n200", etc.
# e o valor será outro dicionário: { thread -> tempo_medio_avanc }
time_data = {}

# Buscar todos os arquivos CSV dentro de Resultats/Amdahl/
csv_files = sorted(glob.glob("../Resultats/Amdahl/*.csv"))

if not csv_files:
    print("❌ Nenhum arquivo CSV encontrado na pasta Resultats/Amdahl/")
    exit(1)

# 1) Ler cada arquivo e extrair n_value, t_value e time_avanc_medio
for csv_filename in csv_files:
    basename = os.path.basename(csv_filename)  # Exemplo: n300_t4_Amdahl.csv
    parts = basename.split("_")
    # Ex: n_value = "n300", t_value = "t4"
    n_value = parts[0]     # "n300"
    t_value = parts[1]     # "t4"

    # Ler o arquivo, ignorando as primeiras linhas de parâmetros
    with open(csv_filename, "r") as f:
        lines = f.readlines()

    # Encontrar a linha onde começam os dados numéricos
    start_line = next(i for i, line in enumerate(lines) if "TimeStep" in line)
    data_lines = lines[start_line:]

    from io import StringIO
    df = pd.read_csv(StringIO("".join(data_lines)))

    # Calcular o valor médio de 'Time_Avancement' (como se fosse o tempo total)
    time_avanc_medio = df["Time_Avancement"].mean()

    # Guardar no dicionário
    if n_value not in time_data:
        time_data[n_value] = {}
    threads = int(t_value[1:])  # remove o 't' => '4'
    time_data[n_value][threads] = time_avanc_medio

# 2) Para cada n_value, calculamos speedup = tempo_seq / tempo_paralelo, usando (t=1) como base
for n_value, thread_time_map in time_data.items():
    # Precisamos do tempo com threads=1
    if 1 not in thread_time_map:
        continue

    tempo_seq = thread_time_map[1]  # tempo médio (t=1)

    # Ordenar as chaves (número de threads)
    threads_list = sorted(thread_time_map.keys())

    # Calcular speedup
    speedups = []
    for t in threads_list:
        tempo_paralelo = thread_time_map[t]
        speedup = tempo_seq / tempo_paralelo
        speedups.append(speedup)

    # 3) Plotar Speedup x Número de Threads
    plt.figure(figsize=(8, 5))
    plt.plot(threads_list, speedups, marker="o", linestyle="-")
    plt.xlabel("Número de Threads")
    plt.ylabel("Speedup (baseado em Time_Avancement médio)")
    plt.title(f"Speedup - Amdahl - Grid: {n_value}")
    plt.grid(True)

    # Salvar o gráfico
    output_filename = f"{output_folder}speedup_{n_value}.png"
    plt.savefig(output_filename)
    print(f"📂 Gráfico salvo: {output_filename}")
    plt.close()
