import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
from io import StringIO

# Pasta onde salvaremos os gráficos
output_folder = "graphs/Amdahl/"
os.makedirs(output_folder, exist_ok=True)

# Procurar arquivos CSV no formato "mpi_*.csv"
# Ex.: "mpi_n100_t4_process2.csv"
csv_files = sorted(glob.glob("../Resultats/mpi_*.csv"))

if not csv_files:
    print("❌ Nenhum arquivo CSV encontrado em ../Resultats/ com padrão 'mpi_*.csv'")
    exit(1)

# time_data[n_value][p_value] = tempo_medio_de_Time_Avancement
time_data = {}

for csv_filename in csv_files:
    base = os.path.basename(csv_filename)    # Ex: "mpi_n100_t4_process2.csv"
    parts = base.split('_')                  # ["mpi", "n100", "t4", "process2.csv"]

    if len(parts) < 4:
        # Não segue o formato esperado
        continue

    # Extrair n_value (ex.: "n100" → "100")
    n_part = parts[1]                       # "n100"
    if not n_part.startswith('n'):
        continue
    n_value = n_part[1:]                    # remove 'n' => "100"

    # Extrair número de processos do quarto elemento: "process2.csv" → "process2"
    # Depois removemos a extensão e "process", restando ex.: "2"
    process_part = parts[3].split('.')[0]   # "process2"
    if not process_part.startswith('process'):
        continue
    p_str = process_part.replace('process', '')  # "2"

    try:
        p_value = int(p_str)
    except ValueError:
        # não é número
        continue

    # Ler o CSV só a partir da linha que contém "TimeStep"
    with open(csv_filename, "r") as f:
        lines = f.readlines()
    try:
        start_line = next(i for i, line in enumerate(lines) if "TimeStep" in line)
    except StopIteration:
        # Não contém a coluna "TimeStep"
        continue

    data_lines = lines[start_line:]
    data_str = "".join(data_lines)
    df = pd.read_csv(StringIO(data_str))

    # Calculamos a média de 'Time_Avancement'
    # (poderia escolher outra métrica, se quiser)
    time_avanc_medio = df["Time_Avancement"].mean()

    # Armazenar no dicionário
    if n_value not in time_data:
        time_data[n_value] = {}
    time_data[n_value][p_value] = time_avanc_medio

# ------------------------------------
# Definimos baseline = 2 processos
baseline_p = 2

# Calcular e plotar Speedup para cada n_value
for n_value, process_time_map in time_data.items():
    if baseline_p not in process_time_map:
        print(f"⚠ Sem baseline (process{baseline_p}) para n={n_value}, pulando.")
        continue

    tempo_seq = process_time_map[baseline_p]  # tempo base (p=2)

    # Ordenar a lista de processos disponíveis
    processes_list = sorted(process_time_map.keys())

    # Calcular Speedup = tempo_seq / tempo_paralelo
    speedups = []
    for p in processes_list:
        tempo_paralelo = process_time_map[p]
        speedup = tempo_seq / tempo_paralelo
        speedups.append(speedup)

    # Plot: Speedup vs número de processos
    plt.figure(figsize=(8, 5))
    plt.plot(processes_list, speedups, marker="o", linestyle="-", color="b")
    plt.xlabel("Número de Processos (MPI)")
    plt.ylabel("Speedup (em relação a p=2)")
    plt.title(f"Speedup (baseline p={baseline_p}) - Grid: {n_value}")
    plt.grid(True)

    # Salvar gráfico
    output_filename = os.path.join(output_folder, f"speedup_n{n_value}_basep{baseline_p}.png")
    plt.savefig(output_filename)
    print(f"📂 Gráfico salvo: {output_filename}")
    plt.close()
