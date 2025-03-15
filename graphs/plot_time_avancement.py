import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
from io import StringIO

# Cria a pasta de saída se não existir
output_folder = "graphs/Amdahl/"
os.makedirs(output_folder, exist_ok=True)

# Localiza todos os arquivos CSV que seguem o padrão "mpi_*.csv"
csv_files = sorted(glob.glob("../Resultats/mpi_*.csv"))

if not csv_files:
    print("❌ Nenhum arquivo CSV encontrado em ../Resultats/ com o padrão 'mpi_*.csv'")
    exit(1)

# Dicionário para agrupar arquivos por "n" (exemplo: n=100, 200, 300...)
categories = {}

# Agrupar os arquivos por valor de n (exemplo: "n100" → "100")
for csv_filename in csv_files:
    base = os.path.basename(csv_filename)  # ex.: "mpi_n100_t1_process2.csv"
    parts = base.split('_')                # ["mpi", "n100", "t1", "process2.csv"]

    if len(parts) < 2:
        continue  # não segue o padrão

    # Extrair "n100" e remover o 'n' => "100"
    n_part = parts[1]  # "n100"
    if not n_part.startswith('n'):
        continue
    n_value = n_part[1:]  # remove 'n'

    # Adiciona ao dicionário: categories["100"] -> lista de arquivos
    if n_value not in categories:
        categories[n_value] = []
    categories[n_value].append(csv_filename)

# Definir algumas cores para cada linha plotada
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange', 'purple', 'brown']

# Gerar um gráfico para cada valor de n
for n_value, files in categories.items():
    if not files:
        continue

    plt.figure(figsize=(10, 5))  # Tamanho do gráfico

    for i, csv_filename in enumerate(files):
        with open(csv_filename, "r") as f:
            lines = f.readlines()

        # Encontrar a linha onde começam de fato os dados (que contém "TimeStep")
        try:
            start_line = next(idx for idx, line in enumerate(lines) if "TimeStep" in line)
        except StopIteration:
            print(f"❗ 'TimeStep' não encontrado em {csv_filename}, arquivo ignorado.")
            continue

        # Montar o CSV somente com a parte numérica
        data_csv = "".join(lines[start_line:])
        df = pd.read_csv(StringIO(data_csv))

        # Extrair threads e processos do nome, para compor a legenda
        base = os.path.basename(csv_filename)  # "mpi_n100_t1_process2.csv"
        parts = base.split('_')                # ["mpi", "n100", "t1", "process2.csv"]

        # threads => parts[2], ex.: "t1"
        if len(parts) > 2 and parts[2].startswith('t'):
            t_value = parts[2][1:]  # remove 't'
        else:
            t_value = "?"

        # processos => parts[3], ex.: "process2.csv"
        proc_value = "?"
        if len(parts) > 3 and parts[3].startswith('process'):
            # remover extensão: "process2.csv" -> "process2"
            proc_str = parts[3].split('.')[0]
            # remover "process" -> "2"
            proc_value = proc_str.replace('process', '')

        # Plotar Time_Avancement (ms) vs TimeStep
        plt.plot(df["TimeStep"],
                 1000 * df["Time_Avancement"],  # multiplicando por 1000 => ms
                 linestyle="-",
                 color=colors[i % len(colors)],
                 label=f"GridSize={n_value}, Threads={t_value}, Processes={proc_value}")

    plt.xlabel("TimeStep")
    plt.ylabel("Time_Avancement (ms)")
    plt.title(f"Comparaison du Temps d'Avancement par TimeStep (GridSize = {n_value})")
    plt.grid(True)
    plt.legend()

    # Salvar o gráfico
    output_filename = os.path.join(output_folder, f"time_avancement_GridSize_{n_value}.png")
    plt.savefig(output_filename)
    print(f"📂 Gráfico salvo em: {output_filename}")

    # Fecha o gráfico para liberar memória
    plt.close()
