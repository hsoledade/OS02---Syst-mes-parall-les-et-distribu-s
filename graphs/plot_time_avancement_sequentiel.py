import pandas as pd
import matplotlib.pyplot as plt
import glob
import os

# Criar a pasta de saída se não existir
output_folder = "graphs/Sequentiel/"
os.makedirs(output_folder, exist_ok=True)

# Criar um dicionário para armazenar arquivos por categoria (n1*, n2*, ..., n5*)
categories = {f"n{i}*": [] for i in range(1, 16)}

# Buscar todos os arquivos CSV dentro de Resultats/Sequentiel/
csv_files = sorted(glob.glob("../Resultats/Sequentiel/*.csv"))

if not csv_files:
    print("❌ Nenhum arquivo CSV encontrado na pasta Resultats/Sequentiel/")
    exit(1)

# Organizar os arquivos por categoria
for csv_filename in csv_files:
    filename = os.path.basename(csv_filename)  # Exemplo: "n300_t4_sequentiel.csv"
    for category in categories:
        if filename.startswith(category[:-1]):  # Remover '*' do nome
            categories[category].append(csv_filename)
            break

# Cores diferentes para cada linha
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange']

# Gerar um gráfico para cada categoria (n1*, n2*, ..., n5*)
for category, files in categories.items():
    if not files:
        continue  # Pular categorias sem arquivos

    plt.figure(figsize=(10, 5))

    # Iterar sobre os arquivos encontrados e plotar
    for i, csv_filename in enumerate(files):
        # Ler o arquivo ignorando as primeiras linhas de parâmetros
        with open(csv_filename, "r") as f:
            lines = f.readlines()

        # Encontrar a linha onde começam os dados numéricos
        start_line = next(i for i, line in enumerate(lines) if "TimeStep" in line)

        # Criar um novo CSV apenas com os dados
        data_lines = lines[start_line:]
        data_csv = "".join(data_lines)

        # Criar um DataFrame Pandas a partir dos dados filtrados
        from io import StringIO
        df = pd.read_csv(StringIO(data_csv))

        # Extrair informações do nome do arquivo
        filename = os.path.basename(csv_filename)  # Exemplo: "n300_t4_sequentiel.csv"
        parts = filename.split("_")
        n_value = parts[0][1:]  # Pega o número após "n"
        t_value = parts[1][1:]  # Pega o número após "t"

        # Plotar os dados (convertendo Time_Avancement para milissegundos)
        plt.plot(df["TimeStep"], 1000 * df["Time_Avancement"], linestyle="-", color=colors[i % len(colors)], label=f"Grid Size={n_value}, Thread={t_value}")

    # Configurar rótulos e título
    plt.xlabel("TimeStep")
    plt.ylabel("Time_Avancement (ms)")
    plt.title(f"Comparaison du Temps d'Avancement par TimeStep - Séquentiel (GridSize = {n_value})")
    plt.legend()  # Adiciona legenda para diferenciar os threads
    plt.grid(True)

    # Salvar o gráfico na pasta graphs/Sequentiel/
    output_filename = f"{output_folder}time_avancement_GridSize_{n_value}.png"
    plt.savefig(output_filename)
    print(f"📂 Gráfico salvo: {output_filename}")

    # Fechar o gráfico para liberar memória
    plt.close()
