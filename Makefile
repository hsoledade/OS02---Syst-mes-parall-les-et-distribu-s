# Inclusão do arquivo de configuração do sistema operacional
include inc/Make_linux.inc
#include inc/Make_msys2.inc
#include inc/Make_osx.inc

# Definição das flags do compilador
CXXFLAGS = -std=c++17 -Wall -Iinc
ifdef DEBUG
CXXFLAGS += -g -O0 -fbounds-check -pedantic -D_GLIBCXX_DEBUG -fsanitize=address
else
CXXFLAGS += -O3 -march=native
endif

# Diretórios
SRC_DIR = src
OBJ_DIR = obj
INC_DIR = inc
RESULTS_DIR = Resultats
SEQUENTIEL_DIR = $(RESULTS_DIR)/Sequentiel
AMDAHL_DIR = $(RESULTS_DIR)/Amdahl
GUSTAFSON_DIR = $(RESULTS_DIR)/Gustafson

# Criar pastas se não existirem
$(shell mkdir -p $(OBJ_DIR) $(RESULTS_DIR) $(SEQUENTIEL_DIR) $(AMDAHL_DIR) $(GUSTAFSON_DIR))

# Encontrar todos os arquivos .cpp dentro de SRC_DIR
SRC_FILES = $(wildcard $(SRC_DIR)/*.cpp)

# Gerar a lista de objetos equivalente em OBJ_DIR
OBJ_FILES = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRC_FILES))

# Executável final
ALL = simulation.exe

# Alvo padrão
default: help

# Compilar todos os executáveis
all: $(ALL)

# Regras de dependências usando apenas .hpp
$(OBJ_DIR)/display.o: $(SRC_DIR)/display.cpp $(INC_DIR)/display.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR)/model.o: $(SRC_DIR)/model.cpp $(INC_DIR)/model.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR)/simulation.o: $(SRC_DIR)/simulation.cpp $(INC_DIR)/model.hpp $(INC_DIR)/display.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Linkagem final do executável
simulation.exe: $(OBJ_FILES)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LIB)	

# Limpeza dos arquivos gerados
clean:
	@rm -fr $(OBJ_DIR)/*.o 
	@rm -f $(RESULTS_DIR)/*.txt $(RESULTS_DIR)/*.log # Remove arquivos de texto, logs, mas mantém as pastas
	@rm -f *.exe *~

# Exibir ajuda
help:
	@echo "Available targets : "
	@echo "    all            : compile all executables"
	@echo "Add DEBUG=yes to compile in debug"
	@echo "Configuration :"
	@echo "    CXX      :    $(CXX)"
	@echo "    CXXFLAGS :    $(CXXFLAGS)"

# Geração de resultados
results: all
	@echo "Executando simulação e armazenando resultados..."
	./simulation.exe > $(RESULTS_DIR)/Resultados_Parte_1.txt

# Converter arquivos Markdown para HTML (opcional)
%.html: %.md
	pandoc -s --toc $< --css=./github-pandoc.css --metadata pagetitle="OS202 - TD1" -o $@
