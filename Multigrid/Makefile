.DEFAULT_GOAL := all
CC = gcc
CFLAGS = -Wall -Wextra -O2
EXTRA_CFLAGS =
LDLIBS = -lm

include config.mk

TARGET = multigrid
COMMON_SRCS = main.c hierarquia.c util.c matriz.c operadores.c

ifeq ($(PROBLEMA),POISSON)
PROBLEM_SRCS = poisson.c poisson_cli.c poisson_plot.c problema.c residuo.c suavizacao.c erro.c vcycle_poisson.c
else ifeq ($(PROBLEMA),CALOR)
PROBLEM_SRCS = eqcalor.c eqcalor_cli.c calor_operacoes.c calor_fmg.c calor_plot.c temporal.c espacial.c calor.c residuo_calor.c suavizacao_calor.c vcycle_calor.c
else
$(error PROBLEMA deve ser POISSON ou CALOR em config.mk)
endif

ifeq ($(DIMENSAO),2D)
CFLAGS += -DDIM2
else ifeq ($(DIMENSAO),3D)
CFLAGS += -DDIM3
else ifeq ($(DIMENSAO),1D)
CFLAGS += -DDIM1
else
$(error DIMENSAO deve ser 1D 2D ou 3D em config.mk)
endif

ifeq ($(PRINT), ON)
CFLAGS += -DPRINT
endif

ifeq ($(DEBUG), ON)
CFLAGS += -DDEBUG
endif

# Cada configuração tem objetos próprios: trocar de problema ou de dimensão
# não exige apagar nem sobrescrever os objetos da outra configuração.
SOURCES = $(COMMON_SRCS) $(PROBLEM_SRCS)
BUILD = build/$(PROBLEMA)_$(DIMENSAO)_$(PRINT)_$(DEBUG)
OBJS = $(addprefix $(BUILD)/normal/,$(SOURCES:.c=.o))
PLOT_OBJS = $(addprefix $(BUILD)/plot/,$(SOURCES:.c=.o))
PROFILE_OBJS = $(addprefix $(BUILD)/profile/,$(SOURCES:.c=.o))
-include $(OBJS:.o=.d) $(PLOT_OBJS:.o=.d) $(PROFILE_OBJS:.o=.d)

.PHONY: all clean poisson zip heat problem config plot unprint print 1d 2d 3d help valid debug undebug video profile plot-build

all: $(TARGET)

$(TARGET): $(OBJS) config.mk
	$(CC) $(CFLAGS) $(EXTRA_CFLAGS) $(OBJS) -o $@ $(LDLIBS)

$(BUILD)/normal/%.o: %.c
	@mkdir -p "$(dir $@)"
	$(CC) $(CFLAGS) $(EXTRA_CFLAGS) -MMD -MP -c $< -o $@

$(BUILD)/plot/%.o: %.c
	@mkdir -p "$(dir $@)"
	$(CC) $(CFLAGS) $(EXTRA_CFLAGS) -DPLOT -MMD -MP -c $< -o $@

$(BUILD)/profile/%.o: %.c
	@mkdir -p "$(dir $@)"
	$(CC) $(CFLAGS) $(EXTRA_CFLAGS) -pg -MMD -MP -c $< -o $@

# Executa tarefas auxiliares em uma cópia temporária do código.
# Nenhum .o ou executável do diretório do usuário é removido.
# O nome do executável auxiliar é exclusivo de cada tarefa.
# Ferramentas auxiliares possuem executáveis distintos e jamais substituem multigrid.
multigridvalid: $(OBJS) config.mk
	$(CC) $(CFLAGS) $(EXTRA_CFLAGS) $(OBJS) -o $@ $(LDLIBS)

multigridplot: $(PLOT_OBJS) config.mk
	$(CC) $(CFLAGS) $(EXTRA_CFLAGS) $(PLOT_OBJS) -o $@ $(LDLIBS)

multigridvideo: $(PLOT_OBJS) config.mk
	$(CC) $(CFLAGS) $(EXTRA_CFLAGS) $(PLOT_OBJS) -o $@ $(LDLIBS)

multigridprofile: $(PROFILE_OBJS) config.mk
	$(CC) $(CFLAGS) $(EXTRA_CFLAGS) -pg $(PROFILE_OBJS) -o $@ $(LDLIBS)

plot-build: multigridplot

poisson:
	@sed -i 's/^PROBLEMA = .*/PROBLEMA = POISSON/' config.mk
	@echo "Problema selecionado: Poisson"

heat:
	@sed -i 's/^PROBLEMA = .*/PROBLEMA = CALOR/' config.mk
	@echo "Problema selecionado: Equacao do calor"

1d:
	@sed -i 's/^DIMENSAO = .*/DIMENSAO = 1D/' config.mk
	@echo "Dimensao selecionada: 1D"

2d:
	@sed -i 's/^DIMENSAO = .*/DIMENSAO = 2D/' config.mk
	@echo "Dimensao selecionada: 2D"

3d:
	@sed -i 's/^DIMENSAO = .*/DIMENSAO = 3D/' config.mk
	@echo "Dimensao selecionada: 3D"

problem:
	@echo "Problema selecionado: $(PROBLEMA)"
	@echo "Dimensao selecionada: $(DIMENSAO)"

clean:
	rm -rf build
	rm -f multigrid multigridvalid multigridplot multigridvideo multigridprofile *.o *.txt *.out *.zip *.vtk *.dat

config:
	@cat config.mk

valid: multigridvalid
ifeq ($(PROBLEMA),POISSON)
	@set -e; trap 'rm -f multigridvalid' EXIT; LC_NUMERIC=C bash ./s.sh
else
	@set -e; trap 'rm -f multigridvalid' EXIT; LC_NUMERIC=C bash ./calor.sh
endif

print:
	@sed -i 's/^PRINT = .*/PRINT = ON/' config.mk
	@$(MAKE) config

unprint:
	@sed -i 's/^PRINT = .*/PRINT = OFF/' config.mk
	@$(MAKE) config

debug:
	@sed -i 's/^DEBUG = .*/DEBUG = ON/' config.mk
	@$(MAKE) config

undebug:
	@sed -i 's/^DEBUG = .*/DEBUG = OFF/' config.mk
	@$(MAKE) config

plot: multigridplot
	@set -e; trap 'rm -f multigridplot' EXIT; LC_NUMERIC=C bash plot.sh

video:
ifeq ($(PROBLEMA),CALOR)
	@set -e; trap 'rm -f multigridvideo' EXIT; $(MAKE) --no-print-directory multigridvideo; LC_NUMERIC=C bash ./video.sh
else
	@echo "Nao e possivel fazer video do problema estacionario de Poisson. Selecione 'make heat' primeiro."
endif

zip:
	@bash zip.sh

ifeq ($(firstword $(MAKECMDGOALS)),profile)
PROFILE_N := $(word 2,$(MAKECMDGOALS))
ifneq ($(PROFILE_N),)
.PHONY: $(PROFILE_N)
$(PROFILE_N):
	@:
endif
endif

profile: multigridprofile
	@set -e; trap 'rm -f multigridprofile' EXIT; bash profile.sh $(PROFILE_N)

help:
	
	@echo "make               compila somente o executavel do usuario: ./multigrid"
	@echo "make problem/config mostra problema e dimensao selecionados"
	@echo "make poisson/heat  escolhe qual implementacao o Makefile vai ligar"
	@echo "make 1d/2d/3d     seleciona a dimensao na compilacao"
	@echo ""
	@echo "ATENÇÃO: o make vai sempre compilar um único executável e para compilar é sempre 'make'. Nunca outro. Se a configuração for POISSON 2D, então 'make' compilará executável para essa configuração. Se quiser trocar para, digamos EQUAÇÃO DO CALOR 1D, fazer primeiro 'make heat'; 'make 1d'; e somente então compilar."
	@echo ""	
	@echo "make valid         imprime a tabela de validação da discretização, com cálculo da ordem de convergência apropriada"
	@echo "make plot          faz os gráficos da solução analítica, da solução computada e do erro absoluto"
	@echo "make video         gera vídeo da evolução da equação do calor"
	@echo "make profile       imprime relatório gprof"
	@echo "make print/unprint e make debug/undebug configuram a proxima compilacao. Aqui, PRINT imprimirá a evolução do resíduo a cada ciclo na tela e DEBUG gerará relatório completo em arquivo"
	@echo "make clean         remove executaveis, objetos e arquivos temporarios"
	@echo "make zip           compacta o diretorio"
	@echo ""
	@echo "Sintaxe para problema de Poisson: ./multigrid <malha> <malha fina> ou ./multigrid v 1 <...>"
	@echo "Sintaxe para Equação do Calor: ./multigrid <malha> <dt> | s <malha1> <malha2> <dt> | t <malha> <dt1> <dt2> | n <malha> <passos1> [passos2]"
