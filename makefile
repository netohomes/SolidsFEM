# VARIABLES DE COMPILACIÓN
# Compilador
CC = gcc
# Opciones de compilación
FLAGS = -Wall -g
# Nombre del programa
EXE = output
# Source files
SOURCE_FILES = cond_cont.c esfuerzos.c grad_conj.c lectura.c matriz_elemental.c matriz_global.c memoria.c oper_matrices.c principal.c result.c

all:
	@echo Compilando...
	$(CC) $(FLAGS) $(SOURCE_FILES) -o $(EXE)
clean:
	@echo Borrando ejecutable...
	rm $(EXE)
