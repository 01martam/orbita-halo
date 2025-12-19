CFLAGS=-g -Wall -fdiagnostics-color=always -Wno-misleading-indentation

ifeq ($(OS),Windows_NT)
    RM = cmd /C del /Q
    EXT = .exe
else
    RM = rm -f
    EXT =
endif

PROGS=PR4ex1$(EXT) PR4ex2$(EXT) PR4ex3$(EXT) 

all : $(PROGS)


QRsolve_b.o : QRsolve_b.c
	gcc -c $(CFLAGS) $<

rk45_b.o : rk45_b.c
	gcc -c $(CFLAGS) $< -lm

flux_b.o : flux_b.c rk45_b.h
	gcc -c $(CFLAGS) flux_b.c -lm

PR4ex1 : PR4ex1.c flux_b.o rk45_b.o
	gcc -o $@ $(CFLAGS) PR4ex1.c flux_b.o rk45_b.o -lm

PR4ex2 : PR4ex2.c flux_b.o rk45_b.o
	gcc -o $@ $(CFLAGS) PR4ex2.c flux_b.o rk45_b.o -lm

PR4ex3 : PR4ex3.c flux_b.o rk45_b.o QRsolve_b.o
	gcc -o $@ $(CFLAGS) PR4ex3.c flux_b.o rk45_b.o QRsolve_b.o -lm
clean :

	$(RM) *.o

realclean : clean
	$(RM) $(PROGS)
