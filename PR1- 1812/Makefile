CFLAGS=-g -Wall -fdiagnostics-color=always -Wno-misleading-indentation

ifeq ($(OS),Windows_NT)
    RM = cmd /C del /Q
    EXT = .exe
else
    RM = rm -f
    EXT =
endif

PROGS=PR1ex2$(EXT) PR1ex3$(EXT) PR1ex4$(EXT) 

all : $(PROGS)

QRsolve_b.o : QRsolve_b.c
	gcc -c $(CFLAGS) $<

PR1ex2$(EXT) : PR1ex2.c QRsolve_b.o
	gcc -o $@ $(CFLAGS) PR1ex2.c QRsolve_b.o -lm
PR1ex3$(EXT) : PR1ex3.c QRsolve_b.o
	gcc -o $@ $(CFLAGS) PR1ex3.c QRsolve_b.o -lm
PR1ex4$(EXT) : PR1ex4.c QRsolve_b.o
	gcc -o $@ $(CFLAGS) PR1ex4.c QRsolve_b.o -lm

clean :
	$(RM) *.o

realclean : clean
	$(RM) $(PROGS)
