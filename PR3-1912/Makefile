CFLAGS=-g -Wall -fdiagnostics-color=always -Wno-misleading-indentation

ifeq ($(OS),Windows_NT)
    RM = cmd /C del /Q
    EXT = .exe
else
    RM = rm -f
    EXT =
endif

PROGS=PR3ex2$(EXT)  PR3ex3$(EXT) 

all : $(PROGS)

rk45_b.o : rk45_b.c
	gcc -c $(CFLAGS) $< -lm

PR3ex2$(EXT) : PR3ex2.c rk45_b.o
	gcc -o $@ $(CFLAGS) PR3ex2.c rk45_b.o -lm

PR3ex3$(EXT) : PR3ex3.c rk45_b.o
	gcc -o $@ $(CFLAGS) PR3ex3.c rk45_b.o -lm

clean :

	$(RM) *.o

realclean : clean
	$(RM) $(PROGS)
