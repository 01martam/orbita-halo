CFLAGS=-g -Wall -fdiagnostics-color=always -Wno-misleading-indentation

ifeq ($(OS),Windows_NT)
    RM = cmd /C del /Q
    EXT = .exe
else
    RM = rm -f
    EXT =
endif

PROGS=PR3ex2$(EXT) PR3ex2v2$(EXT)  PR3ex3$(EXT) Test-RK45$(EXT)

all : $(PROGS)

rk45_b.o : rk45_b.c
	gcc -c $(CFLAGS) $< -lm

PR3ex2 : PR3ex2.c rk45_b.o
	gcc -o $@ $(CFLAGS) PR3ex2.c rk45_b.o -lm

PR3ex2v2 : PR3ex2v2.c rk45_b.o
	gcc -o $@ $(CFLAGS) PR3ex2v2.c rk45_b.o -lm

PR3ex3 : PR3ex3.c rk45_b.o
	gcc -o $@ $(CFLAGS) PR3ex3.c rk45_b.o -lm

Test-RK45 : Test-RK45.c rk45_b.o
	gcc -o $@ $(CFLAGS) Test-RK45.c rk45_b.o -lm

clean :

	$(RM) *.o

realclean : clean
	$(RM) $(PROGS)
