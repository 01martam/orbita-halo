CFLAGS=-g -Wall -fdiagnostics-color=always -Wno-misleading-indentation

ifeq ($(OS),Windows_NT)
    RM = cmd /C del /Q
    EXT = .exe
else
    RM = rm -f
    EXT =
endif

PROGS=PR2ex2$(EXT) PR2ex3$(EXT) 

all : $(PROGS)

QRsolve_b.o : QRsolve_b.c
	gcc -c $(CFLAGS) $<

TestFunction.o : TestFunction.c
	gcc -c $(CFLAGS) $<

Newton_b.o : Newton_b.c QRsolve_b.h
	gcc -c $(CFLAGS) Newton_b.c

PR2ex2 : PR2ex2.c QRsolve_b.o Newton_b.o TestFunction.o
	gcc -o $@ $(CFLAGS) PR2ex2.c QRsolve_b.o Newton_b.o TestFunction.o -lm


PR2ex3 : PR2ex3.c QRsolve_b.o Newton_b.o TestFunction.o
	gcc -o $@ $(CFLAGS) PR2ex3.c QRsolve_b.o Newton_b.o TestFunction.o -lm
clean :
	$(RM) *.o

realclean : clean
	$(RM) $(PROGS)
