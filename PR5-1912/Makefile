CFLAGS=-g -Wall -fdiagnostics-color=always -Wno-misleading-indentation

ifeq ($(OS),Windows_NT)
    RM = cmd /C del /Q
    EXT = .exe
else
    RM = rm -f
    EXT =
endif

PROGS=PR5ex2$(EXT) opham_pendol$(EXT) opham_rtbp$(EXT) rtbp_int$(EXT)

all : $(PROGS)

QRsolve_b.o : QRsolve_b.c
	gcc -c $(CFLAGS) $<

Newton_b.o : Newton_b.c QRsolve_b.h
	gcc -c $(CFLAGS) Newton_b.c

rk45_b.o : rk45_b.c
	gcc -c $(CFLAGS) $< -lm

flux_b.o : flux_b.c rk45_b.h
	gcc -c $(CFLAGS) flux_b.c -lm

opham_b.o : opham_b.c QRsolve_b.h flux_b.h rk45_b.h Newton_b.h
	gcc -c $(CFLAGS) opham_b.c -lm

rtbp.o : rtbp.c
	gcc -c $(CFLAGS) $<

pendol.o : pendol.c
	gcc -c $(CFLAGS) $<

PR5ex2$(EXT) : PR5ex2.c opham_b.o pendol.o QRsolve_b.o flux_b.o rk45_b.o Newton_b.o
	gcc -o $@ $(CFLAGS) PR5ex2.c opham_b.o pendol.o QRsolve_b.o flux_b.o rk45_b.o Newton_b.o -lm

opham_pendol$(EXT) : opham_pendol.c opham_b.o pendol.o flux_b.o rk45_b.o Newton_b.o QRsolve_b.o
	gcc -o $@ $(CFLAGS) opham_pendol.c opham_b.o pendol.o flux_b.o rk45_b.o Newton_b.o QRsolve_b.o -lm

opham_rtbp$(EXT) : opham_rtbp.c opham_b.o rtbp.o flux_b.o rk45_b.o Newton_b.o QRsolve_b.o
	gcc -o $@ $(CFLAGS) opham_rtbp.c opham_b.o rtbp.o flux_b.o rk45_b.o Newton_b.o QRsolve_b.o -lm

rtbp_int$(EXT) : rtbp_int.c opham_b.o rtbp.o flux_b.o rk45_b.o Newton_b.o QRsolve_b.o
	gcc -o $@ $(CFLAGS) rtbp_int.c opham_b.o rtbp.o flux_b.o rk45_b.o Newton_b.o QRsolve_b.o -lm

clean :

	$(RM) *.o

realclean : clean
	$(RM) $(PROGS)
