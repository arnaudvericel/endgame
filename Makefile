FC=gfortran
FFLAGS= -g -Wall -Wextra -fdefault-real-8 -fdefault-double-8
DEBUG=-fcheck=all
#FFLAGS+=${DEBUG}
ifeq ($(THANOS), yes)
    FFLAGS += -DTHANOS
endif
SRC=config.F90 functions.F90 initial.F90 evol.F90 endgame.F90
OBJ=${SRC:.F90=.o}

%.o: %.F90
	$(FC) $(FFLAGS) -o $@ -c $<

endgame: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

all: endgame clean

clean:
	rm *.o *.mod
	rm -r *.dSYM

endgame.o: config.o functions.o evol.o

evol.o: config.o functions.o

functions.o: config.o
