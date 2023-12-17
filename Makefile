# Makefile Example

# *****************************************************
# Variables to control Makefile operation

PGCX = g++
PGCXFLAGSOMP = -O2 -fopenmp

# ****************************************************
# Targets needed to bring the executable up to date
	
mainpg: main.o nearestIBpoint.o assign_flag.o evalMaxError.o solverUVGaussSeidelSOR.o solverPGaussSeidelSOR.o
	$(PGCX) $(PGCXFLAGSOMP) -o run_executable.out main.o nearestIBpoint.o assign_flag.o evalMaxError.o solverUVGaussSeidelSOR.o solverPGaussSeidelSOR.o
	
# The main.o target can be written more simply
	
main.o: main.o nearestIBpoint.h assign_flag.h evalMaxError.h solverUVGaussSeidelSOR.h solverPGaussSeidelSOR.h
	$(PGCX) $(PGCXFLAGSOMP) -c main.cpp
	
nearestIBpoint.o: nearestIBpoint.h
	$(PGCX) $(PGCXFLAGSOMP) -c nearestIBpoint.cpp
    
assign_flag.o: assign_flag.h nearestIBpoint.h
	$(PGCX) $(PGCXFLAGSOMP) -c assign_flag.cpp

evalMaxError.o: evalMaxError.h
	$(PGCX) $(PGCXFLAGSOMP) -c evalMaxError.cpp

solverUVGaussSeidelSOR.o: solverUVGaussSeidelSOR.h evalMaxError.h
	$(PGCX) $(PGCXFLAGSOMP) -c solverUVGaussSeidelSOR.cpp

solverPGaussSeidelSOR.o: solverPGaussSeidelSOR.h evalMaxError.h
	$(PGCX) $(PGCXFLAGSOMP) -c solverPGaussSeidelSOR.cpp
