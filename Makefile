CC        = g++
CFLAGS    = -O3 
OBJECTS =  -llapack -lblas -lgfortran 
LIBS =
INCLUDES =-I . 

all: state  

state: state.cpp exactdiag.h spmatrix.h lanczos.h aux.h green.h
	$(CC) state.cpp $(CFLAGS)  $(OBJECTS)  $(INCLUDES) -o state.out

clean:
	rm -f *.out #*# 

