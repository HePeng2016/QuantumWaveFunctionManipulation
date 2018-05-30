Objects = WavafunctionManipulation.o
CC = g++  -std=c++11   -g   #-Wall
LIB_FLAGS = -lgsl -lgslcblas
install : $(Objects)

WavafunctionManipulation.o : WavefunctionManipulation.cpp WavafunctionManipulation.h 
	$(CC) -c $<
clean : 
	rm  $(objects)





