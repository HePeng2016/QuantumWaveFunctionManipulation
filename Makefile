objects = WFM/WavefunctionManipulation.o
CC = g++  -std=c++11  -g  -I WFM/
LIB_FLAGS = -lgsl -lgslcblas
all:    prepareall
	$(CC)$(CXXFLAGS)   -o   WaveFun    main.cpp  $(objects) $(LIB_FLAGS)
prepareall:    subsystem
subsystem:
	$(MAKE) -C WFM 
clean :  cleansub
	rm     WaveFun
cleansub :
	rm  $(objects)
  
