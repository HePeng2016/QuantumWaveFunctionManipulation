objects = cints.o  chgp.o  
CC = g++  -g   #-Wall
install : $(objects) 

cints.o : cints.c  cints.h   
	$(CC) -c $<
chgp.o : chgp.c	chgp.h cints.h
	$(CC) -c $<
clean : 
	rm  $(objects)


