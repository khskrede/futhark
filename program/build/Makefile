
CC           = g++ 
CC-FLAGS     = -c -Wall -g -fopenmp
LINKER       = g++ -fopenmp
LINKER-FLAGS = -o
RM           = rm -f
MAKE         = make
OBJS         = main.o cgsolver.o diagmatrix.o
TARGET       = main

$(TARGET): $(OBJS) 
	$(LINKER) $(LINKER-FLAGS) $(TARGET) $(OBJS)

main.o : ../source/main.cpp
	$(CC) $(CC-FLAGS) ../source/main.cpp

cgsolver.o : ../source/cgsolver.cpp
	$(CC) $(CC-FLAGS) ../source/cgsolver.cpp

diagmatrix.o : ../source/diagmatrix.cpp
	$(CC) $(CC-FLAGS) ../source/diagmatrix.cpp


clean : 
	rm -f ./*.o

distclean :
	rm -f ./*.o main

rebuild: clean build
