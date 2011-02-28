
CC           = g++ 
CC-FLAGS     = -c -Wall -g 
LINKER       = g++
LINKER-FLAGS = -o
RM           = rm -f
MAKE         = make
OBJS         = main.o cgsolver.o
TARGET       = main

$(TARGET): $(OBJS) 
	$(LINKER) $(LINKER-FLAGS) $(TARGET) $(OBJS)

main.o : ../source/main.cpp
	$(CC) $(CC-FLAGS) ../source/main.cpp

cgsolver.o : ../source/cgsolver.cpp
	$(CC) $(CC-FLAGS) ../source/cgsolver.cpp


clean : 
	rm -f ./*.o

distclean :
	rm -f ./*.o main

rebuild: clean build
