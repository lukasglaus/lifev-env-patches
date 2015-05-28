# define the C compiler to use
CC = g++-4.9

# define any compile-time flags
CFLAGS = -std=c++0x -O3 #-Og -g3

# define any directories containing header files other than /usr/include
INCLUDES = -I /usr/local/Cellar/eigen/3.2.4/include/eigen3

# define the C source files
SRCS = main.cpp Circulation.hpp CirculationCoupling.hpp CirculationIO.hpp CirculationElementFactory.hpp

# Replace the .c with .o for objects
OBJS = main.o Circulation.o CirculationCoupling.o CirculationIO.o CirculationElementFactory.o

# define the executable file 
MAIN = main


default: $(MAIN)

$(MAIN): $(OBJS) 
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $(OBJS)

main.o: $(SRCS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

Circulation.o: Circulation.hpp CirculationCoupling.hpp CirculationIO.hpp CirculationElementFactory.hpp
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

CirculationCoupling.o: CirculationCoupling.hpp
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

CirculationIO.o: CirculationIO.hpp
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

CirculationElementFactory.o: CirculationElementFactory.hpp
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

clean:
	-rm -f *.o core *.core $(MAIN)



# left ventricle treated explicitly

explicit: mainExplicit.o Circulation.o CirculationCoupling.o CirculationIO.o
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ mainExplicit.o Circulation.o CirculationCoupling.o CirculationIO.o

mainExplicit.o: mainExplicit.cpp $(SRCS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<
