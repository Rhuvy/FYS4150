all: compile link run clean

Particle:
	g++ -c -I include src/Particle.cpp

PenningTrap:
	g++ -c -I include src/PenningTrap.cpp

full:
	g++ -c -I include main.cpp

compile: Particle PenningTrap full

link:
	g++ -o main Particle.o PenningTrap.o main.o -larmadillo

run:
	./main

clean:
	@rm -rf *.o
