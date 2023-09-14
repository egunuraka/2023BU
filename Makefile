all:
	g++ -Wall src/methods/explicit_rk.hpp src/methods/explicit_rk.cpp src/parse.hpp src/parse.cpp src/main.cpp -o main
clear:
	rm -f trajectory barycenter_coorinates barycenter_speed energy potential_energy kinetic_energy angular_momentum analyt main
