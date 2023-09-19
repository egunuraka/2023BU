all:
	g++ -Wall src/methods/explicit_rk.hpp src/methods/explicit_rk.cpp src/parse.hpp src/parse.cpp src/methods/matrix.cpp src/methods/matrix.hpp src/methods/vector.cpp src/methods/vector.hpp src/main.cpp -o main
clear:
	rm -f trajectory main
