all: instance

instance: instance.cpp
	g++ -std=c++11 -Wall -g -O3 -DNDEBUG -DUSE_BOOST_RANDOM -funroll-loops -ffast-math -I/usr/include/eigen3 -I/usr/local/include/eigen3 -o instance instance.cpp

.PHONY: clean
clean:
	rm -rf instance instance.dSYM
