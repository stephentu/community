all: instance

instance: instance.cpp
	g++ -std=c++11 -Wall -g -O3 -DNDEBUG -DUSE_BOOST_RANDOM -funroll-loops -ffast-math -I/usr/include/eigen3 -I/usr/local/include/eigen3 -o instance instance.cpp

.PHONY: build_on_amp
build_on_amp: instance.cpp
	/opt/amp/gcc/bin/g++ -std=c++0x -Wall -g -O3 -DNDEBUG -DUSE_BOOST_RANDOM -funroll-loops -ffast-math -I/home/eecs/sltu/eigen -I/home/eecs/sltu/boost/include -o instance instance.cpp -lrt

.PHONY: clean
clean:
	rm -rf instance instance.dSYM
