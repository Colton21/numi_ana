CXXFLAGS=-std=c++11 $(shell root-config --cflags)
LIBS=$(shell root-config --libs)

run: xsec_calc
	./xsec_calc

xsec_calc: xsec_calc.o
	g++ -std=c++11 -ggdb -O0 -o xsec_calc xsec_calc.C ${LIBS}

