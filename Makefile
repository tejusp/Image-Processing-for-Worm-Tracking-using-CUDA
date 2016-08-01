CC = g++

CFLAGS = `pkg-config --cflags opencv`
LIBS = `pkg-config --libs opencv`
CUDA_INCLUDEPATH=/usr/local/cuda-7.5/include

all: cube.o main.o
	nvcc -o main `pkg-config --libs opencv` main.o cube.o 

cube.o: cube.cu
	nvcc -c cube.cu

main.o: 
	gcc -c main.cpp
  
clean:
	rm -rf main main.o cube.o
