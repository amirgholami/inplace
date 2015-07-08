

all: client

client: client.o utils.o
	gcc -o $@ $^ -L/usr/local/cuda/lib64 -lcudart -lcudart -lgomp -lstdc++ ../build/libinplace.a

utils.o: utils.cpp
	gcc -o $@ -c $^
client.o: client.cu
	nvcc -o $@ -c -arch=sm_35 -Xcompiler -O2 -Xcompiler -Wall -I /usr/local/cuda/include -I /home/agholami/inplace/inplace $^
