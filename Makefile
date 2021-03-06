.PHONY: all compile io1 io2 io3

all: compile

compile:
	@mpic++ -o main ./main.cpp
io1:
	@mpirun --oversubscribe -np 6 ./main ./input1.tsv
io2:
	@mpirun --oversubscribe -np 6 ./main ./input2.tsv
io3:
	@mpirun --oversubscribe -np 11 ./main ./input3.tsv
