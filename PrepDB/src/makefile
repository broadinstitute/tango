PROG=../bin/prepdb

COMPILER=/broad/software/free/Linux/redhat_5_x86_64/pkgs/gcc_4.9.0/bin/g++ 
FLAGs=-std=c++0x
OMP=-fopenmp

all:
	$(COMPILER) $(OMP) $(FLAGs) -O3 *.cpp -o $(PROG)
	
clean:
	rm -rf $(PROG)
