# on mac: require export 
# export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/Users/xiaoyang/Desktop/XIAO/program/LIBRARY/boost_1_55_0/stage/lib

PROG=../bin/skeletondb

#BOOSTROOT=/Users/xiaoyang/Desktop/XIAO/program/LIBRARY/boost_1_55_0
#COMPILER=/opt/gcc-4.8.2-build/bin/g++

# server 
BOOSTROOT=/seq/viral/analysis/xyang/programs/Library/boost_1_55_0
COMPILER=/broad/software/free/Linux/redhat_5_x86_64/pkgs/gcc_4.7.2/bin/g++ 
FLAGs=-std=c++0x
OMP=-fopenmp

# mac
#all:
#	$(COMPILER) -I$(BOOSTROOT) -L$(BOOSTROOT)/stage/lib $(OMP) $(FLAGs) -O3 *.cpp -o $(PROG) -lboost_system -lboost_filesystem   	

#linux
all:
	$(COMPILER) -I$(BOOSTROOT) -L$(BOOSTROOT)/stage/lib $(OMP) $(FLAGs) -O3 *.cpp -o $(PROG) -lboost_system -lboost_filesystem -Wl,-rpath,$(BOOSTROOT)/stage/lib

clean:
	rm -rf $(PROG)
