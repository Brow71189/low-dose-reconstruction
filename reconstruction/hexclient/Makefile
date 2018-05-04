#clang has better diagnostics and performance since gcc has degraded with ubuntu 16.04
VER=7
#CC=clang++
CC=g++-${VER}
#useful for benchmarking but doubles size of executable
#CC=mpic++

ifeq ($(CC),mpic++)
	MPI_FLAGS=-DUSE_MPI
	EXECUTABLE=hexclientMPI
else
	MPI_FLAGS=
	EXECUTABLE=hexclient
endif
# -U and -D define and undefine the custom flags
# USE_CCTABLE either cache all likelyhoods or recompute them from black probability //heavy memory cost
# USE_MPFR either use built in floating point or handle sepparate mantissa and exponent //heavy performance cost
# MUL_PTBL either look up old and new pixelprobability or precompute all combinations //slight perfromance gain
# USE_MAX enable use_max(0/1) option// slight performance cost for rarely useful feature
# FANCY_FEATURES enables most general transformations to enable EM and HISTO //hardly useful
# MACRO_MPFR provides macro version instead of inline qfr_add and qfr_exp2Z // about 2% extra performance, but no type safety	
MYFLAGS=-DUSE_CCTABLE -DUSE_MPFR -UMUL_PTBL -UUSE_MAX -UFANCY_FEATURES -DMACRO_MPFR
#We have to see if they do anything, well they dont
#ALIGNFLAGS=-falign-functions=16 -falign-loops=16
#math -fassociative-math -freciprocal-math -fno-signed-zeros -fno-trapping-math 
MATHFLAGS=-ffast-math -funsafe-math-optimizations 
#the following flags are not supported by clang and dont help g++ either
#GCSEFLAGS=-fgcse-sm -fgcse-las -fgcse-after-reload 
#-funsafe-loop-optimizations -Wunsafe-loop-optimizations //there is only one case in the default destructor of command
# -fdelete-null-pointer-checks -finline-functions //no effect 
#LOOPFLAGS2=-ftree-loop-im -ftree-loop-if-convert -ftree-loop-if-convert-stores -floop-strip-mine -fgraphite-identity -floop-block -floop-interchange -floop-nest-optimize
LOOPFLAGS=-funswitch-loops -fprefetch-loop-arrays -ftree-loop-distribution 
TOPFLAGS=-ftoplevel-reorder
ifeq ($(CC),clang++)
	SPECIFIC= -O3 
else
	SPECIFIC= -Ofast -flto -fipa-pta -fmodulo-sched -fmodulo-sched-allow-regmoves ${LOOPFLAGS} ${TOPFLAGS}
 
endif
OFLAGS= ${SPECIFIC}
#${MATHFLAGS}  //They just cancel Kahan summations in initialization, and dont matter for updates
DFLAGS=-g3
#good 
#bad -fwhole-program  -static -funroll-loops -fdelete-null-pointer-checks -ftree-vectorize -fno-crossjumping -fweb ${GCSEFLAGS} 
#new -fvariable-expansion-in-unroller -ftracer   
CFLAGS=-c -Wall -pedantic -Wextra -static -std=c++0x ${OFLAGS} ${MPI_FLAGS} ${MYFLAGS}
#CFLAGS=-c -Wall -O0 -ggdb -std=c++0x ${MYFLAGS}
#-s -static -static
#//usr/lib/x86_64-linux-gnu/crt1.o //usr/lib/x86_64-linux-gnu/gcrt1.o //compiles and crashes
#//usr/lib/x86_64-linux-gnu/crtn.o //usr/lib/x86_64-linux-gnu/Mcrt1.o //cannot fint _start
SYSLIBS=//usr/lib/x86_64-linux-gnu/crti.o //usr/lib/x86_64-linux-gnu/Scrt1.o
LDFLAGS= -static -flto -L //usr/lib/gcc/x86_64-linux-gnu/${VER}
SOURCES=hexworker.cpp communicator.cpp command.cpp globals.cpp hexlikelyhood.cpp transformations.cpp zippedimg.cpp hex_em.cpp histogram.cpp abberation.cpp
OBJECTS=$(SOURCES:.cpp=.o)


all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
	
clean:
	rm -rf *o $(EXECUTABLE)

clear:
	rm -rf *o $(EXECUTABLE)	
