#clang is about 35% faster for hexsampler
#CC=clang++
CC=g++-7

EXECUTABLE=hexsampler

#math -fassociative-math -freciprocal-math -fno-signed-zeros -fno-trapping-math 
MATHFLAGS=-ffast-math -funsafe-math-optimizations 
#the following flags are not supported by clang and dont help g++ either
#GCSEFLAGS=-fgcse-sm -fgcse-las -fgcse-after-reload
#-funsafe-loop-optimizations -Wunsafe-loop-optimizations //ther is only one case in the default destructor of command 
LOOPFLAGS2=-ftree-loop-im -ftree-loop-if-convert -ftree-loop-if-convert-stores -floop-strip-mine -floop-block 
LLOPFLAGS3=-floop-interchange -fgraphite-identity -floop-nest-optimize -fprefetch-loop-arrays -ftree-loop-distribution
# -fprofile-instr-generate -fprofile-instr-use=hexsampler.profdata 
LOOPFLAGS= -funswitch-loops 
ifeq ($(CC),clang++)
	SPECIFIC= -O3 
else
	SPECIFIC=-Ofast -funswitch-loops 	
endif
OFLAGS=-g ${SPECIFIC} ${MATHFLAGS} 

DFLAGS=-g3  
#-rdynamic included in g3? -DNOCATCH
#good -flto -fipa-pta -fmodulo-sched -fmodulo-sched-allow-regmoves
#bad -fwhole-program  -static -funroll-loops -fdelete-null-pointer-checks -ftree-vectorize -fno-crossjumping -fweb ${GCSEFLAGS} 
#new -fvariable-expansion-in-unroller -ftracer   
CFLAGS=-c -Wall -pedantic -static -std=c++0x ${OFLAGS} 
#CFLAGS=-c -O0 -Wall -pedantic -static -std=c++0x ${DFLAGS}
#-s 
#LDFLAGS= -lmcheck -static -flto -flto  
LDFLAGS= -static -s 
SOURCES1=hexsampler.cpp communicator.cpp command.cpp globals.cpp basis.cpp unitcell.cpp offset.cpp optimizer.cpp hex_optimizer.cpp ellipse.cpp
SOURCES2=hexagons.cpp subframes.cpp xorshift1024star.cpp correlate_target.cpp graphene.cpp
SOURCES=${SOURCES1} ${SOURCES2}

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
