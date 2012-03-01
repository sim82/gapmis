SRC = gapmis.c
OBJ = $(SRC:.c=.o)
CC  = gcc

#SRC_VEC = vec_cpu/gapmis_vec.cpp vec_cpu/main.c
OBJ_VEC = vec_cpu/gapmis_vec.o vec_cpu/main.o
OBJ_NTOM = vec_cpu/gapmis_vec.o vec_cpu/gapmis_ntom.o


CFLAGS = -Wall -O2
CXXFLAGS = -O3 -Wall -msse3 -I.
OFLAGS = -msse3 -O2 -fomit-frame-pointer -funroll-loops  

all: gapmis gapmis_vec gapmis_ntom

gapmis: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(OFLAGS)

gapmis_vec: $(OBJ_VEC)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ_VEC)

gapmis_ntom: $(OBJ_NTOM)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ_NTOM) -lpthread -lboost_thread
clean:
	rm -f gapmis $(OBJ) $(OBJ_VEC) $(OBJ_NTOM) *~ 

$(OBJ) : EDNAFULL.h EBLOSUM62.h gapmis.h errors.h makefile
$(OBJ_VEC) : vec_cpu/aligned_buffer.h  vec_cpu/gapmis_vec.h vec_cpu/cycle.h vec_cpu/vec_unit.h vec_cpu/fasta.h makefile
