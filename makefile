SRC = gapmis.c
OBJ = $(SRC:.c=.o)
CC  = gcc

CFLAGS = -Wall -msse3 -O3 -fomit-frame-pointer -funroll-loops  

all: gapmis

gapmis: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ)

clean:
	rm -f gapmis $(OBJ) *~ 

$(OBJ) : EDNAFULL.h EBLOSUM62.h gapmis.h errors.h makefile
