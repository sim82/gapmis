SRC = gapmis.c
OBJ = $(SRC:.c=.o)
CC  = gcc

CFLAGS = -Wall
OFLAGS = -msse3 -O2 -fomit-frame-pointer -funroll-loops  

all: gapmis

gapmis: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(OFLAGS)

clean:
	rm -f gapmis $(OBJ) *~ 

$(OBJ) : EDNAFULL.h EBLOSUM62.h gapmis.h makefile
