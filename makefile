SRC = gapmis.c
OBJ = $(SRC:.c=.o)
CC  = gcc

CFLAGS = -Wall
OFLAGS = 

all: gapmis

gapmis: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(OFLAGS)

clean:
	rm -f gapmis $(OBJ) *~ 

$(OBJ) : EDNAFULL.h EBLOSUM62.h gapmis.h makefile
