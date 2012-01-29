SRC = gapmis.c functions.c output.c
OBJ = $(SRC:.c=.o)
CC  = gcc

CFLAGS = -g -Wall
OFLAGS = 

all: gapmis

gapmis: $(OBJ)
	$(CC) $(CFLAGS)  -o $@  $(OBJ) -lm $(OFLAGS)

clean:
	rm -f gapmis $(OBJ) *~ gapmis.out

$(OBJ) : functions.h EDNAFULL.h EBLOSUM62.h output.h types.h makefile
