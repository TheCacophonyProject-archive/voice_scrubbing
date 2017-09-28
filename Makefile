TARGET = mute_low
SRC    = $(wildcard *.c)
HDR    = $(wildcard *.h)
OBJ    = $(SRC:.c=.o)
LIB    = -lm

CC     = gcc
CFLAGS := -W -Wall -pedantic -pedantic -O4

.PHONY: clean

%.o: %.c $(HDR)
	$(CC) -c $< -o $@

$(TARGET): $(OBJ)
	$(CC) $(OBJ) -o $@ $(LIB)

clean:
	rm -f $(OBJ) $(TARGET)
