CC=gcc
CFLAGS=-Wall -Wextra
LDFLAGS=-lm

cortex_bin_reader: cortex_bin_reader.c
	$(CC) $(CFLAGS) -o cortex_bin_reader cortex_bin_reader.c $(LDFLAGS)

all: cortex_bin_reader

clean:
	rm -rf cortex_bin_reader

.PHONY: all clean
