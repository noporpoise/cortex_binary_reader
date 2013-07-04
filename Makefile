ifndef CC
	CC=gcc
endif

CFLAGS=-Wall -Wextra
LDFLAGS=-lm
DEBUG_ARGS=

ifdef DEBUG
	DEBUG_ARGS=-g -ggdb -DDEBUG=1
endif

cortex_bin_reader: cortex_bin_reader.c
	$(CC) $(CFLAGS) $(DEBUG_ARGS) -o cortex_bin_reader cortex_bin_reader.c $(LDFLAGS)

all: cortex_bin_reader

clean:
	rm -rf cortex_bin_reader

.PHONY: all clean
