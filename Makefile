cortex_bin_reader: cortex_bin_reader.c
	gcc -Wall -Wextra -o cortex_bin_reader cortex_bin_reader.c -lm

all: cortex_bin_reader

clean:
	rm -rf cortex_bin_reader

.PHONY: all clean
