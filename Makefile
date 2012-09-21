all: clean
	gcc -Wall -Wextra -o cortex_bin_reader cortex_bin_reader.c

clean:
	rm -rf cortex_bin_reader
