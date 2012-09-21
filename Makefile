all:
	gcc -Wall -Wextra -o cortex_bin_check cortex_bin_check.c

clean:
	rm -rf cortex_bin_check
