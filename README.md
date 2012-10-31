cortex_binary_reader
====================
    
Reads cortex binary files and prints contents.  

Build
-----

Just type:

    make

and to run:

    ./cortex_bin_reader

To print header and exit

    cortex_bin_reader --print_info in.ctx

To print header and parse kmers (checks for corruption) -- this is the default
behaviour so these two are the same

    cortex_bin_reader --print_info --parse_kmers in.ctx
    cortex_bin_reader in.ctx

To only print kmers:

    cortex_bin_reader --print_kmers in.ctx

To sum coverage:

    cortex_bin_reader --print_kmers in.ctx | awk '{total += $2} END { print total}'

Usage
-----

    usage: cortex_bin_reader [OPTIONS] <binary.ctx>
      Prints out header information and kmers for cortex_var binary files.  Runs
      several checks to test if binary file is valid. 

      OPTIONS:
      --print_info    Print header info and exit. If used on its own kmers are not
                      printed or checked (fast option).

      --print_kmers   Print each kmer. If used on its own, other information
                      (i.e. headers) is not printed out

      --parse_kmers   Print header info, parse but don't print kmers [default]

      If no options are specified '--parse_kmers --print_info' is used.

      Kmers are printed in the order they are listed in the file. 
      For each kmer we print: <kmer_seq> <covg_in_col0 ...> <edges_in_col0 ...>
        e.g. GTAAGTGCCA 6 4 ..g....T .c..A..T
             means col 0: covg 6 [G]GTAAGTGCCA[T]
                   col 1: covg 4 [C]GTAAGTGCCA[A|T]

      Header checks:
        * binary version is 4, 5 or 6
        * Kmer size is an odd number > 1
        * number of bitfields is compatible with kmer size
        * number of colours is > 0
        * Strings are correct length (don't have premature \0)

      Kmer checks:
        * each kmer's top bits are all zeroed (i.e. kmer is not 'oversized')
        * no more than one kmer is all As i.e. no multiple 'AAAAAAAA' kmers
        * each kmer has coverage greater than zero in at least one colour

      Comments/bugs/requests: <turner.isaac@gmail.com>
