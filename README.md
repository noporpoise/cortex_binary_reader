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
      --print_info    Print out header info [default]
      --print_kmers   Print out each kmer
      --parse_kmers   Parse but don't print kmers [default]

      If no options are specified '--print_info --parse_kmers' is used

      Kmers are printed in the order they are listed in the file. 
      For each kmer we print: <kmer_seq> <covg_in_col0 ...> <edges_in_col0 ...>
        e.g. GTAAGTGCCA 6 4 ..g....T .c..A..T
             means col 0: covg 6 [G]GTAAGTGCCA[T]
                   col 1: covg 4 [C]GTAAGTGCCA[A|T]

      Header checks:
        * binary version is 6
        * Kmer size is an odd number > 1
        * number of bitfields is compatible with kmer size
        * number of colours is > 0

      Kmer checks:
        * each kmer's top bits are all zeroed (i.e. kmer is not 'oversized')
        * no more than one kmer is all As i.e. no multiple 'AAAAAAAA' kmers
        * each kmer has covg greater than zero in at least one colour

      Comments/bugs/requests: <turner.isaac@gmail.com>
