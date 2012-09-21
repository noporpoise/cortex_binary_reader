cortex_binary_reader
====================
    
Reads cortex binary files and prints contents.  

Build
-----

Just type:

    make

and to run:

    ./cortex_bin_reader

Usage
-----

    usage: cortex_bin_reader [--print_kmers] <binary.ctx>
      Prints out header information and kmers for cortex_var binary files.  Runs
      several checks to test if binary file is valid. 
    
      --print_kmers    prints kmers in the order they are listed in the file.
                       For each kmer prints:
    
        <kmer_seq> <covg_in_col0 ...> <edges_in_col0 ...>
    
        e.g. GTAAGTGCCA 6 4 ..g....T .c..A..T
             meaning:
                col 0: covg 6 [G]GTAAGTGCCA[T]
                col 1: covg 4 [C]GTAAGTGCCA[A|T]
    
      Current Tests:
        * Checks binary version is 6
        * Checks kmer size is an odd number > 1
        * Checks number of bitfields is compatible with kmer size
        * Checks number of colours is > 0
        * Checks each kmer's top bits are all zeroed (i.e. kmer is not 'oversized')
        * Checks if more than one kmer is all As i.e. multiple 'AAAAAAAA' kmers

      Comments/bugs/requests: <turner.isaac@gmail.com>
