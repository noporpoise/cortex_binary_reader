#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <errno.h>

typedef struct
{
  char tip_cleaning;
  char remove_low_covg_supernodes;
  char remove_low_covg_kmers;
  char cleaned_against_graph;
  int remove_low_covg_supernodes_thresh;
  int remove_low_covg_kmer_thresh;
  char* name_of_graph_clean_against;
} CleaningInfo;

typedef enum
{
  Adenine   = 0,
  Cytosine  = 1,
  Guanine   = 2,
  Thymine   = 3,
  Undefined = 4,
} Nucleotide;

// Are we printing kmers?
char print_info = 1;
char print_kmers = 0;
char parse_kmers = 1;

// Reading stats
int num_of_kmers_read = 0;

// Checks
unsigned long num_of_all_zero_kmers = 0;
unsigned long num_of_oversized_kmers = 0;
unsigned long num_of_zero_covg_kmers = 0;

void print_kmer_stats()
{
  if(num_of_all_zero_kmers > 1)
    fprintf(stderr, "Error: %lu all-zero-kmers seen\n", num_of_all_zero_kmers);

  if(num_of_oversized_kmers > 0)
    fprintf(stderr, "Error: %lu oversized kmers seen\n", num_of_oversized_kmers);

  if(num_of_zero_covg_kmers > 0)
  {
    fprintf(stderr, "Error: %lu kmers have no coverage is any colour\n",
            num_of_zero_covg_kmers);
  }

  if(print_info || parse_kmers)
    printf("kmers read: %lu\n", (unsigned long)num_of_kmers_read);
}

void my_fread(void *ptr, size_t size, size_t nitems, FILE *stream,
              char* entry_name)
{
  size_t read;
  if((read = fread(ptr, size, nitems, stream)) != nitems)
  {
    fprintf(stderr, "Couldn't read '%s': expected %li; recieved: %li; (exiting...)\n",
            entry_name, (long)nitems, (long)read);

    if(print_kmers)
      printf("----\n");

    print_kmer_stats();
    exit(EXIT_FAILURE);
  }

  int err;
  if((err = ferror(stream)) != 0)
  {
    fprintf(stderr, "An error occurred on file reading: %i\n", err);
  }
}

void print_binary(FILE* stream, uint64_t binary)
{
  int i;
  for(i = 63; i >= 0; i--)
    fprintf(stream, "%c", ((binary >> i) & 0x1 ? '1' : '0'));
}

char binary_nucleotide_to_char(Nucleotide n)
{
  switch (n)
  {
    case Adenine:
      return 'A';
    case Cytosine:
      return 'C';
    case Guanine:
      return 'G';
    case Thymine:
      return 'T';
    default:
      fprintf(stderr, "Non existent binary nucleotide %d\n", n);
      exit(EXIT_FAILURE);
  }
}

char char_rev_comp(char c)
{
  switch(c)
  {
    case 'T':
      return 'A';
    case 'G':
      return 'C';
    case 'C':
      return 'G';
    case 'A':
      return 'T';
    case 't':
      return 'a';
    case 'g':
      return 'c';
    case 'c':
      return 'g';
    case 'a':
      return 't';
    default:
      fprintf(stderr, "Non existent char nucleotide %c\n", c);
      exit(EXIT_FAILURE);
  }
}

void binary_kmer_right_shift_one_base(uint64_t* kmer, int num_of_bitfields)
{
  int i;
  for(i = num_of_bitfields-1; i > 0; i--)
  {
    kmer[i] >>= 2;
    kmer[i] |= (kmer[i-1] << 62); // & 0x3
  }

  kmer[0] >>= 2;
}

char* get_edges_str(char edges, char* kmer_colour_edge_str)
{
  int i;

  char *str = "acgtACGT";

  char left = edges >> 4;
  char right = edges & 0xf;

  for(i = 0; i < 4; i++)
    kmer_colour_edge_str[i] = (left & (0x1 << (3-i)) ? str[i] : '.');

  for(i = 0; i < 4; i++)
    kmer_colour_edge_str[i+4] = (right & (0x1 << i) ? str[i+4] : '.');

  kmer_colour_edge_str[8] = '\0';

  return kmer_colour_edge_str;
}

char* binary_kmer_to_seq(uint64_t* bkmer, char * seq,
                         int kmer_size, int num_of_bitfields)
{
  uint64_t local_bkmer[num_of_bitfields];

  int i;

  // Copy over a word at a time
  for(i = 0; i < num_of_bitfields; i++)
  {
    local_bkmer[i] = bkmer[i];
  }

  // Loop backwards over bases
  for(i = kmer_size-1; i >= 0; i--)
  {
    seq[i] = binary_nucleotide_to_char(local_bkmer[num_of_bitfields-1] & 0x3);
    binary_kmer_right_shift_one_base(local_bkmer, num_of_bitfields);
  }

  seq[kmer_size] = '\0';

  return seq;
}

void print_usage()
{
  fprintf(stderr,
"usage: cortex_bin_reader [OPTIONS] <binary.ctx>\n"
"  Prints out header information and kmers for cortex_var binary files.  Runs\n"
"  several checks to test if binary file is valid. \n"
"\n"
"  OPTIONS:\n"
"  --print_info    Print out header info [default]\n"
"  --print_kmers   Print out each kmer\n"
"  --parse_kmers   Parse but don't print kmers [default]\n"
"\n"
"  If no options are specified '--print_info --parse_kmers' is used\n"
"\n"
"  Kmers are printed in the order they are listed in the file. \n"
"  For each kmer we print: <kmer_seq> <covg_in_col0 ...> <edges_in_col0 ...>\n"
"    e.g. GTAAGTGCCA 6 4 ..g....T .c..A..T\n"
"         means col 0: covg 6 [G]GTAAGTGCCA[T]\n"
"               col 1: covg 4 [C]GTAAGTGCCA[A|T]\n"
"\n"
"  Header checks:\n"
"    * binary version is 6\n"
"    * Kmer size is an odd number > 1\n"
"    * number of bitfields is compatible with kmer size\n"
"    * number of colours is > 0\n"
"\n"
"  Kmer checks:\n"
"    * each kmer's top bits are all zeroed (i.e. kmer is not 'oversized')\n"
"    * no more than one kmer is all As i.e. no multiple 'AAAAAAAA' kmers\n"
"    * each kmer has covg greater than zero in at least one colour\n"
"\n"
"  Comments/bugs/requests: <turner.isaac@gmail.com>\n");

  exit(EXIT_FAILURE);
}

int main(int argc, char** argv)
{
  char* filepath;

  if(argc < 2)
  {
    print_usage();
  }
  else if(argc > 2)
  {
    print_info = 0;
    print_kmers = 0;
    parse_kmers = 0;

    int i;

    for(i = 1; i < argc-1; i++)
    {
      if(strcasecmp(argv[i], "--print_info") == 0)
        print_info = 1;
      else if(strcasecmp(argv[i], "--print_kmers") == 0)
        print_kmers = 1;
      else if(strcasecmp(argv[i], "--parse_kmers") == 0)
        parse_kmers = 1;
      else
        print_usage();
    }
  }

  filepath = argv[argc-1];

  if(print_info)
    printf("Loading file: %s\n", filepath);

  FILE* fh = fopen(filepath, "r");

  if(fh == NULL)
  {
    fprintf(stderr, "Error: cannot open file '%s'\n", filepath);
    exit(EXIT_FAILURE);
  }

  /*
  // Check sizes
  printf("-- Datatypes --\n");
  printf("int: %i\n", (int)sizeof(int));
  printf("long: %i\n", (int)sizeof(long));
  printf("long long: %i\n", (int)sizeof(long long));
  printf("double: %i\n", (int)sizeof(double));
  printf("long double: %i\n", (int)sizeof(long double));
  */

  if(print_info)
    printf("----\n");

  int i;

  // Read magic word at the start of header
  char magic_word[7];
  magic_word[6] = '\0';

  my_fread(magic_word, sizeof(char), 6, fh, "Magic word");

  if(strcmp(magic_word, "CORTEX") != 0)
  {
    fprintf(stderr, "Magic word doesn't match 'CORTEX' (start)\n");
    exit(EXIT_FAILURE);
  }

  // Read version number
  int version;
  my_fread(&version, sizeof(int), 1, fh, "binary version");

  int kmer_size;
  my_fread(&kmer_size, sizeof(int), 1, fh, "kmer size");

  int num_of_bitfields;
  my_fread(&num_of_bitfields, sizeof(int), 1, fh, "number of bitfields");

  int num_of_colours;
  my_fread(&num_of_colours, sizeof(int), 1, fh, "number of colours");

  if(print_info)
  {
    printf("binary version: %i\n", version);
    printf("kmer size: %i\n", kmer_size);
    printf("bitfields: %i\n", num_of_bitfields);
    printf("colours: %i\n", num_of_colours);
  }

  // Checks

  if(version != 6)
    fprintf(stderr, "Error: binary version is not '6'\n");

  if(kmer_size % 2 == 0)
    fprintf(stderr, "Error: kmer size is not an odd number\n");

  if(kmer_size < 3)
    fprintf(stderr, "Error: kmer size is less than three\n");

  if(num_of_bitfields * 32 < kmer_size)
    fprintf(stderr, "Error: Not enough bitfields for kmer size\n");

  if((num_of_bitfields-1)*32 > kmer_size)
    fprintf(stderr, "Error: using more than the minimum number of bitfields\n");

  if(num_of_colours == 0)
    fprintf(stderr, "Error: number of colours is zero\n");

  //

  int *mean_read_lens_per_colour = (int*)malloc(num_of_colours*sizeof(int));
  my_fread(mean_read_lens_per_colour, sizeof(int), num_of_colours, fh,
           "mean read length for each colour");

  long *total_seq_loaded_per_colour = (long*)malloc(num_of_colours*sizeof(long));
  my_fread(total_seq_loaded_per_colour, sizeof(long), num_of_colours, fh,
           "total sequance loaded for each colour");


  char **sample_names = (char**)malloc(sizeof(char*) * num_of_colours);

  for(i = 0; i < num_of_colours; i++)
  {
    int str_length;
    my_fread(&str_length, sizeof(int), 1, fh, "sample name length");

    if(str_length == 0)
    {
      sample_names[i] = NULL;
    }
    else
    {
      sample_names[i] = (char*)malloc(str_length * sizeof(char)+1);
      my_fread(sample_names[i], sizeof(char), str_length, fh, "sample name");
      sample_names[i][str_length] = '\0';
    }
  }

  long double *seq_error_rates = malloc(sizeof(long double) * num_of_colours);
  my_fread(seq_error_rates, sizeof(long double), num_of_colours, fh,
           "seq error rates");

  CleaningInfo *cleaning_infos = (CleaningInfo*)malloc(sizeof(CleaningInfo) * num_of_colours);

  for(i = 0; i < num_of_colours; i++)
  {
    my_fread(&(cleaning_infos[i].tip_cleaning), sizeof(char), 1, fh,
             "tip cleaning");
    my_fread(&(cleaning_infos[i].remove_low_covg_supernodes), sizeof(char), 1, fh,
             "remove low covg supernodes");
    my_fread(&(cleaning_infos[i].remove_low_covg_kmers), sizeof(char), 1, fh,
             "remove low covg kmers");
    my_fread(&(cleaning_infos[i].cleaned_against_graph), sizeof(char), 1, fh,
             "cleaned against graph");

    my_fread(&(cleaning_infos[i].remove_low_covg_supernodes_thresh), sizeof(int),
             1, fh, "remove low covg supernode threshold");
  
    my_fread(&(cleaning_infos[i].remove_low_covg_kmer_thresh), sizeof(int),
             1, fh, "remove low covg kmer threshold");

    int name_length;
    my_fread(&name_length, sizeof(int), 1, fh, "graph name length");

    if(name_length == 0)
    {
      cleaning_infos[i].name_of_graph_clean_against = NULL;
    }
    else
    {
      cleaning_infos[i].name_of_graph_clean_against
        = malloc((name_length+1) * sizeof(char));

      my_fread(cleaning_infos + i, sizeof(char), name_length, fh,
               "graph name length");

      cleaning_infos[i].name_of_graph_clean_against[name_length] = '\0';
    }
  }

  // Print colour info

  if(print_info)
  {
    for(i = 0; i < num_of_colours; i++)
    {
      printf("-- Colour %i --\n", i);

      printf("  sample name: '%s'\n", sample_names[i]);
      printf("  mean read length: %i\n", mean_read_lens_per_colour[i]);
      printf("  total sequence loaded: %li\n", total_seq_loaded_per_colour[i]);
      printf("  sequence error rate: %Lf\n", seq_error_rates[i]);

      printf("  tip clipping: %s\n", cleaning_infos[i].tip_cleaning ? "yes" : "no");

      if(cleaning_infos[i].remove_low_covg_supernodes)
      {
        printf("  remove_low_coverage_supernodes: yes [threshold: %i]\n",
               cleaning_infos[i].remove_low_covg_supernodes_thresh);
      }
      else
      {
        printf("  remove_low_coverage_supernodes: no\n");
      }

      if(cleaning_infos[i].remove_low_covg_kmers)
      {
        printf("  remove_low_coverage_kmers: yes [threshold: %i]\n",
               cleaning_infos[i].remove_low_covg_kmer_thresh);
      }
      else
      {
        printf("  remove_low_coverage_kmers: no\n");
      }

      if(cleaning_infos[i].cleaned_against_graph)
      {
        printf("  cleaned against graph: yes [against: '%s']\n",
               cleaning_infos[i].name_of_graph_clean_against);
      }
      else
      {
        printf("  cleaned against graph: no\n");
      }
    }

    printf("----\n");
  }

  // Read magic word at the end of header
  my_fread(magic_word, sizeof(char), 6, fh, "magic word (end)");

  if(strcmp(magic_word, "CORTEX") != 0)
  {
    fprintf(stderr, "Error: magic word doesn't match 'CORTEX' (end): '%s'\n",
            magic_word);
    exit(EXIT_FAILURE);
  }

  // Finished parsing header
  if(!parse_kmers && !print_kmers)
  {
    fclose(fh);
    exit(EXIT_SUCCESS);
  }

  // Kmer data
  uint64_t* kmer = (uint64_t*)malloc(sizeof(uint64_t) * num_of_bitfields);
  char* seq = (char*)malloc(sizeof(char) * kmer_size);
  uint32_t* covgs = (uint32_t*)malloc(sizeof(uint32_t) * num_of_colours);
  char* edges = (char*)malloc(sizeof(char) * kmer_size);
  char kmer_colour_edge_str[9];

  int bits_in_top_word = 2 * (kmer_size % 32);
  uint64_t top_word_mask = (~(uint64_t)0) << bits_in_top_word;

  // Read kmer in bytes so we can see if there are extra bytes at the end of
  // the file
  size_t chars_read;

  while((chars_read = fread(kmer, sizeof(char), 8*num_of_bitfields, fh)) > 0)
  {
    if(chars_read != (size_t)8*num_of_bitfields)
    {
      fprintf(stderr, "Error: unusual extra bytes [%i] at the end of the file\n",
              (int)chars_read);
      break;
    }

    my_fread(covgs, sizeof(uint32_t), num_of_colours, fh, "kmer covg");
    my_fread(edges, sizeof(char), num_of_colours, fh, "kmer edges");

    //
    // Kmer checks
    //

    // Check top bits of kmer
    if(kmer[0] & top_word_mask)
    {
      if(num_of_oversized_kmers == 0)
      {
        fprintf(stderr, "Error: oversized kmer\n");

        for(i = 0; i < num_of_bitfields; i++)
        {
          fprintf(stderr, "  word %i: ", i);
          print_binary(stderr, kmer[i]);
          fprintf(stderr, "\n");
        }
      }

      num_of_oversized_kmers++;
    }

    // Check for all-zeros (i.e. all As kmer: AAAAAA)
    uint64_t kmer_words_or = 0;

    for(i = 0; i < num_of_bitfields; i++)
      kmer_words_or |= kmer[i];

    if(kmer_words_or == 0)
    {
      if(num_of_all_zero_kmers == 1)
        fprintf(stderr, "Error: more than one all 'A's kmers seen\n");

      num_of_all_zero_kmers++;
    }

    // Check covg is == 0 for all colours
    char kmer_has_covg = 0;

    for(i = 0; i < num_of_colours; i++)
      if(covgs[i] > 0)
        kmer_has_covg = 1;

    if(kmer_has_covg == 0)
    {
      if(num_of_zero_covg_kmers == 0)
        fprintf(stderr, "Error: kmer coverage is zero in all colours\n");

      num_of_zero_covg_kmers++;
    }

    // Print?
    if(print_kmers)
    {
      binary_kmer_to_seq(kmer, seq, kmer_size, num_of_bitfields);
      printf("%s", seq);

      // Print coverages
      for(i = 0; i < num_of_colours; i++)
        printf(" %li", (unsigned long)covgs[i]);

      // Print edges
      for(i = 0; i < num_of_colours; i++)
        printf(" %s", get_edges_str(edges[i], kmer_colour_edge_str));

      printf("\n");
    }

    num_of_kmers_read++;
  }

  if(print_kmers && print_info)
    printf("----\n");

  if(errno != 0)
    fprintf(stderr, "Error: errno set [%i]\n", (int)errno);

  int err;
  if((err = ferror(fh)) != 0)
    fprintf(stderr, "Error: occurred after file reading [%i]\n", err);

  print_kmer_stats();

  fclose(fh);

  exit(EXIT_SUCCESS);
}
