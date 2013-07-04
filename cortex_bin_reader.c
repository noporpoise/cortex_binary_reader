#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <sys/stat.h>
#include <string.h>
#include <inttypes.h>
#include <errno.h>
#include <math.h>
#include <ctype.h> // toupper

#include "stream_buffer.h"

// Set buffer to 1MB
#define BUFFER_SIZE (1<<20)

#define MIN2(x,y) ((x) <= (y) ? (x) : (y))
#define MAX2(x,y) ((x) >= (y) ? (x) : (y))


// Calculates log2 of number since log2 is only available in some libc versions
#define Log2(n) (log(n) / log(2))

typedef struct
{
  char tip_cleaning;
  char remove_low_covg_supernodes;
  char remove_low_covg_kmers;
  char cleaned_against_graph;
  uint32_t remove_low_covg_supernodes_thresh;
  uint32_t remove_low_covg_kmers_thresh;
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

//
// What should we do
//
// Are we printing kmers?
char print_info = 1;
char print_kmers = 0;
char parse_kmers = 1;

buffer_t *buffer;

//
// File data
//
uint32_t version;
uint32_t kmer_size;
uint32_t num_of_bitfields;
uint32_t num_of_colours;

// version 7
uint64_t expected_num_of_kmers;
uint32_t num_of_shades, shade_bytes;

// version 6 only below here
char **sample_names = NULL;
long double *seq_error_rates = NULL;
CleaningInfo *cleaning_infos = NULL;

//
// Data about file contents
//
off_t file_size;
size_t num_bytes_read = 0;

// Does this file pass all tests?
uint32_t num_errors = 0, num_warnings = 0;

// Reading stats
unsigned long num_of_kmers_read = 0;
unsigned long sum_of_covgs_read = 0;
unsigned long sum_of_seq_loaded = 0;

// Checks
unsigned long num_of_all_zero_kmers = 0;
unsigned long num_of_oversized_kmers = 0;
unsigned long num_of_zero_covg_kmers = 0;

static void report_warning(const char* fmt, ...)
{
  num_warnings++;

  va_list argptr;
  va_start(argptr, fmt);
  fprintf(stderr, "Warning: ");
  vfprintf(stderr, fmt, argptr);
  va_end(argptr);
}

static void report_error(const char* fmt, ...)
{
  num_errors++;

  va_list argptr;
  va_start(argptr, fmt);
  fprintf(stderr, "Error: ");
  vfprintf(stderr, fmt, argptr);
  va_end(argptr);
}

static unsigned long round_up_ulong(unsigned long num, unsigned long nearest)
{
  return nearest * ((num + nearest - 1) / nearest);
}

static unsigned int num_of_digits(unsigned long num)
{
  unsigned int digits;

  for(digits = 1; num >= 10; digits++)
    num /= 10;

  return digits;
}

// result must be long enough for result + 1 ('\0'). Max length required is:
// strlen('18,446,744,073,709,551,615')+1 = 27
// returns pointer to result
static char* ulong_to_str(unsigned long num, char* result)
{
  int digits = num_of_digits(num);
  int num_commas = (digits-1) / 3;

  int i;
  char *p = result + digits + num_commas;

  *p = '\0';
  p--;

  for(i = 0; i < digits; i++)
  {
    if(i > 0 && i % 3 == 0)
    {
      *p = ',';
      p--;
    }

    *p = '0' + (num % 10);
    p--;
    num /= 10;
  }

  return result;
}

// result must be long enough for result + 1 ('\0').
// Max length required is: 26+1+decimals+1 = 28+decimals bytes
//   strlen('-9,223,372,036,854,775,808') = 27
//   strlen('.') = 1
//   +1 for \0
static char* double_to_str(double num, int decimals, char* str)
{
  unsigned long whole_units = (unsigned long)num;
  num -= whole_units;

  ulong_to_str(whole_units, str);

  if(decimals > 0)
  {
    // Horrible hack to save character being overwritten with a leading zero
    // e.g. 12.121 written as '12' then '0.121', giving '10.121', put back '2'
    // '12.121'
    size_t offset = strlen(str);
    char c = str[offset-1];
    sprintf(str+offset-1, "%.*lf", decimals, num);
    str[offset-1] = c;
  }

  return str;
}

// str must be 26 + 3 + 1 + num decimals + 1 = 31+decimals bytes
// breakdown:
//   strlen('18,446,744,073,709,551,615') = 26
//   strlen(' GB') = 3
//   strlen('.') = 1
//   +1 for '\0'
static char* bytes_to_str(unsigned long num, int decimals, char* str)
{
  const unsigned int num_unit_sizes = 7;
  char *units[] = {"B", "KB", "MB", "GB", "TB", "PB", "EB"};

  unsigned long unit;
  unsigned long num_cpy = num;

  for(unit = 0; num_cpy >= 1024 && unit < num_unit_sizes; unit++)
    num_cpy /= 1024;

  unsigned long bytes_in_unit = 0x1UL << (10 * unit);
  double num_of_units = (double)num / bytes_in_unit;

  double_to_str(num_of_units, decimals, str);
  size_t offset = strlen(str);
  strcpy(str+offset, units[unit]);

  return str;
}


// str must be at least 32 bytes long
// max lenth: strlen '18,446,744,073,709,551,615.0 GB' + 1 = 32 bytes
static void set_memory_required_str(unsigned long num_of_hash_entries, char* str)
{
  // Size of each entry is rounded up to nearest 8 bytes
  unsigned long num_of_bytes
    = num_of_hash_entries *
      round_up_ulong(8*num_of_bitfields + 5*num_of_colours + 1, 8);

  bytes_to_str(num_of_bytes, 1, str);
}

// Returns -1 on failure
static off_t get_file_size(char* filepath)
{
  struct stat st;

  if (stat(filepath, &st) == 0)
      return st.st_size;

  fprintf(stderr, "Error: Cannot determine size of %s: %s\n",
          filepath, strerror(errno));

  return -1;
}

static void print_kmer_stats()
{
  char num_str[50];

  if(num_of_all_zero_kmers > 1)
  {
    report_error("%s all-zero-kmers seen\n",
                 ulong_to_str(num_of_all_zero_kmers, num_str));
  }

  if(num_of_oversized_kmers > 0)
  {
    report_error("%s oversized kmers seen\n",
                 ulong_to_str(num_of_oversized_kmers, num_str));
  }

  if(num_of_zero_covg_kmers > 0)
  {
    report_error("%s kmers have no coverage in any colour\n",
                 ulong_to_str(num_of_zero_covg_kmers, num_str));
  }

  if((print_kmers || parse_kmers) && print_info)
  {
    printf("kmers read: %s\n", ulong_to_str(num_of_kmers_read, num_str));
    printf("covgs read: %s\n", ulong_to_str(sum_of_covgs_read, num_str));
    printf("seq loaded: %s\n", ulong_to_str(sum_of_seq_loaded, num_str));
  }

  if(print_info)
  {
    // Memory calculations
    // use expected number of kmers if we haven't read the whole file
    unsigned long kmer_count
      = (print_kmers || parse_kmers ? num_of_kmers_read : expected_num_of_kmers);

    // Number of hash table entries is 2^mem_height * mem_width
    // Aim for 80% occupancy once loaded
    float extra_space = 10.0/8;
    unsigned long hash_capacity = extra_space * kmer_count;

    // mem_width must be within these boundaries
    unsigned int min_mem_width = 5;
    unsigned int max_mem_width = 50;
    unsigned int min_mem_height = 12;
    // min mem usage = 2^12 * 5 = 20,480 entries = 320.0 KB with k=31,cols=1

    unsigned long mem_height = min_mem_height;
    unsigned long mem_width = max_mem_width;
    unsigned long hash_entries = (0x1UL << mem_height) * mem_width;

    if(hash_capacity > hash_entries)
    {
      // Resize
      mem_height = Log2((double)hash_capacity / (max_mem_width-1))+0.99;
      mem_height = MIN2(mem_height, 32);
      mem_height = MAX2(mem_height, min_mem_height);

      mem_width = hash_capacity / (0x1UL << mem_height) + 1;

      printf("mem_width: %lu; mem_height: %lu;\n", mem_width, mem_height);

      if(mem_width < min_mem_width)
      {
        // re-calculate mem_height
        mem_height = Log2((double)hash_capacity / min_mem_width)+0.99;
        mem_height = MIN2(mem_height, 32);
        mem_height = MAX2(mem_height, min_mem_height);
        mem_width = hash_capacity / (0x1UL << mem_height) + 1;
        mem_width = MAX2(mem_width, min_mem_width);
      }

      hash_entries = (0x1UL << mem_height) * mem_width;
    }

    char min_mem_required[32];
    char rec_mem_required[32];

    set_memory_required_str(kmer_count, min_mem_required);
    set_memory_required_str(hash_entries, rec_mem_required);

    printf("Memory required: %s\n", min_mem_required);
    printf("Memory suggested: --mem_width %lu --mem_height %lu\n",
           mem_width, mem_height);

    char hash_entries_numstr[32];
    ulong_to_str(hash_entries, hash_entries_numstr);

    printf("  [%s entries; %s memory]\n", hash_entries_numstr, rec_mem_required);
  }
}

static void my_fread(FILE *fh, void *ptr, int size, const char* entry_name)
{
  int read = fread_buf(fh, ptr, size, buffer);

  if(read != size)
  {
    report_error("Couldn't read '%s': expected %li; recieved: %li; (fatal)\n",
                 entry_name, (long)size, (long)read);

    if(print_kmers)
      printf("----\n");

    print_kmer_stats();
    exit(EXIT_FAILURE);
  }

  num_bytes_read += read;
}

static void print_binary(FILE* fh, uint64_t binary)
{
  int i;
  for(i = 63; i >= 0; i--)
    fprintf(fh, "%c", ((binary >> i) & 0x1 ? '1' : '0'));
}

static char binary_nucleotide_to_char(Nucleotide n)
{
  switch (n)
  {
    case Adenine:  return 'A';
    case Cytosine: return 'C';
    case Guanine:  return 'G';
    case Thymine:  return 'T';
    default:
      fprintf(stderr, "Non existent binary nucleotide %d\n", n);
      exit(EXIT_FAILURE);
  }
}

static void binary_kmer_right_shift_one_base(uint64_t* kmer, int num_of_bitfields)
{
  int i;
  for(i = num_of_bitfields-1; i > 0; i--)
  {
    kmer[i] >>= 2;
    kmer[i] |= (kmer[i-1] << 62); // & 0x3
  }

  kmer[0] >>= 2;
}


#define rev_nibble(x) (((x&0x1)<<3) | ((x&0x2)<<1) | ((x&0x4)>>1) | ((x&0x8)>>3))

static char* get_edges_str(char edges, char* kmer_colour_edge_str)
{
  int i;

  char str[] = "acgt";

  char left = edges >> 4;
  left = rev_nibble(left);
  char right = edges & 0xf;

  for(i = 0; i < 4; i++)
    kmer_colour_edge_str[i] = (left & (0x1 << i) ? str[i] : '.');

  for(i = 0; i < 4; i++)
    kmer_colour_edge_str[i+4] = toupper(right & (0x1 << i) ? str[i] : '.');

  kmer_colour_edge_str[8] = '\0';

  return kmer_colour_edge_str;
}

static char* binary_kmer_to_seq(uint64_t* bkmer, char * seq,
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

#define has_shade(p,n)   (((p)[(n) >> 3] >> ((n) & 0x7)) & 0x1)

static char get_shade_char(uint8_t *shades, uint8_t *shends, int p)
{
  char shend = has_shade(shends,p);
  char shade = has_shade(shades,p);
  if(shade && shend) return '-';
  else if(shend) return 'A'+(p % 26);
  else if(shade) return 'a'+(p % 26);
  else return '.';
}

static void print_colour_shades(uint8_t *shades, uint8_t *shends)
{
  size_t i;
  for(i = 0; i < num_of_shades; i++)
    putc(get_shade_char(shades, shends, i), stdout);
}

static void print_usage()
{
  fprintf(stderr,
"usage: cortex_bin_reader [OPTIONS] <binary.ctx>\n"
"  Prints out header information and kmers for cortex_var binary files.  Runs\n"
"  several checks to test if binary file is valid. \n"
"\n"
"  OPTIONS:\n"
"  --print_info    Print header info and exit. If used on its own kmers are not\n"
"                  printed or checked (fast option).\n"
"\n"
"  --print_kmers   Print each kmer. If used on its own, other information\n"
"                  (i.e. headers) is not printed out\n"
"\n"
"  --parse_kmers   Print header info, parse but don't print kmers [default]\n"
"\n"
"  If no options are specified '--parse_kmers --print_info' is used.\n"
"\n"
"  Kmers are printed in the order they are listed in the file. \n"
"  For each kmer we print: <kmer_seq> <covg_in_col0 ...> <edges_in_col0 ...>\n"
"    e.g. GTAAGTGCCA 6 4 ..g....T .c..A..T\n"
"         means col 0: covg 6 [G]GTAAGTGCCA[T]\n"
"               col 1: covg 4 [C]GTAAGTGCCA[A|T]\n"
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
      {
        print_info = 1;
      }
      else if(strcasecmp(argv[i], "--print_kmers") == 0)
      {
        print_kmers = 1;
      }
      else if(strcasecmp(argv[i], "--parse_kmers") == 0)
      {
        print_info = 1;
        parse_kmers = 1;
      }
      else
        print_usage();
    }
  }

  filepath = argv[argc-1];

  if(print_info)
    printf("Loading file: %s\n", filepath);

  file_size = get_file_size(filepath);

  FILE* fh = fopen(filepath, "r");

  if(fh == NULL)
  {
    report_error("cannot open file '%s'\n", filepath);
    exit(EXIT_FAILURE);
  }

  if(file_size != -1 && print_info)
  {
    char str[31];
    bytes_to_str(file_size, 0, str);
    printf("File size: %s\n", str);
  }

  buffer = buffer_new(BUFFER_SIZE);

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

  unsigned int i;

  // Read magic word at the start of header
  char magic_word[7];
  magic_word[6] = '\0';

  my_fread(fh, magic_word, strlen("CORTEX"), "Magic word");

  if(strcmp(magic_word, "CORTEX") != 0)
  {
    fprintf(stderr, "Magic word doesn't match 'CORTEX' (start)\n");
    exit(EXIT_FAILURE);
  }

  // Read version number
  my_fread(fh, &version, sizeof(uint32_t), "binary version");
  my_fread(fh, &kmer_size, sizeof(uint32_t), "kmer size");
  my_fread(fh, &num_of_bitfields, sizeof(uint32_t), "number of bitfields");
  my_fread(fh, &num_of_colours, sizeof(uint32_t), "number of colours");

  if(print_info)
  {
    printf("binary version: %i\n", (int)version);
    printf("kmer size: %i\n", (int)kmer_size);
    printf("bitfields: %i\n", (int)num_of_bitfields);
    printf("colours: %i\n", (int)num_of_colours);
  }

  if(version >= 7)
  {
    my_fread(fh, &expected_num_of_kmers, sizeof(uint64_t), "number of kmers");
    my_fread(fh, &num_of_shades, sizeof(uint32_t), "number of shades");

    if(print_info)
    {
      char tmp[256];
      printf("kmers: %s\n", ulong_to_str(expected_num_of_kmers,tmp));
      printf("shades: %i\n", (int)num_of_shades);
    }
  }

  // Checks

  if(version > 7 || version < 4)
    report_error("Sorry, we only support binary versions 4, 5, 6 & 7\n");

  if(kmer_size % 2 == 0)
    report_error("kmer size is not an odd number\n");

  if(kmer_size < 3)
    report_error("kmer size is less than three\n");

  if(num_of_bitfields * 32 < kmer_size)
    report_error("Not enough bitfields for kmer size\n");

  if((num_of_bitfields-1)*32 >= kmer_size)
    report_error("using more than the minimum number of bitfields\n");

  if(num_of_colours == 0)
    report_error("number of colours is zero\n");

  if(num_of_shades != 0 && (num_of_shades & (num_of_shades-1)))
    report_error("number of shades is not a power of 2\n");

  //

  // Read array of mean read lengths per colour
  uint32_t *mean_read_lens_per_colour
    = (uint32_t*)malloc(num_of_colours*sizeof(uint32_t));

  my_fread(fh, mean_read_lens_per_colour, sizeof(uint32_t) * num_of_colours,
           "mean read length for each colour");

  // Read array of total seq loaded per colour
  uint64_t *total_seq_loaded_per_colour
    = (uint64_t*)malloc(num_of_colours*sizeof(uint64_t));

  my_fread(fh, total_seq_loaded_per_colour, sizeof(uint64_t) * num_of_colours,
           "total sequance loaded for each colour");

  for(i = 0; i < num_of_colours; i++)
  {
    sum_of_seq_loaded += total_seq_loaded_per_colour[i];
  }

  if(version >= 6)
  {
    sample_names = (char**)malloc(sizeof(char*) * num_of_colours);

    for(i = 0; i < num_of_colours; i++)
    {
      uint32_t str_length;
      my_fread(fh, &str_length, sizeof(uint32_t), "sample name length");

      if(str_length == 0)
      {
        sample_names[i] = NULL;
      }
      else
      {
        sample_names[i] = (char*)malloc((str_length+1) * sizeof(char));
        my_fread(fh, sample_names[i], str_length, "sample name");
        sample_names[i][str_length] = '\0';

        // Check sample length is as long as we were told
        size_t sample_name_len = strlen(sample_names[i]);

        if(sample_name_len != str_length)
        {
          // Premature \0 in string
          report_warning("Sample %i name has length %lu but is only %lu chars "
                         "long (premature '\\0')\n",
                         i, str_length, sample_name_len);
        }
      }
    }

    seq_error_rates = malloc(sizeof(long double) * num_of_colours);
    my_fread(fh, seq_error_rates, sizeof(long double) * num_of_colours,
             "seq error rates");

    cleaning_infos = (CleaningInfo*)malloc(sizeof(CleaningInfo) * num_of_colours);

    for(i = 0; i < num_of_colours; i++)
    {
      my_fread(fh, &(cleaning_infos[i].tip_cleaning), 1, "tip cleaning");
      my_fread(fh, &(cleaning_infos[i].remove_low_covg_supernodes), 1,
               "remove low covg supernodes");
      my_fread(fh, &(cleaning_infos[i].remove_low_covg_kmers), 1,
               "remove low covg kmers");
      my_fread(fh, &(cleaning_infos[i].cleaned_against_graph), 1,
               "cleaned against graph");

      my_fread(fh, &(cleaning_infos[i].remove_low_covg_supernodes_thresh),
               sizeof(uint32_t), "remove low covg supernode threshold");
    
      my_fread(fh, &(cleaning_infos[i].remove_low_covg_kmers_thresh),
               sizeof(uint32_t), "remove low covg kmer threshold");

      if(!cleaning_infos[i].remove_low_covg_supernodes &&
         cleaning_infos[i].remove_low_covg_supernodes_thresh > 0)
      {
        report_warning("Binary header gives sample %i a cleaning threshold of "
                       "%u for supernodes when no cleaning was performed\n",
                       i, cleaning_infos[i].remove_low_covg_supernodes_thresh);
      }

      if(!cleaning_infos[i].remove_low_covg_kmers &&
         cleaning_infos[i].remove_low_covg_kmers_thresh > 0)
      {
        report_warning("Binary header gives sample %i a cleaning threshold of "
                       "%u for kmers when no cleaning was performed\n",
                       i, cleaning_infos[i].remove_low_covg_kmers_thresh);
      }

      uint32_t name_length;
      my_fread(fh, &name_length, sizeof(uint32_t), "graph name length");

      if(name_length == 0)
      {
        cleaning_infos[i].name_of_graph_clean_against = NULL;
      }
      else
      {
        cleaning_infos[i].name_of_graph_clean_against
          = (char*)malloc((name_length + 1) * sizeof(char));

        my_fread(fh, cleaning_infos[i].name_of_graph_clean_against,
                 name_length, "graph name length");

        cleaning_infos[i].name_of_graph_clean_against[name_length] = '\0';
      
        // Check sample length is as long as we were told
        size_t cleaned_name_len
          = strlen(cleaning_infos[i].name_of_graph_clean_against);

        if(cleaned_name_len != name_length)
        {
          // Premature \0 in string
          report_warning("Sample [%i] cleaned-against-name has length %u but is "
                         "only %u chars long (premature '\\0')\n",
                         i, name_length, cleaned_name_len);
        }
      }
    }
  }

  // Print colour info

  if(print_info)
  {
    for(i = 0; i < num_of_colours; i++)
    {
      printf("-- Colour %i --\n", i);

      if(version >= 6)
      {
        // Version 6 only output
        printf("  sample name: '%s'\n", sample_names[i]);
      }

      char tmp[32];

      printf("  mean read length: %u\n",
             (unsigned int)mean_read_lens_per_colour[i]);
      printf("  total sequence loaded: %s\n",
             ulong_to_str(total_seq_loaded_per_colour[i], tmp));
      
      if(version >= 6)
      {
        // Version 6 only output
        printf("  sequence error rate: %Lf\n", seq_error_rates[i]);

        printf("  tip clipping: %s\n",
               (cleaning_infos[i].tip_cleaning == 0 ? "no" : "yes"));

        printf("  remove low coverage supernodes: %s [threshold: %i]\n",
               cleaning_infos[i].remove_low_covg_supernodes ? "yes" : "no",
               cleaning_infos[i].remove_low_covg_supernodes_thresh);

        printf("  remove low coverage kmers: %s [threshold: %i]\n",
               cleaning_infos[i].remove_low_covg_kmers ? "yes" : "no",
               cleaning_infos[i].remove_low_covg_kmers_thresh);

        printf("  cleaned against graph: %s [against: '%s']\n",
               cleaning_infos[i].cleaned_against_graph ? "yes" : "no",
               (cleaning_infos[i].name_of_graph_clean_against == NULL
                  ? "" : cleaning_infos[i].name_of_graph_clean_against));
      }
    }

    printf("--\n");
  }

  // Read magic word at the end of header
  my_fread(fh, magic_word, strlen("CORTEX"), "magic word (end)");

  if(strcmp(magic_word, "CORTEX") != 0)
  {
    report_error("magic word doesn't match 'CORTEX' (end): '%s'\n", magic_word);
    exit(EXIT_FAILURE);
  }

  // Calculate number of kmers
  if(version < 7 && file_size != -1)
  {
    size_t bytes_remaining = file_size - num_bytes_read;
    size_t num_bytes_per_kmer = sizeof(uint64_t) * num_of_bitfields +
                                sizeof(uint32_t) * num_of_colours +
                                sizeof(uint8_t) * num_of_colours;

    expected_num_of_kmers = bytes_remaining / num_bytes_per_kmer;

    size_t excess = bytes_remaining - (expected_num_of_kmers * num_bytes_per_kmer);

    if(excess > 0)
    {
      report_error("Excess bytes. Bytes:\n  file size: %lu;\n  for kmers: %lu;"
                   "\n  num kmers: %lu;\n  per kmer: %lu;\n  excess: %lu\n",
                   file_size, bytes_remaining, expected_num_of_kmers,
                   num_bytes_per_kmer, excess);
    }
  }

  if(print_info)
  {
    char num_str[50];
    printf("Expected number of kmers: %s\n",
           ulong_to_str(expected_num_of_kmers, num_str));
    printf("----\n");
  }

  // Finished parsing header
  if(!parse_kmers && !print_kmers)
  {
    print_kmer_stats();
    fclose(fh);
    exit(EXIT_SUCCESS);
  }


  shade_bytes = num_of_shades >> 3;
  size_t shade_array_bytes = shade_bytes * num_of_colours;

  // Kmer data
  uint64_t* kmer = (uint64_t*)malloc(sizeof(uint64_t) * num_of_bitfields);
  uint32_t* covgs = (uint32_t*)malloc(sizeof(uint32_t) * num_of_colours);
  uint8_t* edges = (uint8_t*)malloc(sizeof(uint8_t) * kmer_size);
  uint8_t* shade_data = (uint8_t*)malloc(shade_array_bytes);
  uint8_t* shend_data = (uint8_t*)malloc(shade_array_bytes);

  // Convert values to strings
  char* seq = (char*)malloc(sizeof(char) * kmer_size);
  char kmer_colour_edge_str[9];

  // Check top word of each kmer
  int bits_in_top_word = 2 * (kmer_size % 32);
  uint64_t top_word_mask = (~(uint64_t)0) << bits_in_top_word;

  size_t num_bytes_per_bkmer = sizeof(uint64_t)*num_of_bitfields;

  // Read kmer in bytes so we can see if there are extra bytes at the end of
  // the file
  size_t bytes_read;

  while((bytes_read = fread_buf(fh, kmer, num_bytes_per_bkmer, buffer)) > 0)
  {
    if(bytes_read != num_bytes_per_bkmer)
    {
      report_error("unusual extra bytes [%i] at the end of the file\n",
                   (int)bytes_read);
      break;
    }
    num_bytes_read += bytes_read;

    my_fread(fh, covgs, sizeof(uint32_t) * num_of_colours, "kmer covg");
    my_fread(fh, edges, sizeof(uint8_t) * num_of_colours, "kmer edges");

    if(version >= 7)
    {
      uint8_t *shades = shade_data, *shends = shend_data;
      for(i = 0; i < num_of_colours; i++)
      {
        my_fread(fh, shades, sizeof(uint8_t) * shade_bytes, "shades");
        my_fread(fh, shends, sizeof(uint8_t) * shade_bytes, "shade ends");
        shades += shade_bytes;
        shends += shade_bytes;
      }
    }

    //
    // Kmer checks
    //

    // Check top bits of kmer
    if(kmer[0] & top_word_mask)
    {
      if(num_of_oversized_kmers == 0)
      {
        report_error("oversized kmer [index: %lu]\n", num_of_kmers_read);

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
      {
        report_error("more than one all 'A's kmers seen [index: %lu]\n",
                     num_of_kmers_read);
      }

      num_of_all_zero_kmers++;
    }

    // Check covg is 0 for all colours
    char kmer_has_covg = 0;

    for(i = 0; i < num_of_colours; i++)
    {
      if(covgs[i] > 0)
      {
        kmer_has_covg = 1;
        break;
      }
    }

    if(kmer_has_covg == 0)
    {
      if(num_of_zero_covg_kmers == 0)
      {
        report_error("a kmer has zero coverage in all colours [index: %lu]\n",
                     num_of_kmers_read);
      }

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

      if(version >= 7 && num_of_shades > 0)
      {
        for(i = 0; i < num_of_colours; i++)
        {
          putc(' ', stdout);
          print_colour_shades(shade_data + i*shade_bytes, shend_data + i*shade_bytes);
        }
      }

      putc('\n', stdout);
    }

    num_of_kmers_read++;

    for(i = 0; i < num_of_colours; i++)
      sum_of_covgs_read += covgs[i];
  }

  if(num_of_kmers_read != expected_num_of_kmers)
  {
    report_error("Expected %lu kmers, read %lu\n",
                 expected_num_of_kmers, num_of_kmers_read);
  }

  if(print_kmers && print_info)
    printf("----\n");

  // check for various reading errors
  if(errno != 0)
  {
    report_error("errno set [%i]\n", (int)errno);
  }

  int err;
  if((err = ferror(fh)) != 0)
  {
    report_error("occurred after file reading [%i]\n", err);
  }

  // For testing output
  //num_of_bitfields = 2;
  //num_of_kmers_read = 3600000000;
  //num_of_kmers_read = 12345;
  //num_of_kmers_read = 3581787;
  //num_of_kmers_read = 0;

  print_kmer_stats();

  fclose(fh);

  free(kmer);
  free(covgs);
  free(edges);
  free(shade_data);
  free(shend_data);

  buffer_free(buffer);

  if((print_kmers || parse_kmers) && print_info)
  {
    printf("----\n");
    if(num_warnings > 0 || num_errors > 0)
      printf("Warning: %u; Errors: %u", num_warnings, num_errors);
    if(num_errors == 0)
      printf(num_warnings ? "Binary may be ok\n" : "Binary is valid\n");
  }

  exit(EXIT_SUCCESS);
}
