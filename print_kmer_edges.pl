#!/usr/bin/perl

use strict;
use warnings;

use FASTNFile;
use GeneticsModule;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "" .
"Usage: ./print_kmer_edges.pl <kmer> [file1 ..]
  Takes a kmer and a sequence file and print its edges.\n";

  exit;
}

if(@ARGV < 2)
{
  print_usage();
}

my $kmer = shift;
my @seq_files = @ARGV;

if(@ARGV == 0)
{
  # Read from stdin
  push(@ARGV, "-");
}

my ($seq) = read_all_from_files(@ARGV);

my $strand = "+";

for(my $i = 0; $i < 2; $i++)
{
  if($i == 1)
  {
    $kmer = rev_comp($kmer);
    $strand = "-";
  }

  while(my ($name, $read) = each(%$seq))
  {
    if($read =~ /(.?)$kmer(.?)/i)
    {
      my ($left,$right) = (defined($1) ? $1 : "",defined($2) ? $2 : "");

      if($strand eq "-")
      {
        ($left,$right) = (rev_comp($right),rev_comp($left));
      }

      print "$strand $left:$right\n";

      exit;
    }
  }
}

