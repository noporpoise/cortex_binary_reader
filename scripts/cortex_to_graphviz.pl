#!/usr/bin/perl

use strict;
use warnings;

use File::Basename;
use IPC::Open2;

my %complement = ('a' => 't', 'A' => 'T',
                  'c' => 'g', 'C' => 'G',
                  'g' => 'c', 'G' => 'C',
                  't' => 'a', 'T' => 'A');

sub print_usage
{
  for my $err (@_)
  {
    print STDERR "Error: $err\n";
  }
  
  print STDERR "" .
"Usage: ./cortex_to_graphviz.pl [--point|--shaded] <in.ctx>
  Prints graphviz `dot' output.  Not to be used with large graphs!

  --point   Don't print kmer values, only points
  --shaded  Print shades

  Example: ./cortex_to_graphviz.pl small.ctx > small.dot
           dot -Tpng small.dot > small.png\n";

  exit(-1);
}

my $use_points = 0;
my $print_shades = 0;

while(@ARGV > 1) {
  if($ARGV[0] =~ /^-?-p(oints?)?$/i) {
    shift;
    $use_points = 1;
  }
  elsif($ARGV[0] =~ /^-?-s(hade[sd]?)?$/i) {
    shift;
    $print_shades = 1;
  }
  else { print_usage("Unknown option '$ARGV[0]'"); }
}

if(@ARGV != 1)
{
  print_usage();
}

my $file = shift;

if(!(-r $file))
{
  print STDERR "Error: Cannot read file: $file\n";
}

my $cmd = dirname(__FILE__)."/../cortex_bin_reader";

if(!(-e $cmd))
{
  print STDERR "Error: executable cortex_bin_reader doesn't exist -- " .
               "did you compile?\n";
}
elsif(!(-x $cmd))
{
  print STDERR "Error: cortex_bin_reader doesn't appear to be executable\n";
  exit(-1);
}

my $cmdline = "$cmd --print_kmers $file 2>&1";

my ($pid, $in, $out);

print "digraph G {\n";

my @fcols = qw(red green blue orange purple pink brown black);

if($print_shades)
{
  $pid = open2($in, $out, $cmdline) or die("Cannot run cmd: '$cmdline'");

  while(defined(my $line = <$in>))
  {
    my ($kmer, $covgs, $edges, $shades) = parse_ctx_line($line);
    if(defined($kmer))
    {
      my @flav = defined($shades) ? split('', $shades) : ('.');

      print $kmer . ' [shape=none label=<<table ' .
            'border="'.(defined($shades) && $shades =~ /^\-+$/ ? '1' : '0').'" '.
            'cellborder="0">
<tr><td PORT="'.$kmer.'" colspan="'.@flav.'" cellpadding="0" cellspacing="0">
<font face="courier" point-size="9">'.($use_points ? '.' : $kmer).'</font></td>
</tr><tr>';

      for(my $i = 0; $i < @flav; $i++) {
        print '<td fixedsize="true" width="3" height="3" ' .
              'cellpadding="0" cellspacing="0" border="1" ';
        if($flav[$i] ne '.') {
          if($flav[$i] eq '-') { print 'bgcolor="black"'; }
          elsif($flav[$i] eq lc($flav[$i])) { print 'style="rounded"'; }
          else { print 'style="rounded" bgcolor="'.$fcols[$i].'"'; }
          print ' color="'.$fcols[$i].'"';
        }
        else { print 'color="white"'; }
        print '></td>'."\n";
      }

      print "</tr></table>>];\n";
    }
  }

  close($in);
  close($out);
}
else
{
  print "  node [" . ($use_points ? "shape=point label=none" : "shape=none") ."]\n";
}

$pid = open2($in, $out, $cmdline) or die("Cannot run cmd: '$cmdline'");

print "  edge [dir=both arrowhead=none arrowtail=none]\n";

while(defined(my $line = <$in>))
{
  my ($kmer, $covgs, $edges, $shades) = parse_ctx_line($line);
  if(defined($kmer))
  {
    my $num_edges_printed = 0;

    for(my $i = 0; $i < 4; $i++)
    {
      if((my $edge = substr($edges, $i, 1)) ne ".")
      {
        my $prev_kmer = uc($edge) . substr($kmer,0,-1);
        my $right_base = substr($kmer,-1);
        dump_edge($prev_kmer, $right_base, 0);
        $num_edges_printed++;
      }
    }

    for(my $i = 4; $i < 8; $i++)
    {
      if((my $edge = substr($edges, $i, 1)) ne ".")
      {
        dump_edge($kmer, uc($edge), 1);
        $num_edges_printed++;
      }
    }

    if($num_edges_printed == 0)
    {
      print "  ".kmer_key($kmer)."\n";
    }
  }
}

print "}\n";

close($in);
close($out);

waitpid($pid, 1);

sub parse_ctx_line
{
  my ($line) = @_;

  chomp($line);
  if($line =~ /Error/i)
  {
    print STDERR "$line\n";
    return undef;
  }
  elsif($line =~ /^([acgt]+) (\d+ )+([acgt\.]{8})(?: ([a-z\-\.]+))?$/i)
  {
    # return: $kmer, $covgs, $edges, $shades
    return ($1, $2, $3, $4);
  }
  else
  {
    print STDERR "Cannot parse line:\n";
    print STDERR "  $line\n";
    exit(-1);
  }
}

sub kmer_key
{
  my ($kmer) = @_;
  my $kmer_revcmp = revcmp($kmer);
  return $kmer lt $kmer_revcmp ? $kmer : $kmer_revcmp;
}

sub revcmp
{
  my ($seq) = @_;
  for(my $i = 0; $i < length($seq); $i++)
  {
    my $b = substr($seq, $i, 1);
    substr($seq, $i, 1) = $complement{$b};
  }
  return reverse($seq);
}

sub dump_edge
{
  my ($kmer1, $rbase, $going_right) = @_;

  my $kmer2 = substr($kmer1, 1) . $rbase;

  my $key1 = kmer_key($kmer1);
  my $key2 = kmer_key($kmer2);

  my $rev1 = ($kmer1 ne $key1);
  my $rev2 = ($kmer2 ne $key2);

  # When doing right hand edges, do only those that go to left
  #                              or those to a greater key
  # When doing left hand edges, do only those that go to a greater key
  if(($going_right && $rev1 == $rev2) || ($rev1 != $rev2 && $key1 le $key2))
  {
    # Print a coloured edge for each colour that is in both nodes
    $key1 = $print_shades ? "$key1:$key1:" : "$key1:";
    $key2 = $print_shades ? "$key2:$key2:" : "$key2:";
    print "  $key1" . ($rev1 ? 'w' : 'e') . " -> " .
             $key2  . ($rev2 ? 'e' : 'w') . "\n";
  }
}
