#!/usr/bin/perl

use strict;
use warnings;

use File::Basename;
use IPC::Open2;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use CortexGraph;

sub print_usage
{
  for my $err (@_)
  {
    print STDERR "Error: $err\n";
  }
  
  print STDERR "" .
"Usage: ./cortex_to_graphviz.pl [--point|--shaded|--simplify] <in.ctx>
  Prints graphviz `dot' output.  Not to be used with large graphs!

  --point    Don't print kmer values, only points
  --shaded   Print shades
  --simplify Simplify supernodes

  Example: ./cortex_to_graphviz.pl small.ctx > small.dot
           dot -Tpng small.dot > small.png\n";

  exit(-1);
}

my $use_points = 0;
my $print_shades = 0;
my $simplify = 0;

while(@ARGV > 1) {
  if($ARGV[0] =~ /^-?-p(oints?)?$/i) {
    shift;
    $use_points = 1;
  }
  elsif($ARGV[0] =~ /^-?-shade[sd]??$/i) {
    shift;
    $print_shades = 1;
  }
  elsif($ARGV[0] =~ /^-?-simplify$/i) {
    shift;
    $simplify = 1;
  }
  else { print_usage("Unknown option '$ARGV[0]'"); }
}

if(@ARGV != 1) { print_usage(); }

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

# colours for shades
my @cols = qw(red green blue orange purple pink brown black);

# Print a warning only once if shades mismatch
my $shade_mismatch = 0;

# cortex_bin_reader command
my $cmdline = "$cmd --print_kmers $file 2>&1";
my ($pid, $in, $out);

print "digraph G {\n";
print "  edge [dir=both arrowhead=none arrowtail=none]\n";
print "  node [".($use_points ? "shape=point, label=none" : "shape=none").", ".
      "fontname=courier, fontsize=9]\n";

if($simplify)
{
  $pid = open2($in, $out, $cmdline) or die("Cannot run cmd: '$cmdline'");

  my $graph = new CortexGraph();

  # Construct graph
  while(defined(my $line = <$in>))
  {
    my ($kmer, $covgs, $edges, $shades) = parse_ctx_line($line);
    $graph->add_kmer($kmer);
    $graph->{$kmer}->{'shades'} = $shades;

    my @edges_arr = split('', uc($edges));

    for my $prev_edge (grep {$_ ne '.'} @edges_arr[0..3]) {
      $graph->add_edges_between(uc($prev_edge).substr($kmer,0,-1), $kmer);
    }

    for my $next_edge (grep {$_ ne '.'} @edges_arr[4..7]) {
      $graph->add_edges_between($kmer, substr($kmer,1).$next_edge);
    }
  }

  # $graph->dump();
  # exit;

  # Get kmer size
  my $kmer_size = $graph->get_kmer_size();

  # print "kmer size: $kmer_size\n";

  # Simplify graph into supernodes
  # Hash of edge kmers -> supernodes
  my %super_graph = ();
  my @supernodes = ();

  for my $key (keys %$graph) {
    if(!defined($graph->{$key}->{'visited'})) {
      my $contig = $graph->get_supernode($key);
      $graph->mark_kmers_visited($contig);
      my $supernode = {'seq' => $contig};
      push(@supernodes, $supernode);
      my $key0 = kmer_key(substr($contig, 0, $kmer_size));
      my $key1 = kmer_key(substr($contig, -$kmer_size));
      $super_graph{$key0} = $supernode;
      $super_graph{$key1} = $supernode;
    }
  }

  # Print nodes
  if($print_shades)
  {
    for my $supernode (@supernodes)
    {
      my $seq = $supernode->{'seq'};
      my $kmer0 = substr($seq,0,$kmer_size);
      my $shades0 = $graph->{kmer_key($kmer0)}->{'shades'};


      my $kmer_len = length($seq)+1-$kmer_size;
      my $num_shaded_nodes = $kmer_len < 3 ? $kmer_len : 3;
      my $num_of_shades = length($shades0);
      my $num_of_cols = $num_shaded_nodes*$num_of_shades+($num_shaded_nodes-1);

      print $seq . ' [shape=none label=<<table ' .
  'border="'.(defined($shades0) && $shades0 =~ /^\-+$/ ? '1' : '0').'" '.
  'cellborder="0" cellpadding="0" cellspacing="0">
  <tr><td PORT="top" colspan="'.$num_of_cols.'" cellpadding="0" cellspacing="0" border="0">
  '.($use_points ? '.' : $seq).'</td>
  </tr><tr>';

      # print first kmer shades
      print_kmer_shades($shades0);

      if($num_shaded_nodes == 3) {
        # Get middle shades
        my $kmer1 = substr($seq,1,$kmer_size);
        my $shades1 = $graph->{kmer_key($kmer1)}->{'shades'};

        for(my $i = 2; $i < length($seq)-$kmer_size; $i++)
        {
          $kmer1 = substr($seq,$i,$kmer_size);
          my $tmp_shades = $graph->{kmer_key($kmer1)}->{'shades'};

          # Merge if not equal
          if($tmp_shades ne $shades1)
          {
            if(!$shade_mismatch) {
              warn("Shades mismatch within a supernode");
              $shade_mismatch = 1;
            }

            # Merge
            for(my $j = 0; $j < length($shades1); $i++) {
              my ($a,$b) = map {substr($_,$j,1)} ($shades1,$tmp_shades);
              my ($uc,$lc) = (0,0);
              if($a ne '.') {
                if($a eq uc($a)) {$uc = $a;}
                if($a eq lc($a)) {$lc = $a;}
              }
              if($b ne '.') {
                if($b eq uc($b)) {$uc = $b;}
                if($b eq lc($b)) {$lc = $b;}
              }
              my $c = '-';
              if($uc && $lc) { $c = '-'; }
              elsif($uc) { $c = $uc; }
              elsif($lc) { $c = $lc; }
              substr($shades1,$j,1) = $c;
            }
          }
        }

        # Print middle kmer shades
        print '<td>|</td>'."\n";
        print_kmer_shades($shades1);
      }

      if($num_shaded_nodes > 1) {
        # Print last kmer shades
        my $kmer2 = substr($seq,-$kmer_size);
        my $shades2 = $graph->{kmer_key($kmer2)}->{'shades'};
        print '<td>|</td>'."\n";
        print_kmer_shades($shades2);
      }

      print "</tr></table>>];\n";
    }
  }
  else
  {
    for my $supernode (@supernodes) {
      print "  $supernode->{'seq'}\n";
    }
  }

  # Print edges
  for my $supernode (@supernodes)
  {
    my $kmer0 = substr($supernode->{'seq'}, 0, $kmer_size);
    my $kmer1 = substr($supernode->{'seq'}, -$kmer_size);
    my ($key0, $key1) = map {kmer_key($_)} ($kmer0, $kmer1);
    my $reverse0 = get_orientation($kmer0, $key0);
    my $reverse1 = get_orientation($kmer1, $key1);

    my @prev_edges = $graph->get_edges($key0,!$reverse0);
    my @next_edges = $graph->get_edges($key1,$reverse1);

    # print "@prev_edges:$kmer0  $kmer1:@next_edges\n";

    for my $next (@next_edges) {
      my $kmer = substr($supernode->{'seq'},-$kmer_size+1).$next;
      my $next_supernode = $super_graph{kmer_key($kmer)};
      print_supernode($supernode, $next, $next_supernode, $kmer, 1);
    }

    for my $prev (@prev_edges) {
      my $kmer = revcmp($prev).substr($supernode->{'seq'},0,$kmer_size-1);
      my $prev_supernode = $super_graph{kmer_key($kmer)};
      print_supernode($supernode, $prev, $prev_supernode, $kmer, 0);
    }
  }
}
else
{
  # Not 'simplifying' contigs
  if($print_shades)
  {
    $pid = open2($in, $out, $cmdline) or die("Cannot run cmd: '$cmdline'");

    while(defined(my $line = <$in>))
    {
      my ($kmer, $covgs, $edges, $shades) = parse_ctx_line($line);
      if(defined($kmer))
      {
        my $num_shades = defined($shades) ? length($shades) : 1;
        print $kmer . ' [shape=none label=<<table ' .
              'border="'.(defined($shades) && $shades =~ /^\-+$/ ? '1' : '0').'" '.
              'cellborder="0">
  <tr><td PORT="top" colspan="'.$num_shades.'" cellpadding="0" cellspacing="0">
  '.($use_points ? '.' : $kmer).'</td>
  </tr><tr>';

        print_kmer_shades($shades);

        print "</tr></table>>];\n";
      }
    }

    close($in);
    close($out);
  }

  $pid = open2($in, $out, $cmdline) or die("Cannot run cmd: '$cmdline'");

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
}

close($in);
close($out);

print "}\n";

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

sub print_kmer_shades
{
  my ($shades_txt) = @_;

  my @sharr = defined($shades_txt) ? split('', $shades_txt) : ('.');

  for(my $i = 0; $i < @sharr; $i++) {
    print '<td fixedsize="true" width="3" height="3" ' .
          'cellpadding="0" cellspacing="0" border="1" ';
    if($sharr[$i] ne '.') {
      if($sharr[$i] eq '-') { print 'bgcolor="black"'; }
      elsif($sharr[$i] eq lc($sharr[$i])) { print 'style="rounded"'; }
      else { print 'style="rounded" bgcolor="'.$cols[$i].'"'; }
      print ' color="'.$cols[$i].'"';
    }
    else { print 'color="white"'; }
    print '></td>'."\n";
  }
}

sub print_supernode
{
  my ($supernode,$rbase,$next_supernode,$kmer0,$going_right) = @_;

  my $kmer_size = length($kmer0);

  my $kmer1a = substr($next_supernode->{'seq'}, 0, $kmer_size);
  my $kmer1b = substr($next_supernode->{'seq'}, -$kmer_size);
  my ($key0,$key1a,$key1b) = map {kmer_key($_)} ($kmer0, $kmer1a, $kmer1b);

  my ($kmer1, $key1);
  if($key0 eq $key1a) { $kmer1 = $kmer1a; $key1 = $key1a; }
  elsif($key0 eq $key1b) { $kmer1 = $kmer1b; $key1 = $key1b; }
  else { die("Mismatch in supernode edges"); }

  my $arrive_left = $kmer0 eq ($going_right ? $kmer1a : revcmp($kmer1a));

  if(($supernode->{'seq'} lt $next_supernode->{'seq'}) ||
     ($supernode->{'seq'} le $next_supernode->{'seq'} &&
      ($arrive_left != $going_right || $arrive_left && $going_right)))
  {
    my $from = $supernode->{'seq'};
    my $to = $next_supernode->{'seq'};
    if($print_shades) { ($from,$to) = map {"$_:top"} ($from,$to); }
    print "  $from:" . ($going_right ? 'e' : 'w') . " -> " .
          "$to:"  . ($arrive_left ? 'w' : 'e') . "\n";
  }
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
    $key1 = $print_shades ? "$key1:top:" : "$key1:";
    $key2 = $print_shades ? "$key2:top:" : "$key2:";
    print "  $key1" . ($rev1 ? 'w' : 'e') . " -> " .
             $key2  . ($rev2 ? 'e' : 'w') . "\n";
  }
}
