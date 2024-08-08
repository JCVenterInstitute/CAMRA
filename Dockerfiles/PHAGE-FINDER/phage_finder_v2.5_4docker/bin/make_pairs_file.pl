#!/usr/bin/env perl

# Perl script to auto-generate the pairwise BLAST results file, pairs.file, for display using LinearDisplay.pl

my $prog = $0;
$prog =~ s/.*\///;
my $invocation = $prog . " @ARGV";

use strict;
use Getopt::Std;
getopts ('b:de:f:g:hi:l:p:V');

############## Declare global variables #################
our ($opt_h,$opt_b,$opt_d,$opt_e,$opt_f,$opt_g,$opt_i,$opt_l,$opt_p,$opt_V);
my ($btab_file, $frag_file, $gene_att, $evalue, $percentid, $min_hit_length, $pep_file, $DEBUG);
my $version = "1.0";

############## parse options ###############
if ($opt_h) { &option_help; }
if ($opt_V) {die "$version\n";}
if (-s $opt_b) {
    $btab_file = $opt_b;
}
else {
    print "You must provide a WUBLAST or NCBI-formated btab file as input\n";
    &option_help;
}
if (-s "$opt_f") {
    $frag_file = $opt_f;
} 
else {
    print "You must provide a frag.file as input\n";
    &option_help;
}
if (-s "$opt_g") {
    $gene_att = $opt_g;
}
else {
    print "You must provide a gene attribute file as input\n";
    &option_help;
}
if (-s "$opt_p") {
    $pep_file = $opt_p;
} 
else { # if no value for option p (pepfile), quit with help menu
    print "You must provide a pep file as input\n";
    &option_help;
}

if ($opt_e) {$evalue = $opt_e;} else { $evalue = 0.00001; } # if no E-value cut-off given, make default 0.00001
if ($opt_i) {$percentid = $opt_i;} else { $percentid = 35.0; } # the minimum cutoff to use blast scores for possible matches
if ($opt_l) {$min_hit_length = $opt_l;} else { $min_hit_length = 1; }
if ($opt_d) {$DEBUG = 1;} else {$DEBUG = 0;} # Debug mode.

my @order = ();  # array to store the frag_id order
my %Locus_lookup = (); # Key = fragment id, value = locus name
my %feat_hash = (); # only used to store protein lengths
my %relationship_hash = (); # Key1 = Query ($a[0]), Key2 = Subject (hit), Key3 = % id

###############
# Subroutines #
###############

sub option_help {
    print <<_EOB_;
    $prog - To auto-generate a pairs.file for making Linear figures with BLAST results
   Usage: $prog <options>
  Switch: -h for help\
  Option:
    -b: btab file (WUBLAST-format or NCBI -m 8 or 9 or BLAST+ -outfmt 6) [REQUIRED]
    -e: e-value [DEFAULT = 0.00001]
    -f: frag.file used by LinearDisplay.pl [REQUIRED]
    -g: gene_att.file used by LinearDisplay.pl [REQUIRED]
    -p: name of concatinated .pep file [REQUIRED to calc protein lengths]
    -i: aa % identity cut-off [DEFAULT = 35]
    -l: minimum % match length [DEFAULT = 1]
    -V: print version information
    -d: DEBUG MODE (default = 0)
  Format: 
  Output:
_EOB_
     exit;
}

sub get_frag_ids { # from frag.file

    my ($infile) = shift;
    my $prefix = "";
    my @a = ();
    my @b = ();
    my %Prefixes = ();

    open (INFILE, "<$infile") || die ("ERROR: can't open file $infile: $!\n");
    print "Getting the genome plot order from $infile ...\n";
    while (<INFILE>)  {
	chomp;
	if (/^\d+$/) {next}; # skip blank lines
	@a = split(/\t/);
	## code to skip over same genome matches (ie, when we have multiple frag_ids from the same genome)
	#@b = split(/_/,$a[1]);
	$prefix = $a[1]; # not sure what this prefix nonsense if all about... to avoid duplications I guess
	#$prefix = join("_", @b[1..$#b-1]);
	print "prefix = <$prefix>\n" if ($DEBUG);
	if (!$Prefixes{$prefix}) {
	    push (@order, $a[1]);
	    $Prefixes{$prefix} = 1; # added this new prefix to the list of prefixes.
	}
    }
    close (INFILE);
}

sub locus_lookup { # from gene_att.file

    my ($infile) = shift;
    my @a = ();
    my @b = ();
    my $locus = "";
    open (INFILE, "<$infile") || die ("ERROR: can't open file $infile: $!\n");
    print "Obtaining genome tag and locus information from $infile ...\n";
    while (<INFILE>)  {
	chomp;
	if (/^$/) {next}; # skip blank lines
	@a = split(/\t/);
	#@b = split(/_/,$a[7]);
	#$locus = "$a[4]-$b[0]";
	$locus = "$a[2]"; # not sure why I didn't pull from column 2 of the gene_att for for this
        $Locus_lookup{$a[7]}{$locus} = 1; # key1 = fragment id, key2 = locus name, value = 1
	print "Locus_lookup:  tag = <$a[7]>, locus = <$locus>\n" if ($DEBUG);
    }
    close (INFILE);
}

sub get_protein_info { #modified from PanOCT

  my ($infile) = shift;
  my $id;
  my $tag;
  my $sequence = "";
  my $length = "";

  unless (open (PEPFILE, "<$infile"))  {
    die ("can't open file $infile.\n");
  }
  print "Retrieving protein length information\n";
  my ($save_input_separator) = $/;
  $/="\n>";
  while (<PEPFILE>) {
    ($id,$tag,$sequence) = /^>?\s*([\w|\.]+)-(\w+).*\n([^>]+)>?/; # split the header line and sequence (very cool)
    $id =~ s/>//;
    $id = $id . "-" . $tag;
    $sequence =~ s/[^a-zA-Z]//g; # remove any non-alphabet letters
    $length = length($sequence);
    $feat_hash{$id}->{'length'}= $length; # only store the length, not the sequence
    print "$id = $feat_hash{$id}->{'length'}\n" if ($DEBUG);
    $sequence = ""; #clear out the sequence for the next round.
  }
  $/ = $save_input_separator; # restore the input separator
  close (PEPFILE);
  return;
}

sub select_data_from_btab { # get tab-delimited BLAST results in either WUBLAST or NCBI blastall (-m8/-m9) or blast+ (-outfmt 6) formats; modified from PanOCT

    my ($infile) = shift; 
    my @btab_line = (); # array variable to store split btab lines
    my $btab = ""; # 1 = WU, 0 = NCBI
    my $qmatch_length = ""; # stores the query length
    my $pct_qmatch_len = ""; # stores the percent query match length
    my $smatch_length = ""; # stores the subject length
    my $pct_smatch_len = ""; # stores the percent subject match length
    my $qid = ""; # query id
    my $sid = ""; # subject id (from database)
    my $query; # $btab_line[0]
    my $subject; # $btab_line[5]
    my $qbegin = ""; # start query
    my $qend = ""; # end query
    my $sbegin = ""; # start subject
    my $send = ""; # end subject
    my $pid = ""; # percent identity
    my $evlu = ""; # e-value
    my $qlength = ""; # size of query protein sequence
    my $slength = ""; # size of subject (database match) protein sequence
    open (INFILE, "<$infile") || die ("ERROR: can't open file $infile: $!\n");

    # this while loop is to interrogate the file to determine which tabular output style it is
    while (<INFILE>)  {
	chomp;
	if (!/^#/) { # don't look at REMs
	    @btab_line = split(/\t/);
	    last; # we have our info, break the loop
	}
    }
    close(INFILE);

    if ($#btab_line >= "19") { # adjusted because Badgers btab.pl script only goes to e-value or perl column 19
	print "Detected WUBLAST-style btab file ...\n";
	$btab = 1; # WU
    }
    elsif ($#btab_line == "11") {
	print "NCBI blastall (-m 8 or -m 9) or BLAST+ (-outfmt 6) option btab file detected ...\n";
	$btab = 0; # NCBI
    }
    else {
	die ("ERROR:  BLAST data must be either WUBLAST btab or NCBI blastall -m8 or -m9, or BLAST+ (-outfmt 6) formats.\n");
    }
    ### process BLAST results ###
    open (INFILE, "<$infile");
    while (<INFILE>)  {
	chomp;
	@btab_line = split(/\t/);
        # same variables for both btab styles
	$qbegin = $btab_line[6];
	$qend = $btab_line[7];
	$sbegin = $btab_line[8];
	$send = $btab_line[9];

	if ($btab) { # WU

            # adjusted because Badgers WU btab.pl script only goes to e-value or perl column 19
	    # ========================================================
	    # btab output for WUBLAST output
	    # column number Description (for Perl), add 1 for Unix
	    # 0       Query Sequence Name
	    # 1       Date of the Analysis
	    # 2       Query Sequence Length
	    # 3       Search Method  --  Blast family application name
	    # 4       Database Name
	    # 5       Subject Sequence Name  --  Database entry name
	    # 6       Start of alignment on query (5' nucleotide match in query)
	    # 7       End of alignment on query (3' nucleotide match in query)
	    # 8       Start of alignment on subject (5' nucleotide match in db hit)
	    # 9       End of alignment on subject (3' nucleotide match in db hit)
	    # 10      % Identity 
	    # 11      % Similarity 
	    # 12      Score (bits)
	    # 13      File Offset for Beginning of Alignment
	    # 14      File Offset for End of Alignment
	    # 15      Description (annotatioon)
	    # 16      Frame  --  1 through 6, or NULL
	    # 17      Query Strand  --  Plus, Minus or NULL
	    # 18      DB sequence length
	    # 19      Expect -- expected value
	    # 20      P-Value  --  Poisson ratio
	    # ========================================================

	    # propigate variables for WUBLAST
	    $query = $btab_line[0];
	    $subject = $btab_line[5];
	    $pid = $btab_line[10];
	    $evlu = $btab_line[19];
	}
	else { # NCBI
	    # ========================================================
	    # btab output from NCBI blastall (-m 8 or -m 9) or blast+ (-outfmt 6) option:
	    # column number Description (for Perl), add 1 for Unix
	    # 0      Query_id
	    # 1	     subject_id (Hit from db)
	    # 2	     % Identity
	    # 3	     length of alignment
	    # 4	     number or mismatches
	    # 5	     number of gaps
	    # 6	     start of alignment on query (5' nucleotide match in query)
	    # 7	     end of alignment on query (3' nucleotide match in query)
	    # 8	     start of alignment on subject (5' nucleotide match in db hit)
	    # 9	     end of alignment on subject (3' nucleotide match in db hit)
	    # 10     e-value
	    # 11     score (bits)
	    # ========================================================
	    
	    # propigate variables for NCBI m8/m9
	    if ($btab_line[0] =~ /^#/) { next;} # skip the lines beginning with #
	    $query = $btab_line[0];
	    $subject = $btab_line[1];
	    $pid = $btab_line[2];
	    $evlu = $btab_line[10];
	}
	### Generic processing ###
	#if ($query =~ /-/) { # if a tag in the feat_name, parse the feat_name and the tag depending on how many dashes there are
	#    my @tmpSplitQ = split(/-/, $query);
	#    $qid = $tmpSplitQ[0];
	#}
	#else { # no tag provided with feat_name, so use whole query name
	    $qid = $query;
	#}
	#if ($subject =~ /-/) {
	#    my @tmpSplitS = split(/-/, $subject);
	#    $sid = $tmpSplitS[0];
	#}
	#else { # no tag provided with feat_name, so use whole subject name
	    $sid = $subject;
	#}
	print "qid = $qid, sid = $sid\n" if ($DEBUG);
	$qlength = $feat_hash{$qid}->{'length'};
	$slength = $feat_hash{$sid}->{'length'};
	$qmatch_length = abs(($qbegin - $qend) + 1);
	$pct_qmatch_len = $qmatch_length/$qlength*100;
	print "$qid length = $qlength, qmatch_length = $qmatch_length, %query_match = $pct_qmatch_len, min_hit_length = $min_hit_length\n" if ($DEBUG);
        $smatch_length = abs(($sbegin - $send) + 1);
	$pct_smatch_len = $smatch_length/$slength*100;
	if ($pid < $percentid)  {next;}
	if ($evlu > $evalue)  {next;}
	if ($pct_qmatch_len < $min_hit_length)  {next;}
	if ($pct_smatch_len < $min_hit_length)  {next;}
	if (!defined $relationship_hash{$qid}{$sid})  {
	    $relationship_hash{$qid}{$sid}->{'id'} = $pid;
	    print "Query: $qid X Subject $sid = $relationship_hash{$qid}{$sid}->{'id'}\n" if ($DEBUG);
	}
    }
    close (INFILE);
}

################
#  M  A  I  N  #
################

&get_frag_ids($frag_file);
&locus_lookup($gene_att);
&get_protein_info($pep_file);
&select_data_from_btab($btab_file);
# for loop to loop through @order array, take $i and $i + 1 until $#a - 1, match frag and print pairs.file
# may need to play around with varibles and how I will actually accomplish this task
my $i = "";
my $j = "";
my $key1 = ""; #qid
my $key2 = ""; #sid
my @K1 = ();
my @K2 = ();
print "Writing pairs information to autopairs.file ...\n";
unless (open (OUTFILE, ">autopairs.file")) {
    die ("ERROR: can write to file autopairs.file.\n");
}
for ($i=0; $i<($#order); $i++){
    $j = $i+1;
    print "$i = $order[$i]\n" if ($DEBUG);
    foreach $key1 (sort {$a <=> $b} keys %{$Locus_lookup{$order[$i]}}) {
	#@K1 = split(/-/,$key1); # don't want to print the tags in this case
	#print "\t$key1\n" if ($DEBUG);
	foreach $key2 (keys %{$Locus_lookup{$order[$j]}}) {
	    #@K2 = split(/-/,$key2);
	    #print "\t\t$key2\n" if ($DEBUG);
	    if ($relationship_hash{$key1}{$key2}->{'id'}) {
	      #print OUTFILE "$order[$i]\t@K1[0]\t$order[$j]\t@K2[0]\t$relationship_hash{$key1}{$key2}->{'id'}\n";
	      print OUTFILE "$order[$i]\t$key1\t$order[$j]\t$key2\t$relationship_hash{$key1}{$key2}->{'id'}\n";
	    }
        }
    }
}
close (OUTFILE);
exit(0);
