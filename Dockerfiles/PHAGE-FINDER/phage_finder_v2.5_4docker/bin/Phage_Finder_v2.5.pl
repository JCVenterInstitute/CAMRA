#!/usr/bin/env perl
#
#Copy (C) 2022-2023  The J. Craig Venter Institute (JCVI).  All rights reserved
#Written by Derrick E. Fouts, Ph.D.

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

# ver 2.5 added the ability to test ncRNAs as targets of integration. Cleaned up redundant code in find_att_sites.  Fixed an issue with house-keeping genes with bad com_names being filled in to include a t/tm/ncRNA to test, but were not removed if no att sites were found in the t/tm/ncRNA genes.  Now, it only addes these genes when an att site has been found. Fixed a bug that inabled tRNAs to be "targets" when they were inside of the integrase (i.e., on the wrong side).
# ver 2.4 added the ability to identify receptor-binding proteins (RBPs) analogous to other "core" genes that are currently identified
# ver 2.4 also adds the ability to generate LinearDisplay.pl frag.file, which is needed along with the gene.att file to generate high-quality linear images of each prophage region.
# ver 2.3 added the ability to generte the LinearDisplay.pl gene.att file, which is needed to generate linear ORF map images of each prophage region.
# ver 2.2 to update code to enable the path to the phage_finder home directory to be provided on the command line (removed home directory requirement).  Also fixed a bug where the end5 of the phage region would not be reported in cases where a tRNA gene was at the 5' end and the tRNA gene was on the negative strand.  Also fixed bugs associated with the xgraph and reports files where negative numbers were being printed.  Reverted back to printing all hits in the xgraph plot rather than those > minnumberhits.  The pico, medio, and regions files were fixed to remove "singleton" att sites that would be printed when an att site was removed after finding house-keeping genes during autofill to a potential tRNA att site.  Fixed a bug where an att site was found, but not reported in the PFPR_tab.txt file because a tRNA gene was not being included (ie, was considered a bad com_name).  Finally, added more comments and cleaned-up the code as recommended by Dr. Phil Goetz also at JCVI that were discovered during integration of Phage_Finder into the JCVI automated annotation pipeline.
# ver 2.1 to fix bug in Phage_subs.pm that prevented the proper creation and usage of the phage_finder_info.txt file when providing a Genbank .ptt file.
# ver 2.0 to fix failure to properly loop through every contig.  Added HMMER3 input capability. Standardized region .pep, .seq, and .con output filenames.  Changed .tab to PFPR_tab.txt to easy opening in MS Excel. Added curated lists of core phage functions to aid in identification of these genes.
# ver 1.9beta to fix # of hits per window that is allowable to be considered as a valid phage region (gdv phage region 8 had 3 hits per window)
# Perl Script to find prophages in bacterial genomes using btab output and the end5 of each orf
# ver 1.8 12/30/02: added capability to search for phage HMM hits and include data in output files
# ver 1.8 12/18/02: restructured subs find_5prime_end and find_3prime_end to look for hits slightly upstream of orf (within window or next window) and use an array to move through end5s of orfs
# ver 1.7 11/14/02: fix problem where 1 window with no BLAST hits within a larger phage region prematurely terminates the region - look at next window as well (sub find_regions) 
# ver 1.6 11/14/02: check for existence of ($phage_hash{1}->{'seed_end5'}) [ie a putative hit region] so that the program will not waist time extending regions that do not exist
# ver 1.6 streamlined counting loop within MAIN to delete hash keys of $searchhash when already counted and also end looping through keys of searchhash when keys >= window size
# ver 1.6 fixed problem where only one role_id was stored per ORF (now $hithash and $rehash can hold multiple role_ids)
# ver 1.2 fills in all featnames into %hithash (in sub select_featnames_from_btab), but only hit=1 if $a[10] != -100
# ver 1.2 use SQL to only pull end5s of orfs in specific role categories (typical of phage regions) and record that role_id
# ver 1.3 will add a max peak height counter and store median of max peak in $fighash{$n}->{'median'} (will later be used for outward scans to find more potential phage genes
# ver 1.1 includes step size separate from window size (like selfsim)
#
# Boolean logic convention True = 1, False = 0
my $prog = $0;
$prog =~ s/.*\///;
my $invocation = $prog . " @ARGV";

use strict;
use FindBin qw($Bin); # Find where the executable is.  Added per Phil Goetz's suggestion 04/16/13
use lib "$Bin/../lib"; #phage library location
use Math::Round qw(nlowmult nhimult round);
use PHAGE::Phage_subs;
#use Bio::SeqFeature::Generic; # added 05/20/2022 to write GenBank files
#use Bio::SeqIO; # added 05/20/2022 to write GenBank files
use Getopt::Std;
getopts ('b:c:df:hi:m:n:p:r:s:t:w:A:E:F:H:P:SV'); # Rearranged alphabetically per Phil Goetz's suggestion 04/16/13 (except I separate lc from uc)

our ($opt_b,$opt_c,$opt_d,$opt_f,$opt_h,$opt_i,$opt_m,$opt_n,$opt_p,$opt_r,$opt_s,$opt_t,$opt_w,$opt_A,$opt_E,$opt_F,$opt_H,$opt_S,$opt_V);
my ($blastbin,$basedir,$mode_basedir,$logfile,$default_hpw,$window,$step,$evalue,$asmbl_id,$btabpath,$btabfile,$infofile,$hmmdir,$hmmfile,$tRNA_data,$tmRNA_data,$ncRNA_data,$asmbly_file_path,$asmbly_file,$search_method,$strict,$mode,$badcomnamesfile,$frag_line_num,$DEBUG);
my $version = "2.5";

if ($opt_d) {$DEBUG = 1;} else {$DEBUG = 0;} # Debug mode.
if ($opt_h) { &option_help; }
if ($opt_V) {die "$version\n";}
if ($opt_S) { # if strict mode (print only regions with core HMM hits or Mu-like and > 10 Kbp)
    $strict = 1;
    $mode = "strict";
}
else {
    $strict = 0;
    $mode = "nonstrict";
}

## if no value for option b (base or working directory), use current working directory
if ($opt_b) {
    $basedir = $opt_b;
}
else { 
    $basedir = $ENV{'PWD'};
}

## make directory structure for output ## moved 03/07/14 dfouts
print "preMAIN:  basedir = <$basedir>\n" if ($DEBUG);
$mode_basedir = &make_basedir($basedir,$mode,$DEBUG); # make new base directory (strict or nonstrict)
$badcomnamesfile = "$mode_basedir/bad_comnames.txt"; # added 09/11/14
print "preMAIN:  mode_basedir = <$mode_basedir>\n" if ($DEBUG);
##

### Create log file in base directory ###
print "preMAIN:  Creating the log file ...\n" if ($DEBUG);
$logfile = &create_log("$mode_basedir/phage_phinder");
&write_log("0", $invocation);

## if no path to btab file provided, make default = base directory
if ($opt_p) {
    $btabpath = $opt_p;
}
else {
    $btabpath = $basedir;
}
## if no value for option t (name of btab file), quit with help menu with error message
if (($opt_t) && (-e "$btabpath/$opt_t" == 1)) {
    $btabfile = "$btabpath/$opt_t";
}
elsif (-z "$btabpath/$opt_t") { # if file is present, but not BLAST results, then don't bother looking for prophages
    &write_log("5", "No BLASTP results, aborting prophage finding"); # added 03/07/14 by dfouts since small contigs may not have BLAST matches
}
else {
    print "FATAL error: Must specify -t <BLASTP output>\n";
    &option_help;
}
## if no value for option t (name of info file), quit with help menu with error message
if (($opt_i)  && (-e "$basedir/$opt_i" == 1) && (-z "$basedir/$opt_i" == 0)) {
    $infofile = "$basedir/$opt_i";
}
else {
    print "FATAL error: Must specify -i <phage_finder_info.txt file>\n";
    &option_help;
}
## if no value for option m (name of hmm file), quit with help menu
$hmmdir = $basedir; # need this variable to capture the location of the original basedir
if ((length($opt_m) >0)  && (-e "$hmmdir/$opt_m" == 1) && (-z "$hmmdir/$opt_m" == 0)) {
    $hmmfile = "$hmmdir/$opt_m";
}
else {
    undef $hmmfile;
}
if ($opt_w) {$window = $opt_w;} else { $window = 10000; } # if no window size given, make default "10000"
if ($opt_s) {$step = $opt_s;} else { $step = 5000; } # if no step size given, make default "5000"
if ($opt_E) {$evalue = $opt_E;} else { $evalue = 0.00001; } # if no E-value given, make default 0.00001
if ($opt_H) {$default_hpw = $opt_H;} else {$default_hpw = 4;} # if the # of hits per window is not defined, default = 4
# Removed defined_asmbl_id function (opt_a) - not necissary - data stored in the info_file 04/17/13 dfouts
## Check whether the tRNAscan output file is present (either full path or relative to $basedir).  Not dies if not provided.
# added 04/17/13 dfouts - suggested by pgoetz
if ($opt_r) {
    if (-s $opt_r) {
	$tRNA_data = $opt_r;
    }
    elsif (-s "$basedir/$opt_r") {
	$tRNA_data = "$basedir/$opt_r";
    }
    else {
	undef $tRNA_data; # set as undefined since no suitable data
    }
}
## Check whether the Aragorn (tmRNA finder) output file is present (either full path or relative to $basedir).  Dies if not provided.
# added 04/17/13 dfouts - suggested by pgoetz
if ($opt_n) {
    if (-s $opt_n) {
	$tmRNA_data = $opt_n;
    }
    elsif (-s "$basedir/$opt_n") {
	$tmRNA_data = "$basedir/$opt_n";
    }
    else {
	undef $tmRNA_data; #set as undefined since no suitable data
    }
}
## Check whether the cmscan (ncRNA finder) output file is present (either full path or relative to $basedir).  Dies if not provided.
# added 01/12/22 dfouts
if ($opt_c) {
    if (-s $opt_c) {
	$ncRNA_data = $opt_c;
    }
    elsif (-s "$basedir/$opt_c") {
	$ncRNA_data = "$basedir/$opt_c";
    }
    else {
	undef $ncRNA_data; #set as undefined since no suitable data
    }
}
## Check whether the Assembly or contig/scaffold fasta file is present (either full path or relative to $basedir).  Dies if not provided.
# added 04/17/13 dfouts - suggested by pgoetz
if ($opt_A) {
    if (-s $opt_A) {
	$asmbly_file = $opt_A;
    }
    elsif (-s "$basedir/$opt_A") {
	$asmbly_file = "$basedir/$opt_A";
    }
    else {
	die "FATAL error: Could not find assembly/contig/scaffold fasta file $opt_A";
    }
}
else {
    print "FATAL error: Must specify -A <assembly/contig/scaffold fasta file>\n";
    &option_help;
}
### att site prediction methods ###
if (($opt_F eq "M") || ($opt_F eq "m"))  {$search_method = "mummer";}
elsif (($opt_F eq "F") || ($opt_F eq "f"))  {$search_method = "fasta33";}
elsif (($opt_F eq "B") || ($opt_F eq "b")) {
    $search_method = "blast";
}
else {
    $search_method = "blast";
}

############## Declare variables #################
my %fighash = (); # gets cleared with each assembly
my %hithash = (); # 1D, key = feat_name
my %rehash = (); # 2D, key1 = asmbl_id, key2 = end5
my %searchhash = (); # 2D, key1 = asmbl_id, key2 = end5, value = feat_name
my %phage_hash = (); # gets cleared with each assembly
my %exclude_hash = (); # list of accessions to exclude
my %ok_comnames = ();
my %DB_info = (); # key is the BLAST DB tag, values are full name of phage and the taxonomy
my %asmbl_idhash = ();
my %assembly_hash = (); # gets cleared with each assembly
my %tRNA_hash = ();
my %serine_HMM_hash = ('PF00239' => 1);
my %tyrosine_HMM_hash = ('PF00589' => 1,
                         'PF02899' => 1);
my %lists_hash = (); # stores the accessions for head morphogenesis Small and Large Terminases, Portal, Protease, Scaffolding, and Major Capsid proteins
my %HMMs_hash = (); # stores the curated HMM accessions
#my $blastbin = "/usr/local/packages/ncbi-blast+-2.4.0/bin"; #created 07/06/2016 [dfouts]
#my $blastbin = "/usr/local/devel/ANNOTATION/dfouts/ncbi-blast-2.8.1+/bin"; # created 11/10/21 [dfouts]
#my $blastbin = "/Users/dfouts/Installers/blast+/ncbi-blast-2.11.0+/bin"; # for OS X
my $blastbin = "/scratch/dfouts/bin/ncbi-blast-2.14.0+/bin"; # new blastbin location
my $phome = "$Bin/.."; # home removed. Using $phome to make relative to what was executed (using FindBin)
my $write_dir = "";
my @genomearray = ();
my $HMMversion = ""; # new variable to store the HMM version (currently 2 or 3)
my @hmm_data = ();
my $comfile = "$phome/phage_com_names_combo_v3.txt";
my @filenames = ();
# future update: pull crib file from config.file set-up
my $DBfile = "$phome/DB/phage_03_25_19.crib";
my $num_contigs = "";
my $phage_regions = "";
my $asmbly_status = 0; # boolean 0 = undefined, missing or inconsistant assembly data, 1 = ok assembly data
my $genome_phage_cnt = 0; # set the genome-wide phage counter to zero.
my $prefix = "";
my $hitsperwindow = $default_hpw; # initialize the hitsperwindow variable

sub find_regions {

  my ($asmbl_id,$hitsperwindow,$figref,$hitref,$phageref,$DEBUG) = @_;
  my $regions = 0;
  my $max_peak = 0;
  my $newcluster = 0;
  my $end5 = "";
  my $fig = "";
  my $hold_fig = "";
  my $prk = ""; # stores the phageref key counts

  foreach $fig (sort {$a <=> $b} keys %{$figref}) {
    #if (($newcluster == 1) && ($figref->{$fig}->{'counts'} < ($hitsperwindow-1)) && (!defined($phageref->{$regions}->{'seed_end5'})))  { # if really the end of the region and not a window with 0 hits in the middle
    # note sure why we were checking if the seed_end5 was set or not
    if (($newcluster == 1) && ($figref->{$fig}->{'counts'} < ($hitsperwindow-1))) {
      $regions++;
      $phageref->{$regions}->{'seed_end5'} = $end5;
      $figref->{$hold_fig}->{'peak_number'} = $regions;
      $figref->{$hold_fig}->{'peak_value'} = $max_peak;
      $max_peak = 0; # reset max peak indicator to zero
      $newcluster = 0;
      print "find_regions: region = $regions, end5 = $end5, seed_end5 = $phageref->{$regions}->{'seed_end5'}\n" if ($DEBUG);
    }
    else {
      $newcluster = 0;
    }
    if (($figref->{$fig}->{'counts'} >= $hitsperwindow) && ($figref->{$fig}->{'counts'} > $max_peak))  {  #### Default = 4, ver 1.9beta changed to user defined ?
      $max_peak = $figref->{$fig}->{'counts'};
      $end5 = $hitref->{${$figref->{$fig}->{'featnames'}}[0]}->{'end5'};
      $hold_fig = $fig; # added 04/16/13 by dfouts because the value of the step with max_peak was not being held in short contigs
      print "find_regions: HITSperWINDOW = $hitsperwindow\n" if ($DEBUG);
      print "find_regions: FOUND asmbl_id ($asmbl_id): $regions, max_peak = $max_peak, and end5 = $end5\n" if ($DEBUG);
    }
    elsif ($max_peak > 0) {
    #elsif (($max_peak > 0) && ($figref->{$fig}->{'counts'} == 0) && (($end5 - $phageref->{$regions}->{'seed_end5'}) > 70000)) {
      # we want to merge regions if we have 2 peaks within 70000 bp (was set to 35000 bp, but Acinetobacter baumannii strain MRSN 3527 has a ~ 80 kbp phage with a 10kb ORF spanning two regions that needed merged, with a distance of
      # ~ 65 kb between peaks, so I bumped this up to 70kb (dfouts 10/23/2013).
      # note this > 70000 condition is a temporary fix.  Really want both peaks processed and then a decision made at the end to determine if 1 or 2 phages (piggy back or same phage).
      # brought the > 70000 bp condition down from above because some regions were getting skipped over by having their region-specific data rewritten in the phageref hash
      # this did not work.  I think we will just let close regions get lumped together for now and write a separate subroutine to tease piggy-back phage apart in the future 02/03/2017
      $newcluster = 1;
      #$hold_fig = $fig;
      print "find_regions(zero count):  maxpeak = $max_peak, newcluster = $newcluster\n" if ($DEBUG);
    }
  }
  if ($max_peak > 0) { #new simplified code for recording bulk phage regions at the end of a contig dfouts 11/26/14
    $regions++;
    $phageref->{$regions}->{'seed_end5'} = $end5;
    $figref->{$hold_fig}->{'peak_number'} = $regions; # added 04/16/13 by dfouts to store values for xgraph drawing in short contigs/scaffolds
    $figref->{$hold_fig}->{'peak_value'} = $max_peak; # added 04/16/13 by dfouts to store values for xgraph drawing in short contigs/scaffolds
    print "find_regions(PHAGE_AT_END):  $end5\t$phageref->{$regions}->{'seed_end5'}\t$regions\thold_fig=$hold_fig\n" if ($DEBUG);
  }
  # if (($regions == 0) && ($max_peak > 0) && (!defined($phageref->{$regions+1}->{'seed_end5'})))  {  # case where we run out of sequence and never got 0 phage hits (found bug with scaffold searches).  Changed from if $newcluster = 0 to checking that no regions were defined dfouts 11/07/2013.  Added check if $regions == 0 07/16/14 because bogus phage regions were being created when the was a merged peak at the end of a contig.  Only want to check this exception if no regions were defined since they will just be merged.  Original code was written before I checked the distance between peaks to be considered a new region.
#    $regions++;
#    $phageref->{$regions}->{'seed_end5'} = $end5;
#    $figref->{$hold_fig}->{'peak_number'} = $regions; # added 04/16/13 by dfouts to store values for xgraph drawing in short contigs/scaffolds
#    $figref->{$hold_fig}->{'peak_value'} = $max_peak; # added 04/16/13 by dfouts to store values for xgraph drawing in short contigs/scaffolds
#    print "FIX:  $end5\t$phageref->{$regions}->{'seed_end5'}\t$regions\thold_fig=$hold_fig\n" if ($DEBUG);
#  }
#  elsif (($newcluster == 0) && ($max_peak > 0) && (!defined($phageref->{$regions+1}->{'seed_end5'})))  { # needed to add this code back for cases where the phage is at the end of the contig/assembly and never got recorded! 
#      $regions++;
#      $phageref->{$regions}->{'seed_end5'} = $end5;
#      $figref->{$hold_fig}->{'peak_number'} = $regions;
#      $figref->{$hold_fig}->{'peak_value'} = $max_peak;
#      print "PHAGE_at_END:  $end5\t$phageref->{$regions}->{'seed_end5'}\t$regions\thold_fig=$hold_fig\n" if ($DEBUG);
#  }
  #delete ($phageref->{"0"}); # delete the zero key
  $prk = keys %{$phageref};
  print "find_regions[end]: Total region count = <$regions>, hash key count = <$prk>\n" if ($DEBUG);
  return ($regions);
}

sub find_5prime_end {

  my ($n,$asmbl_id,$aref,$phageref,$hitref,$tRNAref,$reref,$okref,$figref,$DEBUG) = @_;
  my $i = "";
  my $lastkey = $aref->[$#{$aref}]+1; # start $lastkey as the 3' most coordinate
  my $hold = "";
  my $j = "";
  my $find = "";
  my $rRNAdist = 2000;
  my $feat_name = "";
  my $m = ""; # holds the next 5' region counter (n-1)
  my $regions = $n; # holds the total number of phage regions and returns it if we merge any
  my $prk = ""; # holds the number of phage_ref keys

  $prk = keys %{$phageref};
  print "find_5prime_end[start]: #phage_regions = $regions, #phageref keys = <$prk>\n" if ($DEBUG);

  for ($i = $#{$aref}; $i >= 0; $i--)  {
    $prk = keys %{$phageref};
    print "|====>i = [$i]<==== n = <$n> ---->#phageref keys = <$prk><----|\n" if ($DEBUG);
    if ($n == "0") {
      last; #end for loop if n=0
    }
    elsif ((exists $phageref->{$n}) && ($aref->[$i] <= $phageref->{$n}->{'seed_end5'}) && ($aref->[$i] < $lastkey))  {
      print "------->n = $n, i = $i, aref = <$aref->[$i]>, lastkey = <$lastkey>, #phageref keys = <$prk><-------\n" if ($DEBUG);
      my $lm = nlowmult($window, $aref->[$i])-$window;
      print "find_5prime_end: n = $n\tworking on $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}, clean_name = <$hitref->{$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}}->{'clean_name'}>, lowmult = \[$lm\]...\n" if ($DEBUG);
      print "CURRENT FEAT_NAME:  $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}, clean_name = <$hitref->{$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}}->{'clean_name'}>\n" if ($DEBUG);
      print "ok_comname: <$okref->{$hitref->{$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}}->{'clean_name'}}>\n" if ($DEBUG);
      print "dbHit:  <$hitref->{$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}}->{'hit'}>\n" if ($DEBUG);
      #print "COUNTS: $figref->{nlowmult($window, $aref->[$i])}->{'counts'}\tNAME: $figref->{nlowmult($window, $aref->[$i])}->{'featnames'}[0]\n" if ($DEBUG);
      #print "nextCOUNTS: $figref->{nlowmult($window, $aref->[$i])-$step}->{'counts'}\tnextNAME: $figref->{nlowmult($window, $aref->[$i])-$step}->{'featnames'}[0]\n" if ($DEBUG);
      print "currentSTEPcounts:  $figref->{nlowmult($step, $aref->[$i])}->{'counts'}\tNAME: $figref->{nlowmult($step, $aref->[$i])}->{'featnames'}[0]..$figref->{nlowmult($step, $aref->[$i])}->{'featnames'}[$#{ $figref->{nlowmult($step, $aref->[$i])}->{'featnames'} }]\n" if ($DEBUG);
      print "nextSTEPcounts:  $figref->{nlowmult($step, $aref->[$i]+$step)}->{'counts'}\tNAME: $figref->{nlowmult($step, $aref->[$i]+$step)}->{'featnames'}[0]..$figref->{nlowmult($step, $aref->[$i])+$step}->{'featnames'}[$#{ $figref->{nlowmult($step, $aref->[$i]+$step)}->{'featnames'} }]\n" if ($DEBUG);
      print "nextWINDOWcounts: $figref->{nlowmult($window, $aref->[$i])+$window}->{'counts'}\tNAME: $figref->{nlowmult($window, $aref->[$i])+$window}->{'featnames'}[0]..$figref->{nlowmult($window, $aref->[$i]+$window)}->{'featnames'}[$#{ $figref->{nlowmult($window, $aref->[$i]+$window)}->{'featnames'} }]\n" if ($DEBUG);

      ############################### Check for running past next region ##############################
      #moved from near end of loop
      if ((exists $phageref->{$n-1}) && ($aref->[$i] < $phageref->{$n-1}->{'seed_end5'})) { # case where phage region n has gone past the seed end 5 of n-1 01/16/19
	#note: to avoid perl adding in the reference to "0" when n-1 is equal to 0 (via autovivification) we need to also check $phageref->{}
	# move phage region $n to $n-1 and remove the old region $n
	$m = $n-1;
	print "find_5prime_end_passthru: -----> Region = <$m> will be merged with Region <$n>, Region n-1 <$m> seed <$phageref->{$m}->{'seed_end5'}> will become seed <$phageref->{$n}->{'seed_end5'}>\n" if ($DEBUG);
	$phageref->{$m} = $phageref->{$n}; # make the new phage region equal to the old one
	print "find_5prime_end_passthru: -----> New seed = <$phageref->{$m}->{'seed_end5'}>\n" if ($DEBUG);
	print "find_5prime_end_passthru_delete: n = $n\n" if ($DEBUG);
	delete($phageref->{$n}); # delete phage region $n that is being replaced (and overlapping)
	&rekey_regions($m, $phageref, $DEBUG); # call subroutine to renumber phage regions
	$regions--; #decrement the phage region counter
	$prk = keys %{$phageref};
	print "find_5prime_end_passthru_check_region_counts: #phage_regions = $regions, #phageref keys = <$prk>\n" if ($DEBUG);
	$n = $m;  # set n equal to n-1 (next region
      }
      ################################################################################################
      if (exists $hitref->{$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}}->{'hmm'}) {
	$phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
	$lastkey = $aref->[$i];
	$phageref->{$n}->{'last5ORF'} = $lastkey;
	print "find_5prime_end: HMM5(1a): <$aref->[$i]> = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG);
      }
      # if this ORF has a hit in the Phage db, then include
      elsif ($hitref->{$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}}->{'hit'} == 1) {
	$phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
	$lastkey = $aref->[$i];
	$phageref->{$n}->{'last5ORF'} = $lastkey;
	print "find_5prime_end: HIT5(1a): <$aref->[$i]> = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG);
      }
      # if tRNA, tmRNA or ncRNA is present and the distance between last orf and tRNA is < 2 Kbp (to avoid getting those between rRNAs), include $phageref->{$key}->{'tRNA'}{$t}
      elsif ((exists $tRNAref->{$asmbl_id}{$aref->[$i]}) && (($lastkey - $aref->[$i]) < $rRNAdist))  {
	$phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
	$phageref->{$n}->{'tRNA'}{$aref->[$i]} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
	$lastkey = $aref->[$i];
	print "find_5prime_end: tRNA5: <$aref->[$i]> = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG);
      }
      elsif ((exists $tRNAref->{$asmbl_id}{$aref->[$i-1]}) && (($lastkey - $aref->[$i-1]) < $rRNAdist) && ($i-1 >= 0))  {
	$phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
	$phageref->{$n}->{'memberhash'}{$aref->[$i-1]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i-1]}->{'featname'};
	$phageref->{$n}->{'tRNA'}{$aref->[$i-1]} = $reref->{$asmbl_id}{$aref->[$i-1]}->{'featname'};
	$lastkey = $aref->[$i-1];
	print "find_5prime_end: tRNA5bypass(1a): <$aref->[$i]> = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG);
	print "find_5prime_end: tRNA5bypass(-1a): <$aref->[$i-1]> = $reref->{$asmbl_id}{$aref->[$i-1]}->{'featname'}\n" if ($DEBUG);
      }
      elsif ((exists $tRNAref->{$asmbl_id}{$aref->[$i-2]}) && (($lastkey - $aref->[$i-2]) < $rRNAdist) && ($i-2 >= 0) ) {
	$phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
	$phageref->{$n}->{'memberhash'}{$aref->[$i-1]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i-1]}->{'featname'};
	$phageref->{$n}->{'memberhash'}{$aref->[$i-2]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i-2]}->{'featname'};
	$phageref->{$n}->{'tRNA'}{$aref->[$i-2]} = $reref->{$asmbl_id}{$aref->[$i-2]}->{'featname'};
	$lastkey = $aref->[$i-2];
	print "find_5prime_end: tRNA5bypass(1b): <$aref->[$i]> = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG);
	print "find_5prime_end: tRNA5bypass(-1b): <$aref->[$i-1]> = $reref->{$asmbl_id}{$aref->[$i-1]}->{'featname'}\n" if ($DEBUG);
	print "find_5prime_end: tRNA5bypass(-2b): <$aref->[$i-2]> = $reref->{$asmbl_id}{$aref->[$i-2]}->{'featname'}\n" if ($DEBUG);
      }
      elsif ((exists $tRNAref->{$asmbl_id}{$aref->[$i-3]}) && (($lastkey - $aref->[$i-3]) < $rRNAdist) && ($i-3 >= 0)) {
	$phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
	$phageref->{$n}->{'memberhash'}{$aref->[$i-1]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i-1]}->{'featname'};
	$phageref->{$n}->{'memberhash'}{$aref->[$i-2]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i-2]}->{'featname'};
	$phageref->{$n}->{'memberhash'}{$aref->[$i-3]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i-3]}->{'featname'};
	$phageref->{$n}->{'tRNA'}{$aref->[$i-3]} = $reref->{$asmbl_id}{$aref->[$i-3]}->{'featname'};
	$lastkey = $aref->[$i-3];
	print "find_5prime_end: tRNA5bypass(1c): <$aref->[$i]> = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG);
	print "find_5prime_end: tRNA5bypass(-1c): <$aref->[$i-1]> = $reref->{$asmbl_id}{$aref->[$i-1]}->{'featname'}\n" if ($DEBUG);
	print "find_5prime_end: tRNA5bypass(-2c): <$aref->[$i-2]> = $reref->{$asmbl_id}{$aref->[$i-2]}->{'featname'}\n" if ($DEBUG);
	print "find_5prime_end: tRNA5bypass(-3c): <$aref->[$i-3]> = $reref->{$asmbl_id}{$aref->[$i-3]}->{'featname'}\n" if ($DEBUG);
      }
      # if ok common name and there are still hits in the current window, then include
      elsif (($okref->{$hitref->{$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}}->{'clean_name'}} == 1) &&
	     ($figref->{nlowmult($window, $aref->[$i])}->{'counts'} >= 1) && 
	     ($hitref->{$figref->{nlowmult($window, $aref->[$i])}->{'featnames'}[0]}->{'end5'} < $aref->[$i]))
	{
	  $phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
	  $lastkey = $aref->[$i];
	  $phageref->{$n}->{'last5ORF'} = $lastkey;
	  print ">>$figref->{nlowmult($window, $aref->[$i])}->{'counts'}\t$figref->{nlowmult($window, $aref->[$i])}->{'featnames'}[0]<<\n" if ($DEBUG);
	  print "find_5prime_end: OKcom_name5: <$aref->[$i]> = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG);
	}

      # if ok common name and there are still hits in the next window, then include
      ### special note:  changed counts from >= 2 to >= 1 for the Campylobacter jejuni 84-25 genome since there were 2 peaks separated with a single hit in one window
      elsif (($okref->{$hitref->{$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}}->{'clean_name'}} == 1) &&
	     ($figref->{nlowmult($window, $aref->[$i])-$window}->{'counts'} >= 1) && 
	     ($hitref->{$figref->{nlowmult($window, $aref->[$i])-$window}->{'featnames'}[0]}->{'end5'} < $aref->[$i]))
	{
	  $phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
	  $lastkey = $aref->[$i];
	  $phageref->{$n}->{'last5ORF'} = $lastkey;
	  print ">>$figref->{nlowmult($window, $aref->[$i])-$step}->{'counts'}\t$figref->{nlowmult($window, $aref->[$i])-$step}->{'featnames'}[0]<<\n" if ($DEBUG);
	  print "find_5prime_end: OKcom_name5b: <$aref->[$i]> = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG);
	}

      # case where >= 3 hits in next window, but no clean name (bad annotation perhaps or something novel

      elsif (($figref->{nlowmult($window, $aref->[$i])-$step}->{'counts'} >= 3) && 
	     ($hitref->{$figref->{nlowmult($window, $aref->[$i])-$step}->{'featnames'}[0]}->{'end5'} < $aref->[$i]))
	{
	  $phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
	  $lastkey = $aref->[$i];
	  $phageref->{$n}->{'last5ORF'} = $lastkey;
	  print ">>$figref->{nlowmult($window, $aref->[$i])}->{'counts'}\t$figref->{nlowmult($window, $aref->[$i])}->{'featnames'}[0]<<\n" if ($DEBUG);
	  print "find_5prime_end: OKcom_name5c: <$aref->[$i]> = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG);
	}

      #        # if there are at least 2 hits in the current step and end5 of first featname having hit in current window is less than end 5 of current orf, include
      elsif (($figref->{nlowmult($step, $aref->[$i])}->{'counts'} >= 2) && 
	     ($hitref->{$figref->{nlowmult($step, $aref->[$i])}->{'featnames'}[0]}->{'end5'} < $aref->[$i]))  {
	$hold = $hitref->{$figref->{nlowmult($step, $aref->[$i])}->{'featnames'}[0]}->{'end5'}; # store coordinate of 5' most hit in current window
	$phageref->{$n}->{'memberhash'}{$hold}->{'featname'} = $reref->{$asmbl_id}{$hold}->{'featname'};
	$j = $i;
	until ($hold == $aref->[$j])  { # fill in the other featnames until hit reached
	  $phageref->{$n}->{'memberhash'}{$aref->[$j]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$j]}->{'featname'};
	  print "find_5prime_end: fill (j = <$j>): <$aref->[$j]> = $phageref->{$n}->{'memberhash'}{$aref->[$j]}->{'featname'}\n" if ($DEBUG);
	  if ((exists $tRNAref->{$asmbl_id}{$aref->[$j]}) && (($aref->[$j-1] - $aref->[$j]) < $rRNAdist)) {
	    $phageref->{$n}->{'tRNA'}{$aref->[$j]} = $reref->{$asmbl_id}{$aref->[$j]}->{'featname'};
	  }
	  elsif ((exists $tRNAref->{$asmbl_id}{$aref->[$j]}) && (($aref->[$j-1] - $aref->[$j]) >= $rRNAdist)) { # likely within a rRNA operon, remove
	    $hold = $lastkey;
	    print "find_5prime_end: DELETING possible rRNA operon <$aref->[$j]> = $reref->{$asmbl_id}{$aref->[$j]}->{'featname'}\n" if ($DEBUG);
	    until ($hold == $aref->[$j])  {
	      print "find_5prime_end: deleting: <$aref->[$j]> = $phageref->{$n}->{'memberhash'}{$aref->[$j]} and $phageref->{$n}->{'tRNA'}{$aref->[$j]}\n" if ($DEBUG);
	      delete($phageref->{$n}->{'memberhash'}{$aref->[$j]});
	      delete($phageref->{$n}->{'tRNA'}{$aref->[$j]});
	      $j++;
	    }
	    $feat_name = $reref->{$asmbl_id}{$lastkey}->{'featname'};
	    if (($feat_name =~ /tm?RNA/) || ($feat_name =~ /ncRNA/)) { # added 04/15/2013 and mod 01/13/22 by dfouts for cases where tRNA/tmRNA gene is at end
	      #RNA featnames stored differently in hithash (to make unique added asmbl_id to freatname)than rehash
	      $feat_name = $asmbl_id . "_" . $feat_name;
	    }
	    if ($lastkey < $hitref->{$feat_name}->{'end3'})  {  # if 5 prime ORF ----->
	      if (exists $phageref->{$n}->{'tRNA'}{$lastkey}) { # case where tRNA is at 5' end
		$phageref->{$n}->{'5prime'} = $phageref->{$n}->{'last5ORF'};
	      }
	      else {
		$phageref->{$n}->{'5prime'} = $lastkey; # if find "housekeeping gene" quit the search and record the 5' end as the last good orf
	      }
	    }
	    else {   # <---- make 3' end of 5' most ORF the 5' boundary of phage region
	      $phageref->{$n}->{'5prime'} = $hitref->{$feat_name}->{'end3'};
	    }
	    $phageref->{$n}->{'ORF5'} = $reref->{$asmbl_id}{$lastkey}->{'featname'}; # store feat_name of 5' most ORF
	    print "find_5prime_end_1: ==> Region $n 5prime end = $phageref->{$n}->{'5prime'}\tORF = $phageref->{$n}->{'ORF5'}\n" if ($DEBUG);
	    $n = $n - 1; # decrement phage region key counter when finished with phage region
	    #if ($n == 0)  {return;} #not sure now if this is really necissary
	    next;
	  }
	  $j--;
	}
	$lastkey = $hold;
	$phageref->{$n}->{'last5ORF'} = $lastkey if (!exists $phageref->{$n}->{'tRNA'}{$lastkey});
      }

      # specifically check for the word "transposase" in the com_name - these are ok and are often not the first word
      elsif ($hitref->{$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}}->{'com_name'} =~ /transposase/) {
	$phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
	$lastkey = $aref->[$i];
	$phageref->{$n}->{'last5ORF'} = $lastkey;
	print "find_5prime_end: OKtransposase5: <$aref->[$i]> = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG);
      }

      else {
	###### added 10/14/05 to include a tRNA that was target when there are no db hits leading to it up to 10 kb
	my $pos = $i;
	my $end5 = $aref->[$pos];
	## added ($pos < 0) condition 05/14/07 because if phage region was at 5prime end of a contig, this would loop indefinately
	## changed to ($pos < 1) 01/10/19 because $pos can become -1, which caused the array position to go to the other end of the array and pick up a tRNA on the opposite (3prime) end
	until ((($lastkey-$end5) > 10000) || ($pos < 1)){ # only check within 10 kb of lastkey
	  $pos--;
	  $end5 = $aref->[$pos];
	  if (exists $tRNAref->{$asmbl_id}{$end5}) {
	    $phageref->{$n}->{'memberhash'}{$end5}->{'featname'} = $reref->{$asmbl_id}{$end5}->{'featname'};
	    $phageref->{$n}->{'tRNA'}{$end5} = $reref->{$asmbl_id}{$end5}->{'featname'};
	    $lastkey = $end5;
	    print "find_5prime_end: tRNAextend: <i = $i, pos = $pos, aref = $aref->[$i], end5 = $end5> = $reref->{$asmbl_id}{$end5}->{'featname'}\n" if ($DEBUG);
	  }
	}
	$feat_name = $reref->{$asmbl_id}{$phageref->{$n}->{'last5ORF'}}->{'featname'};
	if (($feat_name =~ /tm?RNA/) || ($feat_name =~ /ncRNA/)) { # added 04/15/2013 and mod 01/13/22 by dfouts for cases where tRNA/tmRNA gene is at end
	  #RNA featnames stored differently in hithash (to make unique added asmbl_id to featname) than rehash
	  $feat_name = $asmbl_id . "_" . $feat_name;
	}

	#else { # DEFAULT case where we are not at the 5prime end of the contig
	print "find_5prime_end_2c: -----> Region = <$n>, Seed = <$phageref->{$n}->{'seed_end5'}>\n" if ($DEBUG);
	if ($lastkey < $hitref->{$feat_name}->{'end3'})  {  # if 5 prime ORF ----->
	  $phageref->{$n}->{'5prime'} = $phageref->{$n}->{'last5ORF'}; # if find "housekeeping gene" quit the search and record the 5' end as the last good orf
	  print "find_5prime_end_2c: -----> feat_name = $feat_name, last5ORF = $phageref->{$n}->{'last5ORF'}, 5prime end of region = $phageref->{$n}->{'5prime'}\n" if ($DEBUG);
	}
	else {   # <---- make 3' end of 5' most ORF the 5' boundary of phage region
	  $phageref->{$n}->{'5prime'} = $hitref->{$feat_name}->{'end3'};
	  print "find_5prime_end_2c: <----- feat_name = $feat_name, last5ORF = $phageref->{$n}->{'last5ORF'}, 5prime end of region = $phageref->{$n}->{'5prime'}\n" if ($DEBUG);
	}
	$phageref->{$n}->{'ORF5'} = $reref->{$asmbl_id}{$phageref->{$n}->{'last5ORF'}}->{'featname'};
	print "find_5prime_end_2c: ==> Region $n 5prime end = $phageref->{$n}->{'5prime'}\tORF = $phageref->{$n}->{'ORF5'}\n" if ($DEBUG);
	$prk = keys %{$phageref};
	print "find_5prime_end_2c_check_region_counts: #phage_regions = $regions, #phageref keys = <$prk>\n" if ($DEBUG);
	#### found a case where a bad com_name near int was in the region but the int was skipped because of this.  Need to check if next ORF has a hit or an a phage HMM before decrementing the region counter?
	$n = $n - 1; # decrement phage region key counter when finished with phage region
#	}
      }
    }
  }
  if ($n > 0)  {  # weird case when phage is at the end of a molecule (when plasmid in genome project is actually the replicative form of an induced prophage or draft genome contig)
    $feat_name = $reref->{$asmbl_id}{$lastkey}->{'featname'};
    if (($feat_name =~ /tm?RNA/) || ($feat_name =~ /ncRNA/)) { # added 04/15/2013 and od 01/13/22 by dfouts for cases where tRNA/tmRNA gene is at end
      #RNA featnames stored differently in hithash (to make unique added asmbl_id to freatname)than rehash
      $feat_name = $asmbl_id . "_" . $feat_name;
    }
    # if we are truely at the 5prime end of the contig or genome, then set region to 1 #
    if ($phageref->{$n}->{'last5ORF'} < $hitref->{$feat_name}->{'end3'})  {  # if 5 prime ORF ----->
      $phageref->{"1"}->{'5prime'} = $phageref->{$n}->{'last5ORF'}; # if find "housekeeping gene" quit the search and record the 5' end as the last good orf
    }
    else {   # <---- make 3' end of 5' most ORF the 5' boundary of phage region
      $phageref->{"1"}->{'5prime'} = $hitref->{$feat_name}->{'end3'};
    }
    $phageref->{"1"}->{'ORF5'} = $reref->{$asmbl_id}{$phageref->{$n}->{'last5ORF'}}->{'featname'};
    print "find_5prime_end_3: ==> Region $n 5prime end = $phageref->{$n}->{'5prime'}\tORF = $phageref->{$n}->{'ORF5'}, lastkey = <$lastkey>\n" if ($DEBUG);
    print "find_5prime_end_3: ==> Region 1 5prime end = $phageref->{\"1\"}->{'5prime'}\tORF = $phageref->{\"1\"}->{'ORF5'}, lastkey = <$lastkey>\n" if ($DEBUG);
    until  ($n == "1") { # 01/10/19
      print "find_5prime_end_3_delete: n = $n\n" if ($DEBUG);
      delete($phageref->{$n}); # delete phage region that is being replaced (and overlapping) from the nth region until the first region
      $regions--; #decrement the phage region counter
      $n = $n - 1; # decrement phage region key counter when finished with phage region
      &rekey_regions($n, $phageref, $DEBUG); # call subroutine to renumber phage regions
    }
    $n = 0;
  }
  $prk = keys %{$phageref};
  print "find_5prime_end:[END] ==> region counter n = <$n> phage regions = <$regions>, #phageref keys = <$prk>\n" if ($DEBUG);
  print ".....................................................[ END of find_5prime_end ] ...................................................\n" if ($DEBUG);
  return($regions);
}

sub find_3prime_end { # relative to contig 5' --->

  my ($max,$asmbl_id,$aref,$phageref,$hitref,$tRNAref,$reref,$okref,$figref,$DEBUG) = @_;
  my $i = "";
  my $j = "";
  my $m = ""; # placeholder for the value of $n, the region counter
  my $n = 1; # region counter
  my $hold = "";
  my $lastkey = "";
  my $find = "";
  my $rRNAdist = 2000;
  my $feat_name = "";
  my $regions = $max; # holds the total number of phage regions and returns it if we merge any
  my $prk = ""; # holds the number of phage_ref keys
  my $pos = "";
  my $end5 = "";
  my $clean = ""; # boolean if clean com_name observed in region 0 = no, 1 = yes;

  $prk = keys %{$phageref};
  print "find_3prime_end[start]: #phage_regions = <$regions>, #phage_ref keys = <$prk>\n" if ($DEBUG);

  for $i ( 0 .. $#{$aref} )  {
    $prk = keys %{$phageref};
    print "|====>i = [$i]<==== n = <$n> ---->#phageref keys = <$prk><----|\n" if ($DEBUG);
    if ((exists $phageref->{$n}) && ($aref->[$i] >= $phageref->{$n}->{'seed_end5'}) && ($n <= $max) && ($aref->[$i] > $lastkey))  {
      #print "------->n = $n, i = $i, aref = $aref->[$i], next seed <$phageref->{$n+1}->{'seed_end5'}>, #phageref keys = <$prk><-------\n" if ($DEBUG);
      print "CURRENT FEAT_NAME:  $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}, clean_name = <$hitref->{$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}}->{'clean_name'}>\n" if ($DEBUG);
      print "COUNTS: $figref->{nlowmult($window, $aref->[$i])}->{'counts'}\tNAME: $figref->{nlowmult($window, $aref->[$i])}->{'featnames'}[0]..$figref->{nlowmult($window, $aref->[$i])}->{'featnames'}[$#{ $figref->{nlowmult($window, $aref->[$i])+$step}->{'featnames'} }]\n" if ($DEBUG);
      print "currentSTEPcounts:  $figref->{nlowmult($step, $aref->[$i])}->{'counts'}\tNAME: $figref->{nlowmult($step, $aref->[$i])}->{'featnames'}[0]..$figref->{nlowmult($step, $aref->[$i])}->{'featnames'}[$#{ $figref->{nlowmult($step, $aref->[$i])}->{'featnames'} }]\n" if ($DEBUG);
      print "nextSTEPcounts:  $figref->{nlowmult($step, $aref->[$i]+$step)}->{'counts'}\tNAME: $figref->{nlowmult($step, $aref->[$i]+$step)}->{'featnames'}[0]..$figref->{nlowmult($step, $aref->[$i])+$step}->{'featnames'}[$#{ $figref->{nlowmult($step, $aref->[$i]+$step)}->{'featnames'} }]\n" if ($DEBUG);
      print "nextWINDOWcounts: $figref->{nlowmult($window, $aref->[$i])+$window}->{'counts'}\tNAME: $figref->{nlowmult($window, $aref->[$i])+$window}->{'featnames'}[0]..$figref->{nlowmult($window, $aref->[$i]+$window)}->{'featnames'}[$#{ $figref->{nlowmult($window, $aref->[$i]+$window)}->{'featnames'} }]\n" if ($DEBUG);
      print "find_3prime_end: last hit in next step = <$figref->{nlowmult($step, $aref->[$i])+$step}->{'featnames'}[$#{ $figref->{nlowmult($step, $aref->[$i])+$step}->{'featnames'} }]>, end5 = <$hitref->{$figref->{nlowmult($step, $aref->[$i])+$step}->{'featnames'}[$#{ $figref->{nlowmult($step, $aref->[$i])+$step}->{'featnames'} }]}->{'end5'}>, current end5 = <$aref->[$i]>\n" if ($DEBUG);
      ######################### check if we have gone past the seed end 5 of the next region #####################
      #moved from just above the n increment at the end of the loop, but this resulted in multiple regions being passed before checked
      if ((exists $phageref->{$n+1}) && ($aref->[$i] > $phageref->{$n+1}->{'seed_end5'})) { # case where phage region n has gone past the seed end 5 of n+1 01/18/19
	#$m = $n; # hold the value of $n
	#$n = $n + 1; # increment phage region key counter
	print "find_3prime_end_passthru_delete: n = $n, i = <$i>, aref = <$aref->[$i]>, next seed end 5 = <$phageref->{$n+1}->{'seed_end5'}>\n" if ($DEBUG);
	delete($phageref->{$n+1}); # delete phage region that is being replaced (and overlapping) from the nth region until the first region
	&rekey_regions($n, $phageref, $DEBUG); # call subroutine to renumber phage regions
	$regions--; #decrement the total number of phages counter
      }
      ############################################################################################################
      if (exists $hitref->{$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}}->{'hmm'}) {
	print "find_3prime_end:HMM3: <$aref->[$i]> = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG);
	$phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
	$lastkey = $aref->[$i];
	$phageref->{$n}->{'last3ORF'} = $lastkey;
      }
      # if this ORF has a hit in the Phage db, then include
      elsif ($hitref->{$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}}->{'hit'} == 1) {
	print "find_3prime_end: HIT3(1a): <$aref->[$i]> = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG);
	$phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
	$lastkey = $aref->[$i];
	$phageref->{$n}->{'last3ORF'} = $lastkey;
      }
      # if tRNA, tmRNA or ncRNA is present, include
      elsif ((exists $tRNAref->{$asmbl_id}{$aref->[$i]}) && (($aref->[$i] - $lastkey) < $rRNAdist)) {
	print "find_3prime_end: tRNA3: <$aref->[$i]> = $aref->[$i]\t$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG);
	$phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
	$phageref->{$n}->{'tRNA'}{$aref->[$i]} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
	$lastkey = $aref->[$i];
      }
      elsif ((exists $tRNAref->{$asmbl_id}{$aref->[$i+1]}) && (($aref->[$i+1] - $lastkey) < $rRNAdist) && ($i+1 <= $#{$aref})) {
	print "find_3prime_end: tRNA3bypass(1a): <$aref->[$i]> = $aref->[$i]\t$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG);
	print "find_3prime_end: tRNA3bypass(+1a): <$aref->[$i+1]> = $aref->[$i+1]\t$reref->{$asmbl_id}{$aref->[$i+1]}->{'featname'}\n" if ($DEBUG);
	$phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
	$phageref->{$n}->{'memberhash'}{$aref->[$i+1]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i+1]}->{'featname'};
	$phageref->{$n}->{'tRNA'}{$aref->[$i+1]} = $reref->{$asmbl_id}{$aref->[$i+1]}->{'featname'};
	$lastkey = $aref->[$i+1];
      }
      elsif ((exists $tRNAref->{$asmbl_id}{$aref->[$i+2]}) && (($aref->[$i+2] - $lastkey) < $rRNAdist) && ($i+1 <= $#{$aref})) {
	print "find_3prime_end: tRNA3bypass(1b): $aref->[$i]\t$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG);
	print "find_3prime_end: tRNA3bypass(+1b): $aref->[$i+1]\t$reref->{$asmbl_id}{$aref->[$i+1]}->{'featname'}\n" if ($DEBUG);
	print "find_3prime_end: tRNA3bypass(+2b): $aref->[$i+2]\t$reref->{$asmbl_id}{$aref->[$i+2]}->{'featname'}\n" if ($DEBUG);
	$phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
	$phageref->{$n}->{'memberhash'}{$aref->[$i+1]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i+1]}->{'featname'};
	$phageref->{$n}->{'memberhash'}{$aref->[$i+2]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i+2]}->{'featname'};
	$phageref->{$n}->{'tRNA'}{$aref->[$i+2]} = $reref->{$asmbl_id}{$aref->[$i+2]}->{'featname'};
	$lastkey = $aref->[$i+2];
      }
      elsif ((exists $tRNAref->{$asmbl_id}{$aref->[$i+3]}) && (($aref->[$i+3] - $lastkey) < $rRNAdist) && ($i+1 <= $#{$aref})) {
	print "find_3prime_end: tRNA3bypass(1c): <$aref->[$i]> = $aref->[$i]\t$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG);
	print "find_3prime_end: tRNA3bypass(+1c): <$aref->[$i+1]> = $aref->[$i+1]\t$reref->{$asmbl_id}{$aref->[$i+1]}->{'featname'}\n" if ($DEBUG);
	print "find_3prime_end: tRNA3bypass(+2c): <$aref->[$i+2]> = $aref->[$i+2]\t$reref->{$asmbl_id}{$aref->[$i+2]}->{'featname'}\n" if ($DEBUG);
	print "find_3prime_end: tRNA3bypass(+3c): <$aref->[$i+3]> = $aref->[$i+3]\t$reref->{$asmbl_id}{$aref->[$i+3]}->{'featname'}\n" if ($DEBUG);
	$phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
	$phageref->{$n}->{'memberhash'}{$aref->[$i+1]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i+1]}->{'featname'};
	$phageref->{$n}->{'memberhash'}{$aref->[$i+2]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i+2]}->{'featname'};
	$phageref->{$n}->{'memberhash'}{$aref->[$i+3]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i+3]}->{'featname'};
	$phageref->{$n}->{'tRNA'}{$aref->[$i+3]} = $reref->{$asmbl_id}{$aref->[$i+3]}->{'featname'};
	$lastkey = $aref->[$i+3];
      }
      # if ok common name and there are still hits in the current window, then include
      elsif (($okref->{$hitref->{$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}}->{'clean_name'}} == 1) &&
	     (($figref->{nlowmult($window, $aref->[$i])}->{'counts'}) >= 1) &&
	     ($hitref->{$figref->{nlowmult($window, $aref->[$i])}->{'featnames'}[$#{ $figref->{nlowmult($window, $aref->[$i])}->{'featnames'} }]}->{'end5'} > $aref->[$i]))
	{
	  print "find_3prime_end: OKcom_name3: <$aref->[$i]> = ($hitref->{$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}}->{'clean_name'}) $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG);
	  $phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
	  $lastkey = $aref->[$i];
	  $phageref->{$n}->{'last3ORF'} = $lastkey;
	}

      # if ok common name and there are still hits in the next window, then include
      elsif (($okref->{$hitref->{$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}}->{'clean_name'}} == 1) &&
	     (($figref->{nlowmult($window, $aref->[$i])+$window}->{'counts'}) >= 1) &&
	     ($hitref->{$figref->{nlowmult($window, $aref->[$i])+$window}->{'featnames'}[$#{ $figref->{nlowmult($window, $aref->[$i])+$window}->{'featnames'} }]}->{'end5'} > $aref->[$i]))
	{
	  print "find_3prime_end: OKcom_name3b: <$aref->[$i]> = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG); 
	  $phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
	  $lastkey = $aref->[$i];
	  $phageref->{$n}->{'last3ORF'} = $lastkey;
	}

      # if not good com_name and there are still >=2 (changed from 3) hits in the next step, then include
      # changed nlowmult($window, to nlowmult($step, dfouts 11/07/2013
      elsif ((($figref->{nlowmult($step, $aref->[$i])+$step}->{'counts'}) >= 2) &&
	     ($hitref->{$figref->{nlowmult($step, $aref->[$i])+$step}->{'featnames'}[$#{ $figref->{nlowmult($step, $aref->[$i])+$step}->{'featnames'} }]}->{'end5'} > $aref->[$i]))
	{
	  print "find_3prime_end: OKcom_name3c: <$aref->[$i]> = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG); 
	  $phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
	  $lastkey = $aref->[$i];
	  $phageref->{$n}->{'last3ORF'} = $lastkey;
	}

      # if there is a hit in the current step and end3 of the last featname having hit in current window is less than end 5 of current orf, include
      elsif ((($figref->{nlowmult($step, $aref->[$i])}->{'counts'}) >= 2) && 
	     ($hitref->{$figref->{nlowmult($step, $aref->[$i])}->{'featnames'}[$#{ $figref->{nlowmult($step, $aref->[$i])}->{'featnames'} }]}->{'end5'} > $aref->[$i])) {
	$hold = $hitref->{$figref->{nlowmult($step, $aref->[$i])}->{'featnames'}[$#{ $figref->{nlowmult($step, $aref->[$i])}->{'featnames'} }]}->{'end5'};
	# if there are more than 3 ORFs until the next hit, then bail
	$phageref->{$n}->{'memberhash'}{$hold}->{'featname'} = $reref->{$asmbl_id}{$hold}->{'featname'};
	$j = $i;
	until ($hold == $aref->[$j])  {
	  print "find_3prime_end:  fill (j = <$j>): <$aref->[$j]> = to hit in next step [ $reref->{$asmbl_id}{$aref->[$j]}->{'featname'} ]\n" if ($DEBUG);
	  $phageref->{$n}->{'memberhash'}{$aref->[$j]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$j]}->{'featname'};
	  if ((exists $tRNAref->{$asmbl_id}{$aref->[$j]}) && (($aref->[$j] - $aref->[$j-1]) < $rRNAdist)) { 
	    $phageref->{$n}->{'tRNA'}{$aref->[$j]} = $reref->{$asmbl_id}{$aref->[$j]}->{'featname'};
	    print "find_3prime_end: STEP: <$aref->[$j]> = $phageref->{$n}->{'memberhash'}{$aref->[$j]}->{'featname'}\n" if ($DEBUG);
	  }
	  elsif ((exists $tRNAref->{$asmbl_id}{$aref->[$j]}) && (($aref->[$j] - $aref->[$j-1]) >= $rRNAdist)) { # likely within a rRNA operon, remove
	    $hold = $lastkey;
	    until ($hold == $aref->[$j])  {
	      print "find_3prime_end: deleting: <$aref->[$j]> = $phageref->{$n}->{'memberhash'}{$aref->[$j]} and $phageref->{$n}->{'tRNA'}{$aref->[$j]}\n" if ($DEBUG);
	      delete($phageref->{$n}->{'memberhash'}{$aref->[$j]});
	      delete($phageref->{$n}->{'tRNA'}{$aref->[$j]});
	      $j--;
	    }
	    $feat_name = $reref->{$asmbl_id}{$lastkey}->{'featname'};
	    if (($feat_name =~ /tm?RNA/) || ($feat_name =~ /ncRNA/)){ # added 04/15/2013 and mod 01/13/22 by dfouts for cases where tRNA/tmRNA gene is at end
	      #RNA featnames stored differently in hithash (to make unique added asmbl_id to freatname)than rehash
	      $feat_name = $asmbl_id . "_" . $feat_name;
	    }
	    ### come back to this to fix?
#	    if ($lastkey > $hitref->{$feat_name}->{'end3'})  {  # if 3 prime ORF <-----
#	      if (exists $phageref->{$n}->{'tRNA'}{$lastkey}) { # case where tRNA is at 3' end
#		$phageref->{$n}->{'3prime'} = $phageref->{$n}->{'last3ORF'}; #fix this so that the ORF 3' end is the outer end
#	      }
#	      else {
#		$phageref->{$n}->{'3prime'} = $lastkey; # if find "housekeeping gene" quit the search and record the 3' end as the last good orf
#	      }
#	    }
#	    else {   # ----> make 3' end of 5' most ORF the 3' boundary of phage region
#	      $phageref->{$n}->{'3prime'} = $hitref->{$feat_name}->{'end3'};
#	    }
	    #### new code ####
	    if ($lastkey > $hitref->{$feat_name}->{'end3'}) { # gene is <------5'
		print "find_3prime_end: making last ORF 5' end <----* the 3'prime end of the phage region\n" if ($DEBUG);
		$phageref->{$n}->{'3prime'} = $hitref->{$feat_name}->{'end5'};
	    }
	    else { # gene is ------>
		print "find_3prime_end: making last ORF 3' end <----* the 3'prime end of the phage region\n" if ($DEBUG);
		$phageref->{$n}->{'3prime'} = $hitref->{$feat_name}->{'end3'};
	    }
	    #### end new code ####
	    $phageref->{$n}->{'ORF3'} = $reref->{$asmbl_id}{$lastkey}->{'featname'};
	    print "find_3prime_end: <$aref->[$i]> = STEP ==> Region $n 3prime end = $phageref->{$n}->{'3prime'}\n" if ($DEBUG);
	    $n = $n + 1; # increment phage region key counter when finished with phage region
	    #if ($n == $max+1) {return;} #not sure this is really necissary
	    next;
	  }
	  $j++;
	}
	$lastkey = $hold;
	$phageref->{$n}->{'last3ORF'} = $lastkey if (!exists $phageref->{$n}->{'tRNA'}{$lastkey});
      }
      # specifically check for the word "transposase" in the com_name - these are ok and are often not the first word
      elsif ($hitref->{$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}}->{'com_name'} =~ /transposase/) {
	print "find_3prime_end: OKtransposase3: <$aref->[$i]> = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG);
	$phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
	$lastkey = $aref->[$i];
      }
      else {
	########### added 05/14/07 to include a tRNA that was target when there are no db hits leading to it up to 10 kb
	$pos = $i;
	$hold = $pos; # save starting position
	$end5 = $aref->[$pos];
	## added ($pos > $#{$aref}) condition 05/14/07 because if phage region was at 3prime end of a contig, this would loop indefinately
	print "find_3prime_end: else: i = <$i>, pos = <$pos>, aref = <$aref->[$i]>, end5 = <$end5> = <$reref->{$asmbl_id}{$end5}->{'featname'}>\n" if ($DEBUG);
	until ((($end5-$lastkey) > 10000) || ($pos > $#{$aref})){ # only check within 10 kb of lastkey
	  if ((exists $tRNAref->{$asmbl_id}{$end5}) && ($clean)) { # if tRNA and genes between our lastkey and the tRNA are acceptable.  Should not have a case where the tRNA is the first gene since that would have been picked up earlier
	    $phageref->{$n}->{'memberhash'}{$end5}->{'featname'} = $reref->{$asmbl_id}{$end5}->{'featname'};
	    $phageref->{$n}->{'tRNA'}{$end5} = $reref->{$asmbl_id}{$end5}->{'featname'};
	    $lastkey = $end5;
	    print "find_3prime_end: tRNAextend: <i = $i, pos = $pos, aref = $aref->[$i], end5 = $end5> = $reref->{$asmbl_id}{$end5}->{'featname'}\n" if ($DEBUG);
	    # need code to fill in memberhash up to this tRNA # -> commented out this because housekeeping genes were sneaking in since there was no check of the com_name.  Genes are added back if an att site is found in sub find_att_sites
	    #if ($hold < $pos) { #only need to backfill if we actually skipped past ORFs to get to the tRNA
	    #  until ($hold == ($pos-1)) {
	#	$phageref->{$n}->{'memberhash'}{$aref->[$hold]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$hold]}->{'featname'};
	#	print "find_3prime_end_fill: hold ($hold) = <$aref->[$hold]>, pos ($pos) = <$aref->[$pos]>\n" if ($DEBUG);
	#	$hold++;
	 #     }
	  #  }
	    #####
	  }
	  elsif ($okref->{$hitref->{$reref->{$asmbl_id}{$end5}->{'featname'}}->{'clean_name'}} == 1) {
	    $clean = 1;
	    print "find_3prime_end_clean:  clean = <$clean>, pos = <$pos>, end5 = <$end5>, locus = <$reref->{$asmbl_id}{$end5}->{'featname'}>, annotation = <$hitref->{$reref->{$asmbl_id}{$end5}->{'featname'}}->{'com_name'}>\n" if ($DEBUG);
	  }
	  else {
	    $clean = ""; # found a bad com_name, so bonk it
	    print "find_3prime_end_not_clean:  clean = <$clean>, pos = <$pos>, end5 = <$end5>, locus = <$reref->{$asmbl_id}{$end5}->{'featname'}>, annotation = <$hitref->{$reref->{$asmbl_id}{$end5}->{'featname'}}->{'com_name'}>\n" if ($DEBUG);
	  }
	  $pos++;
	  $end5 = $aref->[$pos];
	}
	###########
	$feat_name = $reref->{$asmbl_id}{$lastkey}->{'featname'};
	if (($feat_name =~ /tm?RNA/) || ($feat_name =~ /ncRNA/)) { # added 04/15/2013 and mod 01/13/22 by dfouts for cases where tRNA/tmRNA gene is at end
	  #RNA featnames stored differently in hithash (to make unique added asmbl_id to freatname)than rehash
	  $feat_name = $asmbl_id . "_" . $feat_name;
	}
	print "find_3prime_end:  lastkey = <$lastkey>, hitref_end5 = <$hitref->{$feat_name}->{'end5'}>, hitref_end3 = <$hitref->{$feat_name}->{'end3'}>\n" if ($DEBUG);
	######################### crazy stuff - need to rethink #########################
	#############include tRNA here and then remove later?
#	if (($lastkey > $hitref->{$feat_name}->{'end3'}) && ($hitref->{$feat_name}->{'end3'} < $hitref->{$feat_name}->{'end5'}))  {  # if 3 prime ORF 3'<-----'5 (fixed 07/07/22) changed ">" to "<" 07/13/22
#	    if (exists $phageref->{$n}->{'tRNA'}{$lastkey}) { # case where tRNA is at 3' end and facing proper direction (<-----?), pick last 3'ORF instead of the tRNA
#		print "find_3prime_end: tRNA at 3'prime end, making last ORF 5' end <----* the 3'prime end of the phage region\n" if ($DEBUG);
#		#$phageref->{$n}->{'3prime'} = $phageref->{$n}->{'last3ORF'};
#		$phageref->{$n}->{'3prime'} = $hitref->{$feat_name}->{'end5'};
#	    }
#	    else {
#		print "find_3prime_end: making last ORF 5' end <----* the 3'prime end of the phage region\n" if ($DEBUG);
#		$phageref->{$n}->{'3prime'} = $hitref->{$feat_name}->{'end3'};; # if find "housekeeping gene" quit the search and record the 3' end as the last good orf
#	    }
#	    print "find_3prime_end_2: -----> region = $n, feat_name = $feat_name, last3ORF = $phageref->{$n}->{'last3ORF'}, 3prime end of region = $phageref->{$n}->{'3prime'}\n" if ($DEBUG);
#	} #might want to think about this condition as I think it is not quite right
#	else {   # ----> make 3' end of 5' most ORF the 3' boundary of phage region
#	    print "find_3prime_end: making last ORF 3' <featname = $feat_name> end ---->* the 3'prime end of the phage region\n" if ($DEBUG);
#	    $phageref->{$n}->{'3prime'} = $hitref->{$feat_name}->{'end3'};
#	    print "find_3prime_end_2: -----> region = $n, feat_name = $feat_name, last3ORF = $phageref->{$n}->{'last3ORF'}, 3prime end of region = $phageref->{$n}->{'3prime'}\n" if ($DEBUG);
#	}

	##begin new cleaner code:
	# to set proper phage 3' end based on orientation of the last gene in the region
	# note: we will be including tRNA, tmRNA, ncRNA in the phage region
	#       if it is the target of insertion, then it will remain
	#       if not, it will be removed from the region later
	if ($lastkey > $hitref->{$feat_name}->{'end3'}) { # gene is <------5'
	    print "find_3prime_end: making last ORF 5' end <----* the 3'prime end of the phage region\n" if ($DEBUG);
		$phageref->{$n}->{'3prime'} = $hitref->{$feat_name}->{'end5'};
	}
	else { # gene is ------>
	    print "find_3prime_end: making last ORF 3' end <----* the 3'prime end of the phage region\n" if ($DEBUG);
		$phageref->{$n}->{'3prime'} = $hitref->{$feat_name}->{'end3'};
	}
	#################################################################################
	$phageref->{$n}->{'ORF3'} = $reref->{$asmbl_id}{$phageref->{$n}->{'last3ORF'}}->{'featname'};
	print "find_3prime_end_2: ==> Region $n 3prime end = $phageref->{$n}->{'3prime'}\tORF = $phageref->{$n}->{'ORF3'}\n" if ($DEBUG);
	$n = $n + 1; # increment phage region key counter when finished with phage region
      }
    }
    if ($n > "$regions") {
      print "find_3prime_end_max: n <$n> exceeded max number of phage regions <$regions>!\n" if ($DEBUG);
      last; #end for loop if n=$max (the maximum number of phage regions)
    }
  }
  # weird case when phage is at the end of a molecule (when plasmid in genome project is actually the replicative form of an induced prophage
  if ($n <= $regions)  {
    $feat_name = $reref->{$asmbl_id}{$lastkey}->{'featname'};
    if (($feat_name =~ /tm?RNA/) || ($feat_name =~ /ncRNA/)) { # added 04/15/2013 and mod 01/13/22 by dfouts for cases where tRNA/tmRNA gene is at end
      #RNA featnames stored differently in hithash (to make unique added asmbl_id to freatname)than rehash
      $feat_name = $asmbl_id . "_" . $feat_name;
    }

    if ($phageref->{$n}->{'last3ORF'} > $hitref->{$feat_name}->{'end3'})  {  # if 3 prime ORF <-----
      $phageref->{$n}->{'3prime'} = $phageref->{$n}->{'last3ORF'}; # if find "housekeeping gene" quit the search and record the 3' end as the last good orf
    }
    else {   # ----> make 3' end of 5' most ORF the 3' boundary of phage region
      $phageref->{$n}->{'3prime'} = $hitref->{$feat_name}->{'end3'};
    }

    $phageref->{$n}->{'ORF3'} = $reref->{$asmbl_id}{$phageref->{$n}->{'last3ORF'}}->{'featname'};
    print "find_3prime_end_3: ==> Region $n 3prime end = $phageref->{$n}->{'3prime'}\tORF = $phageref->{$n}->{'ORF3'}\n" if ($DEBUG);
    print "find_3prime_end_3: contig length = <$aref->[$#{$aref}]>, last nucleotide position = <$lastkey>\n" if ($DEBUG);
    if ((exists $phageref->{$n+1}) && (exists $phageref->{$n+1}->{'seed_end5'})) {
      print "find_3prime_end_3*: Seed next region = <$phageref->{$n+1}->{'seed_end5'}>, Array position = <$i>, Current nuc. position = <$aref->[$i]>\n" if ($DEBUG);
    }
    if ((exists $phageref->{$n+1}) && ($lastkey > $phageref->{$n+1}->{'seed_end5'})) { # case where phage region n has gone past the seed end 5 of n+1 01/18/19
      # $aref->[$i] should be $lastkey since if we are at the end of the array for loop and $i is wiped after the loop
      $m = $n; # hold the value of $n
      $n = $n + 1; # increment phage region key counter
      print "find_3prime_end_3_delete: n = $n\n" if ($DEBUG);
      delete($phageref->{$n}); # delete phage region that is being replaced (and overlapping)
      &rekey_regions($m, $phageref, $DEBUG); # call subroutine to renumber phage regions
      $regions--; # decrement the total number of phages counter
    }
      $n = $n + 1; # increment phage region key counter when finished with phage region
  }
  $prk = keys %{$phageref};
  print "find_3prime_end:[END] == phage regions = <$regions>, #phage_ref keys = <$prk>\n" if ($DEBUG);
  print ".....................................................[ END of find_3prime_end ] ...................................................\n" if ($DEBUG);
  return($regions);
}

sub find_att_sites { # try to determine the attachment site for each phage region

    my ($prefix,$asmbl_id,$search_method,$ahref,$phageref,$hitref,$reref,$tRNAref,$aref,$figref,$window,$step,$hitsperwindow,$okref,$DEBUG) = @_;
    my @k = ();
    my @grase = ();
    my $key = ""; # stores keys of %phage_hash (phage number)
    my $yek = ""; # stores keys of %memberhash (feat_names)
    my $i = ""; # reused loop var
    my $phage_5prime = "" ; # predicted beginning of phage region
    my $phage_3prime = ""; # predicted end of phage region
    my $halfsize = ""; # middle of phage region coordinate
    my $Bint = ""; # leftmost end of left most integrase (either 5' or 3' of orf, but 5' relative to phage region)
    my $Eint = ""; # rightmost end of right most integrase (either 5' or 3' of orf, but 3' relative to phage region)
    my $graseref = ""; # holds pointer to @grase, containing featnames of integrases within a phage region
    my $end5 = ""; # 5' end of integrase
    my $end3 = ""; # 3' end of integrase
    my $head = ""; # 400 bp region that is outward-facing from Bint or Eint
    my $tail = ""; # region that is the remainder of the phage region plus 5000 bp for slop
    my $size = ""; # size from Bint or Eint until the end of the predicted phage region 
    my $size_head = ""; # how large of head to pull
    my $size_tail = "";  # how large of tail to pull
    my $start_head = ""; # beginning of head region to pull
    my $start_tail = ""; # beginning of tail region to pull
    my $fivetend5 = "";
    my $fivetend3 = "";
    my $threetend5 = "";
    my $threetend3 = "";
    my $tRNA5 = "";
    my $tRNA3 = "";
    my $t = "";
    my $found_att = 0;
    my $int = 0;
    my $feat_name = "";
    my $hold_attB = ""; # placeholder for first attB found (beginning of region with tRNA att site
    my $hold_attE = ""; # placeholder for first attE found (end of region with tRNA att site
    my $start = ""; # contains coordinate in genomearray to start checking for target site
    my $finish = ""; # contains coordinate in genomearray to end the checking for target site
    my $lenB = ""; # stores the length of the att site for beginning int
    my $lenE = ""; # stores the length of the att site for ending int
    my $B_region_length = "";
    my $E_region_length = "";
    my $Bextension = 15000; # set size of extension (additional sequence past initial boundaries to search for TSD)
    my $extension = "";
    my $direction = ""; #added 02/24/22 by dfouts to determine the orientation of the prophage based on which side of the phage the integrase is on (under tRNA_subs)
    my $int = ""; # added 02/24/22 by dfouts to store either $Bint or $Eint depending of which has value to avoid passing one that is empty (under tRNA_subs)
    my $Dist2end = ""; # variable to store the distance of the phage region to the end of the contig (assembly)
    my $LastWindow = ""; # variable to store the last window on the assembly

###########################
    local *beginning_sub = sub {  # create a function local to sub find_att_sites (a nested subroutine) to minimize duplication of identical code and reusing vars
	($key,$Bint,$graseref) = @_;  # passed the $key variable because the local function was unable to see it otherwise
        $extension = $Bextension;
	$size = abs($Bint - $phage_3prime);
	$size_head = 400;
	$start_head = $Bint-$size_head;
        if ($start_head < 0) { # when working with contigs, we don't want to be running off the 5' end of the contig [ 05/15/07 ]
	    $start_head = 1;
            $size_head = $Bint;
	}
	$start_tail = $phage_3prime-(round($size*0.2));
	$size_tail = $extension+(round($size*0.2));
        if (($start_tail+$size_tail) > $ahref->{$asmbl_id}->{'length'}) { # added 05/15/07 to make sure we don't run off the 3' end of a contig
	    $size_tail = $ahref->{$asmbl_id}->{'length'}-$start_tail;
            $extension = 0;
	}
	print "find_att_sites:beginning_sub: phage 5 prime end = $phage_5prime\n" if ($DEBUG);
	print "find_att_sites:beginning_sub: phage 3 prime end = $phage_3prime\n" if ($DEBUG);
	print "find_att_sites:beginning_sub: size = $size\n" if ($DEBUG);
	print "find_att_sites:beginning_sub: start_head = $start_head\n" if ($DEBUG);
	print "find_att_sites:beginning_sub: size_head = $size_head\n" if ($DEBUG);
	print "find_att_sites:beginning_sub: start_tail = $start_tail\n" if ($DEBUG);
	print "find_att_sites:beginning_sub: size_tail = $size_tail\n" if ($DEBUG);
	($found_att, $hold_attB, $hold_attE) = &get_TSD($blastbin,$prefix,$asmbl_id,$hitref,$phageref,$reref,$aref,$DEBUG,"+","phage",$search_method,$key,$ahref,$start_head,$start_tail,$size_head,$size_tail,"0","0",$figref,$window,$step,$hitsperwindow,$halfsize,$Bint,$graseref);
	### could move to get_TSD ###
	if ($found_att)  { # calculate the distance from outermost integrase coord to outermost att coord
	    $phageref->{$key}->{'att_distance'} = $Bint - $phageref->{$key}->{'left'};
	}
	### end potential code move ###
    };
############################
    local *ending_sub = sub {  # create a function local to sub find_att_sites (a nested subroutine) to minimize duplication of identical code and reusing vars
	($key,$Eint,$graseref) = @_; # passed the $key variable because the local function was unable to see it otherwise
	$extension = $Bextension;
	$size = abs($Eint - $phage_5prime);
	$size_head = 400;
        $start_head = $Eint;
        if (($start_head+$size_head) > $ahref->{$asmbl_id}->{'length'}) {  #05/15/07 to make sure we don't run off the 3' end of a contig/scaffold
	    $size_head = $ahref->{$asmbl_id}->{'length'}-$start_head;
        }
	$start_tail = $Eint-1-($size+$extension);
        if ($start_tail < 0) { # added 05/15/07 to make sure we don't run off the 5' end of a contig/scaffold
	    $start_tail = 1;
            $extension = 0;
	}
        $size_tail = $extension+(round($size*0.2));
	print "find_att_sites:ending_sub: phage 5 prime end = $phage_5prime\n" if ($DEBUG);
        print "find_att_sites:ending_sub: extension = $extension\n" if ($DEBUG);
	print "find_att_sites:ending_sub: size = $size\n" if ($DEBUG);
	print "find_att_sites:ending_sub: start_head = $start_head\n" if ($DEBUG);
	print "find_att_sites:ending_sub: size_head = $size_head\n" if ($DEBUG);
	print "find_att_sites:ending_sub: start_tail = $start_tail\n" if ($DEBUG);
	print "find_att_sites:ending_sub: size_tail = $size_tail\n" if ($DEBUG);
	($found_att, $hold_attB, $hold_attE) = &get_TSD($blastbin,$prefix,$asmbl_id,$hitref,$phageref,$reref,$aref,$DEBUG,"-","phage",$search_method,$key,$ahref,$start_head,$start_tail,$size_head,$size_tail,"0","0",$figref,$window,$step,$hitsperwindow,$halfsize,$Eint,$graseref);
	### could move to get_TSD ###
	if ($found_att)  { # calculate the distance from outermost integrase coord to outermost att coord
	    $phageref->{$key}->{'att_distance'} = $phageref->{$key}->{'right'} - $Eint;
	}
	### end potential code move ###
    };
############################ Find tRNA att sites ##############################################################################
    local *tRNA_sub = sub {  # create a function local to sub find_att_sites to find tRNA TSDs
	($key,$Bint,$Eint,$graseref) = @_; # passed the $key variable because the local function was unable to see it otherwise
        $extension = $Bextension;
	$fivetend5 = $tRNA5;
	$fivetend3 = $tRNAref->{$asmbl_id}{$tRNA5}->{'end3'};
	$threetend5 = $tRNA3;
	$threetend3 = $tRNAref->{$asmbl_id}{$tRNA3}->{'end3'};
	## start new code added 02/24/22 by dfouts ##
	if (($Bint > "") && ($Eint == "")) { # phage with Bint----->
	    $int = $Bint;
	    $direction = "+";
	    print "-> tRNA_subB:  int = Bint ($int), direction = $direction\n" if ($DEBUG);
	}
	elsif (($Bint == "") && ($Eint > "")) { # phage with <-----Eint
	    $int = $Eint;
	    $direction = "-";
	    print "-> tRNA_subE:  int = Eint ($int), direction = $direction\n" if ($DEBUG);
	}
	elsif (($Bint > "") && ($Eint > "")) { # phage with 2 ints, choose the one closest to the tRNA
	    if (defined $fivetend5) {
		$int = $Bint;
		$direction = "+";
		print "-> tRNA_subBothB:  int = Bint ($int), direction = $direction\n" if ($DEBUG);
	    }
	    else {
		$int = $Eint;
		$direction = "-";
		print "-> tRNA_subBothE:  int = Eint ($int), direction = $direction\n" if ($DEBUG);
	    }
	}
	## end new code ##    
	print "find_att_sites:tRNA_sub: fivetend5: $fivetend5\n" if ($DEBUG);
	print "find_att_sites:tRNA_sub: fivetend3: $fivetend3\n" if ($DEBUG);
	print "find_att_sites:tRNA_sub: threetend5: $threetend5\n" if ($DEBUG);
	print "find_att_sites:tRNA_sub: threetend3: $threetend3\n" if ($DEBUG);
	print "find_att_sites:tRNA_sub: halfsize: $halfsize\n" if ($DEBUG);
	# only viably possiblities are -----> at beginning or <----- at end of region
	# however, ran into a case where the integrase is at the beginning but the tRNAs are at the end (3') so phage direction was set to (-) and Eint was passed to get_TSD, but with no value!  I need to consider the direction of the integrase here
	if (defined $fivetend5)  { # if tRNA is on the 5' side of the region
	    if ($fivetend5 < $fivetend3) { # direction ---> and on 5' end of phage region
		$start_head = $fivetend5-1;
                if ($start_head < 0) { # 05/15/07
		    $start_head = 1;
		}
		$start_tail = $fivetend3 + 1;
		$size_head = ($fivetend3 - $fivetend5) + 1; # pull complete tRNA and search this first
		$size_tail = ($phage_3prime - $start_tail) + 1 + $extension;
                if (($start_tail+$size_tail) > $ahref->{$asmbl_id}->{'length'}) { # 05/15/07
		    $size_tail = $ahref->{$asmbl_id}->{'length'}-$start_tail;
		    $extension = 0;
		}
		print "find_att_sites:tRNA_sub: tRNA B----->\n" if ($DEBUG);
		print "find_att_sites:tRNA_sub: start_head = $start_head\n" if ($DEBUG);
		print "find_att_sites:tRNA_sub: size_head = $size_head\n" if ($DEBUG);
		print "find_att_sites:tRNA_sub: start_tail = $start_tail\n" if ($DEBUG);
		print "find_att_sites:tRNA_sub: size_tail = $size_tail\n" if ($DEBUG);
		($found_att, $hold_attB, $hold_attE) = &get_TSD($blastbin,$prefix,$asmbl_id,$hitref,$phageref,$reref,$aref,$DEBUG,$direction,"phage\_$reref->{$asmbl_id}{$tRNA5}->{'featname'}",$search_method,$key,$ahref,$start_head,$start_tail,$size_head,$size_tail,"0","0",$figref,$window,$step,$hitsperwindow,$halfsize,$int,$graseref);
		#### NEW CODE:  added passing $Bint to sub get_TSD so that we can check to make sure that any potential tRNA being used as a target is outside of the predicted integrase
		#01/10/2022 by dfouts
		if ($found_att == 1)  { # if match tRNA sequence, then check to see if we can expand on the match on 3' side
		    print "find_att_sites:tRNA_sub: second round expanded tRNA search...\n" if ($DEBUG);
		    $start_head = $fivetend5-1;
                    if ($start_head < 0) { # 05/15/07
			$start_head = 1;
		    }
		    $start_tail = $fivetend3 + 1 + 200;
		    $size_head = ($fivetend3 - $fivetend5) + 1 + 200; # pull complete tRNA + 200 extra nucleotides on 3' end
		    $size_tail = ($phage_3prime - $start_tail) + 1 + $extension + 200;
                    if (($start_tail+$size_tail) > $ahref->{$asmbl_id}->{'length'}) { # 05/15/07
			$size_tail = $ahref->{$asmbl_id}->{'length'}-$start_tail;
			$extension = 0;
		    }
		    print "find_att_sites:tRNA_sub: tRNA B----->\n" if ($DEBUG);
		    print "find_att_sites:tRNA_sub: start_head = $start_head\n" if ($DEBUG);
		    print "find_att_sites:tRNA_sub: size_head = $size_head\n" if ($DEBUG);
		    print "find_att_sites:tRNA_sub: start_tail = $start_tail\n" if ($DEBUG);
		    print "find_att_sites:tRNA_sub: size_tail = $size_tail\n" if ($DEBUG);
		    &get_TSD($blastbin,$prefix,$asmbl_id,$hitref,$phageref,$reref,$aref,$DEBUG,$direction,"phage\_$reref->{$asmbl_id}{$tRNA5}->{'featname'}\_r2",$search_method,$key,$ahref,$start_head,$start_tail,$size_head,$size_tail,$hold_attB,$hold_attE,$figref,$window,$step,$hitsperwindow,$halfsize,$int,$graseref);
		    $phageref->{$key}->{'target'} = $reref->{$asmbl_id}{$tRNA5}->{'featname'};
		}
		else { print "find_att_sites:tRNA_sub: NO att site found for $reref->{$asmbl_id}{$tRNA5}->{'featname'}\n" if ($DEBUG);}
		}
	    else { print "find_att_sites:tRNA_sub: $reref->{$asmbl_id}{$tRNA5}->{'featname'} facing wrong direction <---- at beginning of phage region $key\n" if ($DEBUG);}
	    }
	if ((defined $threetend5) && ($found_att == 0)) { # if tRNA is on the 3' side of the region
	    if (($threetend5 > $threetend3) && ($found_att == 0) ) { # direction <---- and not 1 tRNA and 3' end of phage region and no att site found previously
		$start_head = $threetend3; # pull complete tRNA and search this first
		$size_head = ($threetend5 - $threetend3) + 1;
		$start_tail = $phage_5prime - $extension;
                if ($start_tail < 0) { # added 05/15/07 to make sure we don't run off the 5' end of a contig/scaffold
		    $start_tail = 1;
		    $extension = 0;
		}
		$size_tail = (($start_head-1) - $start_tail) + 1;
		print "find_att_sites:tRNA_sub: tRNA <-----E\n" if ($DEBUG);
		print "find_att_sites:tRNA_sub: start_head = $start_head\n" if ($DEBUG);
		print "find_att_sites:tRNA_sub: size_head = $size_head\n" if ($DEBUG);
		print "find_att_sites:tRNA_sub: start_tail = $start_tail\n" if ($DEBUG);
		print "find_att_sites:tRNA_sub: size_tail = $size_tail\n" if ($DEBUG);
		($found_att, $hold_attB, $hold_attE) = &get_TSD($blastbin,$prefix,$asmbl_id,$hitref,$phageref,$reref,$aref,$DEBUG,$direction,"phage\_$reref->{$asmbl_id}{$tRNA3}->{'featname'}",$search_method,$key,$ahref,$start_head,$start_tail,$size_head,$size_tail,"0","0",$figref,$window,$step,$hitsperwindow,$halfsize,$int,$graseref);
		if ($found_att == 1)  { # if match tRNA sequence, then check to see if we can expand on the match on 3' side
		    print "find_att_sites:tRNA_sub: second round expanded tRNA search...\n" if ($DEBUG);
		    $start_head = $threetend3-201; # pull 200 extra on 3' side to search
		    $size_head = ($threetend5 - $threetend3) + 1 + 200;
		    $start_tail = $phage_5prime - $extension;
		    if ($start_tail < 0) { # added 05/15/07 to make sure we don't run off the 5' end of a contig/scaffold
			$start_tail = 1;
			$extension = 0;
		    }
		    $size_tail = (($start_head-1) - $start_tail) + 1;
		    print "find_att_sites:tRNA_sub: tRNA <-----E\n" if ($DEBUG);
		    print "find_att_sites:tRNA_sub: start_head = $start_head\n" if ($DEBUG);
		    print "find_att_sites:tRNA_sub: size_head = $size_head\n" if ($DEBUG);
		    print "find_att_sites:tRNA_sub: start_tail = $start_tail\n" if ($DEBUG);
		    print "find_att_sites:tRNA_sub: size_tail = $size_tail\n" if ($DEBUG);
		    &get_TSD($blastbin,$prefix,$asmbl_id,$hitref,$phageref,$reref,$aref,$DEBUG,$direction,"phage\_$reref->{$asmbl_id}{$tRNA3}->{'featname'}\_r2",$search_method,$key,$ahref,$start_head,$start_tail,$size_head,$size_tail,$hold_attB,$hold_attE,$figref,$window,$step,$hitsperwindow,$halfsize,$int,$graseref);
		    $phageref->{$key}->{'target'} = $reref->{$asmbl_id}{$tRNA3}->{'featname'};
		}
		else { print "find_att_sites:tRNA_sub: NO att site found for $reref->{$asmbl_id}{$tRNA3}->{'featname'}\n" if ($DEBUG);}
		}
	    else { print "find_att_sites:tRNA_sub: $reref->{$asmbl_id}{$tRNA3}->{'featname'} facing wrong direction ----> at end of phage region $key\n" if ($DEBUG);}
	}
	#### the following code could be moved to get_TSD and consulidated to be used for any case (beginning, end, RNA) but will need to pass both Bint and Eint to get_TSD####
	### the idea here would be A) to reduce redundant code and B) to use the distance fron int to att site as a means to filter out bogus att sites that are far away from the integrase (e.g., 10-15 kb from int)
        if ($found_att)  {
	    if (($Bint > "") && ($Eint > "")) {
		if (($Bint - $phageref->{$key}->{'left'}) == ($phageref->{$key}->{'right'} - $Eint)) { # if for some reason there are 2 integrases and they have the same distance to the att sites
		    if (defined $fivetend5)  {
			$phageref->{$key}->{'att_distance'} = $Bint - $phageref->{$key}->{'left'};
		    }
		    else {
			$phageref->{$key}->{'att_distance'} = $phageref->{$key}->{'right'} - $Eint;
		    }
		}
		elsif (($Bint - $phageref->{$key}->{'left'}) < ($phageref->{$key}->{'right'} - $Eint)) {
		    $phageref->{$key}->{'att_distance'} = $Bint - $phageref->{$key}->{'left'};
		}
		else  {
		    $phageref->{$key}->{'att_distance'} = $phageref->{$key}->{'right'} - $Eint;
		}
	    }
	    elsif (($Bint > "") && ($Eint == "")) {
		$phageref->{$key}->{'att_distance'} = $Bint - $phageref->{$key}->{'left'};
	    }
	    elsif (($Bint == "") && ($Eint > "")) {
		$phageref->{$key}->{'att_distance'} = $phageref->{$key}->{'right'} - $Eint;
	    }
	}
	#### end potential code move ####
    };
############################ End tRNA att-finding ##################################################################################
    foreach $key (sort {$a <=> $b} keys %{$phageref}) {
	print "find_att_sites:MAIN LOOP: REGION = $key\n" if ($DEBUG);
	$phage_5prime = $phageref->{$key}->{'5prime'};
	$phage_3prime = $phageref->{$key}->{'3prime'};
        print "find_att_sites:  The 5prime end before att is $phage_5prime, the 3prime end is $phage_3prime\n" if ($DEBUG);
	$size = $phage_3prime - $phage_5prime;
	print "find_att_sites:MAIN LOOP: SIZE = $size\n" if ($DEBUG);
	$halfsize = ($size*.5) + $phage_5prime;
	print "find_att_sites:MAIN LOOP: Phage_5prime = <$phage_5prime>, HALFSIZE = $halfsize, Phage_3prime = <$phage_3prime>\n" if ($DEBUG);
	$tRNA5 = undef;
	$tRNA3 = undef;
	if ($phageref->{$key}->{'class'} eq "Mu-like")  { next; } # skip if the region is a Mu-like phage (no integrase and no tRNA specificity)
	print "find_att_sites:PREINNER LOOP: END5 = $end5, END3 = $end3\n" if ($DEBUG);
	if ($phageref->{$key}->{'int_HMM'} > 0)  { # of outermost integrases

########## check for valid integrases (first and last ###############
	    $int = 1;

# new logic to find 5' and 3' integrases
            #push (@grase, shift @{$phageref->{$key}->{'integrases'}}); # gets 5' most integrase
	    push (@grase, ${$phageref->{$key}->{'integrases'}}[0]); # new code get 5' int while keeping the integrase in the original array so we can query them later
	    #push (@grase, pop @{$phageref->{$key}->{'integrases'}}); # gets 3' most integrase if present
	    push (@grase, ${$phageref->{$key}->{'integrases'}}[$#{ $phageref->{$key}->{'integrases'} }]); # new code get 3' int while keeping the integrase in the original array so we can query them later
            print "5GRASE>>>>>> $grase[0]\t$grase[1]\n" if ($DEBUG);
            if ($grase[0] eq $grase[1]) { pop @grase; } # if the same integrase listed, remove last occurance
            for $i (@grase) {
		$end5 = $hitref->{$i}->{'end5'};
                $end3 = $hitref->{$i}->{'end3'};
		print "find_att_sites:INNER LOOP: END5 = $end5, END3 = $end3\n" if ($DEBUG);
		if (($end5 < $end3) && ($end5 < $halfsize) && ($Bint == "")) { # int is 5'--->3' and is in 5' half of phage region
		    $Bint = $end5; # define the 5' most part of the int gene
		}
		elsif (($end5 > $end3) && ($end5 < $halfsize) && ($Bint == "")) { # int is 3'<---5' and is in 5' half of phage region
		    $Bint = $end3; # define the 5' most part of the int gene 
		}
		elsif (($end5 > $end3) && ($end5 > $halfsize)) { 
		    $Eint = $end5;
		}
		elsif (($end5 < $end3) && ($end5 > $halfsize)) {
		    $Eint = $end3;
		}
		print ">valid_integrases:  Bint = <$Bint>, Eint = <$Eint>\n" if ($DEBUG);
	    }
	}
        else  {
	    print "find_att_sites: NO integrase HMMs found\n" if ($DEBUG);
	}
	if ((exists $phageref->{$key}->{'tRNA'}) && ($int == 1))  { # tRNA is present and an integrase in the region, find tRNA TSD
	    print "tRNA and integrase detected in region <$key>\n" if ($DEBUG);
	    if (keys %{$phageref->{$key}->{'tRNA'}} == 1) { # if only 1 key (tRNA)
		print "SINGLE-tRNA found\n" if ($DEBUG);
		@k = keys %{$phageref->{$key}->{'tRNA'}};
		print "end5 = <$k[0]>, halfsize = $halfsize\n" if ($DEBUG);
		if ($k[0] < $halfsize) {
		    $tRNA5 = $k[0];
		    print "===> tRNA5 = <$tRNA5>, $reref->{$asmbl_id}{$tRNA5}->{'featname'}\n" if ($DEBUG);
		}
		else {
		    $tRNA3 = $k[0];
		    print "===> tRNA3 = <$tRNA3>, $reref->{$asmbl_id}{$tRNA3}->{'featname'}\n" if ($DEBUG);
		}
	    }
	    elsif (keys %{$phageref->{$key}->{'tRNA'}} > 1) { # multiple tRNAs in region
		print "MULTI-tRNAs found\n" if ($DEBUG);
		foreach $t (sort {$a <=> $b} keys %{$phageref->{$key}->{'tRNA'}}) { # get most 5' tRNA
		    print "t: $t\ttRNA: $phageref->{$key}->{'tRNA'}{$t}\n" if ($DEBUG);
		    last if (defined($tRNA5) && (($t-$tRNA5) > 1000) && ($tRNA5 < $tRNAref->{$asmbl_id}{$tRNA5}->{'end3'})); # stop looking if we have a 5' tRNA and the distance between 2 tRNAs is > 1000 bp and tRNA must be in proper orientation
		    $tRNA5 = $t if ($t < $halfsize);
		    print "tRNA5===> $tRNA5, $reref->{$asmbl_id}{$tRNA5}->{'featname'}\n" if ($DEBUG);
		}
		foreach $t (sort {$b <=> $a} keys %{$phageref->{$key}->{'tRNA'}}) { # get most 3' tRNA
		    print "t: $t\ttRNA: $phageref->{$key}->{'tRNA'}{$t}\n" if ($DEBUG);
		    last if (defined($tRNA3) && (($tRNA3-$t) > 1000) && ($tRNA3 > $tRNAref->{$asmbl_id}{$tRNA3}->{'end3'}));
		    $tRNA3 = $t if ($t > $halfsize);
		    print "tRNA3===> $tRNA3, $reref->{$asmbl_id}{$tRNA3}->{'featname'}\n" if ($DEBUG);
		}
	    }
	    &tRNA_sub($key,$Bint,$Eint,\@{$phageref->{$key}->{'integrases'}}); # added passing of array pointer so that we can check for piggyback insertions in Phage_subs.pm sub get_TSD
	    print "find_att_sites: tRNA found_att -----> tRNA5: <$tRNA5> ($reref->{$asmbl_id}{$tRNA5}->{'featname'})\ttRNA3: <$tRNA3> ($reref->{$asmbl_id}{$tRNA3}->{'featname'})\n"if ($DEBUG);
            print "find_att_sites: Bint = $Bint, Eint = $Eint\n" if ($DEBUG);
            print "find_att_sites: FOUND_att = $found_att, int = $int\n" if ($DEBUG);
	}
	if (($found_att == "0") && ($int))  { # if no tRNA att found and there are integrases in region, proceed

	    ########## if no tRNA found, check for TSDs ###############
	    
	    print "check_terminal_RNA5:  tRNA5 = <$tRNA5>, phage5prime = <$phageref->{$key}->{'5prime'}>\n" if ($DEBUG);
	    print "check_terminal_RNA3:  tRNA3 = <$tRNA3>, phage3prime = <$phageref->{$key}->{'3prime'}>\n" if ($DEBUG);
	    
	    ### first, remove any RNA that is at the ends of the phage region before moving on with checking for other potential att sites 01/19/22 ###
	    if ($tRNA5 == $phageref->{$key}->{'5prime'}) { #checking if the 5' end of the phage is an RNA.  If so, then remove it since it is not the att site 01/19/22 dfouts
		# fallback to nearest good 5' ORF
		print "remove_tRNA_5prime: last5ORF = $phageref->{$key}->{'last5ORF'}, 5prime end = $phageref->{$key}->{'5prime'}, RNA = $reref->{$asmbl_id}{$tRNA5}->{'featname'}\n" if ($DEBUG);
		if ($phageref->{$key}->{'last5ORF'} < $hitref->{$reref->{$asmbl_id}{$phageref->{$key}->{'last3ORF'}}->{'featname'}}->{'end3'}) { #if 5'--->
		    $phageref->{$key}->{'3prime'} = $phageref->{$key}->{'last5ORF'}; #07/13/22
		}
		else { # if 3'<---
		    $phageref->{$key}->{'3prime'} = $hitref->{$reref->{$asmbl_id}{$phageref->{$key}->{'last3ORF'}}->{'featname'}}->{'end3'};
		}
	    }
	    elsif ($tRNA3 == $phageref->{$key}->{'3prime'}) { # 3prime RNA
		print "remove_tRNA_3prime: last3ORF = $phageref->{$key}->{'last3ORF'}, 3prime end = $phageref->{$key}->{'3prime'}, RNA = $reref->{$asmbl_id}{$tRNA3}->{'featname'}\n" if ($DEBUG);
		#$phageref->{$key}->{'3prime'} =  $phageref->{$key}->{'last3ORF'}; # reset the phage 3' end to the last 3' ORF
		if ($phageref->{$key}->{'last3ORF'} < $hitref->{$reref->{$asmbl_id}{$phageref->{$key}->{'last3ORF'}}->{'featname'}}->{'end3'}) { #if --->3'
		    $phageref->{$key}->{'3prime'} = $hitref->{$reref->{$asmbl_id}{$phageref->{$key}->{'last3ORF'}}->{'featname'}}->{'end3'}; #07/13/22
		}
		else { # if <---5'
		    $phageref->{$key}->{'3prime'} = $phageref->{$key}->{'last3ORF'};
		}
	    }
	    ### end removing of terminal RNAs that are not att sites ###
	    
	    if (($Bint > "") && ($Eint > "")) { # If there are more than one putative integrases present, check both and pick best att (1) target, 2)longest length)
		print "find_att_sites: Phage region $key has at least 2 integrases, at $Bint and $Eint\n" if ($DEBUG);
		print "find_att_sites: Checking Beginning integrase (1 of 2) for phage region $key begins at Beginning $Bint A\n" if ($DEBUG);
		&beginning_sub($key,$Bint,\@{$phageref->{$key}->{'integrases'}});
		$lenB = length($hitref->{"attL_$asmbl_id\_$key"}->{'annotation'});
		if (!defined ($phageref->{$key}->{'target'})) {
		    if ($lenB > 0) { delete $hitref->{"attL_$asmbl_id\_$key"}; } # removing old data so new data can be stored
		    print "find_att_sites: Checking Ending integrase (2 of 2) for phage region $key begins at End $Eint B\n" if ($DEBUG);
		    delete $hitref->{"attL_$asmbl_id\_$key"}; # delete data from beginning_sub
		    $phageref->{$key}->{'target'} = "";
		    print "hold_attB = $hold_attB\thold_attE = $hold_attE\n" if ($DEBUG);
		    delete $phageref->{$key}->{'memberhash'}{$hold_attB};  # remove old attB information
		    delete $phageref->{$key}->{'memberhash'}{$hold_attE};  # remove old attE information
                    $B_region_length = abs($hold_attB-$hold_attE);
		    &ending_sub($key,$Eint,\@{$phageref->{$key}->{'integrases'}});
                    $lenE = length($hitref->{"attL_$asmbl_id\_$key"}->{'annotation'});
                    $E_region_length = abs($hold_attB-$hold_attE);
                    print "Blength = $B_region_length\tElength = $E_region_length\n" if ($DEBUG);
		    if (($B_region_length > $E_region_length) || ($lenB > $lenE))  { # if length of phage region with B int is > or if the length of the att site is >
			print "find_att_sites: Going with beginning integrase (lenB = $lenB, lenE = $lenE) when 2 present\n" if ($DEBUG);
			delete $hitref->{"attL_$asmbl_id\_$key"}; # delete data from ending_sub
			$phageref->{$key}->{'target'} = "";
                        print "hold_attB = $hold_attB\thold_attE = $hold_attE\n" if ($DEBUG);
                        delete $phageref->{$key}->{'memberhash'}{$hold_attB};  # remove old attB information
			delete $phageref->{$key}->{'memberhash'}{$hold_attE};  # remove old attE information
			&beginning_sub($key,$Bint,\@{$phageref->{$key}->{'integrases'}});
		    }
		}
	    }
	    elsif (($Bint > "") && ($Eint == "")) { # If only one at the beginning, use it
		print "find_att_sites: The integrase to study for phage region $key begins at Beginning $Bint C\n" if ($DEBUG);
		&beginning_sub($key,$Bint,\@{$phageref->{$key}->{'integrases'}});
	    }
	    elsif (($Bint == "") && ($Eint > "")) { # If only one at the end, use it
		print "find_att_sites: The integrase to study for phage region $key begins at End $Eint D $Bint\n" if ($DEBUG);
		&ending_sub($key,$Eint,\@{$phageref->{$key}->{'integrases'}});
	    }
	    else { print "find_att_sites: No integrases for phage region $key\n" if (($DEBUG) && ($int == 0)); }
	}

	if ($found_att)  { # if an att site found and the new boundaries are beyond the original phage boundaries, then add ORFs
            print "FOUND_att!\n" if ($DEBUG);
	    if (($phageref->{$key}->{'left'} < $phageref->{$key}->{'5prime'}) || ($phageref->{$key}->{'right'} > $phageref->{$key}->{'3prime'}))  {
### Need to check if extended region overlaps previous or next phage region before extending ###
		print "EXTENDING REGION!\n" if ($DEBUG);
		print "5prime = <$phageref->{$key}->{'5prime'}>, 3prime = <$phageref->{$key}->{'3prime'}>\n" if ($DEBUG);
		my $rv = &add_ORFs_to_region($key,$asmbl_id,$phageref,$hitref,$reref,$aref,$okref,*BCNFILE,$DEBUG);
                if ($rv)  { # removing att site information from the phage hash
		    print "find_att_sites: House-keeping genes found, undoing att info...\n" if ($DEBUG);
		    $phageref->{$key}->{'5prime_att'} = "";
                    $phageref->{$key}->{'3prime_att'} = "";
                    $phageref->{$key}->{'left'} = "";
		    $phageref->{$key}->{'right'} = "";
                    $phageref->{$key}->{'att_distance'} = "";
                    $phageref->{$key}->{'target'} = "";
                    print "EXTENDING_ORFS: Direction = <$phageref->{$key}->{'direction'}>\n" if ($DEBUG);
		    $phageref->{$key}->{'direction'} = "+"; #reset to default (+)
		    if ($phageref->{$key}->{'5prime'} > $phageref->{$key}->{'3prime'}) {
			print "changing directions\n" if ($DEBUG);
			my $hold_5prime = $phageref->{$key}->{'5prime'};
			my $hold_3prime = $phageref->{$key}->{'3prime'};
			$phageref->{$key}->{'5prime'} = $hold_3prime;
			$phageref->{$key}->{'3prime'} = $hold_5prime;
		    }
                    $found_att = 0;
		    $int = 0; # added 09/19/2013 dfouts
		    $Bint = ""; # added 09/19/2013 dfouts
		    $Eint = ""; # added 09/19/2013 dfouts
		    $end5 = ""; # added 09/19/2013 dfouts
		    $end3 = ""; # added 09/19/2013 dfouts
		    @grase = (); # added 09/19/2013 dfouts
                    next;
		}
	    }
	    if ($phageref->{$key}->{'att_distance'} < 0)  { # if att_distance is negative, make zero (case where real int is not translated and found one outside of the region
                print "find_att_sites:  att_distance is $phageref->{$key}->{'att_distance'}!, so setting att_distance to zero\n" if ($DEBUG);
		$phageref->{$key}->{'att_distance'} = 0;
	    }
	}
	############### Check if phage region is near the end of a contig ###############
	# Regardless of whether or not an att site was found, check if any phage regions are close to the end of a contig (the window size),
	# has hits in the last window then fill to the begin and/or end with adjusted coords.

	my $hold_5prime = $phageref->{$key}->{'5prime'}; # store the existing 5 and 3 prime ends
	my $hold_3prime = $phageref->{$key}->{'3prime'};
	my $hold_target = $phageref->{$key}->{'target'};
	my $hold_5prime_att = $phageref->{$key}->{'5prime_att'};
	my $hold_3prime_att = $phageref->{$key}->{'3prime_att'};
	my $hold_attdist = $phageref->{$key}->{'att_distance'};
	my $current_step = nlowmult($step, $hold_5prime); # new variable to store the current step position in figref 
	print "hold_5prime = <$hold_5prime>, hold_3prime = <$hold_3prime>, hold_5primeatt = <$hold_5prime_att>, hold_3primeatt = <$hold_3prime_att>, current_step = <$current_step>\n" if ($DEBUG);
	my $rv = ""; #clear variable just in case
	
	if (($phageref->{$key}->{'5prime'} < $window) && ($figref->{$current_step}->{'counts'} >= 1) || ($hold_5prime == "1")) { # if at or near beginning of contig
	    print "  ->-> Region $key is near the beginning of $asmbl_id: [$figref->{$current_step}->{'counts'} hits]\n" if ($DEBUG);
	    $phageref->{$key}->{'left'} = "1"; # set phage start to the beginning of the contig
	    $phageref->{$key}->{'right'} = $phageref->{$key}->{'3prime'}; # set the right end to the existing 3prime end coordinate. Need to change this when the region was extended above and the 3prime is the same as one of the attprimes
	    #$phageref->{$key}->{'ORF5'} = $reref->{$asmbl_id}{$aref->["0"]}->{'featname'}; # adjust ORF5 to first ORF on contig
	    $phageref->{$key}->{'5prime'} = $phageref->{$key}->{'left'}; # go ahead and set the new 5prime end (can reverse later if bad com_names)
	    $phageref->{$key}->{'5prime_att'} = ""; # clear out att coords
	    $phageref->{$key}->{'3prime_att'} = "";
	    $phageref->{$key}->{'target'} = ""; # clear out target
	    $phageref->{$key}->{'att_distance'} = ""; # clear out att_distance since on att site now
	    print "    -> New LEFT end set to <$phageref->{$key}->{'left'}>, new ORF5: <$phageref->{$key}->{'ORF5'}>\n" if ($DEBUG); 
	    $rv = &add_ORFs_to_region($key,$asmbl_id,$phageref,$hitref,$reref,$aref,$okref,*BCNFILE,$DEBUG);
	    # added this here so that a failed 5' end extension can be fixed
	    if ($rv) { #if bad com_names preventing extension of the phage, restore the previous 5 and 3 prime ends
		print "=> (begin_contig) restoring previous ends to phage region $key...\n" if ($DEBUG);
		print "=> (begin_contig) 5prime = <$hold_5prime>\n" if ($DEBUG);
		print "=> (begin_contig) 3prime = <$hold_3prime>\n" if ($DEBUG);
		print "=> (begin_contig) 5prime_att = <$hold_5prime_att>\n" if ($DEBUG);
		print "=> (begin_contig) 3prime_att = <$hold_3prime_att>\n" if ($DEBUG);
		print "=> (begin_contig) target = <$hold_target>\n" if ($DEBUG);
		print "=> (begin_contig) att_dist = <$hold_attdist>\n" if ($DEBUG);
		$phageref->{$key}->{'5prime'} = $hold_5prime;
		$phageref->{$key}->{'5prime_att'} = $hold_5prime_att;
		$phageref->{$key}->{'3prime_att'} = $hold_3prime_att;
		$phageref->{$key}->{'target'} = $hold_target; # not sure if this need to be here
		$phageref->{$key}->{'att_distance'} = $hold_attdist; # not sure if this need to be here
	    }
	}
	# set new hold vars for any new 5' end coords
	$hold_5prime = $phageref->{$key}->{'5prime'};
	$hold_5prime_att = $phageref->{$key}->{'5prime_att'};
	$hold_3prime_att = $phageref->{$key}->{'3prime_att'};
	
	$Dist2end = $assembly_hash{$asmbl_id}->{'length'} - $window; # find the distance of the phage region to the end of the contig (assembly)
	#$LastWindow = nlowmult($window, $assembly_hash{$asmbl_id}->{'length'}); # the last window on the assembly
	$LastWindow = nhimult($step, $Dist2end); # changed 06/01/2022 to not be so close to the end (e.g., 60000 for a contig that is ~61000 bp).  Only one gene in that region.
	print "->->->-> Dist2end = $Dist2end, LastWindow = $LastWindow\n" if ($DEBUG);
	
	if (($phageref->{$key}->{'3prime'} > $Dist2end) && ($figref->{$LastWindow}->{'counts'} >= 1)) { # if near end of contig
	    print "  ->-> Region $key is near a contig end: $figref->{$LastWindow}->{'counts'} hits]\n" if ($DEBUG);
	    $phageref->{$key}->{'left'} = $phageref->{$key}->{'5prime'}; # set the left end to the existing 5prime end coordinate
	    $phageref->{$key}->{'right'} = $assembly_hash{$asmbl_id}->{'length'}; # set phage start to the end of the contig
	    #$phageref->{$key}->{'ORF3'} = $reref->{$asmbl_id}{$aref->[$#{$aref}]}->{'featname'}; # adjust ORF3 to last ORF on contig
	    $phageref->{$key}->{'3prime'} = $phageref->{$key}->{'right'}; # go ahead and set the new 3prime end (can reverse later if bad com_names)
	    $phageref->{$key}->{'5prime_att'} = ""; # clear out att coords
	    $phageref->{$key}->{'3prime_att'} = "";
	    $phageref->{$key}->{'target'} = ""; # clear out target
	    $phageref->{$key}->{'att_distance'} = ""; # clear out att_distance since on att site now
	    print "    -> New RIGHT end set to <$phageref->{$key}->{'right'}>, new ORF3: <$phageref->{$key}->{'ORF3'}>\n" if ($DEBUG); 
	    $rv = &add_ORFs_to_region($key,$asmbl_id,$phageref,$hitref,$reref,$aref,$okref,*BCNFILE,$DEBUG);
	    if ($rv) { #if bad com_names preventing extension of the phage, restore the previous 5 and 3 prime ends
		print "=> (end_contig) restoring previous ends to phage region $key...\n" if ($DEBUG);
		print "=> (end_contig) 5prime = <$hold_5prime>\n" if ($DEBUG);
		print "=> (end_contig) 3prime = <$hold_3prime>\n" if ($DEBUG);
		print "=> (end_contig) 5prime_att = <$hold_5prime_att>\n" if ($DEBUG);
		print "=> (end_contig) 3prime_att = <$hold_3prime_att>\n" if ($DEBUG);
		print "=> (end_contig) target = <$hold_target>\n" if ($DEBUG);
		print "=> (end_contig) att_dist = <$hold_attdist>\n" if ($DEBUG);
		$phageref->{$key}->{'3prime'} = $hold_3prime;
		$phageref->{$key}->{'5prime_att'} = $hold_5prime_att;
		$phageref->{$key}->{'3prime_att'} = $hold_3prime_att;
		$phageref->{$key}->{'target'} = $hold_target;
		$phageref->{$key}->{'att_distance'} = $hold_attdist;
	    }
	}
	print "> resulting phage region is 5prime = <$phageref->{$key}->{'5prime'}>, 3prime = <$phageref->{$key}->{'3prime'}>, 5primeatt = <$phageref->{$key}->{'5prime_att'}> 3primeatt = <$phageref->{$key}->{'3prime_att'}>\n" if ($DEBUG);
	
	#################
	
	$found_att = 0; # reset so not found
	$int = 0;
	$Bint = "";
	$Eint = "";
	$end5 = "";
	$end3 = "";
	$Dist2end = "";
	$LastWindow = "";
	@grase = ();
    }
}

sub print_regions {
    
    my ($write_dir,$asmbly_status,$asmbl_id,$hitsperwindow,$window,$step,$evalue,$max,$phageref,$aref,$ahref,$hitref,$reref,$dbref,$DEBUG) = @_; 
    my $large = 0;
    my $totalphagebp = 0;
    my $percent = 0;
    my $pico = 0;  # size < 2000 bp
    my $nano = 0;  # size >= 2000, < 5000 bp
    my $micro = 0; # size >= 5000 < 10000 bp
    my $small = 0;
    my $medium = 0;
    my $size = 0;
    my $size_att = 0;
    my $sizer = 0;
    my $totalpicobp = 0;
    my $totalnanobp = 0;
    my $totalmicrobp = 0;
    my $totalsmallbp = 0;
    my $totalmediumbp = 0;
    my $percentsmall = 0;
    my $percentmedium = 0;
    my $hmm_hit_names = "";
    my $end5 = "";
    my $end3 = "";
    my $featname = "";
    my $printfeatname = "";
    my $annotation = "";
    my $header = "";
    my $con_header = "";
    my $region_header = "";
    my $region2_header = "";
    my $region_tag = "";
    #my $seq_tag = "";
    my $key = ""; 
    my $yek = ""; 
    my $HMM_name = "";
    my $small_open = 0;
    my $med_open = 0;
    my $large_open = 0;
    my $phage_end5 = "";
    my $phage_end3 = "";
    my $phage_direction = "";
    my $label = "";
    my $is_piggy = "";
    my $bugtag = $ahref->{$asmbl_id}->{'genus'} . "_" . $ahref->{$asmbl_id}->{'species'} . "_" . $ahref->{$asmbl_id}->{'strain'} . "_" . $asmbl_id;
    my $prefix = "$write_dir";
    my $region_id = "";
    #my $seqfileprefix = "$prefix/$bugtag";
    #my $pepfileprefix = "$prefix/$bugtag";
    #my $confileprefix = "$prefix/$bugtag";
    my $print_files = undef;
    my $m = ""; # variable for storing region before deleted regions

    ## start of local subroutine to print regions to files ##
    
    local *print_sub = sub { # universal print subroutine within print_regions
	my ($key, $FH, $label, $is_piggy, $count, $type, $phage_end5, $phage_end3, $phage_direction) = @_;
	my $hold_locus = "";
	my $print5 = 0;
	my $seq = "";
	my $pep = "";
	my $seq_header = "";
	my $i = "";
	my $B_pos = "";
	my $E_pos = "";
	my %Best_hits = (); # reset with each call to print_sub
	my $k = "";
	my $hit = "";
	my $best = "";
	my $shortag = "";
	my $totalhits = 0;
	my $role_id = "";
	my $ftype = ""; #variable to define feat_type for Linear_Display script
	my $RNA_type = ""; # variable to hold the name of the non-coding RNA (tRNA, tmRNA or ncRNA);
	my $flank5 = ""; # variable to hold the end5 coordinate of the 5' flanking gene for printing to PFPR_tab.txt
	my $flank3 = ""; # variable to hold the end5 coordinate of the 3' flanking gene for printing to PFPR_tab.txt
	my $flank5_gene = ""; # variable to store the featname of the 5' flanking gene
	my $flank3_gene = ""; # variable to store the featname of teh 3' flanking gene
	my $contig_end5 = ""; # variable to store the contig-centric phage boundaries (5' to 3')
	my $contig_end3 = ""; # variable to store the contig-centric phage boundaries (5' to 3')
	my $attR = ""; # variable to store the sequence of attR
	my $attL = ""; # variable to store the sequence of attL
	my $target = ""; # variable to store the target feat_name
	my $region_size = ""; # variable to store the size of the phage region (with or without att sizes)
	
	$phageref->{$key}->{'ORFcnt'} = 0; # set the ORF count to zero
	
	if ($count > 1)  {
	    print $FH "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n";
	    select((select($FH), $= = 67, $- = 0)[0]);
	}
	else  {
	    select((select($FH), $= = 67)[0]);
	}

	## set the boundaries of the phage region relative to the assembly/contig so that we don't have to worry about diretionality when looping through contig coords 5- to 3'

	if ($phage_end5 < $phage_end3) { # had to fix this by adding the check that phage_end5 < phage_end3 instead of direction check 7/11/22
	    $contig_end5 = $phage_end5;
	    $contig_end3 = $phage_end3;
	}
	else {
	    $contig_end5 = $phage_end3;
	    $contig_end3 = $phage_end5;
	}
	print "PRINT_SUB: DIRECTION: <$phage_direction>\n" if ($DEBUG);
	print "PRINT_SUB: contig_end5 = <$contig_end5>, contig_end3 = <$contig_end3>\n" if ($DEBUG);
	print "PRINT_SUB: phage_end5 = <$phage_end5>, phage_end3 = <$phage_end3>\n" if ($DEBUG);

	## populate the $Best_hits has so that we can determine the best database hit for this phage region
	
	foreach $k (sort {$a <=> $b} keys %{$phageref->{$key}->{'memberhash'}}) { # pre loop to get the Best Hit
	    $featname = $phageref->{$key}->{'memberhash'}{$k}->{'featname'};
	    print "PRINT_SUB: BEST_HITS: pre_featname = $featname\n" if ($DEBUG);
	    if ((exists $hitref->{$featname}->{'annotation'}) && ($featname !~ /^att/) && ($featname !~ /^tm?RNA/) && ($featname !~ /^ncRNA/)){
		print "PRINT_SUB: ===> Region $key contig_coords: $contig_end5 -> $contig_end3 (phage_coords: $phage_end5 to $phage_end3), annotation = $hitref->{$featname}->{'annotation'}, end5 = $hitref->{$featname}->{'end5'}\n" if ($DEBUG);
		if (($hitref->{$featname}->{'end5'} >= $contig_end5) && ($hitref->{$featname}->{'end5'} <= $contig_end3)) {
		    $Best_hits{$hitref->{$featname}->{'phage'}}->{'counts'}++;
		    $totalhits++;
		    print "PRINT_SUB: BEST_HITS:    post_featname: $featname\t$hitref->{$featname}->{'annotation'}\t$hitref->{$featname}->{'phage'}\n" if ($DEBUG);
		}
	    }
	}
	
	## sort by phage with the greatest number of hits to get the phage name, then end the loop
	
	foreach $hit (sort {$Best_hits{$b}->{'counts'} <=> $Best_hits{$a}->{'counts'}} keys %Best_hits)  {
	    $best = $hit;
	    #print "Region $key Best Hit = $best, ($Best_hits{$best}->{'counts'})\n";
	    $shortag = $best;
	    $shortag =~ s/^PHAGE_//;
	    last;
	}
	### printing out initial phage region info
	# new, more streamlined code! 07/07/22

	if ($size_att > 1) { # att found
	    print "PRINT_SUB: print_regions: att found ...\n" if ($DEBUG);
	    $attR = "$hitref->{\"attR_$asmbl_id\_$key\"}->{'annotation'}";
	    $attL = "$hitref->{\"attL_$asmbl_id\_$key\"}->{'annotation'}";
	    $region_size = $size_att;
	    
	    if (length($phageref->{$key}->{'target'}) > 0) { # if there is a target (assumes that also has an att site)
		print "PRINT_SUB: print_regions: target detected...\n" if ($DEBUG);
		$target = $phageref->{$key}->{'target'};
	    }
	    else { # no target was detected
		print "print_regions: no target detected but att found ...\n" if ($DEBUG);
		$target = "N.D.";
	    }
	}
	else { # no att found
	    $attR = "N.D.";
	    $attL = "N.D.";
	    $region_size = $size;
	    $target = "N.D."; #need to add this since no att = no target
	}

	print "$label $type region $region_id is from $phage_end5 to $phage_end3  and is $region_size bp in size, gc%: $phageref->{$key}->{'gc'}%, TARGET: $target, Piggback: $is_piggy. Best hit: $best ($Best_hits{$best}->{'counts'}/$totalhits)\n";

	$region_header = "$label $type region $region_id is from $phage_end5 to $phage_end3 and is $region_size bp in size, gc%: $phageref->{$key}->{'gc'}%, TARGET: $target, Piggyback insertion: $is_piggy.";

	print TAB "$region_id\t$asmbl_idhash{$asmbl_id}->{'genomesize'}\t$ahref->{$asmbl_id}->{'gc'}\t$phage_end5\t$phage_end3\t$region_size\t$label\t$is_piggy\t$type\t$attR\t$attL\t$target\t$phageref->{$key}->{'gc'}\t$shortag\t";

	print FRAGFILE "$frag_line_num\t$region_id\t$phage_end5\t$phage_end3\n";
	
#	if (length($phageref->{$key}->{'target'}) > 0) { # if there is a target (assumes that also has an att site)
#	    print "print_regions: target detected...\n" if ($DEBUG);
#	    print "$label $type region $region_id is from $phageref->{$key}->{'5prime'} to $phageref->{$key}->{'3prime'} or ($phageref->{$key}->{'5prime_att'} to $phageref->{$key}->{'3prime_att'}) and is $size ($size_att) bp in size, gc%: $phageref->{$key}->{'gc'}%, TARGET: $phageref->{$key}->{'target'}, Piggback: $is_piggy. Best hit: $best ($Best_hits{$best}->{'counts'}/$totalhits)\n";
#	    $region_header = "$label $type region $region_id is from $phageref->{$key}->{'5prime'} to $phageref->{$key}->{'3prime'} or ($phageref->{$key}->{'5prime_att'} to $phageref->{$key}->{'3prime_att'}) and is $size ($size_att) bp in size, gc%: $phageref->{$key}->{'gc'}%, TARGET: $phageref->{$key}->{'target'}, Piggyback insertion: $is_piggy.";
	    ##
#	    print TAB "$region_id\t$asmbl_idhash{$asmbl_id}->{'genomesize'}\t$ahref->{$asmbl_id}->{'gc'}\t$phageref->{$key}->{'5prime_att'}\t$phageref->{$key}->{'3prime_att'}\t$size_att\t$label\t$is_piggy\t$type\t$hitref->{\"attR_$asmbl_id\_$key\"}->{'annotation'}\t$hitref->{\"attL_$asmbl_id\_$key\"}->{'annotation'}\t$phageref->{$key}->{'target'}\t$phageref->{$key}->{'gc'}\t$shortag\t";
	    ##
#	    print FRAGFILE "$frag_line_num\t$region_id\t$phageref->{$key}->{'5prime_att'}\t$phageref->{$key}->{'3prime_att'}\n";
#	}
#	elsif ($size_att > 1)  {
#	    print "print_regions: no target detected but att found ...\n" if ($DEBUG);
#	    print "$label $type region $key is from $phageref->{$key}->{'5prime'} to $phageref->{$key}->{'3prime'} or ($phageref->{$key}->{'5prime_att'} to $phageref->{$key}->{'3prime_att'}) and is $size ($size_att) bp in size, gc%: $phageref->{$key}->{'gc'}%, Piggyback: $is_piggy. Best hit: $best ($Best_hits{$best}->{'counts'}/$totalhits)\n";
#	    $region_header = "$label $type region $region_id is from $phageref->{$key}->{'5prime'} to $phageref->{$key}->{'3prime'} or ($phageref->{$key}->{'5prime_att'} to $phageref->{$key}->{'3prime_att'}) and is $size ($size_att) bp in size, gc%: $phageref->{$key}->{'gc'}%, Piggyback insertion: $is_piggy.";
	    
	    ##
#	    print TAB "$region_id\t$asmbl_idhash{$asmbl_id}->{'genomesize'}\t$ahref->{$asmbl_id}->{'gc'}\t$phageref->{$key}->{'5prime_att'}\t$phageref->{$key}->{'3prime_att'}\t$size_att\t$label\t$is_piggy\t$type\t$hitref->{\"attR_$asmbl_id\_$key\"}->{'annotation'}\t$hitref->{\"attL_$asmbl_id\_$key\"}->{'annotation'}\tN.D.\t$phageref->{$key}->{'gc'}\t$shortag\t";
	    ##
#	    if ($phage_direction eq "+") {
#		print FRAGFILE "$frag_line_num\t$region_id\t$phageref->{$key}->{'5prime_att'}\t$phageref->{$key}->{'3prime_att'}\n";
#	    }
#	    else {
#		print FRAGFILE "$frag_line_num\t$region_id\t$phageref->{$key}->{'3prime_att'}\t$phageref->{$key}->{'5prime_att'}\n";
#	    }
#	}
#	else  {
#	    print "print_regions: no target or att detected...\n" if ($DEBUG);
#	    if ($phage_direction eq "+") { # phage ----> (added 07/06/22)
#		print "$label $type region $region_id is from $phageref->{$key}->{'5prime'} to $phageref->{$key}->{'3prime'} and is $size bp in size, gc%: $phageref->{$key}->{'gc'}%, Piggyback: $is_piggy. Best hit: $best ($Best_hits{$best}->{'counts'}/$totalhits)\n";
#		$region_header = "$label $type region $region_id is from $phageref->{$key}->{'5prime'} to $phageref->{$key}->{'3prime'} and is $size bp in size, gc%: $phageref->{$key}->{'gc'}%, Piggyback insertion: $is_piggy.";
		
#		print TAB "$region_id\t$asmbl_idhash{$asmbl_id}->{'genomesize'}\t$ahref->{$asmbl_id}->{'gc'}\t$phageref->{$key}->{'5prime'}\t$phageref->{$key}->{'3prime'}\t$size\t$label\t$is_piggy\t$type\tN.D.\tN.D.\tN.D.\t$phageref->{$key}->{'gc'}\t$shortag\t";
	
#		print FRAGFILE "$frag_line_num\t$region_id\t$phageref->{$key}->{'5prime'}\t$phageref->{$key}->{'3prime'}\n";
#	    }
#	    else { # phage <----- (added 07/06/22)
#		print "$label $type region $region_id is from $phageref->{$key}->{'3prime'} to $phageref->{$key}->{'5prime'} and is $size bp in size, gc%: $phageref->{$key}->{'gc'}%, Piggyback: $is_piggy. Best hit: $best ($Best_hits{$best}->{'counts'}/$totalhits)\n";
#		$region_header = "$label $type region $region_id is from $phageref->{$key}->{'3prime'} to $phageref->{$key}->{'5prime'} and is $size bp in size, gc%: $phageref->{$key}->{'gc'}%, Piggyback insertion: $is_piggy.";
		
#		print TAB "$region_id\t$asmbl_idhash{$asmbl_id}->{'genomesize'}\t$ahref->{$asmbl_id}->{'gc'}\t$phageref->{$key}->{'3prime'}\t$phageref->{$key}->{'5prime'}\t$size\t$label\t$is_piggy\t$type\tN.D.\tN.D.\tN.D.\t$phageref->{$key}->{'gc'}\t$shortag\t";

#		print FRAGFILE "$frag_line_num\t$region_id\t$phageref->{$key}->{'3prime'}\t$phageref->{$key}->{'5prime'}\n";
#	    }
#	}
	$region2_header = "Best database match ($Best_hits{$best}->{'counts'}/$totalhits) to $dbref->{$best}->{'name'} \[$best\], taxonomy: $dbref->{$best}->{'taxon'}";
	$header = "$ahref->{$asmbl_id}->{'title'} [asmbl_id: $asmbl_id], $ahref->{$asmbl_id}->{'length'}, gc%: $ahref->{$asmbl_id}->{'gc'}%";
	$end5 = "";
	$end3 = "";
	$printfeatname = "";
	$annotation = "";
	write if ($DEBUG);
	write $FH;

	#### add display_id for each prophage region in the genbank file ####
	#change out vars $phage_end5 and $phage_end3 for $contig_end5 and $contig_end3 that is unidirectional (5' to 3') based on $phage_direction
	
	foreach $yek (sort {$a <=> $b} keys %{$phageref->{$key}->{'memberhash'}}) {
	    $featname = $phageref->{$key}->{'memberhash'}{$yek}->{'featname'};
	    print "PRINTSUB_memberhash_loop: $phage_end5 [ $yek ] $phage_end3, featname = <$featname>\n" if ($DEBUG);
	    #$att_featname = $featname . "_" . $asmbl_id . "_" . $region_tag;
	    $role_id = "2000"; #default reset for each new ORF (dfouts) 07/06/2017
	    if ((($yek >= $contig_end5) && ($yek <= $contig_end3)) || ($featname eq $target)) { # if gene is within phage boundaries +/- 1 gene or is the target
		$end5 = $yek;
		$end3 = $hitref->{$phage_hash{$key}->{'memberhash'}{$end5}->{'featname'}}->{'end3'};
		if (($asmbly_status) && ($featname ne $phageref->{$key}->{'target'}) && ($phageref->{$key}->{'type'} eq "prophage"))  { # don't pull seqs for target
		    $seq = &quickcut($ahref->{$asmbl_id}->{'sequence'}, $end5, $end3);
		 
		    ### Brought this code back to fix tags in .seq and .pep files so we know which genome they come from and can use with other in-house scripts that require unique tags
		    if (exists $hitref->{$featname}->{'genome_tag'}) { # new code 4/6/2010 so we have a clean feat_name-tag printed (ie giBLA-ABPI00000000_PFPR01)
			$seq_header = "gi" . $featname . "-" . $hitref->{$featname}->{'genome_tag'} . "_" . $asmbl_id . "_" . $region_tag . " " . $hitref->{$featname}->{'com_name'}; # combodb format
		    }
		    else {
			$seq_header = $featname . "_" . $region_tag . " " . $hitref->{$featname}->{'com_name'}; # removed the $asmbl_id from header to avoid redundancy (e.g., G952_RS11985-NZ_KE702177_NZ_KE702177_PFPR09)
		    }
		    ##### add features to genbank file here? #
		    #
		    #
		    #####
		    if (($phageref->{$key}->{'memberhash'}{$yek}->{'featname'} !~ /^att/) && ($featname !~ /RNA/)) { #skip over non-gene features and targets for integration
			$phageref->{$key}->{'ORFcnt'}++; # count the number of ORFs per phage region
			&printFasta(*SEQFILE, $seq_header, $seq);
			($pep) = &dna2pep($seq);
			&printFasta(*PEPFILE, $seq_header, $pep);
			if (exists $hitref->{$featname}->{'Small Terminase'})  {
			    &printFasta(*SMALLTERMFILE, $seq_header, $pep);
			    $role_id = "1005";
			}
			elsif (exists $hitref->{$featname}->{'Large Terminase'})  {
			    &printFasta(*LARGETERMFILE, $seq_header, $pep);
			    $role_id = "1006";
			}
			elsif (exists $hitref->{$featname}->{'Portal'})  {
			    &printFasta(*PORTALFILE, $seq_header, $pep);
			    $role_id = "1007";
			}
			elsif (exists $hitref->{$featname}->{'Capsid Protease'})  {
			    &printFasta(*PROTEASEFILE, $seq_header, $pep);
			    $role_id = "1008";
			}
			elsif (exists $hitref->{$featname}->{'Major Capsid'})  {
			    &printFasta(*CAPSIDFILE, $seq_header, $pep);
			    $role_id = "1009";
			}
			elsif (exists $hitref->{$featname}->{'serine'})  {
			    &printFasta(*SERINERECOMBFILE, $seq_header, $pep);
			    $role_id = "1002";
			}
			elsif (exists $hitref->{$featname}->{'tyrosine'})  {
			    &printFasta(*TYROSINERECOMBFILE, $seq_header, $pep);
			    $role_id = "1002";
			}
			elsif (exists $hitref->{$featname}->{'Toxins'})  {
			    &printFasta(*TOXINSFILE, $seq_header, $pep);
			    $role_id = "1025";
			}
			elsif (exists $hitref->{$featname}->{'RBP'})  {
			    &printFasta(*RBPSFILE, $seq_header, $pep);
			    $role_id = "1035";
			    # may want to add capacity to annotate Dit, Tal and Whisker proteins in the future.
			}
			if (exists $hitref->{$featname}->{'Scaffold'})  {
			    &printFasta(*SCAFFOLDFILE, , $pep);
			    $role_id = "1008";
			}
			#if end5 or end3 is < beginning of phage region or > the end of the phage region, draw a truncated ORF (ORFNA) and make either end5 or end3 the begin or end coord of the region so it will print
			if ($end5 < $contig_end5) {
			    $end5 = $contig_end5;
			    $ftype = "ORFNA";
			}
			elsif ($end5 > $contig_end3) {
			    $end5 = $contig_end3;
			    $ftype = "ORFNA";
			}
			if ($end3 < $contig_end5) {
			    $end3 = $contig_end5;
			    $ftype = "ORFNA";
			}
			elsif ($end3 > $contig_end3) {
			    $end3 = $contig_end3;
			    $ftype = "ORFNA";
			}
			if ($ftype eq "ORFNA") {
			    &printatt(*ATTFILE, $end5, $end3, $asmbl_id, $featname, $ftype, $role_id, $region_tag, $hitref->{$featname}->{'com_name'});
			    $ftype = ""; #clear the variable for the next gene
			}
			else {
			    &printatt(*ATTFILE, $end5, $end3, $asmbl_id, $featname, "ORF", $role_id, $region_tag, $hitref->{$featname}->{'com_name'});
			}
			print APISFILE "$hitref->{$featname}->{'header'}\t$region_tag\n";
		    }
		}
		if (exists $hitref->{$featname}->{'annotation'}) {
		    if (($featname =~ /^att/) && ($size_att > 1))  {  # don't print E-VALUE if an att site.  Added check for att sites 04/16/13
			$printfeatname = $phageref->{$key}->{'memberhash'}{$yek}->{'featname'};
			$annotation = "$hitref->{$featname}->{'annotation'}";
			write if ($DEBUG);
			write $FH;
			&printatt(*ATTFILE, $end5, $end5, $asmbl_id, $printfeatname, "TSD", $role_id, $region_tag, $annotation); # print att sites on linear fig
			next;
		    }
		    elsif ($featname =~ /^att/) { # added 04/16/13 so it will not print old, incorrect att sites when they were removed from consideration
			next;
		    }
		    else  {
			$annotation = "$hitref->{$featname}->{'annotation'}, TAG = $hitref->{$featname}->{'phage'}, E-VALUE = $hitref->{$featname}->{'evalue'}";
		    }
		    $printfeatname = $featname;
		    write if ($DEBUG);
		    write $FH;
		    $end5 = "";
		    $annotation = "[ANNO] $hitref->{$featname}->{'com_name'}";
		    $printfeatname = "";
		    write if ($DEBUG);
		    write $FH;
		}
		else { # modified 01/13/22 to streamline all tRNA/tmRNA/ncRNAs into one condition [dfouts]
		    if (($featname =~ /tm?RNA/) || ($featname =~ /ncRNA/)) {
			if ($featname =~ /tRNA/) { # check if this is a tRNA
			    $RNA_type = "tRNA";
			}
			elsif ($featname =~ /tmRNA/) { # check if this is a tmRNA
			    $RNA_type = "tmRNA";
			}
			else { # last option is a ncRNA, no need to check
			    $RNA_type = "ncRNA";
			}
			my $com_feat = $asmbl_id . "_" . $featname;
			$annotation = "[$RNA_type] $hitref->{$com_feat}->{'com_name'} $hitref->{$com_feat}->{'organism'}";
			$end3 = $hitref->{$com_feat}->{'end3'}; # RNA feat_names stored differently in hithash to make them unique
			#if end5 or end3 is < beginning of phage region or > the end of the phage region, draw a truncated ORF (ORFNA) and make either end5 or end3 the begin or end coord of the region so it will print
			if ($end5 < $contig_end5) {
			    $end5 = $contig_end5;
			}
			elsif ($end5 > $contig_end3) {
			    $end5 = $contig_end3;
			}
			if ($end3 < $contig_end5) {
			    $end3 = $contig_end5;
			}
			elsif ($end3 > $contig_end3) {
			    $end3 = $contig_end3;
			}
			print "printatt($RNA_type): featname = <$featname>, com_feat = <$com_feat>, end5 = <$end5>, end3 = <$end3>\n" if ($DEBUG);
			&printatt(*ATTFILE, $end5, $end3, $asmbl_id, $featname, $RNA_type, $role_id, $region_tag, $hitref->{$featname}->{'com_name'});
		    }
		    else {
			$annotation = "[ANNO] $hitref->{$featname}->{'com_name'}";
		    }
		    $printfeatname = $featname;
		    write if ($DEBUG);
		    write $FH;
		}
		##
		#if (($contig_end5 <= $yek ) && ($print5 == 0) && (exists $phageref->{$key}->{'memberhash'}{$yek})) {
		if ($print5 == 0) { # print first ORF in the region
		    print TAB "$featname\t";
		    $print5++;
		    
		    ### code below written on 3/25/22 by dfouts to print the 5' flanking ORF ###
		    $flank5 = $genomearray[ $hitref->{ $phageref->{$key}->{'memberhash'}{$yek}->{'featname'} }->{'array_pos'}-1 ];
		    print "=>5prime ORF end5 = $key ($rehash{$asmbl_id}{$key}->{'featname'}), 5flank ORF end5 = $5flank ($rehash{$asmbl_id}{$flank5}->{'featname'})\n" if ($DEBUG);
		    ###
		}
		elsif ($yek <= $contig_end3) { # changed < to <= to include last gene in tabular output (dfouts 09/11/2014)
		    $hold_locus = $yek;
		}
		##
		if ($featname eq $phageref->{$key}->{'target'}) { # if the gene is the target, note this!
		    $end5 = "";
		    $printfeatname = "";
		    $annotation = "[TARGET]";
		    write if ($DEBUG);
		    write $FH;
		}
		if (exists $hitref->{$featname}->{'hmm'}) {
		    foreach $HMM_name (sort {$a <=> $b} keys %{$hitref->{$featname}->{'hmm'}}) {
			$hmm_hit_names = $hmm_hit_names . $HMM_name . "\t";
			$end5 = "";
			$printfeatname = "";
			$annotation = "[HMM-$hitref->{$featname}->{'hmm'}{$HMM_name}->{'model'}] ($HMM_name): $hitref->{$featname}->{'hmm'}{$HMM_name}->{'hmm_com_name'}, SCORE = $hitref->{$featname}->{'hmm'}{$HMM_name}->{'score'}, TRUSTED = $hitref->{$featname}->{'hmm'}{$HMM_name}->{'trusted'}, NOISE = $hitref->{$featname}->{'hmm'}{$HMM_name}->{'noise'}";
			write if ($DEBUG);
			write $FH;
		    }
		}
	    }
	}
	if (!($phageref->{$key}->{'att_distance'})) { #added 01_11_12 by dfouts
	    $phageref->{$key}->{'att_distance'} = "N.D.";
	}
	print "PRINTSUB: flank5 = <$flank5>, contig_end5 = <$contig_end5>\n" if ($DEBUG);
        print "PRINTSUB: flank3 = <$flank3>, contig_end3 = <$contig_end3>\n" if ($DEBUG);		 
	if ($flank5 < 1) { # check if we run off the beginning of the contig.  Maybe in the future check if circular and wrap around to the end gene?
	    $flank5_gene = "N.A.";
	}
	else {
	    $flank5_gene = "$rehash{$asmbl_id}{$flank5}->{'featname'}";
	}
	$flank3 = $genomearray[ $hitref->{ $phageref->{$key}->{'memberhash'}{$hold_locus}->{'featname'} }->{'array_pos'}+1 ];
	if (($flank3 > $asmbl_idhash{$asmbl_id}->{'genomesize'}) || (!defined($flank3))) { #check to see if we ran off the end of the contig.  Maybe in the future check if circular and wrap around to the beginning gene?
	    $flank3_gene = "N.A.";
	}
	else {
	    $flank3_gene = "$rehash{$asmbl_id}{$flank3}->{'featname'}";
	}
	print "=>3prime ORF end5 = $key ($flank5_gene), 3flank ORF end5 = $flank3 ($flank3_gene), contig_size = <$asmbl_idhash{$asmbl_id}->{'genomesize'}>\n" if ($DEBUG);
	print TAB "$phageref->{$key}->{'memberhash'}{$hold_locus}->{'featname'}\t$flank5_gene\t$flank3_gene\t$phageref->{$key}->{'int_HMM'}\t$phageref->{$key}->{'core_HMM'}\t$phageref->{$key}->{'above_noise_core_HMM'}\t$phageref->{$key}->{'lytic_HMM'}\t$phageref->{$key}->{'tail_HMM'}\t$phageref->{$key}->{'RBP_HMM'}\t$phageref->{$key}->{'Mu_HMM'}\t$phageref->{$key}->{'mucin_HMM'}\t$phage_direction\t$phageref->{$key}->{'att_distance'}\t$phageref->{$key}->{'ORFcnt'}\t$phageref->{$key}->{'serine_recomb'}\t$phageref->{$key}->{'toxin_HMM'}\t$phageref->{$key}->{'RBP'}\n";
	##
	$frag_line_num = $frag_line_num + 2; # increment the frag_line_num global line counter variable used to establish unique line numbers for each phage in a genome.
    };
    ## end of local print subroutine
    
    if (($max > 1) || ($max == 0))  {
	print "There are $max putative large/medium/small phages in $asmbl_id!\n";
	print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    }
    elsif ($max == 1)  {
	print "There is $max putative large/medium/small phage in $asmbl_id!\n";
	print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    }
    
format OUTPUT_TOP =
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$header
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$region_header
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$region2_header

END5     FEATNAME           ANNOTATION OF BEST HIT FROM PHAGE DB, [ANNO]TATION, OR [HMM] HIT                                                                       PAGE: @<<
                                                                                                                                                                         $%
............................................................................................................................................................................
.
format OUTPUT =
@<<<<<<< @<<<<<<<<<<<<<<<<< @<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<...
$end5,   $printfeatname,         $annotation
.
format OUTmedium_TOP =
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$header
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$region_header
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$region2_header

END5     FEATNAME           ANNOTATION OF BEST HIT FROM PHAGE DB, [ANNO]TATION, OR [HMM] HIT                                                                       PAGE: @<<
                                                                                                                                                                         $%
............................................................................................................................................................................
.

format OUTmedium =
@<<<<<<< @<<<<<<<<<<<<<<<<< @<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<...
$end5,   $printfeatname,         $annotation
.

format OUTsmall_TOP =
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$header
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$region_header
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$region2_header

END5     FEATNAME           ANNOTATION OF BEST HIT FROM PHAGE DB, [ANNO]TATION, OR [HMM] HIT                                                                       PAGE: @<<
                                                                                                                                                                         $%
............................................................................................................................................................................
.
format OUTsmall =
@<<<<<<< @<<<<<<<<<<<<<<<<< @<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<...
$end5,   $printfeatname,         $annotation
.

format STDOUT_TOP =
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$header
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$region_header
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$region2_header

END5     FEATNAME           ANNOTATION OF BEST HIT FROM PHAGE DB, [ANNO]TATION, OR [HMM] HIT                                                                       PAGE: @<<
                                                                                                                                                                         $%
............................................................................................................................................................................
.
format STDOUT =
@<<<<<<< @<<<<<<<<<<<<<<<<< @<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<...
$end5,   $printfeatname,         $annotation
.
	foreach $key (sort {$a <=> $b} keys %{$phageref}) { #may want to make this a for loop because some keys are being renumbered on the fly, causing some to not get examined and printed when I actually delete and rekey
	    $label = "";
	    $print_files = undef;
	    $size = abs($phageref->{$key}->{'5prime'} - $phageref->{$key}->{'3prime'})+1;
	    $size_att = abs($phageref->{$key}->{'5prime_att'} - $phageref->{$key}->{'3prime_att'})+1;
	    print "print_regions: phage $key, size_att = $size_att\n" if ($DEBUG);
	    if ($size_att > $size)  {  # if there is an att site and the size of the phage with att sites is bigger, use that to group small, med, large
		$sizer = $size_att;
	    }
	    else  {
		$sizer = $size;
	    }
	    ##
	    if (exists $phageref->{$key}->{'direction'})  {
		$phage_direction = $phageref->{$key}->{'direction'};
		print "Print_regions Main loop: phage direction = <$phage_direction>\n" if ($DEBUG);
	    }
	    #this is necissary since there are some scenarios whereby the direction is not set yet 06/05/2019 dfouts
	    else  {
		$phage_direction = "+"; # set default direction as ------->
		print "Print_regions Main loop: phage direction set to <$phage_direction>\n" if ($DEBUG);
	    }
	    if ($size_att > 1)  { # 5prime_att and 3prime_att are already oriented properly!
	#	if ($phage_direction eq "+") { # phage is ----->
		    $phage_end5 = $phageref->{$key}->{'5prime_att'};
		    $phage_end3 = $phageref->{$key}->{'3prime_att'};
	#	    print "Print_regions: making (+) phage $key\n" if ($DEBUG);
	#	}
	#	else { # phage is <-----, so we want to switch phage_end5 to be the 5' most coord relative to the assembly
	#	    $phage_end5 = $phageref->{$key}->{'3prime_att'};
	#	    $phage_end3 = $phageref->{$key}->{'5prime_att'};
	#	    print "Print_regions: making (-) phage $key.  phage_end5 = $phage_end5, phage_end3 = $phage_end3, 5prime_att = ($phageref->{$key}->{'5prime_att'}), 3prime_att = ($phageref->{$key}->{'3prime_att'})\n" if ($DEBUG);
	#	}
		#$phageref->{$key}->{'seq'} = &quickcut($ahref->{$asmbl_id}->{'sequence'}, $phageref->{$key}->{'5prime_att'}, $phageref->{$key}->{'3prime_att'});
	    }
	    else  { # no att site found, use the coords in 5prime and 3prime
		if ($phage_direction eq "+") { # phage is -----> (brought this code back so phage will be in correct ori if int on one side, but contig end on other side)
		    $phage_end5 = $phageref->{$key}->{'5prime'};
		    $phage_end3 = $phageref->{$key}->{'3prime'};
		    print "Print_regions+:  >NO ATT< phage_end5 = ($phage_end5), phage_end3 = ($phage_end3), direction = ($phage_direction)\n" if ($DEBUG);
		}
		else { # phage is <-----, so we want to switch phage_end5 to be the 5' most coord relative to the assembly
		    $phage_end5 = $phageref->{$key}->{'3prime'};
		    $phage_end3 = $phageref->{$key}->{'5prime'};
		    print "Print_regions-:  >NO ATT< phage_end5 = ($phage_end5), phage_end3 = ($phage_end3), direction = ($phage_direction)\n" if ($DEBUG);
		}
		#$phageref->{$key}->{'seq'} = &quickcut($ahref->{$asmbl_id}->{'sequence'}, $phageref->{$key}->{'5prime'}, $phageref->{$key}->{'3prime'});
	    }
	    $phageref->{$key}->{'seq'} = &quickcut($ahref->{$asmbl_id}->{'sequence'}, $phage_end5, $phage_end3); # consulidated cutting based on phage orientation (07/06/22)
	    
	    if ($asmbly_status) {
		$phageref->{$key}->{'gc'} = &get_gc_content($phageref->{$key}->{'seq'});
	    }
	    else {
		$phageref->{$key}->{'gc'} = "NA";
	    }
	    ##
	    
	    ###
	    # check for type and whether has att or not
	    #
	    if ($phageref->{$key}->{'piggyback'} == "1") { #set variable for pretty printing of is_piggyback
		$is_piggy = "Y";}
	    else {
		$is_piggy = "N";
	    }
	    print "PIGGY_BACK:  <$is_piggy>\n" if ($DEBUG);
	    if (($size_att > 1) && (!defined ($phageref->{$key}->{'type'})))  {  # probable integrated element (plasmid, transposon, other)
		if ($phageref->{$key}->{'above_noise_core_HMM'} > 1) { #need 2 >= noise core HMM hits
		    ## maybe make a more strict core set of portal, major capsid, terminase, capsid protease >= noise definition ##
		    $phageref->{$key}->{'type'} = "prophage";
		}
		else  {
		    $phageref->{$key}->{'type'} = "integrated element";
		}
	    }
	    elsif (($size_att > 1) && ($phageref->{$key}->{'type'} eq "prophage") && (!defined ($phageref->{$key}->{'class'})) && ($phageref->{$key}->{'core_HMM'} == 0) && ($phageref->{$key}->{'core_BLAST'} == 0))  {
		$phageref->{$key}->{'class'} = "satellite"; # if no core HMMs or BLAST hits are found, is a prophage with no class and an att site, call satellite (I know there are always exceptions)
	    }
	    elsif (($size_att == 1) && (!defined ($phageref->{$key}->{'type'}))) {
		$phageref->{$key}->{'type'} = "degenerate";
	    }
	    ###
	    if ($key < 10) {
		$region_tag = "PFPR0$key";
	    }
	    else {
		$region_tag = "PFPR$key";
	    }
	    #$seq_tag = "PFPROPHAGE" . "_" . $bugtag . "_" . "$region_tag";
	    ###
	    $region_id = $asmbl_id . "_" . $region_tag; # added 03/08/18 by dfouts
	    
	    if (($sizer < 10000) && ($strict == 0))  {
		if ($small_open == 0)  { # only open each file if needed and only once per invocation of the script so we can get multiple phages per file
		    open (OUTsmall, ">$write_dir/$asmbl_id\_$hitsperwindow\_phpico\_$window\_$step\_$evalue.out") || die "can't open file $write_dir/$asmbl_id\_$hitsperwindow\_phpico\_$window\_$step\_$evalue.out\n";
		    $small_open = 1;
		}
		if ($sizer >= 5000) { $micro++; $totalmicrobp = $totalmicrobp + $sizer;}
		elsif ($sizer >= 2000) { $nano++; $totalnanobp = $totalnanobp + $sizer;}
		else { $pico++; $totalpicobp = $totalpicobp + $sizer;}
		$small++;
		$print_files = 1;
		$totalsmallbp = $totalsmallbp + $sizer;
		if (!defined ($phageref->{$key}->{'class'}))  {
		    $label = "Small";
		}
		else {
		    $label = "Small $phageref->{$key}->{'class'}";
		}
		&print_sub($key, \*OUTsmall, $label, $is_piggy, $small, $phageref->{$key}->{'type'}, $phage_end5, $phage_end3, $phage_direction);
	    }
	    elsif (($sizer >= 10000) && ($sizer <= 18000))  {
		if (($strict == 0) || (($strict == 1) && (($phageref->{$key}->{'type'} eq "prophage") || ($phageref->{$key}->{'type'} eq "bacteriocin"))))  {
		    if ($med_open == 0)  {
			open (OUTmedium, ">$write_dir/$asmbl_id\_$hitsperwindow\_phmedio\_$window\_$step\_$evalue.out") || die "can't open file $write_dir/$asmbl_id\_$hitsperwindow\_phmedio\_$window\_$step\_$evalue.out\n";
			$med_open = 1;
		    }
		    $medium++;
		    $print_files = 1;
		    $totalmediumbp = $totalmediumbp + $sizer;
		    if (!defined ($phageref->{$key}->{'class'}))  {
			$label = "Medium";
		    }
		    else {
			$label = "Medium $phageref->{$key}->{'class'}";
		    }
		    &print_sub($key, \*OUTmedium, $label, $is_piggy, $medium, $phageref->{$key}->{'type'}, $phage_end5, $phage_end3, $phage_direction);
		}
		else {
		    print "L1: Region $key (type $phageref->{$key}->{'type'}) is unlikely a prophage, deleting ...\n" if ($DEBUG);
		    #delete($phageref->{$key}); # added 01/24/2019 by Derrick Fouts to actually delete the region keys and renumber the existing prophage regions
		    #$m = $key - 1;
		    #&rekey_regions($m, $phageref, $DEBUG); # call subroutine to renumber phage regions
		    $phage_regions--;  # deprecate phage_region counter
		}
	    }
	    elsif ($sizer > 18000) { # Print the data for phage regions > 18 Kbp in size to a file, bug fix was reporting small regions as large
		if (($strict == 0) || (($strict == 1) && (($phageref->{$key}->{'type'} eq "prophage") || ($phageref->{$key}->{'type'} eq "bacteriocin"))))  {
		    if ($large_open == 0)  {
			open (OUTPUT, ">$write_dir/$asmbl_id\_$hitsperwindow\_phregions\_$window\_$step\_$evalue.out") || die "can't open file $write_dir/$asmbl_id\_$hitsperwindow\_phregions\_$window\_$step\_$evalue.out\n";
			$large_open = 1;
		    }
		    $print_files = 1;
		    $large++;
		    $totalphagebp = $totalphagebp + $sizer;
		    if (!defined ($phageref->{$key}->{'class'}))  {
			$label = "Large";
		    }
		    else {
			$label = "Large $phageref->{$key}->{'class'}";
		    }
		    &print_sub($key, \*OUTPUT, $label, $is_piggy, $large, $phageref->{$key}->{'type'}, $phage_end5, $phage_end3, $phage_direction);
		}
		else {
		    print "L2: Region $key (type $phageref->{$key}->{'type'}) is unlikely a prophage, deleting ...\n" if ($DEBUG);
		    #delete($phageref->{$key}); # added 01/24/2019 by Derrick Fouts to actually delete the region keys and renumber the existing prophage regions
		    #$m = $key - 1;
		    #&rekey_regions($m, $phageref, $DEBUG); # call subroutine to renumber phage regions
		    $phage_regions--; # deprecate phage_region counter
		}
	    }
	    else {
		print "L3: Region $key is unlikely a prophage, deleting ...\n" if ($DEBUG);
		#delete($phageref->{$key}); # added 01/24/2019 by Derrick Fouts to actually delete the region keys and renumber the existing prophage regions
		#$m = $key - 1;
		#&rekey_regions($m, $phageref, $DEBUG); # call subroutine to renumber phage regions
		$phage_regions--;  # deprecate phage_region counter
	    }
	    if (($print_files) && (($phageref->{$key}->{'type'} eq "prophage") || ($phageref->{$key}->{'type'} eq "bacteriocin")))  { # make con file only if a prophage or bacteriocin
		$con_header = $ahref->{$asmbl_id}->{'title'};
		$con_header =~ s/\,.*$//; # remove everything past comma
		if (length($label) > 1) { $label = "$label "; }
		$con_header = $region_id . " " . $con_header . " " . $label . $phageref->{$key}->{'type'} . " region " . $key . " (" . $phage_end5 . "-" . $phage_end3 . " bp)";
		&printFasta(*CONFILE, $con_header, $phageref->{$key}->{'seq'});
	    }
    }
    
    $percent = ($totalphagebp/$asmbl_idhash{$asmbl_id}->{'genomesize'})*100;
    $percentsmall = ($totalsmallbp/$asmbl_idhash{$asmbl_id}->{'genomesize'})*100;
    $percentmedium = ($totalmediumbp/$asmbl_idhash{$asmbl_id}->{'genomesize'})*100;
    
    print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    
    if ($phage_regions > "0")  {
	
	printf "There are $large regions > 18 Kb and summing $totalphagebp bp of sequence (%5.2f", $percent;
	print "% of the genome)\n";
	
	printf "There are $medium regions between 10 and 18 Kb and summing $totalmediumbp bp of sequence (%5.2f", $percentmedium;
	print "% of the genome)\n";
	
	printf "There are $small regions < 10 Kb and summing $totalsmallbp bp of sequence (%5.2f", $percentsmall;
	print "% of the genome)\n";
	
	if ($large > 0) {
	    print OUTPUT "............................................................................................................................................................................\n";
	    printf OUTPUT "There are $large regions > 18 Kb and summing $totalphagebp bp of sequence (%5.2f", $percent;
	    print OUTPUT "% of the genome)\n";
	    close (OUTPUT);
	}
	
	if ($medium > 0) {
	    print OUTmedium "............................................................................................................................................................................\n";
	    printf OUTmedium "There are $medium regions between 10 and 18 Kb and summing $totalmediumbp bp of sequence (%5.2f", $percentmedium;
	    print OUTmedium "% of the genome)\n";
	    close (OUTmedium);
	}
	
	if ($small > 0) {
	    print OUTsmall "............................................................................................................................................................................\n";
	    printf OUTsmall "There are $small regions < 10 Kb and summing $totalsmallbp bp of sequence (%5.2f", $percentsmall;
	    print OUTsmall "% of the genome)\n";
	    close (OUTsmall);
	}
    }
}

sub write_output {  # write xgraph file and report file
    
    my ($write_dir,$asmbl_id,$hitsperwindow,$ahref,$figref) = @_;
    my $minimum_hits = $hitsperwindow-1;
    my $n = 0;
    my $max_peak = "";
    my $counts = 0;
    my $title = $ahref->{$asmbl_id}->{'title'};
    $title =~ s/$asmbl_id//;
    my $fig = "";
    my @feat_names = ();
    
    open (OUTPUT, ">$write_dir/$asmbl_id\_$hitsperwindow\_phgraph\_$window\_$step\_$evalue.out") || die "can't open file phgraph$window.out.\n";
    print OUTPUT "title = $title id $asmbl_id (>= $minimum_hits hits per window)\n";
    print OUTPUT "title_x = Position (bp)\n";
    print OUTPUT "title_y = # Hits per $window bp window/$step bp step\n";
    print OUTPUT "color = red\n";
    print OUTPUT "thickness = 1\n";
    open (REPORT, ">$write_dir/$asmbl_id\_$hitsperwindow\_phreport\_$window\_$step\_$evalue.out") || die "can't open file phreport$window.out.\n";
    foreach $fig (sort {$a <=> $b} keys %{$figref}) {
	if (defined($figref->{$fig}->{'counts'})) { # added 04/16/13 by dfouts to remove negative $fig numbers and undef counts (from if statements adding junk to the hash???)
	    $counts = $figref->{$fig}->{'counts'};
	}
	else {
	    next; # added 04/16/13 by dfouts to remove negative $fig numbers and undef counts (from if statements adding junk to the hash???)
	}
	if (defined($figref->{$fig}->{'peak_number'})) { # changed from exists to defined 04/16/13 by dfouts
	    $max_peak = $figref->{$fig}->{'peak_value'} + .1;
	    print OUTPUT "anno\t$fig\t$max_peak\t$figref->{$fig}->{'peak_number'}\n";
	}
	print OUTPUT "$fig\t$counts\n";
	# strange syntax because a runtime exception occurs under use strict because you're dereferencing an undefined reference where autovivification won't occur.
	@feat_names = exists( $figref->{$fig}->{'featnames'} ) 
	    ? @{ $figref->{$fig}->{'featnames'} }
	: ();
	print REPORT "$fig\t$counts\t@feat_names\n";
    }
    close (OUTPUT);
    close (REPORT);
}

sub option_help {
   print <<_EOB_;
$prog - Find prophage regions in DNA sequence from btab (blastp) and hmm data
  Usage: $prog <options>
 Switch: -h for help\
 Option:
     -b: (b)ase directory path [default = PWD]
     -d: (d)ebug mode (default = 0)
     -i: tab-delimitied flat (i)nfofile containing scaffold/contig/assembly_ID size_of_molecule feat_name end5 end3 com_name [REQUIRED]
     -m: htab file containing h(m)m data (REQUIRED for finding integrases and att sites)
     -n: Arago(n) tmRNA-finding output file (-m option in aragon)
     -p: (p)ath to btab file (default = base directory)
     -r: t(r)nascan-SE output file
     -c: (c)mscan output file
     -s: (s)tep size (default = 5000 nucleotides)
     -t: name of wu-blas(t) or NCBI blas(t) (-m 8 option) btab input file [REQUIRED]
     -w: Scanning (w)indow size (default = 10000 nucleotides)
     -A: File name of (A)ssembly/contig/scaffold fasta file if in basedir or supply full path to file
     -E: (E)-value (default = 0.00001)
     -F: att site search (F)unction (B or b for NCBI BLAST, M or m for MUMmer, F or f for FASTA33) (default = BLAST)
     -H: Number of allowable hits per window to mark a region (default = 4)
     -S: (S)trict mode:  print only regions that have core HMM hits or Mu-like and are > 10 Kbp (default = 0)
     -V: print (V)ersion information
 Format: 
 Output:
_EOB_
    exit;
}
################################## M A I N ####################################

### gather all information ###

print "MAIN:  Checking the format of $infofile ...\n" if ($DEBUG);
&write_log("1","Checking the format of $infofile");
$infofile = &check_infofile($basedir,$infofile,$asmbly_file,$asmbl_id,$DEBUG);

print "MAIN:  Retrieving information from flat file $infofile . . .\n" if ($DEBUG);
&write_log("1","Retrieving information from flat file $infofile");
$num_contigs = &get_gene_info($infofile,\%asmbl_idhash,\%hithash,\%rehash,$DEBUG); # removed defined_asmbl_id 04/17/13 dfouts

print "MAIN:  There are $num_contigs contigs/assemblies in $infofile\n" if ($DEBUG);
&write_log("1","There are $num_contigs contigs/assemblies in $infofile");

if (defined($tRNA_data))  { # is there a tRNAscan-SE file?
  print "MAIN:  Looking for tRNA information ...\n" if ($DEBUG);
  &write_log("1","Looking for tRNA information");
  &get_tRNAs($tRNA_data,\%tRNA_hash,\%hithash,\%rehash,$DEBUG);
}
else {
  print "MAIN:  User did not specifiy a tRNAScan file or it contains no data, skipping any searches involving tRNA ...\n" if ($DEBUG);
  &write_log("1","User did not specifiy a tRNAScan file or it contains no data, skipping any searches involving tRNA");
}

if (defined($tmRNA_data))  { # is there an Aragorn file?
  print "MAIN:  Looking for tmRNA information ...\n" if ($DEBUG);
  &write_log("1","Looking for tmRNA information");
  &get_tmRNAs($tmRNA_data,\%tRNA_hash,\%hithash,\%rehash,$DEBUG);
}
else {
  print "MAIN:  User did not specifiy an Aragorn file or it contains no data, skipping any searches involving tmRNA ...\n" if ($DEBUG);
  &write_log("1","User did not specifiy an Aragorn file or it contains no data, skipping any searches involving tmRNA");
}

if (defined($ncRNA_data))  { # is there a cmscan output file?
  print "MAIN:  Looking for ncRNA information ...\n" if ($DEBUG);
  &write_log("1","Looking for ncRNA information");
  &get_ncRNAs($ncRNA_data,\%tRNA_hash,\%hithash,\%rehash,$DEBUG);
}
else {
  print "MAIN:  User did not specifiy an cmscan file or it contains no data, skipping any searches involving ncRNAs ...\n" if ($DEBUG);
  &write_log("1","User did not specifiy an cmscan file or it contains no data, skipping any searches involving ncRNAs");
}

print "MAIN:  Getting list of OK phage common names ($comfile) . . .\n" if ($DEBUG);
&write_log("1","Getting list of OK phage common names ($comfile)");
&get_ok_comnames($comfile,\%ok_comnames,$DEBUG);

print "MAIN:  Getting tags, names and taxon info from BLAST DB file ($DBfile) . . .\n" if ($DEBUG);
&write_log("1","Getting tags, names and taxon info from BLAST DB file ($DBfile)");
&get_DB_info($DBfile,\%DB_info,$DEBUG);

if (defined($hmmfile))  { # is there an hmm file?
    open (HMMFILE, "<$hmmfile") || &write_log("4","can't open file HMM datafile $hmmfile: $!\n");
       while (<HMMFILE>)  {
           if (/^# HMMER 3/) {
              $HMMversion = "3";
              last;
           }
           else {
              $HMMversion = "2";
              last;
           }
       }
    close (HMMFILE);
    push (@hmm_data, $hmmfile);
}
elsif (-s "$hmmdir/combined.hmm3") {
    $HMMversion = "3";
    push (@hmm_data, "$hmmdir/combined.hmm3");
}
else {
    $HMMversion = "2";
    push (@hmm_data, "$hmmdir/combined.hmm_GLOCAL");
    push (@hmm_data, "$hmmdir/combined.hmm_FRAG");
}

print "MAIN:  HMMer $HMMversion output detected ...\n" if ($DEBUG);

#### loop through each filename in the @filenames array in the load_hash_from_file subroutine

## future update:  place all of these input files into a config.ini file set-up.

if ($HMMversion == "2") {
   @filenames = ("phage_exclude_v6.lst",
		 "core_hmm.lst",
		 "lysin_holin_hmm.lst",
		 "tails_hmm.lst",
		 "large_term_hmm.lst",
		 "small_term_hmm.lst",
		 "portal_hmm.lst",
		 "protease_hmm.lst",
		 "capsid_hmm.lst",
		 "Small_term_v2.lst",
		 "Large_term_v2.lst",
		 "Portal_v2.lst",
		 "Protease_v2.lst",
		 "Scaffold_v2.lst",
		 "Major_capsid_v2.lst",
		 "RBP.lst");
}
elsif ($HMMversion == "3") {
   @filenames = ("phage_exclude_v6.lst",
		 "core_hmm3.lst",
		 "lysin_holin_hmm3.lst",
		 "tails_hmm3.lst",
		 "large_term_hmm3.lst",
		 "small_term_hmm3.lst",
		 "portal_hmm3.lst",
		 "protease_hmm3.lst",
		 "capsid_hmm3.lst",
		 "toxins_hmm3.lst",
		 "rbp_hmm3.lst",
		 "mucin_binding_hmm3.lst",
		 "Small_term_v2.lst",
		 "Large_term_v2.lst",
		 "Portal_v2.lst",
		 "Protease_v2.lst",
		 "Scaffold_v2.lst",
		 "Major_capsid_v2.lst",
		 "mucin_binding_hmm3.lst",
		 "RBP.lst");
}

&load_hash_from_file($phome,\@filenames,\%exclude_hash,\%lists_hash,\%HMMs_hash,$DEBUG);

### moved below load_hash_from_file on 09_13_13 by dfouts because proteins were not being excluded ###
print "MAIN:  Obtaining hits to phage database from $btabfile . . .\n" if ($DEBUG);
&write_log("1","Obtaining hits to phage database from $btabfile");
&select_featnames_from_btab($btabfile,$evalue,\%hithash,\%searchhash,\%exclude_hash,\%lists_hash,$DEBUG);
### ---------------------------- ###

if (! &find_hmms($phome,$HMMversion,\@hmm_data,$infofile,\%hithash,\%searchhash,\%HMMs_hash,\%serine_HMM_hash,\%tyrosine_HMM_hash,$DEBUG)) {
    print "MAIN:  User did not specify an hmm file used in finding integrases, skipping this analysis . . .\n" if ($DEBUG);
    &write_log("1","User did not specify an hmm file used in finding integrases, skipping this analysis");
}

# all these datafiles were pulled from the contig directories into the main genome directory to facilitate access to this data when there were many contigs (ie in draft genomes)
# their systematic names will make it easier for a bash script to concatenate all of them from a listfile

my $apisfile = "$mode_basedir/PFPR.apis";
my $tabfile = "$mode_basedir/PFPR_tab.txt"; # added 08/22/2008 by Derrick Fouts
my $smalltermpep = "$mode_basedir/PFPR_small_terminase.pep";
my $largetermpep = "$mode_basedir/PFPR_large_terminase.pep";
my $portalpep = "$mode_basedir/PFPR_portal.pep";
my $proteasepep = "$mode_basedir/PFPR_protease.pep";
my $scaffpep = "$mode_basedir/PFPR_scaffold.pep";
my $capsidpep = "$mode_basedir/PFPR_capsid.pep";
my $serinerecombpep = "$mode_basedir/PFPR_serine_recombinase.pep"; # only the one we think is used for integration is pulled out
my $tyrosinerecombpep = "$mode_basedir/PFPR_tyrosine_recombinase.pep"; # only the one we think is used for integration is pulled out
my $toxinspep = "$mode_basedir/PFPR_toxins.pep"; #added 01/09/2012
my $rbpspep = "$mode_basedir/PFPR_rbps.pep"; #added 03/07/2018
my $confile = "$mode_basedir/PFPR.con"; # added 03/07/14
my $seqfile = "$mode_basedir/PFPR.seq"; # added 03/07/14
my $pepfile = "$mode_basedir/PFPR.pep"; # added 03/07/14
my $attfile = "$mode_basedir/PFPR.att"; # added 07/06/17 by dfouts - for making the Linear_Display gene.att file for all prophage regions
my $fragfile = "$mode_basedir/PFPR.frag"; # added 03/08/18 by dfouts - for making the Linear_Display frag.file for all prophage regions
$frag_line_num = 1; # set frag file line counter to one to start with

open (APISFILE, ">$apisfile") || die "can not open APIS lookup file $apisfile for writing\n";
open (TAB, ">$tabfile") || die "can not open file $tabfile for writing\n";
# modified tabular output by changing asmbl_id to region_id so that a unique identifier will be associated with each prophage region, consisting of asmbl_id-region_tag on 03/08/18 by dfouts
# added RBP counts to tabular output file 03/08/18 dfouts
print TAB "#region_id\tgenome_size\tgenome_gc\tbegin_region\tend_region\tsize_region\tlabel\tpiggy_back\ttype\tattR\tattL\ttarget\tregion_gc\tbest_db_match\tbegin_gene\tend_gene\t5flank\t3flank\t#integrase_HMMs\t#core_HMMs\t#>noise_HMMs\t#lytic_HMMs\t#tail_HMMs\t#RBP_HMMs\t#Mu_HMMs\t#mucin_HMMs\tregion_orientation\tdistance_int_to_att\t#genes\t#serine_recombinases\t#toxins\t#RBPs\n"; # print header line for PFPR.tab file
open (SMALLTERMFILE, ">$smalltermpep") || die "can not open small terminase pepfile $smalltermpep for writing\n";
open (LARGETERMFILE, ">$largetermpep") || die "can not open large terminase pepfile $largetermpep for writing\n";
open (PORTALFILE, ">$portalpep") || die "can not open portal pepfile $portalpep for writing\n";
open (PROTEASEFILE, ">$proteasepep") || die "can not open protease pepfile $proteasepep for writing\n";
open (SCAFFOLDFILE, ">$scaffpep") || die "can not open scaffold pepfile $scaffpep for writing\n";
open (CAPSIDFILE, ">$capsidpep") || die "can not open major capsid pepfile $capsidpep for writing\n";
open (SERINERECOMBFILE, ">$serinerecombpep") || die "can not open serine recombinase pepfile $serinerecombpep for writing\n";
open (TYROSINERECOMBFILE, ">$tyrosinerecombpep") || die "can not open tyrosine recombinase pepfile $tyrosinerecombpep for writing\n";
open (TOXINSFILE, ">$toxinspep") || die "can not open toxins pepfile $toxinspep for writing\n";
open (RBPSFILE, ">$rbpspep") || die "can not open toxins pepfile $rbpspep for writing\n";
open (CONFILE, ">$confile") || die "can't open file $confile\n";
open (SEQFILE, ">$seqfile") || die "can't open file $seqfile\n";
open (PEPFILE, ">$pepfile") || die "can't open file $pepfile\n";
open (BCNFILE, ">$badcomnamesfile") || &write_log("4","can't create file $badcomnamesfile: $!\n");
open (FRAGFILE, ">$fragfile") || die "can't open file $fragfile\n"; # LinearDisplay.pl frag.file
open (ATTFILE, ">$attfile") || die "can't open file $attfile\n"; # LinearDisplay.pl gene.att file

### now, loop through assemblies, looking for phages ###
foreach $asmbl_id (sort {$a cmp $b} keys %asmbl_idhash)  { # changed <=> to cmp because most of the time working with NCBI accessions and wanted the output PFPR_tab.txt file to be sorted correctly
  $num_contigs--;
  print "MAIN_LOOP:  Working on asmbl_id $asmbl_id, which is $asmbl_idhash{$asmbl_id}->{'genomesize'} bp in length ...\n" if ($DEBUG);
  $write_dir = &create_dir($mode_basedir,"$asmbl_id");
  ### main loop to count hits within window, incrementing by step size ### 
  %fighash = (); # clear out %fighash hash for each asmbl_id/contig
  %phage_hash = (); # clear our %phage_hash for each asmbl_id/contig
  @genomearray = (); # clear our @genomearray for each asmbl_id/contig
  for (my $n = 0; $n <= $asmbl_idhash{$asmbl_id}->{'genomesize'}; $n = $n + $step)  {  
    $fighash{$n}->{'counts'} = 0;
    foreach my $keys (sort {$a <=> $b} keys %{$searchhash{$asmbl_id}}) { # key is end5 of each BLAST or HMM hit
       if (($keys > $n) && ($keys < ($n+$window)))  { # if end5 is greater than current position, but less than the window size, count it 
         $fighash{$n}->{'counts'}++;
         $fighash{$n}->{'featnames'}[($fighash{$n}->{'counts'}-1)] = $searchhash{$asmbl_id}{$keys}->{'featname'};  # for each step, store the orf names with hits 
       }
       elsif ($keys >= ($n+$window)) { last; } # end the "foreach" loop when keys (end5 is greater than or equal to window size
       if ($keys < ($n+$step)) { delete $searchhash{$asmbl_id}{$keys}; }  # Remove keys that have already been counted and are less than step size so not resorted and counted again
    }
  }
  $hitsperwindow = $default_hpw;  # reset the user-defined or default hitsperwindow for contigs >= 10000 bp
  if (($asmbl_idhash{$asmbl_id}->{'genomesize'} < 50000) && ($asmbl_idhash{$asmbl_id}->{'genomesize'} >= 7500))  {$hitsperwindow = 3;}

  print "MAIN_LOOP:  Finding prophage region rough boundaries ...\n" if ($DEBUG);
  $phage_regions = &find_regions($asmbl_id,$hitsperwindow,\%fighash,\%hithash,\%phage_hash,$DEBUG);

  print "MAIN_LOOP:  Checking for seed end5 ...\n" if ($DEBUG);
  if (exists ($phage_hash{1}->{'seed_end5'})) { # added 11/14/02 to not to attempt extension of phage regions when there are none! 
    &populate_genomearray(\%rehash,\@genomearray,\%hithash,$asmbl_id);
    $phage_regions = &find_5prime_end($phage_regions,$asmbl_id,\@genomearray,\%phage_hash,\%hithash,\%tRNA_hash,\%rehash,\%ok_comnames,\%fighash,$DEBUG);  # find 5' boundary of phage regions
    my $phk = keys %phage_hash;
    print "MAIN_LOOP: post find_5prime: region counter = <$phage_regions>, phage_hash keys = <$phk>\n" if ($DEBUG);
    $phage_regions = &find_3prime_end($phage_regions,$asmbl_id,\@genomearray,\%phage_hash,\%hithash,\%tRNA_hash,\%rehash,\%ok_comnames,\%fighash,$DEBUG);  # find 3' boundary of phage regions
    $phk = keys %phage_hash;
    print "MAIN_LOOP: post find_3prime: region counter = <$phage_regions>, phage_hash keys = <$phk>\n" if ($DEBUG);
    &write_log("1","Writing output files");

    print "  MAIN_LOOP:  Checking for the presence of the contig file <$asmbly_file>...\n" if ($DEBUG);
    if (!defined($asmbly_file)) {
      print "WARNING: User did not specify a .1con file used in finding att sites, skipping this analysis . . .\n";
      &write_log("1","WARNING: User did not specify a .1con file used in finding att sites, skipping this analysis");
      $asmbly_status = 0;
    }
    else { $asmbly_status = 1 };

    $asmbly_status = &get_assemblies($asmbly_file,$infofile,\%assembly_hash,\%asmbl_idhash,$asmbl_id,$DEBUG) if ($asmbly_status == 1); # if we have a defined asmbl_id, then get the sequence...

    if ($asmbly_status == 1)  {
      print "..................................................................................................\n";
      print "Working on $assembly_hash{$asmbl_id}->{'title'} assembly|contig|scaffold id $asmbl_id which is $assembly_hash{$asmbl_id}->{'length'} bp in size, gc%: $assembly_hash{$asmbl_id}->{'gc'}%\n";
      &write_log("1","Working on $assembly_hash{$asmbl_id}->{'title'} assembly|contig|scaffold id $asmbl_id which is $assembly_hash{$asmbl_id}->{'length'} bp in size, gc%: $assembly_hash{$asmbl_id}->{'gc'}%");
    }
    print "  MAIN_LOOP:  Determining region type ...\n" if ($DEBUG);
    $phage_regions = &determine_region_type($phage_regions,\%phage_hash,\%hithash,\%HMMs_hash,\%serine_HMM_hash,\%tyrosine_HMM_hash,$DEBUG); #added return of phage_regions on 01/18/2023 to prevent empty PFPR_tab.txt files when deleting regions in this subroutine
    if (((defined($hmmfile)) || (defined($tRNA_data))) && ($asmbly_status == 1))  { # must have either hmm data or tRNA data AND an assembly to work with to look for att sites
      print "===> ASMBLY_STATUS: $asmbly_status\n" if ($DEBUG);
      print "Looking for putative phage attachment (att) sites using $search_method . . .\n" if ($DEBUG);
      &write_log("1","Looking for putative phage attachment (att) sites using $search_method");
      $prefix = "$write_dir/$asmbl_id\_$hitsperwindow\_phtest\_$window\_$step\_$evalue";
      print "  MAIN_LOOP:  Finding att sites ...\n" if ($DEBUG);
      &find_att_sites($prefix,$asmbl_id,$search_method,\%assembly_hash,\%phage_hash,\%hithash,\%rehash,\%tRNA_hash,\@genomearray,\%fighash,$window,$step,$hitsperwindow,\%ok_comnames,$DEBUG);
    }
    print "  MAIN_LOOP:  Writing output files ...\n" if ($DEBUG);
    &print_regions($write_dir,$asmbly_status,$asmbl_id,$hitsperwindow,$window,$step,$evalue,$phage_regions,\%phage_hash,\@genomearray,\%assembly_hash,\%hithash,\%rehash,\%DB_info,$DEBUG);    # print output of analysis
    &write_output($write_dir,$asmbl_id,$hitsperwindow,\%assembly_hash,\%fighash);     # print file for xgraph input
    &write_log("2", $num_contigs, $asmbl_id);
  }
  ##############
  print "phage_regions = $phage_regions\n" if ($DEBUG);
  if ($phage_regions == "0")  {  # added 11/14/02 to give default output of NO phages if there are no predicted prophages (ie ntcj01 has none)
    #print "Sorry, no phage regions were found in contig <$asmbl_id> :( . . .\n";
    #unlink glob "$write_dir/*";  # remove the directory if no prophages identified (to keep it clean for metagenomic studies)
    #rmdir $write_dir;
    &write_output($write_dir,$asmbl_id,$hitsperwindow,\%assembly_hash,\%fighash);     # print file for xgraph input anyway so we can see where any matches may be.
    &write_log("3", $num_contigs, $asmbl_id);
  }
  else {
    $genome_phage_cnt = $genome_phage_cnt + $phage_regions; 
  }
  print "assembly: $asmbl_id, phage_regions = $phage_regions, genome_phage_cnt = $genome_phage_cnt\n" if ($DEBUG); 
delete $assembly_hash{$asmbl_id}; # remove old assembly information (sequence, title, ...)
##############
}
close (APISFILE);
close (TAB);
print "END: phage_regions = $phage_regions, genome_phage_cnt = $genome_phage_cnt\n" if ($DEBUG);
if ($genome_phage_cnt == "0") {
  open (TAB, ">$tabfile") || die "can not open file $tabfile for writing\n";
  print TAB "No phages were detected in this genome\n"; # print message indicating no phages found
  close (TAB);
}
close (SMALLTERMFILE);
close (LARGETERMFILE);
close (PORTALFILE);
close (PROTEASEFILE);
close (SCAFFOLDFILE);
close (CAPSIDFILE);
close (SERINERECOMBFILE);
close (TYROSINERECOMBFILE);
close (TOXINSFILE);
close (RBPSFILE);
close (CONFILE);
close (SEQFILE);
close (PEPFILE);
close (BCNFILE);
close (FRAGFILE);
close (ATTFILE);

if (-z "$smalltermpep" == 1)  { # if zero contents, delete file
    print "removing $smalltermpep ...\n" if ($DEBUG);
    unlink $smalltermpep;
}

if (-z "$largetermpep" == 1)  { # if zero contents, delete file
    print "removing $largetermpep ...\n" if ($DEBUG);
    unlink $largetermpep;
}

if (-z "$portalpep" == 1)  {
    print "removing $portalpep ...\n" if ($DEBUG);
    unlink $portalpep;
}

if (-z "$proteasepep" == 1)  {
    print "removing $proteasepep ...\n" if ($DEBUG);
    unlink $proteasepep;
}

if (-z "$scaffpep" == 1)  {
    print "removing $scaffpep ...\n" if ($DEBUG);
    unlink $scaffpep;
}

if (-z "$capsidpep" == 1)  {
    print "removing $capsidpep ...\n" if ($DEBUG);
    unlink $capsidpep;
}

if (-z "$serinerecombpep" == 1)  {
    print "removing $serinerecombpep ...\n" if ($DEBUG);
    unlink $serinerecombpep;
}

if (-z "$tyrosinerecombpep" == 1)  {
    print "removing $tyrosinerecombpep ...\n" if ($DEBUG);
    unlink $tyrosinerecombpep;
}

if (-z "$toxinspep" == 1)  {
    print "removing $toxinspep ...\n" if ($DEBUG);
    unlink $toxinspep;
}

if (-z "$apisfile" == 1)  { # if zero contents, remove
    print "removing $apisfile ...\n" if ($DEBUG);
    unlink $apisfile;
}

if (-z "$confile" == 1)  { # remove empty files
    print "removing $confile...\n" if ($DEBUG);
    unlink $confile;
}
if (-z "$seqfile" == 1)  {
    print "removing $seqfile ...\n" if ($DEBUG);
    unlink $seqfile;
} 
if (-z "$pepfile" == 1)  {
    print "removing $pepfile ...\n" if ($DEBUG);
    unlink $pepfile;
}

exit(0);
