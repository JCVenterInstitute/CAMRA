package PHAGE::Phage_subs;

=head1 NAME

PHAGE::Phage_subs - Phage_finder module

=head1 LICENSE

Copy (C) 2023-2024  The J. Craig Venter Institute (JCVI).  All rights reserved

Written by Derrick E. Fouts, Ph.D.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 CITATION

Fouts, D. E. (2006) "Phage_Finder: automated identification and classification
of prophage regions in complete bacterial genome sequences." 
Nucleic Acids Res 34:5839-51. PMCID: PMC1635311.

=head1 SYNOPSIS

  use PHAGE::Phage_subs;

=head1 DESCRIPTION

This module provides many core subroutines for the phage_finder.pl program that handle filehandling, data processing, and log files.

=cut

   BEGIN {
      require 5.006_00; # error if using Perl < v5.6.0
   }

   use strict;
   use Math::Round qw(nlowmult nhimult);
   require Exporter;

   our @ISA = ('Exporter');
   our @EXPORT = qw(
                    get_gc_content
                    load_hash_from_file
                    get_HMM_size_trusted_noise
                    get_ok_comnames
                    get_DB_info
                    make_basedir
                    create_dir
                    populate_genomearray
                    create_log
                    write_log
                    dna2pep
                    quickcut
		    printatt
                    printFasta
                    check_infofile
                    get_gene_info
                    get_tRNAs
                    get_tmRNAs
                    get_ncRNAs
                    select_featnames_from_btab
                    find_hmms
                    get_assemblies
                    get_TSD
                    add_ORFs_to_region
                    adjust_featname
                    rekey_regions
                    determine_region_type
                   );

    our @EXPORT_OK = qw(
                        get_gc_content
                        load_hash_from_file
                        get_HMM_size_trusted_noise
                        get_ok_comnames
                        get_DB_info
                        time_date
                        make_basedir
                        create_dir
                        populate_genomearray
                        create_log
                        write_log
                        revcomp
                        codon2aa
                        dna2pep
                        quickcut
			printatt
                        print_sequence
                        printFasta
                        pttcoords_to_info
                        check_infofile
                        get_gene_info
                        get_tRNAs
                        get_tmRNAs
                        get_ncRNAs
                        select_featnames_from_btab
                        find_hmm
                        get_fasta_record
                        get_assemblies
                        populate_atthash
                        pick_longest_length
                        pick_best_att
                        get_TSD
                        adjust_featname
                        add_ORFs_to_region
                        rekey_regions
                        determine_region_type
                       );

    our $VERSION = '2.5';

    our %EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

sub get_gc_content {
	my $seq = shift;
	$seq =~ s/^>[^\n]*\n//;
	$seq =~ s/\n//g;	
	my $gc=($seq=~tr/[GgCc]/G/);
	$gc=($gc/length($seq))*100;
	$gc=sprintf("%.2f",$gc);
	return $gc;
}

sub load_hash_from_file {

    my ($path,$fileref,$exclude_ref,$lists_ref,$HMM_ref,$DEBUG) = @_;
    my $file = "";
    my $descriptor = "";
    for $file (@{$fileref})  {
	print "LOAD_FILES:  filename = $file\n" if ($DEBUG);
        open (FILE, "<$path/$file") || &write_log("4","can't read file $file: $!\n");
	while (<FILE>)  {
	    chomp;
            if (/^#(.*)$/) {  #capture the file descriptor line beginning with #
		$descriptor = $1;
                print "Getting list of $descriptor accessions from $file . . .\n" if ($DEBUG);
                &write_log("1","Getting list of $descriptor accessions from $file . . .");
	    }
            elsif ($descriptor =~ /EXCLUDE/)  {
		$exclude_ref->{$_} = 1;  # populate exclude hash (one dimensional)
	    }
            elsif ($descriptor =~ /HMM/)  { # only HMM accessions here
		$HMM_ref->{$descriptor}{$_} = 1;  # populate HMMs_hash two-dimensional hash
            }
            else  { # other listfile data (curated Genbank accessions, etc)
		$lists_ref->{$descriptor}{$_} = 1;  # populate lists_hash two-dimensional hash
            }
	    print "<$descriptor>, !$_!\n" if ($DEBUG);
	}
	close (FILE);
    }
}

sub get_HMM_size_trusted_noise  {
    my ($phome,$HMMversion,$hash) = @_;
    my @a = ();
    my @files = ();
    my @model = ();
    my $key = "";
    my $list = "";
    my $i = "";
    
    if ($HMMversion == "2") {
	push (@files, "$phome/HMM_master.lst");
	push (@files, "$phome/HMM_master_FRAG.lst");
	
	push (@model, "GLOCAL");
	push (@model, "FRAG");
    }
    elsif ($HMMversion == "3") {
	push (@files, "$phome/HMM3_master.lst");
	push (@model, "GLOCAL");
    }
    for ($i=0; $i<=$#files; $i++) {
	open (IN, "<$files[$i]") || &write_log("4","can't open file $files[$i]: $!\n");
	while (<IN>)  {
	    chomp;
	    @a=split(/\t+/);
	    $key = $a[0];
	    $hash->{$model[$i]}{$key}->{'hmm_len'} = $a[1];
	    $hash->{$model[$i]}{$key}->{'trusted'} = $a[2];
	    $hash->{$model[$i]}{$key}->{'noise'} = $a[3];
	    if ($HMMversion == "2")  {
		if (($key eq "PF00589") && ($model[$i] eq "GLOCAL")) { $hash->{$model[$i]}{$key}->{'noise'} = "-44.2"; } # analysis of hmm searches from phage db enables lowering noise (XF0480 from NC_002488 with a score of -43.2
		if (($key eq "PF02316") && ($model[$i] eq "GLOCAL")) { $hash->{$model[$i]}{$key}->{'noise'} = "-45.4"; } # analysis of hmm searches from phage db enables lowering noise (next non-Mu phage is PROPHAGE_pseudo_gpp4_UNPUB-ORF06021 with a score of -50.0
		if (($key eq "PF02914") && ($model[$i] eq "GLOCAL")) { $hash->{$model[$i]}{$key}->{'noise'} = "-139.7"; } # analysis of hmm searches from phage db enables lowering noise (next non-Mu phage is PHAGE_Burkho_BcepNazgul-gi|34610172|ref|NP_918997.1| with a score of -140.7
		if (($key eq "PF06074") && ($model[$i] eq "GLOCAL")) { $hash->{$model[$i]}{$key}->{'noise'} = "-370.9"; } # analysis of hmm searches from phage db enables lowering noise (next non-Mu phage is PHAGE_lactob_A2-gi|22296534|ref|NP_680494.1| with a score of -380.9
		if (($key eq "PF07030") && ($model[$i] eq "GLOCAL")) { $hash->{$model[$i]}{$key}->{'noise'} = "-58.0"; } # analysis of hmm searches from phage db enables lowering noise (next non-Mu phage is PHAGE_strepm_VWB-gi|41057255|ref|NP_958281.1| with a score of -69.1 then found a phage in Shigella that was clearly NOT Mu but had a score of -58.1
	    }
	}
	close (IN);
    }
}

sub get_ok_comnames {

    my ($comfile,$okref,$DEBUG) = @_;
    my @temp = ();
    my $name = "";
    open (COMFILE, "<$comfile") || &write_log("4","can't read file $comfile: $!\n");
    while (<COMFILE>)  {
	chomp;
        @temp = split(/\s+/);
	#$name = lc($temp[0]); # added 02/21/2014 to be same as hitref clean_names so ok com_names are not considered bad because of case
	$name = lc($temp[0] . $temp[1] . $temp[2]); #edited 04/13/2022 by dfouts to add up to 3 words to each com_name to make them more specific (i.e., type ii toxin vs. type ii secretion)
	$name =~ s/\(.*$//; # added 04/13/2022 to remove (locus_id) present on some accessions
	$name =~ s/\W//g; # remove all non-word characters
        if (!exists $okref->{$name})  {
	    print "ok_com_name: >$name<\n" if ($DEBUG);
	    $okref->{$name} = 1;
        }
    }
    close (COMFILE);
}

sub get_DB_info {

    my ($DBfile,$DBref,$DEBUG) = @_;
    my @temp = ();
    my $tag = "";
    open (DBFILE, "<$DBfile") || &write_log("4","can't read file $DBfile: $!\n");
    while (<DBFILE>)  {
	chomp;
        @temp = split(/\t+/);
	$tag = $temp[0];
        if (!exists $DBref->{$tag})  {
	    $DBref->{$tag}->{'name'} = $temp[1];
	    $DBref->{$tag}->{'taxon'} = $temp[2];
	    print "get_DB_info: >$tag<\t>$DBref->{$tag}->{'name'}<\t>$DBref->{$tag}->{'taxon'}<\n" if ($DEBUG);
        }
    }
    close (DBFILE);
}

sub time_date {
    my($sec, $min, $hr, $mday, $m, $year, $wday, $time, $date);
    my(%month, %day);
    ($sec, $min, $hr, $mday, $m, $year, $wday) = localtime(time());
    $month{0}="January"; $month{1}="February"; $month{2}="March";
    $month{3}="April"; $month{4}="May"; $month{5}="June";
    $month{6}="July"; $month{7}="August"; $month{8}="September";
    $month{9}="October"; $month{10}="November"; $month{11}="December";
    $day{0}="Sunday"; $day{1}="Monday"; $day{2}="Tuesday";
    $day{3}="Wednesday"; $day{4}="Thursday"; $day{5}="Friday";
    $day{6}="Saturday";
    
    $time = "$hr:$min:$sec";
    $date = "$day{$wday}, $month{$m} $mday";
    return ($time, $date);
}

sub make_basedir { # new subroutine (to ver. 2.3) to create a mode-based master base directory to store all Phage_Finder-generated files (i.e., "strict" vs. "non"strict")
    my ($basedir, $mode, $DEBUG) = @_;
    my $new_basedir = $basedir . "/" . $mode . "_dir"; # set the write directory
    print "MAKE_BASEDIR:  new_basedir = <$new_basedir> ...\n" if ($DEBUG);
    if (-e "$new_basedir" == 0) {  # check if dir present
	print "MAKE_BASEDIR:  Making new <$mode> directory ...\n" if ($DEBUG);
	system("mkdir $new_basedir");
    }
    elsif (-s "$new_basedir" == 0){  # if present and has contents, then remove old data)
	print "MAKE_BASEDIR:  Directory <$new_basedir> has contents, so deleting for fresh run ...\n" if ($DEBUG);
	system("rm $new_basedir/$mode\.*");
    }
    return ($new_basedir);
}

sub create_dir { # subroutine to create a directory to hold phage_finder data (if not already present)
  my ($basedir,$new_path) = @_;
  my $write_dir = $basedir . "/" . $new_path . "_dir"; # set the write directory
  if (-e "$write_dir" == 0) {  # check for asmbl_id dir present
      system("mkdir $write_dir");
  }
  elsif (-s "$write_dir" == 0){  # if present and has contents, then remove old data)
      system("rm $write_dir/$new_path\.*");
  }
  return ($write_dir);
}

sub populate_genomearray {

  my ($reref,$aref,$hitref,$asmbl_id) = @_;
  my $keys;
  my $pos = 0;
  my $feat_name = "";
  foreach $keys (sort {$a <=> $b} keys %{$reref->{$asmbl_id}}) {
    push(@{$aref}, $keys); # inserting the 5' end of each ORF into array
    $feat_name = &adjust_featname($reref->{$asmbl_id}{$keys}->{'featname'},$asmbl_id);
    $hitref->{$feat_name}->{'array_pos'} = $pos; # a map of the array position of the feat_name so we can get either the next 5' or 3' ORF of a particular feat_name
    $pos++;
  }
}

sub create_log {
  my ($logfile) = @_;
  my $i = 1;
  if (-e "$logfile.log" == 1) {  # check for phage_phinder.log present
    until (-e "$logfile\_$i.log" == 0)  {  #add 1 to the name until not present
      $i++;
    }
    open (LOG, ">$logfile\_$i.log") || &write_log("4","can't open file $logfile\_$i.log: $!\n");
    return ("$logfile\_$i.log");
  }
  else  {  # else, make phage_phinder.log for the first time
    open (LOG, ">$logfile.log") || &write_log("4","can't open file $logfile.log: $!\n");
    return ("$logfile.log");
  }
}

sub write_log {

    my ($state, $text, $asmbl_id) = @_;
    my ($time, $date);
    #
    # state crib sheet #
    # 0 = begin or start and return
    # 1 = write general message to log file and return
    # 2 = end or finish under good status, phages found, return
    # 3 = end or finish under good status, no phages found, return
    # 4 = end or finish with error and exit(1)
    # 5 = end or finish under good status, no blast results = no phages, exit(0)
    #
    ($time,$date) = &time_date;
    if ($state == "0")  {
        print LOG "..................................................................................................\n";
        print LOG "$text\n";
    }
    elsif ($state == "1")  {
	print LOG "$time, $date | $text ...\n";
    }
    elsif ($state == "2")  {
        print LOG "$time, $date | FINISH: $asmbl_id\n";
        close(LOG) if ($text == "0");
    }
    elsif ($state == "3")  {
        print LOG "$time, $date | Sorry, no phages in contig $asmbl_id :( ...\n";
        print LOG "$time, $date | FINISH: $asmbl_id\n";
        close(LOG) if ($text == "0");
    }
    elsif ($state == "4")  {
        print LOG "$time, $date | FINISH: ERROR - $text!\n";
        close(LOG);
        exit(1);
    }
    elsif ($state == "5")  { # added 03/07/14 so we can exit nicely when no BLAST results found
	print LOG "$time, $date | FINISH: $text.\n";
	close(LOG);
	exit(0);
    }
    return;
}

sub revcomp {

    my($dna) = @_;
    # First reverse the sequence
    my $revcom = join("", reverse(split(/ */, $dna)));  # reverse att site sequence (from Bob Deboy's SOS.dbi script :))
    # Next, complement the sequence, dealing with upper and lower case and some IUB codes
    # A->T, T->A, C->G, G->C, y->r, r->y, k->m, m->k
    $revcom =~ tr/ACGTacgtyrkm/TGCAtgcarymk/; # complement att site to generate the reverse complement (from Bob's script again)
    return $revcom;
}

# From Chapter 8 
# Beginning Perl for Bioinformatics
# by James Tisdall

#
# codon2aa
#
# A subroutine to translate a DNA 3-character codon to an amino acid
#   Version 3, using hash lookup

sub codon2aa {
    my ($codon) = @_;
    my %genetic_code = ();

    $codon = uc $codon;
 
    (%genetic_code) = (
    
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TCR' => 'S',    # Serine (Wobble)
    'TCY' => 'S',    # Serine (Wobble)
    'TCM' => 'S',    # Serine (Wobble)
    'TCK' => 'S',    # Serine (Wobble)
    'TCS' => 'S',    # Serine (Wobble)
    'TCW' => 'S',    # Serine (Wobble)
    'TCH' => 'S',    # Serine (Wobble)
    'TCB' => 'S',    # Serine (Wobble)
    'TCV' => 'S',    # Serine (Wobble)
    'TCD' => 'S',    # Serine (Wobble)
    'TCN' => 'S',    # Serine (Wobble)
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTY' => 'F',    # Phenylalanine (Wobble)
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TTR' => 'L',    # Leucine (Wobble)
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAY' => 'Y',    # Tyrosine (Wobble)
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGY' => 'C',    # Cysteine (Wobble)
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CTR' => 'L',    # Leucine (Wobble)
    'CTY' => 'L',    # Leucine (Wobble)
    'CTM' => 'L',    # Leucine (Wobble)
    'CTK' => 'L',    # Leucine (Wobble)
    'CTS' => 'L',    # Leucine (Wobble)
    'CTW' => 'L',    # Leucine (Wobble)
    'CTH' => 'L',    # Leucine (Wobble)
    'CTB' => 'L',    # Leucine (Wobble)
    'CTV' => 'L',    # Leucine (Wobble)
    'CTD' => 'L',    # Leucine (Wobble)
    'CTN' => 'L',    # Leucine (Wobble)
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CCR' => 'P',    # Proline (Wobble)
    'CCY' => 'P',    # Proline (Wobble)
    'CCM' => 'P',    # Proline (Wobble)
    'CCK' => 'P',    # Proline (Wobble)
    'CCS' => 'P',    # Proline (Wobble)
    'CCW' => 'P',    # Proline (Wobble)
    'CCH' => 'P',    # Proline (Wobble)
    'CCB' => 'P',    # Proline (Wobble)
    'CCV' => 'P',    # Proline (Wobble)
    'CCD' => 'P',    # Proline (Wobble)
    'CCN' => 'P',    # Proline (Wobble)
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAY' => 'H',    # Histidine (Wobble)
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CAR' => 'Q',    # Glutamine (Wobble)
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'CGR' => 'R',    # Arginine (Wobble)
    'CGY' => 'R',    # Arginine (Wobble)
    'CGM' => 'R',    # Arginine (Wobble)
    'CGK' => 'R',    # Arginine (Wobble)
    'CGS' => 'R',    # Arginine (Wobble)
    'CGW' => 'R',    # Arginine (Wobble)
    'CGH' => 'R',    # Arginine (Wobble)
    'CGB' => 'R',    # Arginine (Wobble)
    'CGV' => 'R',    # Arginine (Wobble)
    'CGD' => 'R',    # Arginine (Wobble)
    'CGN' => 'R',    # Arginine (Wobble)
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATY' => 'I',    # Isoleucine (Wobble)
    'ATM' => 'I',    # Isoleucine (Wobble)
    'ATW' => 'I',    # Isoleucine (Wobble)
    'ATH' => 'I',    # Isoleucine (Wobble)
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'ACR' => 'T',    # Threonine (Wobble)
    'ACY' => 'T',    # Threonine (Wobble)
    'ACM' => 'T',    # Threonine (Wobble)
    'ACK' => 'T',    # Threonine (Wobble)
    'ACS' => 'T',    # Threonine (Wobble)
    'ACW' => 'T',    # Threonine (Wobble)
    'ACH' => 'T',    # Threonine (Wobble)
    'ACB' => 'T',    # Threonine (Wobble)
    'ACV' => 'T',    # Threonine (Wobble)
    'ACD' => 'T',    # Threonine (Wobble)
    'ACN' => 'T',    # Threonine (Wobble)
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAY' => 'N',    # Asparagine (Wobble)
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AAR' => 'K',    # Lysine (Wobble)
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGY' => 'S',    # Serine (Wobble)
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'AGR' => 'R',    # Arginine (Wobble)
    'MGA' => 'R',    # Arginine (First position)
    'MGG' => 'R',    # Arginine (First position)
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GTR' => 'V',    # Valine (Wobble)
    'GTY' => 'V',    # Valine (Wobble)
    'GTM' => 'V',    # Valine (Wobble)
    'GTK' => 'V',    # Valine (Wobble)
    'GTS' => 'V',    # Valine (Wobble)
    'GTW' => 'V',    # Valine (Wobble)
    'GTH' => 'V',    # Valine (Wobble)
    'GTB' => 'V',    # Valine (Wobble)
    'GTV' => 'V',    # Valine (Wobble)
    'GTD' => 'V',    # Valine (Wobble)
    'GTN' => 'V',    # Valine (Wobble)
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GCR' => 'A',    # Alanine (Wobble)
    'GCY' => 'A',    # Alanine (Wobble)
    'GCM' => 'A',    # Alanine (Wobble)
    'GCK' => 'A',    # Alanine (Wobble)
    'GCS' => 'A',    # Alanine (Wobble)
    'GCW' => 'A',    # Alanine (Wobble)
    'GCH' => 'A',    # Alanine (Wobble)
    'GCB' => 'A',    # Alanine (Wobble)
    'GCV' => 'A',    # Alanine (Wobble)
    'GCD' => 'A',    # Alanine (Wobble)
    'GCN' => 'A',    # Alanine (Wobble)
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAY' => 'D',    # Aspartic Acid (Wobble)
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GAR' => 'E',    # Glutamic Acid (Wobble)
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    'GGR' => 'G',    # Glycine (Wobble)
    'GGY' => 'G',    # Glycine (Wobble)
    'GGM' => 'G',    # Glycine (Wobble)
    'GGK' => 'G',    # Glycine (Wobble)
    'GGS' => 'G',    # Glycine (Wobble)
    'GGW' => 'G',    # Glycine (Wobble)
    'GGH' => 'G',    # Glycine (Wobble)
    'GGB' => 'G',    # Glycine (Wobble)
    'GGV' => 'G',    # Glycine (Wobble)
    'GGD' => 'G',    # Glycine (Wobble)
    'GGN' => 'G',    # Glycine (Wobble)
    'RAC' => 'B',    # Asparagine or Aspartic acid (IUPAC)
    'RAT' => 'B',    # Asparagine or Aspartic acid (IUPAC)
    'SAA' => 'Z',    # Glutamine or Glutamic acid (IUPAC)
    'SAG' => 'Z',    # Glutamine or Glutamic acid (IUPAC)
    'TAA' => '_',    # Stop
    'TAG' => '_',    # Stop
    'TGA' => '_',    # Stop
    );

    if (!exists $genetic_code{$codon}) {
        $genetic_code{$codon} = "X";
        print STDERR "Bad codon $codon => $genetic_code{$codon}!!\n";
    }
    return $genetic_code{$codon};
}

sub dna2pep {

    my($dna) = @_;

    # Initialize variables
    my $protein = '';

    # Translate each three-base codon to an amino acid, and append to a protein 
    for(my $i=0; $i < (length($dna) - 2) ; $i += 3) {
        $protein .= &codon2aa( substr($dna,$i,3) );
    }

    $protein =~ s/^./M/; # make the first codon Methionine
    $protein =~ s/_$//; # remove the stop codon
 
    return $protein;
}

sub quickcut {

  my($sequence, $pos1, $pos2) = @_;
  my $start;
  my $length;
  my $outseq;
  if ($pos1 > $pos2)  { # <--- so revcomp seq
      $start = $pos2;
      $length = ($pos1-$pos2)+1;
      $outseq = &revcomp(substr($sequence, $start-1, $length));
  }
  else {  # ---> 
      $start = $pos1;
      $length = ($pos2-$pos1)+1;
      $outseq = substr($sequence, $start-1, $length);
  }
  return $outseq;
}

sub printatt {
  #Gene attribute file format:
  #End5   End3    feat_name       Gene_sym        Locus   feat_type       Role    Fragment_id
  # example output:
  #271281	270118	G952_RS22620-NZ_AWCV01000000_NZ_KE702177_PFPR02		G952_RS22620	ORF	1002	NZ_KE702177
  # need featname, annotation, end5, end3, asmbl_id, role_id passed to this subroutine
  #*ATTFILE, $end5, $end3, $att_featname, $featname, $role_id, $region_tag, $hitref->{$featname}->{'com_name'} passed from Phage_finder.pl
  my ($file,$end5,$end3,$asmbl_id,$locus,$feat_type,$role_id,$region_tag,$anno) = @_;
  my $feat_name = $locus;
  my $frag_id = $asmbl_id . "_" . $region_tag;

  if (($locus =~ /tm?RNA/) || ($locus =~ /ncRNA/)) { # don't remove dash and beyond if tRNA, tmRNA or ncRNA and keep orignal locus name (tRNA-Ile-1);
    $feat_name =~ s/-//g; # just remove the dashes (e.g., tRNA-Ile-1 will become tRNAIle1)
  }
  else {
    $feat_name =~ s/-.*//; # remove everything including and after the dash to get just the locus name without the asmbl_id
    $locus = $feat_name;
  }

  $feat_name = $feat_name . "-" . $asmbl_id . "_" . $region_tag;
  if ($role_id == "2000") { # if role_id not assigned from HMMs and mini BLAST db results, then assign based on annotation
     if (($anno =~ /small terminase/) || ($anno =~ /ters/) || ($anno =~ /small subunit/)) {$role_id = "1005";}
     elsif (($anno =~ /dna packaging/) || ($anno =~ /large terminase/) || ($anno =~ /terl/) || ($anno =~ /large subunit/)) {$role_id = "1006";}
     elsif ($anno =~ /portal/) {$role_id = "1007";}
     elsif (($anno =~ /scaffold/) || ($anno =~ /protease/)) {$role_id = "1008";}
     elsif (($anno =~ /major capsid/) || ($anno =~ /major head/)) {$role_id = "1009";}
     elsif (($anno =~ /minor capsid/) || ($anno =~ /minor head/)){$role_id = "1017";}
     elsif ($anno =~ /head-tail/) {$role_id = "1026";}
     elsif (($anno =~ /tape measure/) || ($anno =~ /tmp/) || ($anno =~ /pbla/)) {$role_id = "1010";}
     elsif (($anno =~ /minor structural/) || ($anno =~ /spike/) || ($anno =~ /fiber/)) {$role_id = "1033";}
     elsif ($anno =~ /tail/) {$role_id = "1011";}
     elsif (($anno =~ /host specificity/) || ($anno =~ /antireceptor/) || ($anno =~ /anti-receptor/) || ($anno =~ /receptor-binding/) || ($anno =~ /rbp/)) {$role_id = "1035";}
     elsif (($anno =~ /amidase/) || ($anno =~ /hydrolase/) || ($anno =~ /lysin/) || ($anno =~ /ply/) || ($anno =~ /lysozyme/)) {$role_id = "1012";}
     elsif ($anno =~ /holin/) {$role_id = "1013";}
     elsif (($anno =~ /antirepressor/) || ($anno =~ /anti-repressor/)) {$role_id = "1015";}
     elsif (($anno =~ /antiterm/) || ($anno =~ /anti-term/)) {$role_id = "1034";}
     elsif (($anno =~ /resistance/) || ($anno =~ /carbapenemase/) || ($anno =~ /toxin/)) {$role_id = "1025";}
     elsif (($anno =~ /nucleotidohydrolase/) || ($anno =~ /dutpase/) || ($anno =~ /thymidylate/) || ($anno =~ /ribonucleotide reductase/) || ($anno =~ /thioredoxin/) || ($anno =~ /glutaredoxin/) || ($anno =~ /deaminase/) || ($anno =~ /nrdh/)) {$role_id = "1027";}
     elsif (($anno =~ /helix-turn-helix/) || ($anno =~ /repressor/) || ($anno =~ /transcriptional/) || ($anno =~ /regulatory/)) {$role_id = "1000";}
     elsif (($anno =~ /spoIIIE/) || ($anno =~ /ftsk/) || ($anno =~ /segregation/)) {$role_id = "1034";}
     elsif (($anno =~ /sigma/) || ($anno =~ /rnap/) || ($anno =~ /rna polymerase/)) {$role_id = "1029"}
     elsif (($anno =~ /replication/) || ($anno =~ /replisome/) || ($anno =~ /primase/) || ($anno =~ /helicase/) || ($anno =~ /sak/) || ($anno =~ /mcm/) || ($anno =~ /polymerase/) || ($anno =~ /dna. /) || ($anno =~ /dna pol/) || ($anno =~ /topoisomerase/) || ($anno =~ /rII./)) {$role_id = "1003";}
     elsif (($anno =~ /virulence/) || ($anno =~ /leukocidin/) || ($anno =~ /vire/) || ($anno =~ /luk/) || ($anno =~ /dismutase/) || ($anno =~ /yopx/) || ($anno =~ /cupin/) || ($anno =~ /collagen/) || ($anno =~ /6-phosphate/)) {$role_id = "1014";}
     elsif (($anno =~ /site-specific/) && (($anno =~ /recombinase/) || $anno =~ /integrase/) || ($anno =~ /integrase/)) {$role_id = "1002";}
     elsif ($anno =~ /excis/) {$role_id = "1032";}
     elsif (($anno =~ /reca/) || ($anno =~ /invertion/) || ($anno =~ /recombinase/) || ($anno =~ /rec. /)) {$role_id = "1001";}
     elsif (($anno =~ /hnh endonuclease/) || ($anno =~ /homing/)) {$role_id = "1024";}
     elsif (($anno =~ /methyl/) || ($anno =~ /restriction/)) {$role_id = "1004";}
     elsif (($anno =~ /single/) || ($anno =~ /dna binding/) || ($anno =~ /dna-binding/)) {$role_id = "1016";}
     elsif ($anno =~ /ligase/) {$role_id = "1022";}
     elsif (($anno =~ /transposon/) || ($anno =~ /transposase/)) {$role_id = "1023";}
     elsif ($anno =~ /conserved/) {$role_id = "1020";}
     elsif (($anno =~ /hypothetical/) || ($anno =~ /unknown/))  {$role_id = "1028";}
  }
  print $file "$end5\t$end3\t$feat_name\t\t$locus\t$feat_type\t$role_id\t$frag_id\n";
  # to do list:  generate the PFPR.frag file with the beginning and ending of each phage region to draw
}

sub print_sequence  {
    my $file = shift;
    my $seqs = shift;

    for (my $j = 0; $j < length($seqs); $j += 60){
        print $file substr($seqs, $j, 60), "\n";
    }
}

sub printFasta
{
    my($file) = $_[0];
    my($header) = $_[1];
    my($seqs) = $_[2];

    print $file ">$header\n";
    &print_sequence($file, $seqs);
}

sub pttcoords_to_info { # convert Genbank .ptt (.coords) file to phage_finder_info.txt format
    my ($basedir,$seqfile,$pttfile,$default_asmbl_id) = @_;
    my $infofile = "$basedir/phage_finder_info.txt";
    my @a = ();
    my @id = ();
    my @seq = ();
    my @title = ();
    my @length = ();
    my $t = "";
    my $sequence = "";
    my $record = "";
    my $line = "";
    my @b = ();
    my $c = "";
    my $size = "";
    my $end5 = "";
    my $end3 = "";
    my $locus = "";
    my $COG = "";
    my $syn = "";
    my $anno = "";
    my $accession = "";
    my $asmbl_id = "";

### first determine the largest sequence in the supplied seqfile and the asmbl_id for it
# do not store any data #
    if ((-e "$seqfile" == 1) && (-z "$seqfile" == 0)) {
	open (GENOME, "<$seqfile") || &write_log("4","can't find the sequence file $seqfile:$!\n");
	while ($record = &get_fasta_record(*GENOME)) {
	    ($t, @seq) = split(/\n/, $record);
	    $sequence = join('', @seq);
	    @id = split(" ", $t);
            if ($id[0] =~ /.*\|(NC_\d+)/){
		$asmbl_id = $1;
            }
            elsif ($id[0] =~ /.*gb\|(\w+)/){ # 04/09/12 v2.1 fix dfouts
		$asmbl_id = $1;
            }
	    else {
		$asmbl_id = $id[0];
	    }
	    if (length($sequence) > $length[$#length])  {
		push(@title, $asmbl_id);
		push(@length, length($sequence));
	    }
	}
	close(GENOME);
        $asmbl_id = pop(@title);
	$size = pop(@length);
    }
    else  { # case where no sequence is provided
        if (defined($default_asmbl_id))  {
	    $asmbl_id = $default_asmbl_id;
	}
        else  {
            $asmbl_id = "1"; # make up an asmbl_id if none provided and no sequence is provided
        }
        open (INFILE, "<$pttfile") || &write_log("4","can't open the ptt file $pttfile: $!\n");
	$line = <INFILE>; #read first line from .ptt file to get the size
	close (INFILE);
	if ($line =~ /.*\d+\..(\d+)/)  { #04/09/2012 fixed bug preventing .ptt file usage dfouts 
	    $size = $1;
	}
    }
    open (OUTFILE, ">$infofile") || &write_log("4","can't create file $infofile: $!\n");

### now parse the .ptt coords file ###
    open (INFILE, "<$pttfile") || &write_log("4","can't open the ptt file $pttfile: $!\n");
    while ($line = <INFILE>)  {
        chomp($line);
        @a = split(/\s+/, $line);
        if ($a[0] !~ /\d+\..\d+/) {next;} # if not number..number then skip
        @b=split(/\../,$a[0]);
        $accession = $a[3];
        if ($a[4] ne "-") {
	    $syn = $a[4];
        }
        else { $syn = undef; }
        $locus = $a[5];
        if ($a[7] ne "-") {
	    $COG = $a[7];
	}
	else { $COG = undef; }
	$anno = join(" ", @a[8..$#a]) . "; " . $locus;
        if ( $a[1] eq "+" ) {
             $end5 = $b[0];
	     $end3 = $b[1];
        }
        else  {
	    $end5 = $b[1];
	    $end3 = $b[0];
        }
        print OUTFILE "$asmbl_id\t$size\t$accession\t$end5\t$end3\t$anno";
        if ($syn) {
	    print OUTFILE "; $syn";
	}
	if ($COG) {
	    print OUTFILE "; $COG";
	}
	print OUTFILE "\n";
    }
    close (OUTFILE);
    close (INFILE);
}

sub check_infofile {
    my ($basedir,$infofile,$seqfile,$default_asmbl_id,$DEBUG) = @_;
    my @a = ();
    my $rvalue = "";
    open (INFILE, "<$infofile") || &write_log("4","can't open infofile $infofile: $!\n");
    chomp;
    @a = split(/\s+/, <INFILE>);
    if ($a[1] =~ /^\D/) {
        #print "Looks like a Genbank .ptt file, checking for phage_finder.info file ...\n";
        #&write_log ("1","Looks like a Genbank .ptt file, converting to phage_finder_info");
        print "$infofile\n" if ($DEBUG);
	if (($infofile =~ /phage_finder_info.txt/) && (-e "$infofile") && (-z "$infofile" == 0)) { # if file exists and value is > 0, then use it instead of making a new one
            close(INFILE);
            print "$infofile found, using this file ...\n";
            &write_log ("1","$infofile found, using this file ...");
	}
        else  {
            print "Creating phage_finder_info.txt from $infofile ...\n";
            &write_log ("1","Creating phage_finder_info.txt from $infofile");
            print "!$infofile!\t!$seqfile!\n" if ($DEBUG);
	    &pttcoords_to_info($basedir,$seqfile,$infofile,$default_asmbl_id);
	    $infofile = "phage_finder_info.txt";  # added 04/09/2012 by dfouts to fix a bug reading from .ptt file
        }
    }
    else  {
        print "Assuming file to have correct phage_finder_info format ...\n" if ($DEBUG);
	&write_log ("1","Assuming file to have correct phage_finder_info format");
    }
    close(INFILE);
    return($infofile);
}

sub get_gene_info {
  my ($infofile,$asmblref,$hitref,$reref,$DEBUG) = @_;
  my @d = ();
  my @a = ();
  my @temp = ();
  my $count = "";
  my $name = "";
  # hash "rehash" has asmbl_id and end5 as keys, featname as value
  # can link $rehash to $hithash using featname value
  # rehash will be used to go back and pick out more potential phage hits within a region and on the flanks if their annotations are "phage-like"
  open (INFOFILE, "<$infofile") || &write_log("4","can't open file $infofile: $!\n");
  while (<INFOFILE>)  {
      chomp;
      @d=split(/\t/); # removed the "+" after the tab delimiter because sometimes the size of the contig (perl column 1) is missing
      if (($d[3] eq "") || ($d[4] eq "")) {
	next; # skip if gene with missing coordinates added 01/22/2019 by Derrick Fouts
      }
      if (!defined $asmblref->{$d[0]}) { # check only that asmbl_id is novel and count only new asmbl_ids
	  $asmblref->{$d[0]}->{'genomesize'} = $d[1];
	  $count++;
      }
      if ($d[2] =~ /^ncRNA/) { # added 01/16/22 to avoid having multiple ncRNAs (from the GenBank file if present and from the cmscan output that is now required)
	  next;
      }
      ## code not functioning properly due to \.* instead of .*.  Thanks Phil Goetz for finding this. On 07/17/14, discovered this is needed so that genome tag can be captured in .seq and .pep files.  Changed \.* to [\w_]+.
      if (($d[2] =~ /^\bgi\b?\|?(\d+)\-([\w_]+)/) && !exists ($hitref->{$d[2]}->{'genome_tag'}))  { # new code 4/6/2010 to include the genome tags if present, chaged from \w+ to \.* 04/16/13, changed to [\w_]+ and also from gi to [gi]? to cover cases of locus names or feat_names being used
	  $d[2] = $1;
	  $hitref->{$d[2]}->{'genome_tag'} = $2;
      }
      elsif ($d[2] =~ /^gi\|?(\d+)/)  {  # if gi numbers are present, we need to only look at the numbers, not the gi| or the :end5..end3 garbage. changed 07/17/14 to not loose the genome identifier after the dash in the pep files
	  $d[2] = $1;
      }
      elsif ($d[2] =~ /^(\S+_gi\d+)-/) {  # if (gene_sym/locus_tag)_gi123456-NC_bla
          $d[2] = $1;
      }
      if (!exists $hitref->{$d[2]}->{'featname'})  {
	$hitref->{$d[2]}->{'asmbl_id'} = $d[0];
        $hitref->{$d[2]}->{'end5'} = $d[3]; # hash "hithash" has featname as key, end5 as value
        $hitref->{$d[2]}->{'end3'} = $d[4];
        $d[5] =~ s/\s\s+/ /g;  # remove stretch of multiple spaces between comname and [organism]
        $d[5] =~ s/  / /g;  # remove any double spaces
        $d[5] =~ s/\{/\[/g; # if we have {Pseudomonas bla}, make [Pseudomonas bla]
        $d[5] =~ s/\}/\]/g;
        if ($d[5] =~ /\[/)  {
          @a=split(/\[/, $d[5]);
          $a[1] =~ s/^ *//; # remove any spaces at the beginning
          $a[1] =~ s/\s*$//; # remove any spaces at the end of the line
          $a[1] = lc($a[1]); # shift to lower case
          $hitref->{$d[2]}->{'organism'} = "[" . $a[1];
        }
        else {
          $hitref->{$d[2]}->{'organism'} = "[No organism specified]";
        }
        $d[5] =~ s/^ *//; # remove any spaces at the beginning
        $d[5] =~ s/\s*$//; # remove any spaces at the end of the line
        $d[5] = lc($d[5]); # shift to lower case
        $hitref->{$d[2]}->{'com_name'} = $d[5];

        # truncate com_name to first 3 words (separated by spaces)
	@temp = split(/\s+/, $d[5]);
	#$name = $temp[0];
	$name = $temp[0] . $temp[1] . $temp[2]; #added by dfouts 04/13/2022 to add more words to make the com_name identifiers more unique and specific
	$name =~ s/\(.*$//; # added 04/13/2022 to remove (locus_id) that is in some accessions
	$name =~ s/\W//g; # remove all non-word characters
        print "get_gene_info: CLEANname = >$name<\n" if ($DEBUG);
	$hitref->{$d[2]}->{'clean_name'} = $name;

        $reref->{$d[0]}{$d[3]}->{'featname'} = $d[2];
        $hitref->{$d[2]}->{'hit'} = 0; # set default hit values to 0
      }
      else  {
        print "ERROR:  non-unique feat_names detected ($d[2]) - each asmbl_id must have unique names for the ORFs so they can be distinguished from the BLAST data!\n";
        close (INFOFILE);
        &write_log("4","non-unique feat_names detected - each asmbl_id must have unique names for the ORFs so they can be distinguished from the BLAST data");
      }

  }
  close (INFOFILE);
  #### check to make sure user is not specifing an asmbl_id that has no information (not in info file) ###
#  if ((defined($asmbl_id)) && (!exists $asmblref->{$asmbl_id}))  {
#      print "ERROR: There is no information for user-specified asmbl_id $asmbl_id in $infofile - terminating analysis\n";
#      &write_log("4","There is no information for user-specified asmbl_id $asmbl_id in $infofile - terminating analysis");
#  }
  return ($count);
}

sub get_tRNAs  {

    my ($tRNA_file,$tRNAref,$hitref,$reref,$DEBUG) = @_;
    my @a = ();
    my %count = ();
    my $line = "";
    my $asmbl_id = "";
    my $codon = "";
    my $featname = "";
    my $type = "";
    my $key = "";
    my $yek = "";

      open (TRNA, "<$tRNA_file") || &write_log("4","can't read file $tRNA_file: $!\n");
      while($line = <TRNA>) {
        chomp($line);
        if (($line =~ "^Sequence") || ($line =~ "^Name") || ($line =~ "^--------")) {next;}
	@a = split(/\t+/, $line);

        $a[0] =~ s/\s*$//; # remove any spaces at the end of the line
        $a[2] =~ s/\s*$//; # remove any spaces at the end of the line
	$a[3] =~ s/\s*$//; # remove any spaces at the end of the line
	$a[4] =~ s/\s*$//; # remove any spaces at the end of the line
	$a[5] =~ s/\s*$//; # remove any spaces at the end of the line
	$a[8] =~ s/\s*$//; # remove any spaces at the end of the line

	if ($a[0] =~ /.*\|(NC_\d+)/) { # get NC refseq number 
	    $asmbl_id = $1;
        }
        elsif ($a[0] =~ /[:|](\d+)/) { # or get only the number
	    $asmbl_id = $1;
        }
        else  {
          $asmbl_id = $a[0];
        }
        $tRNAref->{$asmbl_id}{$a[2]}->{'end3'} = $a[3];
        $tRNAref->{$asmbl_id}{$a[2]}->{'type'} = $a[4];
        $tRNAref->{$asmbl_id}{$a[2]}->{'anticodon'} = $a[5];
        $tRNAref->{$asmbl_id}{$a[2]}->{'cove'} = $a[8];
      }
      close (TRNA);
      foreach $key (sort {$a <=> $b} keys %{$tRNAref}) { # sort by asmbl_id
        foreach $yek (sort {$a <=> $b} keys %{$tRNAref->{$key}}) { # then sort the tRNAs by end 5
	    $type = $tRNAref->{$key}{$yek}->{'type'};
            if ($type =~ /SeC\(p\)/) { # The fasta33 does not like having the (p) in the filename (probably too long)
		$type =~ s/\(p\)//;
	    }
	    $count{$type}++;
	    $featname = "tRNA" . "-" . $type . "-" . $count{$type};
	    $reref->{$key}{$yek}->{'featname'} = $featname;
	    $codon = &revcomp($tRNAref->{$key}{$yek}->{'anticodon'});
	    $featname = $key . "_" . $featname;
            $hitref->{$featname}->{'end5'} = $yek; # added 05/18/05
            $hitref->{$featname}->{'end3'} = $tRNAref->{$key}{$yek}->{'end3'}; # added 05/18/05
	    $hitref->{$featname}->{'com_name'} = "anticodon = $tRNAref->{$key}{$yek}->{'anticodon'}, codon = $codon, cove score = $tRNAref->{$key}{$yek}->{'cove'}";
	    print "get_tRNAs:  end5 = $hitref->{$featname}->{'end5'}, end3 = $hitref->{$featname}->{'end3'}, key = >$key<, yek = >$yek<, Feat_name = >$featname<, Com_name = >$hitref->{$featname}->{'com_name'}<\n" if ($DEBUG);
        }
      }
    return;
}

sub get_tmRNAs  {
    
    my ($tmRNA_file,$tRNAref,$hitref,$reref,$DEBUG) = @_;
    my @a = ();
    my $count = "";
    my $line = "";
    my $B = "";
    my $E = "";
    my $end5 = "";
    my $end3 = "";
    my $asmbl_id = undef;
    my $featname = "";
    my $hold = "";

      open (TMRNA, "<$tmRNA_file") || &write_log("4","can't read file $tmRNA_file: $!\n");
      while($line = <TMRNA>) {
        chomp($line);
	@a = split(/\s+/, $line); # split on whitespace
        if (($a[1] =~ /nucleotides/) && (!defined($asmbl_id))) {
            @a = split(/\s+/, $hold);
	    if ($a[0] =~ /.*\|(NC_\d+)/) { # get NC refseq number 
		$asmbl_id = $1;
	    }
	    elsif ($a[0] =~ /[:|](\d+)/)  { # or get only the number
		$asmbl_id = $1;
	    }
	    else  {
		$asmbl_id = $a[0];
	    }
        print ">>$asmbl_id<<\n" if ($DEBUG);   
	}
	if (($a[0] =~ /Location/) && ($a[1] =~ /(\d+),(\d+)/)) {
	    $count++;
            $B = $1;
            $E = $2;
            if ($a[1] =~ /^c\[/) { # tmRNA is on other strand <----
		$end5 = $E;
		$end3 = $B;
            }
	    else {
		$end5 = $B;
                $end3 = $E;
	    }
	    $tRNAref->{$asmbl_id}{$end5}->{'end3'} = $end3;
	    $tRNAref->{$asmbl_id}{$end5}->{'type'} = "NA";
	    $tRNAref->{$asmbl_id}{$end5}->{'cove'} = "NA";
            $featname = "tmRNA" . "-" . $count;
	    $reref->{$asmbl_id}{$end5}->{'featname'} = $featname;
            $featname = $asmbl_id . "_" . $featname;
            $hitref->{$featname}->{'end5'} = $end5;
            $hitref->{$featname}->{'end3'} = $end3;
            print "get_tmRNAs:  Asmbl_id = >$asmbl_id<, Feat_name = >$featname<, End5 = >$end5<, End3 = >$end3<, " if ($DEBUG);
        }
	elsif (($a[1] =~ /peptide:/) && ($a[2] =~ /(\w+)\*/)) {
	    $tRNAref->{$asmbl_id}{$end5}->{'anticodon'} = $1;
            $hitref->{$featname}->{'com_name'} = "Tag peptide = $tRNAref->{$asmbl_id}{$end5}->{'anticodon'}";
            print "Tag = >$hitref->{$featname}->{'com_name'}<\n" if ($DEBUG);
            $asmbl_id = undef; # clear this if more than one asmbl_id
	}
        $hold = $line;
    }
    close (TMRNA);
    return;
}

sub get_ncRNAs  { # obtain data from the output of cmscan (ncRNA_cmscan.out)
    
    my ($ncRNA_file,$tRNAref,$hitref,$reref,$DEBUG) = @_;
    my @a = ();
    my %count = ();
    my $line = "";
    my $end5 = "";
    my $end3 = "";
    my $asmbl_id = undef;
    my $type = "";
    my $featname = "";
    my $comname = "";

    #data key:
    #$a[1] = target name
    #$a[2] = accession
    #$a[3] = query name (asmbl_id)
    #$a[9] = seq from
    #$a[10] = seq to
    #$a[16] = score
    #$a[26..$#a] = description of target

    open (NCRNA, "<$ncRNA_file") || &write_log("4","can't read file $ncRNA_file: $!\n");
    while($line = <NCRNA>) {
	chomp($line);
	if ($line =~ /^#/) {next;} # skip all lines beginningn with "#"
	@a = split(/\s+/, $line); # split on whitespace
	if ($a[3] =~ /.*\|(NC_\d+)/) { # get NC refseq number 
	    $asmbl_id = $1;
	}
	elsif ($a[3] =~ /[:|](\d+)/)  { # or get only the number
	    $asmbl_id = $1;
	}
	else  {
	    $asmbl_id = $a[3];
	}
	print "get_ncRNAs:  >>$asmbl_id<<\n" if ($DEBUG);   

	#note that Sean Eddy's group already oriented the genes 5prime to 3prime and tell us the direction in $a[11] so we don't have to figure the orientation
	$end5 = $a[9];
	$end3 = $a[10];

	if ($a[1] =~ /6S/) {
	    $type = "6S_SsrS"
	}
	elsif ($a[1] =~ /SRP/) {
	    $type = "ffs_4.5S"
	}
	elsif ($a[1] =~ /RNaseP/) {
	    $type = "RNaseP"
	}
	else {
	    $type = $a[1];
	}

	$comname = join(" ",@a[26..$#a]);
	$tRNAref->{$asmbl_id}{$end5}->{'end3'} = $end3;
	$tRNAref->{$asmbl_id}{$end5}->{'type'} = $type;
	$tRNAref->{$asmbl_id}{$end5}->{'cove'} = $a[16];
	$tRNAref->{$asmbl_id}{$end5}->{'anticodon'} = "NA";

	$count{$type}++;
	$featname = "ncRNA" . "-" . $tRNAref->{$asmbl_id}{$end5}->{'type'} . "-" . $count{$type};
	$reref->{$asmbl_id}{$end5}->{'featname'} = $featname;
	$featname = $asmbl_id . "_" . $featname;
	$hitref->{$featname}->{'end5'} = $end5;
	$hitref->{$featname}->{'end3'} = $end3;
	$hitref->{$featname}->{'com_name'} = $comname;
	print "get_ncRNAs:  Asmbl_id = >$asmbl_id<, Feat_name = >$featname<, End5 = >$hitref->{$featname}->{'end5'}<, End3 = >$hitref->{$featname}->{'end3'}<, com_name = >$hitref->{$featname}->{'com_name'}<\n" if ($DEBUG);
	$asmbl_id = undef; # clear this if more than one asmbl_id
	$end5 = undef; # clear out end5 and end3
	$end3 = undef;
    }
    close (NCRNA);
    return;
}

sub select_featnames_from_btab { # populate a hash (index) with featnames (key) and 1 as value
    my ($btabfile,$evalue,$hitref,$searchref,$excluderef,$lists_ref,$DEBUG) = @_;
    my @a = ();
    my @b = ();
    my @c = ();
    my @accession = ();
    my $featname = "";
    my $WUBLAST = "";
    my $match = "";
    my $altevalue = 0.0001;
    my $descriptor = "";
    my $linecounter = "";
    
    open (INFILE, "<$btabfile") || &write_log("4","can't open file $btabfile: $!\n");
    chomp;
    @a = split(/\t+/, <INFILE>);
    if ($a[3] =~ /washu/) {
        print "Detected WU-BLAST btab file ...\n";
        $WUBLAST = 1;
    }
    else  {
        print "Assuming NCBI BLAST (-m 8) option btab file ...\n" if ($DEBUG);
        $WUBLAST = 0;
    }
    close(INFILE);
    open (INFILE, "<$btabfile") || &write_log("4","can't open file $btabfile: $!\n");
    while (<INFILE>)  {
      chomp;
      $linecounter++;
      @a = split(/\t+/);
      if ($a[0] =~ /^gi\|?(\d+)/)  {  # if gi numbers are present, we need to only look at the numbers, not the gi| or the :end5..end3 garbage
        $a[0] = $1;  # reset $id to be the actual gi number and nothing else
      }
      elsif ($a[0] =~ /^(\S+_gi\d+)-/) {
	  $a[0] = $1;
      }
      if (exists $hitref->{$a[0]}) {  # add featname ONLY if already present in db
        $featname = $a[0];
      }
      else  {
	next;  # skip if feat_name not in the db
      }
      if ($WUBLAST == 1)  { # WU-BLAST btab file

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
        if (exists $excluderef->{$a[5]})  { next; } # skip if in exclude list
        
	foreach $descriptor (keys %{$lists_ref})  {
	    if (($a[19] <= $altevalue) && (exists $lists_ref->{$descriptor}{$a[5]}) && (!$hitref->{$featname}->{"$descriptor"}))  { # check for match to list of key accessions (small/large terminase, portal, capsid
		print "FEAT_NAME <$featname> - HIT <$descriptor> !!!!!!!!!!!!!!!\n" if ($DEBUG);
		$hitref->{$featname}->{"$descriptor"} = 1;
		$hitref->{$featname}->{'com_name'} = $hitref->{$featname}->{'com_name'} . " " . "{$descriptor BLAST $a[5]}" if ($descriptor !~ /Core/);
	    }
        }
        if (($a[10] != -100.000000) && ($a[19] <= $evalue)){ ### add featname ONLY if already present in db and a good hit####
	   if ($hitref->{$featname}->{'hit'} != 1) {  # NOTE: only best hit annotation and evalue are recorded!
            $hitref->{$featname}->{'hit'} = 1;
            $a[15] =~ s/\s\s+/ /g;  # remove stretch of multiple spaces between comname and [organism]
            $a[15] =~ s/  / /g;  # remove any double spaces
            $a[15] =~ s/\s*$//; # remove any spaces at the end of the line
            @c = split(/-/, $a[5]);
            # replaced old TAG-giaccession by giaccession-TAG scheme in db
            $hitref->{$featname}->{'phage'} = $c[1];  # store the phage-tag of the BLAST hit
            $hitref->{$featname}->{'annotation'} = $c[0] . "," . $a[15]; # store feat_name and annoation of the hit
            $hitref->{$featname}->{'evalue'} = $a[19];  # store the e-value of the hit
            $searchref->{$hitref->{$featname}->{'asmbl_id'}}{$hitref->{$featname}->{'end5'}}->{'featname'} = $featname; # hash "searchhash" has asmbl_id and end5 as keys and featname as value ONLY when there is a hit!
          }
        }
      }
      else  { # must be NCBI BLAST file

# ========================================================
# btab output from NCBI blastall (-m 8) or BLAST+ blastn (-outfmt 6) option:

# column number Description (for Perl), add 1 for Unix
# 0	Query_id
# 1	subject_id (Hit from db)
# 2	% Identity
# 3	length of alignment
# 4	number or mismatches
# 5	number of gaps
# 6	start of alignment on query (5' nucleotide match in query)
# 7	end of alignment on query (3' nucleotide match in query)
# 8	start of alignment on subject (5' nucleotide match in db hit)
# 9	end of alignment on subject (3' nucleotide match in db hit)
# 10	e-value
# 11	score (bits)
# ========================================================
	 
	if ($excluderef->{$a[1]})  { # skip if in exclude list
	    print "$linecounter) $a[0]: $a[1] is excluded\n" if ($DEBUG); 
	    next; 
	}
        # new to prevent very weak housekeeping genes from slipping into a phage region
        if (($a[2] < 35) && ($a[11] < 75)) { # if the percent identity is less than 35% and the BITS score is < 75 then skip
            print "$linecounter) $a[0]: $a[1], pid: $a[2], BITS: $a[11] is excluded due to match below pid and BITS\n" if ($DEBUG);
            next;
        }
        print "SELECT_BTAB:  $linecounter) |$a[0]|\t|$a[1]|\n" if ($DEBUG);
        if ($a[10] <= $evalue){ ### add featname ONLY if already present in db and a good hit####
          if ($hitref->{$featname}->{'hit'} != 1) {  # NOTE: only best hit annotation and evalue are recorded!
            $hitref->{$featname}->{'hit'} = 1;
            @c = split(/-/, $a[1]);
            # replaced TAG-giaccession by giaccession-TAG scheme in the db
            $hitref->{$featname}->{'phage'} = $c[1];  # store the phage-tag of the BLAST hit
            $hitref->{$featname}->{'annotation'} = $c[0]; # since NCBI BLAST does not give annotation, use ORF identifier (feat_name)
            $hitref->{$featname}->{'evalue'} = $a[10];  # store the e-value of the hit
            $searchref->{$hitref->{$featname}->{'asmbl_id'}}{$hitref->{$featname}->{'end5'}}->{'featname'} = $featname; # hash "searchhash" has asmbl_id and end5 as keys and featname as value ONLY when there is a hit!
            # Now check each descriptor for a match
            foreach $descriptor (keys %{$lists_ref})  {
		if ((exists $lists_ref->{$descriptor}{$a[1]}) && (!$hitref->{$featname}->{"$descriptor"}))  { # check for match to list of key accessions (small/large terminase, portal, capsid, etc)
		    print "FEAT_NAME <$featname> - HIT <$descriptor> !!!!!!!!!!!!!!!\n" if ($DEBUG);
		    $hitref->{$featname}->{"$descriptor"} = 1;
		    $hitref->{$featname}->{'com_name'} = $hitref->{$featname}->{'com_name'} . " " . "{$descriptor BLAST $a[1], $a[10]}" if ($descriptor !~ /Core/);
	      }
            }
	  }
        }  
      }
    }
    close(INFILE);
}

sub find_hmms {  # find HMM hits to bacteriophage proteins
    # only record the GLOCAL model results if present, else, record the fragment model results (PFAM only)
    # updated to input the new HMMer3 output file (no GLOCAL alignments with HMMer3)
    my ($phome,$HMMversion,$data,$infofile,$hitref,$searchref,$HMMref,$serineref,$tyrosineref,$DEBUG) = @_;
    my %HMM_master_hash = (); # stores the HMM accession (key) and length, trusted and noise cut-offs
    my ($time,$date,$model);
    my @h = ();
    my $file = "";
    my ($query_HMM,$description,$feat_name,$total_score,$noise,$trusted);
    my $scores = 0;
    my $descriptor = "";
    &get_HMM_size_trusted_noise($phome,$HMMversion,\%HMM_master_hash); # get the size, trusted and noise cut-off values

    for $file (@{$data}) { # loop through list of filenames, either user-specified (GLOCAL) or default (both GLOCAL and FRAG)
	if (-s "$file") { # changed 03/07/14 to simplify code dfouts
	    if ($file =~ /FRAG/) {
		$model = "FRAG";
	    }
	    else {
		$model = "GLOCAL";
	    } 
	    &write_log("1","Retrieving HMM hits information from $file");
	    open (HMMFILE, "<$file") || &write_log("4","can't open file HMM datafile $file: $!\n");
	    while (<HMMFILE>)  {
		if ((/^Query HMM:\s+(\S+)/)
                   || (/^Query:\s+(\S+)/)) {
		    $query_HMM = $1;
		}
		elsif (/^Description:\s+(.+)/) {
		    $description = $1;
		    $description =~ s/\t/ /g; # converst tabs into spaces
		    $description =~ s/^\s+//; # remove beginning spaces
		    $description =~ s/\s+$//; # remove trailing spaces
		}
		elsif ((/^(Sequence\s+)Description\s+Score\s+E-value\s+(N|\#D)/) 
                      || (/^\s+E-value\s+score\s+bias\s+E-value\s+score\s+bias\s+exp\s+N\s+Sequence\s+Description/)){
		    $scores = 1;
		    next;
		}
		elsif (($scores) && (/^--------\s+-----------\s+-----\s+-------\s+-{2,3}/)) {
                    next;
		}
		elsif (($scores) && (/^\s+-------\s------\s-----\s+-------\s------\s-----\s+----\s--\s+--------\s+----------/)) {
		    next;
		}
		elsif (($scores) && (/------ inclusion threshold ------/))  {
		    next;
		}
		elsif (($scores) && (/^\s*$/)) { # if found scores and blank line, reset value
		    $scores = 0;
		}
		elsif ($scores){
                    chomp;
		    @h = split(/\s+/);
		    if ($HMMversion == "2") {
			$feat_name = shift @h;
		    }
		    elsif ($HMMversion == "3") {
			$feat_name = $h[9];
		    }
		    if ($feat_name =~ /^gi\|?(\d+)/)  {  # if gi numbers are present, we need to only look at the numbers, not the gi| or the :end5..end3 garbage
			$feat_name = $1;  # reset $id to be the actual gi number and nothing else
		    }
		    elsif ($feat_name =~ /^(\S+_gi\d+)-/) {
			$feat_name = $1;
		    }
                    print "Model: $model, feat_name: $feat_name, query_HMM:  <$query_HMM>, description:  $description\n" if ($DEBUG);
		    if (!exists $hitref->{$feat_name}) {
			&write_log("1","ERROR:  >$feat_name< There seems to be a descrepency between ORF names in $file and $infofile - they should have the same names!");
			next;
		    }		    
		    if ( !( (exists $hitref->{$feat_name}) && (exists $hitref->{$feat_name}->{'hmm'}) && (exists $hitref->{$feat_name}->{'hmm'}{$query_HMM}) && (exists $hitref->{$feat_name}->{'hmm'}{$query_HMM}->{'score'}) ) ) {
			if ($HMMversion == "2")  {
			    pop @h; # remove #D
			    pop @h; # remove E-value
			    $total_score = pop @h;
			}
			elsif ($HMMversion == "3")  {
			    $total_score = $h[2];
			}
			$noise = $HMM_master_hash{$model}{$query_HMM}->{'noise'};
			$trusted = $HMM_master_hash{$model}{$query_HMM}->{'trusted'};
			if ((($total_score > $noise) && ($model eq "GLOCAL")) || (($total_score > $trusted) && ($model eq "FRAG")))  {  # if score is not already recorded and is above the noise cut-off for GLCAL model or trusted for FRAG model
			    $hitref->{$feat_name}->{'hmm'}{$query_HMM}->{'model'} = $model;
			    $hitref->{$feat_name}->{'hmm'}{$query_HMM}->{'score'} = $total_score;
			    $hitref->{$feat_name}->{'hmm'}{$query_HMM}->{'trusted'} = $trusted;
			    $hitref->{$feat_name}->{'hmm'}{$query_HMM}->{'noise'} = $noise;
			    $hitref->{$feat_name}->{'hmm'}{$query_HMM}->{'hmm_com_name'} = $description;
                            print "find_HMMs: $feat_name = $query_HMM, TRUSTED = $trusted > score = $total_score > noise = $noise\n" if ($DEBUG); 
			    foreach $descriptor (keys %{$HMMref})  {
				if ($HMMref->{$descriptor}{$query_HMM} == 1) {
				    $hitref->{$feat_name}->{'com_name'} = $hitref->{$feat_name}->{'com_name'} . " " . "{$descriptor $query_HMM}" if ($descriptor !~ /Core/);
				    $descriptor =~ s/\sHMM//; # remove ( HMM) from descriptor so that if there is only an HMM match, we still get these put into descritor-bases pep files
				    $hitref->{$feat_name}->{"$descriptor"} = 1;
				}
			    }
			    if ($serineref->{$query_HMM} == 1) {
				$hitref->{$feat_name}->{'serine'} = 1;
				$hitref->{$feat_name}->{'com_name'} = $hitref->{$feat_name}->{'com_name'} . " " . "{Serine Recombinase HMM $query_HMM}";
			    }
			    elsif ($tyrosineref->{$query_HMM} == 1) {
				$hitref->{$feat_name}->{'tyrosine'} = 1;
				$hitref->{$feat_name}->{'com_name'} = $hitref->{$feat_name}->{'com_name'} . " " . "{Tyrosine Recombinase HMM $query_HMM}";
			    }
			    
			    if (($hitref->{$feat_name}->{'hit'} == 0) && (($query_HMM ne "TIGR02224") || ($query_HMM ne "TIGR02225") || ($query_HMM ne "TIGR02249"))) {  # added 02/12/04 to add hmm hit to searchhash only if no blast hit
				$searchref->{$hitref->{$feat_name}->{'asmbl_id'}}{$hitref->{$feat_name}->{'end5'}}->{'featname'} = $feat_name;
				#$hitref->{$feat_name}->{'hit'} = 1; # removed 05/20/05
			    }  
			}
		    }   
		}
	    }
	    close (HMMFILE);
	}
	else { return("0"); } # no hmm files present
    }
    return("1");
}

sub get_fasta_record {

    my ($fh) = @_;
    my $record = '';
    my ($save_input_separator) = $/;
    $/ = "^>"; # fixed bug in generating phage_finder_info.txt from .ptt file, was "\n". dfouts 04/09/2012
    $record = <$fh>;
    $/ = $save_input_separator;
    return $record;
}


sub get_assemblies {

  my ($file,$infofile,$ahref,$asmblref,$asmbl_id,$DEBUG) = @_;
  my $title = "";
  my $id_key = "";
  my @id = ();
  my @cropped_id = ();
  my @seq = ();
  my @temp = ();
  my $sequence = "";
  my $found_asmbly = 0;
  my $record;
  my $status;

  #Open and read file with fasta format sequences.
  open (GENOME, "<$file") || &write_log("4","can't find the file $file: $!\n");
  my ($save_input_separator) = $/;
  $/="\n>";
  while (<GENOME>) {
      ($title,$sequence) = /^>?\s*(.*)\n([^>]+)>?/; # split the header line and sequence (very cool)
      $title =~ s/\t/ /g;  # remove tabs
      $sequence =~ s/[^a-zA-Z]//g; # remove any non-alphabet letters
      if ($title =~ /^([^:]+:)?(\b$asmbl_id\b)\s*(.*)$/)  { # 04/17/13 code by Phil Goetz to fix asmbl_id bug where some asmbl_id ids were missed
	  # fixed a bug where asmbl_id "NODE_53" was being called "NODE_5", so I added "\b" before and after asmbl_id regex above 10/03/2022
	  print "get_assemblies: title:  [ $asmbl_id ] ( $title )\n" if ($DEBUG);
	  $id_key = $2; # changed from $1 to $2 04/17/13 to fix asmbl_id bug
	  print "get_assemblies: id_key:  <$id_key>\n" if ($DEBUG);
	  $found_asmbly = 1;
	  $title = $3; # changed from direct removal of asmbl_id to a more elegant method 04/17/13
	  #$title =~ s/^$asmbl_id\s//; # remove asmbl_id from recorded title
	  $ahref->{$id_key}->{'title'} = $title;
	  @temp = split(/\s+/, $title);
	  $ahref->{$id_key}->{'genus'} = $temp[0];
	  $ahref->{$id_key}->{'species'} = $temp[1];
	  $ahref->{$id_key}->{'strain'} = join("_", @temp[2..$#temp]); #use the remainder of the line as the initial strain designation.
	  $ahref->{$id_key}->{'species'} =~ s/,//g; # remove the comma
	  $ahref->{$id_key}->{'strain'} =~ s/,//g; # remove the comma
	  $ahref->{$id_key}->{'strain'} =~ s/\'//g; # remove the single quotes
	  $ahref->{$id_key}->{'strain'} =~ s/-//g; # remove any dashes
	  $ahref->{$id_key}->{'strain'} =~ s/\///g; # remove any forward slashes
	  $ahref->{$id_key}->{'strain'} =~ s/_complete//g; # remove complete
	  $ahref->{$id_key}->{'strain'} =~ s/_genome$//g; # remove genome from end
	  $ahref->{$id_key}->{'strain'} =~ s/_complete//g; # remove complete
	  $ahref->{$id_key}->{'strain'} =~ s/_sequence$//g; # remove sequence from end
	  $ahref->{$id_key}->{'strain'} =~ s/_chromosomal//g; # remove chromosomal
	  $ahref->{$id_key}->{'strain'} =~ s/_pseudomolecule$//g; # remove pseudomolecule from end
	  $ahref->{$id_key}->{'strain'} =~ s/_chromosome//g; # remove chromosome
	  $ahref->{$id_key}->{'strain'} =~ s/_circular//g; # remove circular
	  $ahref->{$id_key}->{'strain'} =~ s/str\._//g; # remove str. abbreviation
	  $ahref->{$id_key}->{'strain'} =~ s/\.//g; # remove period 
	  $ahref->{$id_key}->{'sequence'} = $sequence;
	  $ahref->{$id_key}->{'gc'} = &get_gc_content($sequence);
          #print "gc content is ($ahref->{$id_key}->{'gc'})\n";
	  $ahref->{$id_key}->{'length'}= length($sequence);
	  $status = 1;
	  last;
      }
  }
  $/ = $save_input_separator;
  close (GENOME);
  if ($found_asmbly == 0)  { # no specified asmblys found, exit with error
    print "WARNING: The file $file does not contain the correct sequence for asmbl_id $asmbl_id, skipping att site analysis\n";
    &write_log("1","WARNING: The file $file does not contain the correct sequence for asmbl_id $asmbl_id, skipping att site analysis");
    $ahref->{$id_key}->{'gc'} = "NA";
    $status = 0;
  }
  elsif ($asmblref->{$asmbl_id}->{'genomesize'} != $ahref->{$id_key}->{'length'})  {
    print "WARNING: The asmbl_id <$id_key> ($ahref->{$id_key}->{'length'} bp) is not the same size as <$asmbl_id> reported in $infofile ($asmblref->{$asmbl_id}->{'genomesize'} bp), skipping att site analysis\n";
    &write_log("1","WARNING: The asmbl_id $asmbl_id ($ahref->{$id_key}->{'length'} bp) is not the same size as reported in $infofile ($asmblref->{$asmbl_id}->{'genomesize'} bp), skipping att site analysis");
    $status = 0;
  }
  print "get_assemblies: DONE READING SEQUENCE\n" if ($DEBUG);   
  return($status);
}


sub populate_atthash {

    my ($attref, $ahref, $asmbl_id, $n, $direction, $DEBUG) = @_;
    my $size_h = $attref->{$n}->{'attL_coordE'}-$attref->{$n}->{'attL_coordB'}+1;
    my $size_t = $attref->{$n}->{'attR_coordE'}-$attref->{$n}->{'attR_coordB'}+1;

################## Crib sheet ###################
# fresh from att search not case of force flip  #
# note: integrase is usually on head side       #
#    head                               tail    #
#    attL            phage (+)          attR    #
#    B  E............------->...........B  E    #
#    |  |                               |  |    #
#   L5  L3                             R5  R3   #
# attB                               attE       #
#    left                               right   #
#   pleft                              pright   #
#                                               #
#                                               #
#                                               #
#    tail                               head    #
#    attR            phage (-)          attL    #
#    B  E............<-------...........B  E    #
#    |  |                               |  |    #
#   R3  R5                             L3  L5   #
#       attE                               attB #
#    left                               right   #
#   pright                              pleft   #
#                                               #
#################################################
 
    if (exists $attref->{$n}->{'ori'})  { # flipping a phage region (direction opposite to integrase)
	my $attL5 = $attref->{$n}->{'attL5'};
	my $attL3 = $attref->{$n}->{'attL3'};
	my $attR5 = $attref->{$n}->{'attR5'};
	my $attR3 = $attref->{$n}->{'attR3'};
	my $attR = $attref->{$n}->{'attR'};
	my $attL = $attref->{$n}->{'attL'};
	my $pleft = $attref->{$n}->{'pleft'};
	my $pright = $attref->{$n}->{'pright'};
	$attref->{$n}->{'attL5'} = $attR3;  # compute genome coord for attL 5' end relative to phage
	$attref->{$n}->{'attL3'} = $attR5;  # compute genome coord for attL 3' end relative to phage
	$attref->{$n}->{'attR5'} = $attL3;  # compute genome coord for attR 5' end relative to phage
	$attref->{$n}->{'attR3'} = $attL5;  # compute genome coord for attR 3' end relative to phage
	$attref->{$n}->{'attL'} = &revcomp($attR); # get reverse complement
	$attref->{$n}->{'attR'} = &revcomp($attL);
	$attref->{$n}->{'pleft'} = $pright; # position of leftmost phage coordinate relative to phage direction
	$attref->{$n}->{'pright'} = $pleft; # position of rightmost phage coordinate relative to phage direction
	if ($direction eq "+") { # set new 5' end of each half att site and direction
	    $attref->{$n}->{'attB'} = $attref->{$n}->{'attL5'};;
	    $attref->{$n}->{'attE'} = $attref->{$n}->{'attR5'};;
	    $attref->{$n}->{'ori'} = "+";
	}
	else  {
	    $attref->{$n}->{'attB'} = $attref->{$n}->{'attL5'};;
	    $attref->{$n}->{'attE'} = $attref->{$n}->{'attR5'};;
	    $attref->{$n}->{'ori'} = "-";
	}
    }
    else {
	if ($direction eq "+")  {
	    $attref->{$n}->{'ori'} = "+";
	    $attref->{$n}->{'attL5'} = $attref->{$n}->{'attL_coordB'};  # compute genome coord for attL 5' end relative to phage
	    $attref->{$n}->{'attL3'} = $attref->{$n}->{'attL_coordE'};  # compute genome coord for attL 3' end relative to phage
	    $attref->{$n}->{'attR5'} = $attref->{$n}->{'attR_coordB'};  # compute genome coord for attR 5' end relative to phage
	    $attref->{$n}->{'attR3'} = $attref->{$n}->{'attR_coordE'};  # compute genome coord for attR 3' end relative to phage
	    $attref->{$n}->{'attB'} = $attref->{$n}->{'attL5'}; # changed from attR5 05/20/05
	    $attref->{$n}->{'attE'} = $attref->{$n}->{'attR5'}; # changed from attL5 05/20/05
	    $attref->{$n}->{'attL'} = substr($ahref->{$asmbl_id}->{'sequence'},$attref->{$n}->{'attL5'}-1,$size_h);
	    $attref->{$n}->{'attR'} = substr($ahref->{$asmbl_id}->{'sequence'},$attref->{$n}->{'attR5'}-1,$size_t);
	    $attref->{$n}->{'left'} = $attref->{$n}->{'attL5'}; # position of leftmost phage coordinate relative to genome
	    $attref->{$n}->{'right'} = $attref->{$n}->{'attR3'}; # position of rightmost phage coordinate relative to genome
	    $attref->{$n}->{'pleft'} = $attref->{$n}->{'attL5'}; # position of leftmost phage coordinate relative to phage direction
	    $attref->{$n}->{'pright'} = $attref->{$n}->{'attR3'}; # position of rightmost phage coordinate relative to phage direction
	}
	else  {
	    $attref->{$n}->{'ori'} = "-";
	    $attref->{$n}->{'attL5'} = $attref->{$n}->{'attL_coordE'};  # compute genome coord for attL 5' end relative to phage
	    $attref->{$n}->{'attL3'} = $attref->{$n}->{'attL_coordB'};  # compute genome coord for attL 3' end relative to phage
	    $attref->{$n}->{'attR5'} = $attref->{$n}->{'attR_coordE'};  # compute genome coord for attR 5' end relative to phage
	    $attref->{$n}->{'attR3'} = $attref->{$n}->{'attR_coordB'};  # compute genome coord for attR 3' end relative to phage
	    $attref->{$n}->{'attB'} = $attref->{$n}->{'attL5'}; # changed from attR3 05/20/05
	    $attref->{$n}->{'attE'} = $attref->{$n}->{'attR5'}; # changed from attL3 05/20/05
	    $attref->{$n}->{'attL'} = &revcomp(substr($ahref->{$asmbl_id}->{'sequence'},$attref->{$n}->{'attL3'}-1,$size_h)); # get reverse complement
	    $attref->{$n}->{'attR'} = &revcomp(substr($ahref->{$asmbl_id}->{'sequence'},$attref->{$n}->{'attR3'}-1,$size_t));
	    $attref->{$n}->{'left'} = $attref->{$n}->{'attR3'}; # position of leftmost phage coordinate relative to genome
	    $attref->{$n}->{'right'} = $attref->{$n}->{'attL5'}; # position of rightmost phage coordinate relative to genome
	    $attref->{$n}->{'pleft'} = $attref->{$n}->{'attL5'}; # position of leftmost phage coordinate relative to phage direction
	    $attref->{$n}->{'pright'} = $attref->{$n}->{'attR3'}; # position of rightmost phage coordinate relative to phage direction
	    print ">>->>->>POPULATE_ATTHASH: <rc>\n" if ($DEBUG);
	}
	$attref->{$n}->{'size'} = abs($attref->{$n}->{'right'} - $attref->{$n}->{'left'}); # compute size of region
    }
}

sub pick_longest_length {

    my ($length1, $length2) = @_;
    my $attL = "";
    my $attR = "";
    my $best = "";

    if ($length1 > $length2) {
	$best = 1;
    }
    elsif ($length2 > $length1)  {
	$best = 2;
    }
    else { # both are 0 length or same length, probably hit a phage db or HMM match
	$best = 0;
    }
    return($best);
}
sub pick_best_att {

    my ($name, $direction, $attref, $hitref, $phageref, $reref, $asmbl_id, $aref, $ahref, $DEBUG) = @_;
    my $n = 1;
    my $i = "";
    my $start = "";
    my $finish = "";
    my $feat_name = "";
    my $TargetLength = "";
    my $oneThirdSize = "";
    my $oneThirdEnd = "";
    my $IDtarget = 0; # false initially
    my $inAgene = 0; # false initially
    #my $direction = ""; # why did I redeclare this variable?
    my $hold_5prime = "";
    my $hold_3prime = "";
    my $ORF5_featname = &adjust_featname($phageref->{$name}->{'ORF5'},$asmbl_id);
    my $ORF3_featname = &adjust_featname($phageref->{$name}->{'ORF3'},$asmbl_id);
    my $return_val = "";
          
  ATT:    until (($n > 2) || (!defined ($attref->{$n}->{'left'})))  { # check to see if TSD is within a gene and at 3 prime end of the gene
      $start = $hitref->{$ORF5_featname}->{'array_pos'} - 5; # start 5 ORFs upstream of initial 5' end (not att coord)
      $finish = $hitref->{$ORF3_featname}->{'array_pos'} + 6; # finish 5?ORFs downstream of initial 3' end (not att coord)
      print ">>->>PICK_BEST_ATT: Checking att-site # $n:  start = $start\tfinish = $finish\n" if ($DEBUG);
      print ">>->>PICK_BEST_ATT: attL: $attref->{$n}->{'attL'}\n" if ($DEBUG);
      print ">>->>PICK_BEST_ATT: attR: $attref->{$n}->{'attR'}\n" if ($DEBUG);
      for ( $i=$start; $i<$finish; $i++ )  {
          $inAgene = 0; # reset value for each round
	  $IDtarget = 0; # reset value for each round
	  $feat_name = &adjust_featname($reref->{$asmbl_id}{$aref->[$i]}->{'featname'},$asmbl_id);
	  $TargetLength = abs($hitref->{$feat_name}->{'end5'} - $hitref->{$feat_name}->{'end3'});
	  $oneThirdSize = $TargetLength*0.67;
	  print ">>->>->>PICK_BEST_ATT: $feat_name-$hitref->{$feat_name}->{'end5'}::$attref->{$n}->{'left'}::$attref->{$n}->{'right'}-$hitref->{$feat_name}->{'end3'}\n" if ($DEBUG);
	  if (($hitref->{$feat_name}->{'end5'} <= $attref->{$n}->{'left'}) && # ------>[ phage ]  (changed from $phageref->{$name}->{'5prime_att'} = $attL5)
	      ($hitref->{$feat_name}->{'end3'} >= $attref->{$n}->{'left'}))
	  {
	      $inAgene = 1; # within a gene is true
              print ">>->>->>->>PICK_BEST_ATT: ----->[ phage ] inAgene = 1\n" if ($DEBUG);
	      $oneThirdEnd = $hitref->{$feat_name}->{'end5'}+$oneThirdSize;
	      if ($attref->{$n}->{'left'} >= $oneThirdEnd) { # if within last 1/3 of gene
		  $IDtarget = 1; # found a target = true, 1
                  $direction = "+";
		  print ">>->>->>->>PICK_BEST_ATT: ----->[ phage ] IDtarget = 1, oneThirdEnd = $oneThirdEnd, direction = $direction\n" if ($DEBUG);
	      }
	      elsif ($feat_name =~ /RNA/)  { $direction = "+"; }
	  }
	  elsif (($hitref->{$feat_name}->{'end5'} >= $attref->{$n}->{'right'}) &&  # [ phage ]<------
		 ($hitref->{$feat_name}->{'end3'} <= $attref->{$n}->{'right'}))
	  {   
	      $inAgene = 1; # within a gene is true
              print ">>->>->>->>PICK_BEST_ATT: [ phage ]<----- inAgene = 1\n" if ($DEBUG);
	      $oneThirdEnd = $hitref->{$feat_name}->{'end5'}-$oneThirdSize;         # 3<--|-----5
	      if ($attref->{$n}->{'right'} <= $oneThirdEnd) { # if within last 1/3 of gene
	          $IDtarget = 1; # found a target = true, 1
                  $direction = "-";
	          print ">>->>->>->>PICK_BEST_ATT: [ phage ]<----- IDtarget = 1, oneThirdEnd = $oneThirdEnd, direction = $direction\n" if ($DEBUG);
	      }
	      elsif ($feat_name =~ /RNA/)  { $direction = "-"; }
	  }
###### new code looking at insertion into 5' end of ORF ######

	  elsif (($hitref->{$feat_name}->{'end3'} <= $attref->{$n}->{'left'}) && # <------[ phage ]
	      ($hitref->{$feat_name}->{'end5'} >= $attref->{$n}->{'left'}))
	  {
	      $inAgene = 1; # within a gene is true
              print ">>->>->>->>PICK_BEST_ATT: <-----[ phage ] inAgene = 1\n" if ($DEBUG);
              $oneThirdSize = $TargetLength*0.90;
	      $oneThirdEnd = $hitref->{$feat_name}->{'end3'}+$oneThirdSize;
	      if ($attref->{$n}->{'left'} >= $oneThirdEnd) { # if within last 1/3 of gene
		  $IDtarget = 1; # found a target = true, 1
                  $direction = "-";
		  print ">>->>->>->>PICK_BEST_ATT: <-----[ phage ] IDtarget = 1, oneThirdEnd = $oneThirdEnd, direction = $direction\n" if ($DEBUG);
	      }
	      elsif ($feat_name =~ /RNA/)  { $direction = "-"; }
	  }
	  elsif (($hitref->{$feat_name}->{'end3'} >= $attref->{$n}->{'right'}) &&  # [ phage ]------>
		 ($hitref->{$feat_name}->{'end5'} <= $attref->{$n}->{'right'}))
	  {   
	      $inAgene = 1; # within a gene is true
              print ">>->>->>->>PICK_BEST_ATT: [ phage ]-----> inAgene = 1\n" if ($DEBUG);
              $oneThirdSize = $TargetLength*0.90;
	      $oneThirdEnd = $hitref->{$feat_name}->{'end3'}-$oneThirdSize;        
	      if ($attref->{$n}->{'right'} <= $oneThirdEnd) { # if within last 1/3 of gene
	          $IDtarget = 1; # found a target = true, 1
                  $direction = "+";
	          print ">>->>->>->>PICK_BEST_ATT: [ phage ]-----> IDtarget = 1, oneThirdEnd = $oneThirdEnd, direction = $direction\n" if ($DEBUG);
	      }
	      elsif ($feat_name =~ /RNA/)  { $direction = "+"; }
	  }


	  # if phage hit or HMM hit, skip and check next att site - do NOT want this
	  if (($inAgene) && ($hitref->{$feat_name}->{'hit'} == 1)) {
	      $attref->{$n}->{'length'} = 0; # remove the length if bogus
	      print ">>->>->>->>PICK_BEST_ATT: Remove bogus hit:  $feat_name\n" if ($DEBUG);
	      $n++;
              $return_val = 0; # new so that atts with bogus hits don't slip through
	      next ATT;
	  }
	  elsif (($IDtarget) || (($inAgene) && ($feat_name =~ /RNA/))) # tRNA, tmRNA or ncRNA does not have to be in 3' end (JBact V184, p859, 2002)
	  {
              if (($feat_name =~ /_tRNA/) || ($feat_name =~ /_tmRNA/) || ($feat_name =~ /_ncRNA/))  { #ncRNA added 01/14/22 by dfouts
		  $phageref->{$name}->{'target'} = $reref->{$asmbl_id}{$hitref->{$feat_name}->{'end5'}}->{'featname'};
	      }
	      else  {
		  $phageref->{$name}->{'target'} = $feat_name;
	      }
	      print ">>->>->>->>->>PICK_BEST_ATT: FINAL TARGET = $phageref->{$name}->{'target'} <<-<<-<<\n" if ($DEBUG);
              print "int ori = $attref->{$n}->{'ori'}\ttarget ori = $direction\n" if ($DEBUG);
              if ($direction ne $attref->{$n}->{'ori'}) { 
		  print "changing directions\n" if ($DEBUG);
		  $hold_5prime = $phageref->{$name}->{'5prime'};
		  $hold_3prime = $phageref->{$name}->{'3prime'};
		  $phageref->{$name}->{'5prime'} = $hold_3prime;
		  $phageref->{$name}->{'3prime'} = $hold_5prime;
		  print "old 5prime = $hold_5prime, new 5prime = $phageref->{$name}->{'5prime'}\n" if ($DEBUG);
                  print "old 3prime = $hold_3prime, new 3prime = $phageref->{$name}->{'3prime'}\n" if ($DEBUG);
		  &populate_atthash($attref, $ahref, $asmbl_id, $n, $direction, $DEBUG);
	      }
              $return_val = $n;
	      return($return_val);
	  }
      }
      $n++;
    }
    print ">>->>->>->>->>PICK_BEST_ATT: NO TARGET FOUND (return value = <$return_val> <<-<<-<<\n" if ($DEBUG);
    return($return_val);
}
sub get_TSD {

    my ($bbin, $prefix, $asmbl_id, $hitref, $phageref, $reref, $aref, $DEBUG, $direction, $type, $method, $name, $ahref, $start_head, $start_tail, $size_head, $size_tail, $old_attB, $old_attE, $figref, $window, $step, $hitsperwindow, $halfsize, $int, $graseref) = @_;  # needed to add int and check if an integrase was detected so that the 5' and 3' end RNAs can be properly checked 01/16/22.
    # variable $bbin added 07/06/2016 in order to pass the path to BLAST bin in order to run BLAST+
    # added $Bint or $Eint being passed from main so that we can check that tRNA is outside of integrase
    my $attL = "";
    my $attR = "";
    my @a = ();
    my $size_h = "";
    my $size_t = "";
    my $half_tRNA = "";
    my $temp_region_size = "";
    my $integrase_coord = "";
    my $head = "";
    my $tail = "";
    my $headl = "";
    my $taill = "";
    my $n = 0;
    my $best = 0;
    my %att_hash = ();
    my $size = $phageref->{$name}->{'3prime'} - $phageref->{$name}->{'5prime'};
    my $minsize = $size * .3; # did have .667 and 0.5, but missed B. anthracis ames LambdaBa03 att size and hacked and E.coli O157H7 tRNA-Arg
    my $data_prefix = "$prefix\_$type\_$name";
    my $head_filename = "$data_prefix\_head.seq";
    my $tail_filename = "$data_prefix\_tail.seq";
    my $featname = "";
    my $piggyback = ""; # variable to flag a region as potentially a piggyback insertion (04/25/2022)
    my $i = ""; # variable to store the featname of integrases while looping through @grase
    my $pright = ""; # to store the right end coord relative to the phage
    my $pleft = ""; # to store the left end coord relative to the phage
    my $right = ""; # to store the right end coord relative to the contig
    my $left = ""; # to store the left end coord relative to the contig
    my $int_end5 = "";
    my $int_end3 = "";
    my $int2use = "";
    my $diffRI = "";

    if ($type =~ /phage/)  {
	print "=======> $type <=======\n" if ($DEBUG);
	open (HEADFILE, ">$head_filename") || &write_log("4","1 can\'t open file $head_filename: $!\n");
	open (TAILFILE, ">$tail_filename") || &write_log("4","2 can\'t open file $tail_filename: $!\n");
    }
    $head = substr($ahref->{$asmbl_id}->{'sequence'},$start_head,$size_head);  # pull the sequence 400 bp upstream of the integrase (att sites no more than 200 bp usually)
    $tail = substr($ahref->{$asmbl_id}->{'sequence'},$start_tail,$size_tail);
    print "get_TSD: start_head = $start_head\n" if ($DEBUG);
    print "get_TSD: start_tial = $start_tail\n" if ($DEBUG);
    $headl = length($head);
    print "get_TSD: head_length = $headl\n" if ($DEBUG);
    &printFasta(*HEADFILE, "Bhead", $head);
    $taill = length($tail);
    print "get_TSD: tail_length = $taill\n" if ($DEBUG);
    &printFasta(*TAILFILE, "Btail", $tail);
    close (HEADFILE);
    close (TAILFILE);

########### Do searches based on $method ############

    if ($method eq "mummer")  {
	system ("mummer -mum -l 5 $tail_filename $head_filename | sort -nrk3 | head -2 > $data_prefix\_mummer.out");
	open (MUM, "<$data_prefix\_mummer.out") || &write_log("4","3 can\'t open file $data_prefix\_mummer.out: $!\n");
	$n = 1;
	while (<MUM>) {
	    chomp;
	    @a = split(/\s+/);
	    if ($a[3] ne " ")  {  # There was a mummer match, process data
		my $head_mummerE = $a[2]-1+$a[3]; # calculate end of Mummer query match (head)
		my $tail_mummerE = $a[1]-1+$a[3]; # calculate end of Mummer subject match (tail)
		
		$att_hash{$n}->{'attL_coordB'} = $start_head+$a[2];  # compute genome coord for attL 5' end
		$att_hash{$n}->{'attL_coordE'} = $start_head+$head_mummerE;  # compute genome coord for attL 3' end
		$att_hash{$n}->{'attR_coordB'} = $start_tail+$a[1];  # compute genome coord for attR 5' end
		$att_hash{$n}->{'attR_coordE'} = $start_tail+$tail_mummerE;  # compute genome coord for attR 3' end
		$att_hash{$n}->{'length'} = $a[3];
                print "get_TSD_mummer: $a[2] -> $head_mummerE = $a[1] -> $tail_mummerE\n" if ($DEBUG);
                &populate_atthash(\%att_hash, $ahref, $asmbl_id, $n, $direction, $DEBUG);
		$n++;
                if ($n == 3) { last; }
	    }
	}
        close (MUM);
    }
    elsif ($method eq "fasta33") {
        system ("fasta33 -H -m 9 -w 60 -O $data_prefix\_fasta.out -n -b 4 -q -d 4 -r +4/-5 -3 $head_filename $tail_filename > /dev/null");
        open (FASTA, "<$data_prefix\_fasta.out") || &write_log("4","3 can\'t open file $data_prefix\_fasta.out: $!\n");
        $n = 1;    
        # =================================================================================================================
        # space-delimited output from FASTA33 (-m 9) option:                                                
        # The best scores are:                                       opt  %_id  %_gid   fa alen  mn0  mx0  mn1  mx1 g_q g_l
        # Btail                                          (7116) [f]  417  0.973 0.973  417  111  139  249 3047 3157   0   0
        #
        # column number Description (for Perl), add 1 for Unix
        # 0     subject_id
        # 1     length of subject
        # 2     [f] = forward, [r] = reverse orientation
        # 3     opt score
        # 4     % id
        # 5     % gid
        # 6     fa
        # 7     length of alignment
        # 8     start of alignment on query (5' nucleotide match in query)
        # 9     end of alignment on query (3' nucleotide match in query)
        # 10    start of alignment on subject (5' nucleotide match in db hit)
        # 11    end of alignment on subject (3' nucleotide match in db hit)
        # 12    g_q
        # 13    g_l
        # =================================================================================================================
        
        while (<FASTA>)  {
	    chomp;
            @a = split(/\s+/);
            if ($a[0] =~ /Btail/)  {
		$att_hash{$n}->{'attL_coordB'} = $start_head+$a[8];
                $att_hash{$n}->{'attL_coordE'} = $start_head+$a[9];
                $att_hash{$n}->{'attR_coordB'} = $start_tail+$a[10];
                $att_hash{$n}->{'attR_coordE'} = $start_tail+$a[11];
                $att_hash{$n}->{'length'} = $a[7];
                print "get_TSD_fasta33: $a[0]\t$a[8]\t$a[9]\t$a[10]\t$a[11]\n" if ($DEBUG);
                &populate_atthash(\%att_hash, $ahref, $asmbl_id, $n, $direction, $DEBUG);
		$n++;
                if ($n == 3) { last; }
            }
	}
        close (FASTA);
    }
    elsif ($method eq "blast")  {
      # code to generate mini BLAST-formatted db and run BLAST+ provided by Dr. Camille Hankel at Dartmouth
      # I had to add "-task blastn" because the default is "megablast", which did not pick up any hits
      # -task blastn-short gave too many short matches [dfouts]
        system ("$bbin/makeblastdb -in $tail_filename -dbtype nucl > /dev/null 2>&1"); # make the tail sequence BLAST searchable
	system ("$bbin/blastn -task blastn -db $tail_filename -query $head_filename -out $data_prefix\_blast.out -outfmt 6 -strand plus -max_target_seqs 4 -num_threads 2 -dust no > /dev/null 2>&1");
        #system ("formatdb -i $tail_filename -p F"); # make the tail sequence BLAST searchable
	#system ("blastall -p blastn -d $tail_filename -i $head_filename -o $data_prefix\_blast.out -m 8 -S 1 -v 4 -b 4 -a 2 -F F > /dev/null 2>&1");
        open (BLAST, "<$data_prefix\_blast.out") || &write_log("4","4 can\'t open file $data_prefix\_blast.out: $!\n");
        $n = 1;
        # ========================================================
        # tab-delimited output from NCBI blastn (-m 8) option:

        # column number Description (for Perl), add 1 for Unix
        # 0	Query_id
        # 1	subject_id (Hit from db)
        # 2	% Identity
        # 3	length of alignment
        # 4	number or mismatches
        # 5	number of gaps
        # 6	start of alignment on query (5' nucleotide match in query)
        # 7	end of alignment on query (3' nucleotide match in query)
        # 8	start of alignment on subject (5' nucleotide match in db hit)
        # 9	end of alignment on subject (3' nucleotide match in db hit)
        # 10	e-value
        # 11	score (bits)
        # ========================================================
	while (<BLAST>) {#Bhead   Btail   98.15   108     2       0       142     249     3050    3157    4e-54    198
	    chomp;
            @a = split(/\t+/);
            $att_hash{$n}->{'attL_coordB'} = $start_head+$a[6];
            $att_hash{$n}->{'attL_coordE'} = $start_head+$a[7];
            $att_hash{$n}->{'attR_coordB'} = $start_tail+$a[8];
            $att_hash{$n}->{'attR_coordE'} = $start_tail+$a[9];
            $att_hash{$n}->{'length'} = $a[3];
            print "get_TSD_blast: $a[0]\t$a[6]\t$a[7]\t$a[8]\t$a[9]\n" if ($DEBUG);
            &populate_atthash(\%att_hash, $ahref, $asmbl_id, $n, $direction, $DEBUG);
            $n++;
            if ($n == 3) { last; } # end loop if $n is 3, only want top 2 hits max
        }
        close (BLAST);
    }

#################### CHOOSE ATT SITE ##################
    print ">get_TSD:  type = $type, region = $name\n" if ($DEBUG);
    if (!defined ($att_hash{'1'}->{'attL_coordB'}))  {
	print ">> get_TSD_1: $method did not find any putative att sites.\n" if ($DEBUG);
    }
    elsif (($type =~ /tm?RNA/) || ($type =~ /ncRNA/) || ($type =~ /_r2/)) { # if first round looking for t/tm/ncRNA att specifically
	# check for longest match
	$featname = $asmbl_id . "_" . $type;
	print "======>$featname<======\n" if ($DEBUG);
	$featname =~ s/phage_//; # remove the leading "phage_" that came from $type
	$featname =~ s/_r2//; # remove trailing "_r2" designation
	print "======>$featname<======\n" if ($DEBUG);
	print ">>get_TSD:  looking at a $type, so pick longest length\n" if ($DEBUG);
	($best) = &pick_longest_length($att_hash{'1'}->{'length'}, $att_hash{'2'}->{'length'});
	#### position of code to check for breaking of a phage region ####
	#### moved down to also work for non-tRNA att sites ####
    }
    else  { # find best att site not in a tRNA.  Need to check if breaking up the phage region like done above for tRNA insertions
	($best) = &pick_best_att($name, $direction, \%att_hash, $hitref, $phageref, $reref, $asmbl_id, $aref, $ahref, $DEBUG);
	print "BEST = <$best>\n" if ($DEBUG);
	if ($best == "")  { # if no target found, then search for the longest length, assuming not bogus
	    ($best) = &pick_longest_length($att_hash{'1'}->{'length'}, $att_hash{'2'}->{'length'});
	    if ($best == "")  { # if same length, then choose longest region
		($best) = &pick_longest_length($att_hash{'1'}->{'size'}, $att_hash{'2'}->{'size'});
	    }  
	}
    }
    print "BEST2 = <$best>\n" if ($DEBUG);

    #### relocated coded from above ####
    ### NEW CODE  if there are still hits >= 2 and phage end is < Eint or midpoint (with no int) when + direction, or phage end is > Bint or midpoint (with no int) when - direction then incorrect att site chosen
    # midpoint mod because Methylococcus bath had this problem with a tmRNA TSD that cut the phage region in half
    # the midpoint did not work for E. faecalis prophage phiEf11 if int present.
    if ($best > 0){
	$pright = $att_hash{$best}->{'pright'}; # store the region right end relative to the phage
	$pleft = $att_hash{$best}->{'pleft'}; # store the region left end relative to the phage, needed to determine if int is on the wrong side of the closest att site
	$right = $att_hash{$best}->{'right'}; # store the region right end relative to the contig
	$left = $att_hash{$best}->{'left'}; # store the region left end relative to the contig
	#ori+:   {pleft}--[int]---------------------------------------------------{pright}>
	#ori-:   <{pright}---------------------------------------------------[int]--{pleft}
	#contig: {left}-------------------------------------------------------------{right}
	
	print ">>>get_TSD:  right = $att_hash{$best}->{'right'}\tleft = $att_hash{$best}->{'left'}\n" if ($DEBUG);
	print ">>>get_TSD:  pright = $att_hash{$best}->{'pright'}\tpleft = $att_hash{$best}->{'pleft'}\n" if ($DEBUG);
	print ">>>get_TSD:  processing best match, phage ori = $att_hash{$best}->{'ori'}\n" if ($DEBUG);
	print ">>>get_TSD:  -hits = $figref->{nlowmult($step, $left)-$step}->{'counts'}, +hits = $figref->{nhimult($step, $right)+$step}->{'counts'}\n" if ($DEBUG);
	print ">>>get_TSD:  phage region $name ORF5 = $phageref->{$name}->{'ORF5'}, ORF3 = $phageref->{$name}->{'ORF3'}\n" if ($DEBUG);
	
	my $pluscoords = nlowmult($step, $right)+$step; # changed from pright to make it relative to the contig (ie, the 3' end of the contig)
	print ">>>get_TSD:  pluscoords(lowmult+) = $pluscoords\n" if ($DEBUG); # checking to see if next step has any phage genes when phage is --->(phage genes?)
	print ">>>get_TSD:  pluscoords_featname = $figref->{$pluscoords}->{'featnames'}[$#{$figref->{$pluscoords}->{'featnames'}}]\n" if ($DEBUG);
	print ">>>get_TSD:  currentplus_featname = $figref->{($pluscoords-$step)}->{'featnames'}[$#{$figref->{$pluscoords-$step}->{'featnames'}}]\n" if ($DEBUG);
	print ">>>get_TSD:  index of last figref array = <$#{$figref->{$pluscoords}->{'featnames'}}>\n" if ($DEBUG);
	
	if (($figref->{$pluscoords}->{'featnames'}[$#{$figref->{$pluscoords}->{'featnames'}}] eq $figref->{($pluscoords-$step)}->{'featnames'}[$#{$figref->{$pluscoords-$step}->{'featnames'}}]) && ($figref->{$pluscoords}->{'featnames'}[$#{$figref->{$pluscoords}->{'featnames'}}])) { #if the next step has the same featname at the first array position (AND has value), then bump up to next step so we don't count it as a hit (should be the last array element, not the first checked for duplication from step to step)
	    print ">>>get_TSD:  $figref->{$pluscoords}->{'featnames'}[$#{$figref->{$pluscoords}->{'featnames'}}] and $figref->{($pluscoords-$step)}->{'featnames'}[$#{$figref->{$pluscoords-$step}->{'featnames'}}] are the same\n" if ($DEBUG);
	    $pluscoords = $pluscoords+$step; # changed from pright
	    print ">>>get_TSD:  new pluscoords = $pluscoords\n" if ($DEBUG);
	}
	print ">>>get_TSD:  aref last window = <$aref->[$#{$aref}]>\n" if ($DEBUG); # check to see what last window is
	if ($pluscoords > $aref->[$#{$aref}]) { # if the pluscoords extends beyond the end of the contig, then set pluscoords to the last step of the contig
	    $pluscoords = nlowmult($step, $right); # trying to correct this by taking the lowmult instead of the himult. Changed from pright
	    print ">>>get_TSD:  off end of contig, new pluscoords = <$pluscoords>\n" if ($DEBUG);
	}
	
	my $minuscoords = nlowmult($step, $left)-$step; # changed from $pright  since if phage ori is + then that is not the 5' end of the contig
	print ">>>get_TSD:  minuscoords(lowmult-) = $minuscoords\n" if ($DEBUG); # (phage genes? <---)
	print ">>>get_TSD:   minuscoords_featname = <$figref->{$minuscoords}->{'featnames'}[0]>\n" if ($DEBUG);
	print ">>>get_TSD:  currentminus_featname = <$figref->{($minuscoords+$step)}->{'featnames'}[0]>\n" if ($DEBUG);
	
	if (($figref->{$minuscoords}->{'featnames'}[0] eq $figref->{($minuscoords+$step)}->{'featnames'}[0]) && ($figref->{$minuscoords}->{'featnames'}[0])) { #if the next step has the same featname at the first array position (AND has value), then bump down to next step
	    print ">>>get_TSD:  $figref->{$minuscoords}->{'featnames'}[0] and $figref->{($minuscoords+$step)}->{'featnames'}[0] are the same\n" if ($DEBUG);
	    $minuscoords = $minuscoords-$step; # changed from pright
	    print ">>>get_TSD:  new minuscoords = $minuscoords\n" if ($DEBUG);
	}
	if ($minuscoords < 0) { # check to make sure minuscoords does not go negative!
	    $minuscoords = 0; # set to 0 if negative
	     print ">>>get_TSD:  off beginning of contig, new minuscoords = <$minuscoords>\n" if ($DEBUG);
	}

	print ">>>get_TSD(check_tRNA_relative_to_int):  int = $int\n" if ($DEBUG);

	## check to see if the att site is OUTSIDE of the integrase and piggyback insertion ##
	if ($int > "") { # int predicted, check where tRNA is relative to int - we want t/tm/ncRNA on the outside of the int to be a target for integration
	    print ">>>> has integrase\n" if ($DEBUG);
	    print ">>>>> plushit counts = $figref->{$pluscoords}->{'counts'}\n" if ($DEBUG);
	    print ">>>>> minushit counts = $figref->{$minuscoords}->{'counts'}\n" if ($DEBUG);
	    print ">>>>> pleft = $pleft, pright = $pright, int = $int\n" if ($DEBUG);
	    my $graselength = scalar(@{$graseref}); # get the length of @grase (the number of integrases in the current phage region)
	    print ">>>>> number of integrases in phage region $name = $graselength\n" if ($DEBUG);

	    ## is the att site outside of the integrase? ##
	    if ((($att_hash{$best}->{'ori'} eq "+") && ($pleft > $int)) || (($att_hash{$best}->{'ori'} eq "-") && ($pleft < $int))) {
		print ">>>>> $type is on the wrong side of the integrase, no att picked\n" if ($DEBUG);
		$best = ""; # no att if incorrect RNA
	    }
	    ## if good hits outside of region, is it because it is a piggyback insertion or did we pick a bogus att side that breaks up our phage region? ##
	    elsif ((($att_hash{$best}->{'ori'} eq "+") && ($figref->{$pluscoords}->{'counts'} >= 2)) ||
		   (($att_hash{$best}->{'ori'} eq "-") && ($figref->{$minuscoords}->{'counts'} >= 2) && ($minuscoords > 0))) { # added 03/15/2022 dfouts
		#added check that first array element of figref->featnames was not the same as the current step
		#need to loop through @grase and check if abs($pright-nearest int end) <= 1000
		print ">>>> checking for piggyback insertion...\n" if ($DEBUG);
		print ">>>> pright = $pright\n" if ($DEBUG);
		print ">>>> piggyback = <$piggyback>\n" if ($DEBUG);
		for $i (@{$graseref}) {
		    # testing for:  [attL|--1st_int-->......phage_region......|attR..up_to_1kb..--2nd_int-->...]
		    # or:           [...<--2nd_int--..up_to_1kb..attR|......phage_region......<--1st_int--|attL]
		    $int_end5 = $hitref->{$i}->{'end5'};
		    $int_end3 = $hitref->{$i}->{'end3'};
		    # determine integrase end nearest the phage region
		    if ((($att_hash{$best}->{'ori'} eq "+") && ($int_end5 < $int_end3)) ||
			(($att_hash{$best}->{'ori'} eq "-") && ($int_end5 > $int_end3))) {
			$int2use = $int_end5;
		    }
		    elsif ((($att_hash{$best}->{'ori'} eq "+") && ($int_end3 < $int_end5)) ||
			   (($att_hash{$best}->{'ori'} eq "-") && ($int_end3 > $int_end5))) {
			$int2use = $int_end3;
		    }
		    
		    $diffRI = abs($pright-$int2use); # determine the distance between pright and the nearest integrase end
		    print ">>>> int ($i) (end5: $int_end5 .. end3: $int_end3) pright ($pright) - int2use ($int2use) = $diffRI \n" if ($DEBUG);
		    if ($diffRI <= 1000) {
			$piggyback = 1;
			print ">>>>> piggyback detected, ending loop\n" if ($DEBUG);
			next; #end loop
		    }
		}
		if ($piggyback) { # potential piggyback insertion, flag it and keep current att site (hopefully of the first phage insertion)
		    print ">>>>>> piggyback insertion predicted, keeping existing att sites...\n" if ($DEBUG);
		    $phageref->{$name}->{'piggyback'} = 1;
		}
	        elsif ((($att_hash{$best}->{'ori'} eq "-") && ($figref->{($minuscoords+$step)}->{'featnames'}[0] eq $figref->{$minuscoords}->{'featnames'}[$#{$figref->{$minuscoords}->{'featnames'}}])) ||
		    (($att_hash{$best}->{'ori'} eq "+") && ($figref->{($pluscoords)}->{'featnames'}[0] eq $figref->{$pluscoords-$step}->{'featnames'}[$#{$figref->{($pluscoords-$step)}->{'featnames'}}])))
		{
		    #if beginning of the next step is the same as the current last ORF, then subtract one from the hit count and retest for 5' area
		    #or  first ORF same as the last ORF in the previous step for the 3' area
		    
		    if ((($att_hash{$best}->{'ori'} eq "+") && (($figref->{$pluscoords}->{'counts'}-1) >= 2)) ||
			(($att_hash{$best}->{'ori'} eq "-") && (($figref->{$minuscoords}->{'counts'}-1) >= 2) && ($minuscoords > 0))) { # added 03/15/2022 dfouts
			print ">>>>> $type predicted att breaks the phage region after considering overlapping hits at the ends, no att picked\n" if ($DEBUG);
			$best = ""; # no att if there are still phage hits outside of the att opposite of the integrase
		    }
		    else {
			print ">>>>> this att site does not break the phage region after all!\n" if ($DEBUG);
		    }
		}
		else { # no piggyback phage, so probably breaking a single phage region
		    print ">>>>> $type predicted att breaks the phage region, no att picked\n" if ($DEBUG);
		    $best = ""; # no att if there are still phage hits outside of the att opposite of the integrase
		}
	    }
	}
	else  { # truncated region with no int so no att site
	    print ">>>> no integrase\n" if ($DEBUG);
	    if (($att_hash{$best}->{'ori'} eq "+") && ($figref->{nlowmult($window, $pright)+$step}->{'counts'} >= 2) && ($pleft < $halfsize)|| 
		(($att_hash{$best}->{'ori'} eq "-") && ($figref->{nlowmult($window, $pright)-$step}->{'counts'} >= 2) && ($pleft > $halfsize))) {
		print ">>>>> truncated region, no att picked\n" if ($DEBUG);
		$best = ""; # no att if truncates region
	    }
	}
    }
    
    #### END relocated code ####
    
    ### check if the att size is too far from the integrase ###
    ## need to check to see if this was implemented elsewhere, like in the main program
    # may want to also do this in pick_best_att subroutine so that alternative att sites can be considered
    # added 02/25/22 by dfouts to prevent the rare cases where distances greater than 10kb were being erroneously chosen
    #if ((abs($att_hash{$best}->{'left'}-$int) > 2000) && (abs($att_hash{$best}->{'right'}-$int) > 2000)) {
    #	$best = ""; # set to zero to unconsider this att site
    #}
    
    ### write results to disk if att found ###
    if ($best > 0) {
        print ">>>>get_TSD: att picked\n" if ($DEBUG);
        $temp_region_size = abs($att_hash{$best}->{'attR_coordE'} - $att_hash{$best}->{'attL_coordB'}) + 1;
	#$phageref->{$name}->{'original_direction'} = $phageref->{$name}->{'direction'}; # added 10/7/2008 to enable recovery of the original direction when deleting an att site because house-keeping genes were added to region.
        $phageref->{$name}->{'direction'} = $att_hash{$best}->{'ori'}; # moved 05/22/2007 when region < 9kb, phage direction was not populated
        print "get_TSD:TEMP_REGION_SIZE: $temp_region_size\t$minsize, $size\n" if ($DEBUG);
	#print "get_TSD:Original phage direction = $phageref->{$name}->{'original_direction'}\n" if ($DEBUG);
        print "get_TSD:Phage direction = $phageref->{$name}->{'direction'}\n" if ($DEBUG);
	if (($temp_region_size >= $minsize) && ($temp_region_size >= 6000)){ # if the region defined by putative att sites at least 9 kbp originally (to elminate short junk regions with tRNA homology), but changed to 6kb on 10/24/2013 for Acinetobacter non-strict region on a plasmid (not phage, but may be interesting)
	#if (($temp_region_size >= $halfsize) || ($att_hash{$best}->{'length'} >= 15)) { # if the region defined by putative att sites is at least 50% (1/2) the size of the BLAST/HMM defined region or the size of the att size is >= 15 bp, ok
	    ## NOTE: this is a quick fix for 2 problems
	    #  1) "Piggy-back" phages: multiphages in tandem that insert into the same target
	    #  2) tRNAs that lack an integrase between the attL and attR
	    $attL = $att_hash{$best}->{'attL'}; # relative to genome
	    $attR = $att_hash{$best}->{'attR'}; # relative to genome
            print "get_TSD: name = $name, attL = $attL, attR = $attR\n" if ($DEBUG);
            open (ATTFILE, ">$data_prefix\_att.out") || &write_log("4","can\'t open file $data_prefix\_att.out: $!\n");

	    # need to also check the distance from integrase to att (should not be larger than 2 kb one would think) - may do this in find_att_sites in main program
	    if ((!exists $hitref->{"attL_$asmbl_id\_$name"}->{'annotation'}) || # if no previous att for this region or if a second round contains an expanded first att
		($attL =~ /$hitref->{"attL_$asmbl_id\_$name"}->{'annotation'}/))  {
		print ATTFILE ">$asmbl_id Phage $name putative half att sites\n";
		print ATTFILE "attL:\t$attL\t$att_hash{$best}->{'attL5'} to $att_hash{$best}->{'attL3'}\n";
		print ATTFILE "attR:\t$attR\t$att_hash{$best}->{'attR5'} to $att_hash{$best}->{'attR3'}\n";
		print "attL:\t$attL\t$att_hash{$best}->{'attL5'} to $att_hash{$best}->{'attL3'}\n" if ($DEBUG);
		print "attR:\t$attR\t$att_hash{$best}->{'attR5'} to $att_hash{$best}->{'attR3'}\n" if ($DEBUG);
                print "5prime_att = $att_hash{$best}->{'pleft'}, 3prime_att = $att_hash{$best}->{'pright'}\n" if ($DEBUG);
		### Adjust the phage 5' and 3' boundarys to the new values
		$phageref->{$name}->{'5prime_att'} = $att_hash{$best}->{'pleft'}; #attL5 before relative to phage
		$phageref->{$name}->{'3prime_att'} = $att_hash{$best}->{'pright'}; #attR3 before relative to phage
		$phageref->{$name}->{'left'} = $att_hash{$best}->{'left'}; # genome coords no matter the phage direction
		$phageref->{$name}->{'right'} = $att_hash{$best}->{'right'}; # genome coords no matter the phage direction
		if ($old_attB > 0)  {
		    delete $phageref->{$name}->{'memberhash'}{$old_attB};  # remove old attB information
		    delete $phageref->{$name}->{'memberhash'}{$old_attE};  # remove old attE information
		}
		$phageref->{$name}->{'memberhash'}{$att_hash{$best}->{'attB'}}->{'featname'} = "attL_$asmbl_id\_$name";
		$phageref->{$name}->{'memberhash'}{$att_hash{$best}->{'attE'}}->{'featname'} = "attR_$asmbl_id\_$name";
		$hitref->{"attR_$asmbl_id\_$name"}->{'annotation'} = "$attR";
		$hitref->{"attL_$asmbl_id\_$name"}->{'annotation'} = "$attL";
		print "$att_hash{$best}->{'attB'} ($attL), $att_hash{$best}->{'attE'} ($attR)\n" if ($DEBUG);
		close (ATTFILE);
		return("1", $att_hash{$best}->{'attB'}, $att_hash{$best}->{'attE'});
	    }
	}
	else {
	    print ">> get_TSD_2: $method had a  putative att site for $type region $name, but the size of the region was < 6 Kb.\n" if ($DEBUG);
            $phageref->{$name}->{'target'} = ""; # clear out target
        }
    }
    else { # no att sites predicted
	print ">> get_TSD_3: $method did not find any putative att sites for $type region $name.\n" if ($DEBUG);
	$phageref->{$name}->{'target'} = ""; # clear out target
    }
    return("0");
}

sub adjust_featname {
    my ($feat_name,$asmbl_id) = @_;

    if ($feat_name =~ /RNA/) {
	$feat_name = $asmbl_id . "_" . $feat_name;
    }
    return ($feat_name);
}

sub add_ORFs_to_region  {

    my ($key,$asmbl_id,$phageref,$hitref,$reref,$aref,$okref,$BCNFH,$DEBUG) = @_;
    my $i = "";
    my $j = "";
    my $new_end5 = $phageref->{$key}->{'left'};
    my $new_end3 = $phageref->{$key}->{'right'};
    my $ORF5_featname = &adjust_featname($phageref->{$key}->{'ORF5'},$asmbl_id);
    my $ORF3_featname = &adjust_featname($phageref->{$key}->{'ORF3'},$asmbl_id);
    my $bad_name = "";
    my $B_pos = "";
    my $B_end5 = "";
    my $B_end3 = "";
    my $B_featname = "";
    my $E_pos = "";
    my $E_end5 = "";
    my $E_featname = "";
    my $E_end3 = "";

#### work on 5' end of phage region first ####

    print "add_ORFs_to_region5:  ORF5 = $ORF5_featname, ORF3 = $ORF3_featname\n" if ($DEBUG);
    $B_pos = $hitref->{$ORF5_featname}->{'array_pos'}; # get position of the first ORF in phage region (before att search). Removed "- 1"
## do I need the next 3 lines?
    $B_end5 = $aref->[$B_pos]; # get end5 of new ORF
    $B_featname = &adjust_featname($reref->{$asmbl_id}{$B_end5}->{'featname'},$asmbl_id); # get cleaned feat_name of new ORF
    $B_end3 = $hitref->{$B_featname}->{'end3'};
    print "add_ORFs_to_region5:  B_pos = $B_pos, B_end5init = $B_end5, new_end5 = $new_end5, B_end3 = $B_end3\n" if ($DEBUG);

## added 10/06/2008 to fix a problem where house-keeping genes were getting added because of a bogus att find
# pre-loop and check if ok_com_name, if >= 2, then adjust to before att coords ### NOTE, we are only checking to make sure the com_name is ok, the next "if" section will add the ORFs to the region
    for ($i = $B_pos; (($B_end5 >= $new_end5) && ($B_end3 >= $new_end5) && ($i >= 0)); $i--)  { #$i > 0 added 09/05/07
	if (($B_featname =~ /tm?RNA/) || ($B_featname =~ /ncRNA/)) { # added 04/15/2013 by dfouts for cases where tRNA/tmRNA gene is at end
	    #ncRNA added 01/14/22 by dfouts
	    #tRNA/tmRNA featnames stored differently in hithash (to make unique added asmbl_id to featname) and rehash
	    $B_featname = $asmbl_id . "_" . $B_featname;
	}
        print "add_ORFs_to_region5_preloop:  Bend5 = $B_end5, Bend3 = $B_end3, new_end5 = $new_end5, Bfeatname = $B_featname\n" if ($DEBUG);
	if ($hitref->{$reref->{$asmbl_id}{$B_end5}->{'featname'}}->{'hit'} == 1) { # if has it hit, it's good (added this 7/11/22)
	    print "add_ORFs_to_region5_preloop:  OK dbHIT Bfeatname = $B_featname\n" if ($DEBUG);
	    if (($okref->{$hitref->{$B_featname}->{'clean_name'}} == 1) || ($B_featname =~ /tm?RNA/) || ($B_featname =~ /ncRNA/)){
		print "add_ORFs_to_region5_preloop:  ->OK feat_name Bfeatname = $B_featname\n" if ($DEBUG);
	    }
	    else {
		print $BCNFH "<$hitref->{$B_featname}->{'clean_name'}>\t<$hitref->{$B_featname}->{'clean_name'}>\t$B_featname\tdb_hit\n"; # record bad com_names anyway so we can fix.
		print "add_ORFs_to_region5_preloop:  ->BAD ($hitref->{$B_featname}->{'clean_name'}), feat_name Bfeatname = $B_featname, bad_name = <$bad_name>\n" if ($DEBUG);
	    }
	}
        elsif (($okref->{$hitref->{$B_featname}->{'clean_name'}} == 1) || ($B_featname =~ /tm?RNA/) || ($B_featname =~ /ncRNA/)){ # added tm?RNA check 04/16/13 because some att sites were being missed, thinking a tRNA is bad to have in a phage region
	#ncRNA check added 01/14/22 by dfouts
            print "add_ORFs_to_region5_preloop:  OK feat_name Bfeatname = $B_featname\n" if ($DEBUG);
	}
        else {
	    $bad_name++;
	    print $BCNFH "<$hitref->{$B_featname}->{'clean_name'}>\t<$hitref->{$B_featname}->{'clean_name'}>\t$B_featname\n"; # added 09/11/14 by dfouts to keep track of bad com_names.
	    print "add_ORFs_to_region5_preloop:  BAD ($hitref->{$B_featname}->{'clean_name'}), feat_name Bfeatname = $B_featname, bad_name = $bad_name\n" if ($DEBUG);
	}
	### moved to end of loop because we were skipping current gene dfouts 09/11/2014
        # problem was that targets of insertion were being added when not supposed to and therefore the att sites were ignored because these targets had "bad com_names"
	# also added some missing ok com names to the com_file
	$B_end5 = $aref->[$i-1] if ($i >= 0); # added -1 05/06/22 because it was not storing the next ORF in B_end5 and counting the same bad com_name twice
	$B_featname = &adjust_featname($reref->{$asmbl_id}{$B_end5}->{'featname'},$asmbl_id); # get cleaned feat_name of new ORF
	$B_end3 = $hitref->{$B_featname}->{'end3'};
    }
##
    if ($bad_name < 1)  {  # changed from $bad_name < 2 09/19/13 dfouts - i.e., not allowing any bogus com_names to slip through
	$B_pos = $hitref->{$ORF5_featname}->{'array_pos'}; # get position of the first ORF in phage region (before att search).  Removed "- 1"
	$B_end5 = $aref->[$B_pos]; # get end5 of new ORF
	$B_featname = &adjust_featname($reref->{$asmbl_id}{$B_end5}->{'featname'},$asmbl_id); # get cleaned feat_name of new ORF
	if (($B_featname =~ /tm?RNA/) || ($B_featname =~ /ncRNA/)) { # added 04/15/2013 by dfouts for cases where tRNA/tmRNA gene is at end
	    #tRNA/tmRNA featnames stored differently in hithash (to make unique added asmbl_id to freatname)than rehash
	    #ncRNA added 01/14/22 by dfouts
	    $B_featname = $asmbl_id . "_" . $B_featname;
	}
	$B_end3 = $hitref->{$B_featname}->{'end3'};
	for ($i = $B_pos; (($B_end5 >= $new_end5) && ($B_end3 >= $new_end5) && ($i >= 0)); $i--)  { #$i > 0 added 09/05/07
	    $B_end5 = $aref->[$i];
	    $B_featname = &adjust_featname($reref->{$asmbl_id}{$B_end5}->{'featname'},$asmbl_id); # get cleaned feat_name of new ORF
	    $B_end3 = $hitref->{$B_featname}->{'end3'};
	    print "add_ORFs_to_region5:  Bend5 = $B_end5, Bend3 = $B_end3, new_end5 = $new_end5, Bfeatname = $B_featname\n" if ($DEBUG);
	    if ((($B_end5 >= $new_end5) && ($B_end3 >= $new_end5)) && ($B_featname ne $phageref->{$key}->{'target'})){
		if (!exists $phageref->{$key}->{'memberhash'}{$B_end5}->{'featname'}) { # added this check to make sure it is not already in the memberhash
		    $phageref->{$key}->{'memberhash'}{$B_end5}->{'featname'} = $B_featname;
		    $phageref->{$key}->{'ORF5'} = $B_featname;
		    print "add_ORFs_to_region5: added orf $B_featname to memberhash for phage region $key\n" if ($DEBUG);
		}
		else { # make it the 5' ORF anyway
		    $phageref->{$key}->{'ORF5'} = $B_featname;
		    print "add_ORFs_to_region5: orf $B_featname is in memberhash already for phage region $key, making ORF5 anyway\n" if ($DEBUG);
		}
	    }
	}
    }
    else {
        print "add_ORFs_to_region5_return:  bad_name = $bad_name.\n" if ($DEBUG);
	return(1); #flag as needed the ends reset and $found = 0
    }

#### now, work on 3' end of phage region ####
    ## changed array_pos + 1 to array_pos - 1 to removing it because this caused an infinate loop situation 05/14/07
    $E_pos = $hitref->{$ORF3_featname}->{'array_pos'}; # get position of the last ORF in the phage region (before att search)
## do I need the next 3 lines?
    $E_end5 = $aref->[$E_pos]; # get end5 of new ORF
    $E_featname = &adjust_featname($reref->{$asmbl_id}{$E_end5}->{'featname'},$asmbl_id); # get cleaned feat_name of new ORF
    $E_end3 = $hitref->{$E_featname}->{'end3'};
    print "add_ORFs_to_region3: arraysize = $#{$aref}, arraypos = $E_pos, feat_name = $E_featname, end5 = $E_end5, end3 = $E_end3, new_end5 = $new_end5, new_end3 = $new_end3\n" if ($DEBUG);

## added 10/06/2008 to fix a problem where house-keeping genes were getting added because of a bogus att find
# pre-loop and check if ok_com_name, if >= 2, then adjust to before att coords
    $bad_name = ""; # reset this variable

    for ($j = $E_pos; (($E_end5 <= $new_end3) && ($E_end3 <= $new_end3) && ($j <= $#{$aref})); $j++)  { # array size condition added 09/05/07
	if (($E_featname =~ /tm?RNA/) || ($E_featname =~ /ncRNA/)) { # added 04/15/2013 by dfouts for cases where tRNA/tmRNA gene is at end
	    #tRNA/tmRNA featnames stored differently in hithash (to make unique added asmbl_id to freatname)than rehash
	    #ncRNA added 01/14/22 by dfouts
	    $E_featname = $asmbl_id . "_" . $E_featname;
	}
        print "add_ORFs_to_region3_preloop:  j = $j, Eend5 = $E_end5, Eend3 = $E_end3, new_end3 = $new_end3, Efeatname = $E_featname\n" if ($DEBUG);
	if ($hitref->{$reref->{$asmbl_id}{$E_end5}->{'featname'}}->{'hit'} == 1) { # if has it hit, it's good (added this 7/11/22)
	    print "add_ORFs_to_region3_preloop:  OK dbHIT Efeatname = $E_featname\n" if ($DEBUG);
	    if (($okref->{$hitref->{$E_featname}->{'clean_name'}} == 1) || ($E_featname =~ /tm?RNA/) || ($E_featname =~ /ncRNA/)){
		print "add_ORFs_to_region3_preloop:  ->OK feat_name Efeatname = $E_featname\n" if ($DEBUG);
	    }
	    else {
		print $BCNFH "<$hitref->{$E_featname}->{'clean_name'}>\t<$hitref->{$E_featname}->{'clean_name'}>\t$E_featname\tdb_hit\n"; # record bad com_names anyway so we can fix.
		print "add_ORFs_to_region3_preloop:  ->BAD ($hitref->{$E_featname}->{'clean_name'}), feat_name Efeatname = $E_featname, bad_name = <$bad_name>\n" if ($DEBUG);
	    }
	}
        elsif (($okref->{$hitref->{$E_featname}->{'clean_name'}} == 1) || ($E_featname =~ /tm?RNA/) || ($E_featname =~ /ncRNA/)) { # added tm?RNA check 04/16/13 because some att sites were being missed, thinking a tRNA is bad to have in a phage region
	    #ncRNA added 01/14/22 by dfouts
            print "add_ORFs_to_region3_preloop:  OK feat_name Efeatname = $E_featname\n" if ($DEBUG);
	}
        else {
	    $bad_name++;
	    print $BCNFH "<$hitref->{$E_featname}->{'clean_name'}>\t<$hitref->{$E_featname}->{'clean_name'}>\t$E_featname\n"; # added 09/11/14 by dfouts to keep track of bad com_names.
	    print "add_ORFs_to_region3_preloop:  BAD ($hitref->{$E_featname}->{'clean_name'}), feat_name Efeatname = $E_featname, bad_name = $bad_name\n" if ($DEBUG);
	}
	### moved to end of loop because we were skipping current gene dfouts 09/11/2014
        # problem was that targets of insertion were being added when not supposed to and therefore the att sites were ignored because these targets had "bad com_names"
	$E_end5 = $aref->[$j+1] if ($j+1 <= $#{$aref}); # added +1 05/06/22 because it was not storing the next ORF in E_end5 and counting the same bad com_name twice
	print "-->next E_end5 = $E_end5, j = $j\n" if ($DEBUG);
	$E_featname = &adjust_featname($reref->{$asmbl_id}{$E_end5}->{'featname'},$asmbl_id); # get cleaned feat_name of new ORF
	$E_end3 = $hitref->{$E_featname}->{'end3'};
    }
    if ($bad_name < 1)  {  # changed from $bad_name < 2 09/19/13 dfouts - i.e., not allowing any bogus com_names to slip through
        $E_pos = $hitref->{$ORF3_featname}->{'array_pos'}; # get position of the last ORF in the phage region (before att search). Removed " -1"
        $E_end5 = $aref->[$E_pos]; # get end5 of new ORF
        $E_featname = &adjust_featname($reref->{$asmbl_id}{$E_end5}->{'featname'},$asmbl_id); # get cleaned feat_name of new ORF
	if (($E_featname =~ /tm?RNA/) || ($E_featname =~ /ncRNA/)) { # added 04/15/2013 by dfouts for cases where tRNA/tmRNA gene is at end
	    #tRNA/tmRNA featnames stored differently in hithash (to make unique added asmbl_id to freatname)than rehash
	    #ncRNA added 01/14/22 by dfouts
	    $E_featname = $asmbl_id . "_" . $E_featname;
	}
        $E_end3 = $hitref->{$E_featname}->{'end3'};
        for ($j = $E_pos; (($E_end5 <= $new_end3) && ($E_end3 <= $new_end3) && ($j <= $#{$aref})); $j++)  { # array size condition added 09/05/07
	    $E_end5 = $aref->[$j];
	    $E_featname = &adjust_featname($reref->{$asmbl_id}{$E_end5}->{'featname'},$asmbl_id); # get cleaned feat_name of new ORF
	    $E_end3 = $hitref->{$E_featname}->{'end3'};
	    print "add_ORFs_to_region3:  Eend5 = $E_end5, Eend3 = $E_end3, new_end3 = $new_end3, Efeatname = $E_featname\n" if ($DEBUG);
	    if ((($E_end5 <= $new_end3) && ($E_end3 <= $new_end3)) && ($E_featname ne $phageref->{$key}->{'target'})){
		if (!exists $phageref->{$key}->{'memberhash'}{$E_end5}->{'featname'}) { # added this check to make sure it is not already in the memberhash
		    $phageref->{$key}->{'memberhash'}{$E_end5}->{'featname'} = $E_featname;
		    $phageref->{$key}->{'ORF3'} = $E_featname;
		    print "add_ORFs_to_region3: added orf $E_featname to memberhash for phage region $key\n" if ($DEBUG);
		}
		else { # make it the 5' ORF anyway
		    $phageref->{$key}->{'ORF3'} = $E_featname;
		    print "add_ORFs_to_region3: orf $E_featname is in memberhash already for phage region $key, making ORF3 anyway\n" if ($DEBUG);
		}
	    }
	}
    }
    else {
        print "add_ORFs_to_region3_return:  bad_name = $bad_name, returning value of 1\n" if ($DEBUG);
	return(1); #flag as needed the ends reset and $found = 0
    }
return(0); #default return (ok_com_names)
}

sub rekey_regions { # subroutine to loop through and renumber phage regions when we merge or delete them to keep consecutive numbering from 1..max
  my ($r,$phages,$DEBUG) = @_;
  my $maxr = keys %{$phages}; # get max number of phage regions currently in %phage_hash
  $maxr = $maxr + 1; # increment this by one since we already deleted the r+1 phage region during the merge before this subroutine is called
  print "Rekey_regions: Starting region = <$r>, max phage regions = <$maxr>\n" if ($DEBUG);
  my $ctr = 0;
  until ($r == ($maxr-1)) {
    $r++;
    $ctr++;
    print "$ctr)*******Rekey before********-> region = $r, maxr = $maxr, Seed_end5 before = $phages->{$r}->{'seed_end5'}\n" if ($DEBUG);
    $phages->{$r} = $phages->{$r+1}; #actual renumbering by swapping pointers
    print "*******Rekey after********-> region = $r, Seed_end5 after = $phages->{$r}->{'seed_end5'}\n" if ($DEBUG);
  }
  delete($phages->{$maxr}); # delete the old max key since we copied the pointer location to a new hash key ($max -1)
  $maxr = keys %{$phages};
  print "Rekey_regions: Ending region = <$r>, new max phage regions = <$maxr>\n" if ($DEBUG);
}

sub determine_region_type  {  # scan through regions and HMM hits to determine region type

  my ($regions,$phageref,$hitref,$HMMref,$serineref,$tyrosineref,$DEBUG) = @_;
  my $key = "";
  my $yek = "";
  my $hmm = "";
  my $phage_5prime = "";
  my $phage_3prime = "";
  my $score = "";
  my $trusted = "";
  my $noise = "";
  my $size = "";
  my $feat_name = "";
  my %mu_hash = ('PF02316' => 1,
                 'PF02914' => 1,
                 'PF06074' => 1,
                 'PF07030' => 1);

  my %trash_HMM = ('TIGR02224' => 1, # XerC
                   'TIGR02225' => 1, # XerD
                   'TIGR02249' => 1);# integron

  my %core_cnts = (); # new (10/05/22) local hash to store core hmms counts in order to flag potential piggyback phage regions
  my $regions = keys %{$phageref};
  my $m = ""; # variable for holding phage region just before a deleted region
  print "determine_region_type: #Phage_hash keys  = <$regions>\n" if ($DEBUG);

LOOP:
  foreach $key (sort {$a <=> $b} keys %{$phageref}) { #added sorting of the regions to make sequencial for renumbering when a regions is deleted
      $phageref->{$key}->{'core_HMM'} = 0;
      $phageref->{$key}->{'core_BLAST'} = 0; # added 08/22/2008 by Derrick Fouts
      $phageref->{$key}->{'above_noise_core_HMM'} = 0;
      $phageref->{$key}->{'lytic_HMM'} = 0;
      $phageref->{$key}->{'tail_HMM'} = 0;
      $phageref->{$key}->{'RBP_HMM'} = 0;
      $phageref->{$key}->{'RBP'} = 0; # added to make the default be zero in the PFPR_tab.txt file
      $phageref->{$key}->{'Mu_HMM'} = 0;
      $phageref->{$key}->{'int_HMM'} = 0;
      $phageref->{$key}->{'serine_recomb'} = 0;
      $phageref->{$key}->{'toxin_HMM'} = 0;
      $phageref->{$key}->{'mucin_HMM'} = 0;
      $phageref->{$key}->{'trash'} = 0;
      $phageref->{$key}->{'PF00239'} = 0;
      $phageref->{$key}->{'att_distance'} = 0; # moved from find_att_sites so that every phage region will have this set to zero
      $phageref->{$key}->{'piggyback'} = 0; # set new piggyback variable to zero 04/25/2022

      foreach $yek (sort {$a <=> $b} keys %{$phageref->{$key}->{'memberhash'}}) {
	$feat_name = $phageref->{$key}->{'memberhash'}{$yek}->{'featname'};
	print "====>determine_region_type: Region <$key> featname = <$feat_name>\n" if ($DEBUG);
	if (exists $hitref->{$feat_name}->{'Large Terminase'}) {
	  $phageref->{$key}->{'core_BLAST'}++;
	  print "====>determine_region_type: Large Terminase core_BLAST = <$hitref->{$feat_name}->{'Large Terminase'}>\n" if ($DEBUG);
	  $core_cnts{'Large Terminase'}++; #advance Large Term core gene counts
	  if ($core_cnts{'Large Terminase'} > 1) { # if more than one Large Terminase gene in this region, flag as a piggyback
	      $phageref->{$key}->{'piggyback'} = 1;
	  }
	}
	elsif (exists $hitref->{$feat_name}->{'Portal'}) {
	  $phageref->{$key}->{'core_BLAST'}++;
	  print "====>determine_region_type: Portal core_BLAST = <$hitref->{$feat_name}->{'Portal'}>\n" if ($DEBUG);
	  $core_cnts{'Portal'}++; #advance Portal core gene counts
	  if ($core_cnts{'Portal'} > 1) { # if more than one Portal gene in this region, flag as a piggyback
	      $phageref->{$key}->{'piggyback'} = 1;
	  }
	}
	elsif (exists $hitref->{$feat_name}->{'Major Capsid'})  {
	  $phageref->{$key}->{'core_BLAST'}++;
	  print "====>determine_region_type: Major Capsid core_BLAST = <$hitref->{$feat_name}->{'Major Capsid'}>\n" if ($DEBUG);
	  $core_cnts{'Major Capsid'}++; #advance Large Term core gene counts
	  if ($core_cnts{'Major Capsid'} > 1) { # if more than one Major Capsid gene in this region, flag as a piggyback
	      $phageref->{$key}->{'piggyback'} = 1;
	  }
	}
	elsif (exists $hitref->{$feat_name}->{'RBP'})  {
	  $phageref->{$key}->{'RBP'}++;
	  print "====>determine_region_type: RBP BLAST = <$hitref->{$feat_name}->{'RBP'}>\n" if ($DEBUG);
	} # added 03/28/18 by Derick Fouts to count # receptor-binding proteins identified per prophage region

	foreach $hmm (keys %{$hitref->{$feat_name}->{'hmm'}}) {
	  $score = $hitref->{$feat_name}->{'hmm'}{$hmm}->{'score'};
	  $trusted = $hitref->{$feat_name}->{'hmm'}{$hmm}->{'trusted'};
	  $noise = $hitref->{$feat_name}->{'hmm'}{$hmm}->{'noise'};
	  if (($tyrosineref->{$hmm} == 1) || ($serineref->{$hmm} == 1)) {
	    $phageref->{$key}->{'int_HMM'}++;
	    push (@{$phageref->{$key}->{'integrases'}}, $feat_name);
	    print "====>determine_region_type: INT: phage# = $key, 5' = $yek, ORF = $feat_name, hmm = >$hmm<\n" if($DEBUG);
	  }
	  # if HMM hit to the XerC, XerD or integron integrase specific HMMs TIGR02224, TIGR02225 or TIGR02249, then forget it, don't look at this region
	  if ($trash_HMM{$hmm} == 1) { 
	    $phageref->{$key}->{'trash'}++;
	    print "====>determine_region_type: hmm = <$hmm>, key = <$key>, feat_name = <$feat_name>, noise = <$noise>, trusted = <$trusted>, score = <$score>\n" if ($DEBUG);
	  }
	  # if int hit first, then find out if it also hits a XerC/D or integron integrase, delete
	  elsif (($phageref->{$key}->{'int_HMM'} > 0) && ($trash_HMM{$hmm} == 1)) {
	    $phageref->{$key}->{'trash'}++
	  }

	  if (($HMMref->{'Core HMM'}{$hmm} == 1) && ($score >= $trusted)){ # core hit must be >= trusted (getting some real junk hits in oral microbes)
	    $phageref->{$key}->{'core_HMM'}++;
	    print ">>>>>>>>CORE: phage = $key, hmm = >$hmm<\n" if ($DEBUG);
	    # new 10/05/22
	    $core_cnts{$hmm}++; #advance specific core gene counts
	    if ($core_cnts{$hmm} > 1) { # if more than one specific core gene hmm in this region, flag as a piggyback
		$phageref->{$key}->{'piggyback'} = 1;
	    }
	  }
	  elsif (($HMMref->{'Core HMM'}{$hmm} == 1) && ($score >= $noise)){
	    $phageref->{$key}->{'above_noise_core_HMM'}++;
	    print ">>>>>>>>CORE_ABOVE_NOISE: phage = $key, hmm = >$hmm<\n" if ($DEBUG);
	  }
	  elsif ($HMMref->{'Tail and Baseplate HMM'}{$hmm} == 1) {
	    $phageref->{$key}->{'tail_HMM'}++;
	    print ">>>>>>>>TAIL: phage = $key, hmm = >$hmm<\n" if ($DEBUG);
	  }
	  elsif ($HMMref->{'RBP HMM'}{$hmm} == 1) {
	    $phageref->{$key}->{'RBP_HMM'}++;
	    print ">>>>>>>>RBP: phage = $key, hmm = >$hmm<\n" if ($DEBUG);
	  }
	  elsif ($HMMref->{'Lytic enzyme HMM'}{$hmm} == 1) {
	    $phageref->{$key}->{'lytic_HMM'}++;
	    print ">>>>>>>>LYTIC: phage = $key, hmm = >$hmm<\n" if ($DEBUG);
	  }
	  elsif ($HMMref->{'Toxins HMM'}{$hmm} == 1) {
	    $phageref->{$key}->{'toxin_HMM'}++;
	    print ">>>>>>>>TOXIN: phage = $key, hmm = >$hmm<\n" if ($DEBUG);
	  }
	  elsif ($HMMref->{'Mucin-binding HMM'}{$hmm} == 1) {
	    $phageref->{$key}->{'mucin_HMM'}++;
	    print ">>>>>>>>MUCIN: phage = $key, hmm = >$hmm<\n" if ($DEBUG);
	  }
	  elsif ($hmm eq "PF04606") {
	    $phageref->{$key}->{'PF04606'} = 1;
	  }
	  elsif ($hmm eq "TIGR01613") {
	    $phageref->{$key}->{'TIGR01613'} = 1;
	  }
	  elsif ($hmm eq "PF00078") {
	    $phageref->{$key}->{'PF00078'} = 1;
	  }
	  elsif ($serineref->{$hmm} == 1) { # Serine recombinase HMM
	    $phageref->{$key}->{'serine_recomb'}++;
	  }
	  if ($mu_hash{$hmm} == 1) {
	    $phageref->{$key}->{'Mu_HMM'}++;
	  }
	}
      }
      # moved this outside the memberhash loop so that it has a chance to record core matches on 01/26/2019
      if (($phageref->{$key}->{'trash'} > 0) && ($phageref->{$key}->{'core_HMM'} == 0) && ($phageref->{$key}->{'core_BLAST'} == 0)) { # if a trash code flagged and no core HMMs, then junk the region: new 05/23/2006. Loosened up on criteria in cases where we hae core BLAST hits, but the proteins are too distant to be picked up by HMMs 01/26/2019 by Derrick Fouts
	print "Removing phage region $key from further consideration - a XerC/D or integron HMM!\n" if ($DEBUG);
	&write_log("1","Removing phage region $key from further consideration - a XerC/D or integron HMM");
	delete($phageref->{$key});
	$regions--; # added 01/18/2023 dfouts
	#$m = $key - 1;
	#&rekey_regions($m, $phageref, $DEBUG); # call subroutine to renumber phage regions
	next LOOP;
      }
      if (!defined ($phageref->{$key}->{'type'})) { # do I really need to check for this since it should not be defined yet?
	if (($phageref->{$key}->{'core_HMM'} > 0) || ($phageref->{$key}->{'core_BLAST'} > 0)) {$phageref->{$key}->{'type'} = "prophage";}
	elsif (($phageref->{$key}->{'lytic_HMM'} > 0) && ($phageref->{$key}->{'tail_HMM'} > 0) && ($phageref->{$key}->{'int_HMM'} == 0)) {
	  $phageref->{$key}->{'type'} = "bacteriocin";
	}
	# This needs rethought - pulling up plasmids with no core phage genes as a prophage (not good)
	#elsif ((($phageref->{$key}->{'lytic_HMM'} > 0) || ($phageref->{$key}->{'tail_HMM'} > 0)) && ($phageref->{$key}->{'int_HMM'} > 0)) { 
	#    $phageref->{$key}->{'type'} = "prophage"; 
	#}
	print "REGION <$key> is type <$phageref->{$key}->{'type'}>\n" if ($DEBUG);
      }
      if ($phageref->{$key}->{'type'} ne "bacteriocin")  { # changed 12/01/04 when P1-like phage with >=noise core HMM and att site
	$phage_5prime = $phageref->{$key}->{'5prime'};
	$phage_3prime = $phageref->{$key}->{'3prime'};
	$size = ($phage_3prime - $phage_5prime) + 1;
	if ($phageref->{$key}->{'Mu_HMM'} > 0) { $phageref->{$key}->{'class'} = "Mu-like"; }
	elsif (($phageref->{$key}->{'PF04606'}) && ($phageref->{$key}->{'TIGR01613'}) && ($phageref->{$key}->{'int_HMM'} > 0)) { # if has Ogr_Delta & primase & integrase
	  if ($phageref->{$key}->{'PF00078'}) { $phageref->{$key}->{'class'} = "Retron"; } # if has reverse transcriptase
	  elsif ($size > 25000) { $phageref->{$key}->{'class'} = "P2-like"; }
	  elsif ($size < 15000) { $phageref->{$key}->{'class'} = "P4-like"; }
	}
      }
}
  return($regions); # added 01/18/2023 to fix a problem with empty PFPR_tab.txt files when a region was deleted in this subroutine
}
1;
