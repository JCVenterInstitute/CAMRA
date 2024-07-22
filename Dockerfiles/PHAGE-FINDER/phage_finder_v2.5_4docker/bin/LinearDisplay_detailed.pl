#!/usr/bin/env perl

use strict;
use FindBin qw($Bin); # Find where the executable is.  Added per Phil Goetz's suggestion
use lib "$Bin/../lib"; #phage library location
#use DBI;
use TIGR::Foundation;
use Math::Round qw(:all);
use Graphics::ColorObject;

require "$Bin/my4complete.ph";
require "$Bin/functions4complete.pl";
#require "/Users/dfouts/bin/LinearDisplay_dir/define_blocks_flatfile.pl";
 
############################################################
# GLOBALS these structures will be considered to be available 
#  to all functions and points in this program.
############################################################

use vars qw($global_struct $rgb $color_lookup $tfobject %gene %chi $color_start $total_colors $draw_ident_key $color_match);

$tfobject = new TIGR::Foundation;
$global_struct = load_global_struct();
($rgb, $color_lookup) = load_color_struct(); 
$color_start = 72; # start value for spectrum
$total_colors = "";
handle_invocation(0);

my $MY_HELPTEXT = qq~
 
 example:
 region_display -a 6330 -L > thing.fig

 -D [D]ebug

 -i display the last [i]nvocation to the user.

 -A fr[A]gment file

 -F gene attribute [F]ile 

 -P Gene [P]air File. A list of gene pairs.

 -R [R]egional matches
 
 -C [C]hi-square file

 -I [I]ndependent genes. Do not display blocks, just display gene pairs individually.

 -g supply a list of [g]ene names to display

 -p [p]age test. Just displays lines only so user can quickly set up white space. 

 -s size. Line [s]pacing, where size is distance from each assembly displayed in pic. 

 -f size. [f]ont size, where size will be size of font.

 -L use [L]ocus names from ident

 -y use gene_s[y]ms from ident, where possible

 -S size. [S]hort genes. Where genes below size cutoff will displayed as empty.
 
 -T [T]oggle to shorten names from BB0123 to 123

 -M file_name. [M]erge in xfig file, usually a legend created by the user.

 -m [m]olecular size

 -t [t]oggle to use tRNA name
 
 -nl [n]o horizontal [l]ine

 -nn [n]o [n]ame.  Do not print a feat or locus name

 -nt [n]o [t]ickmarks

 -oa [o]ld [a]rrowhead stype

 -cm {option} [c]olor of [m]atches {[p]astel or [s]pectrum or [m]onochrome (eg. -cm p) [default = p]
 
 Fragment file format:
 Level	Fragment_id	Left_end	Right_end	Left_margin
 
 Gene attribute file format:
 End5	End3	feat_name	Gene_sym	Locus	feat_type	Role	Fragment_id
 
   -> Acceptable feat_types to draw are:

      CLS - cluster (used to draw clusters as circles)
      ORF - Open Reading Frame with new arrowhead -> or <-
    ORFNA - Open Reading Frame with NO arrowhead -
     PROM - promoter
     PRI  - primer
     PRI2 - primer2
      ENZ - restriction enzyme
     ATTP - phage attachment site (attP)
     TERM - transcriptional terminator
      TSD - Target Site Duplication
    PMARK - contig/assembly separator
   REGION - for noting DNA similarities (shaded region)
     rRNA - ribosomal RNAs (must specify type 16S, 23S, 5S in feat_name)
       16S in feat_name
       23S in feat_name
        5S in feat_name
     sRNA - structural RNAs (must specify type (tmRNA or 7S/SRP or RNaseP/rnpB in feat_name)
       tmRNA in feat_name
       7S/SRP in feat_name
       RNaseP/rnpB in feat_name
     tRNA - transfer RNAs

 Gene pair file format:
 Fragement_id	feat_name       Fragment_id	feat_name     %identity

 Regional match file format:
 Fragment_id    end5     end3     Fragment_id     end5      end3     %identity
 
 Chi-square file format:
 # Fragment_id
 position	chi-square
 position	chi-square
~;

$tfobject->setHelpInfo($MY_HELPTEXT);

############################################################
#               process_args
############################################################

my($fragmt_list, $pairs_list) = &process_args($tfobject);

&read_genes();
&read_chi_square();

############################################################
#                         main
############################################################

if ($color_match eq "p") {$total_colors = 10; $color_start = $color_start - 7;}
elsif ($color_match eq "s") {$total_colors = 255;}

&fig_header() if ($fragmt_list->{'count'} > 0);

for(my $i=0;$i<$fragmt_list->{'count'};$i++) {
    printf STDERR ("DRAWING: %s\n", $fragmt_list->{$i}->{'fragmt_id'});    
    &handle_level($fragmt_list->{$i});
}

if ($pairs_list->{'count'} != 0) {
    print STDERR "drawing matches between levels...\n";
    &handle_links ();
    print STDERR "printing the legend...\n";
    my $line_count = $fragmt_list->{'count'}+1;
    #my $line_count = $fragmt_list->{$count+1}->{'level'};
    my $max_x = &x_transform(100, &create_level_struct($fragmt_list->{$line_count}));
    #my $max_x = 1;
    #my $max_y = line2coord($line_count);
    my $max_y = line2coord($fragmt_list->{$line_count}->{'level'});
    my $cell_width = $global_struct->{'gene_height'};
    my $x_angle = 0.0000; #1.5708
    my $y_angle = 0.0000;
    my $cell_height = $global_struct->{'gene_height'};
    my $x_pad = 1000;
    my $y_pad = 1000;
    if (($draw_ident_key) && ($color_match ne "m"))  { #if we need to draw the key
        my $n = "";
        if ($color_match eq "s")  { $n = 100;} # for full spectrum (the old way), draw 100
        else { $n = 10;} # for pastel colors, only draw 10
	for(my $i = 1; $i <= $n; $i++) { # change $n here to 9 if you want to print to 100..20%
	    my $ratio = $i/$n;
	    my $color = &position_to_color($ratio, $total_colors)-1;
	    if ($i == 1) {
		print &text(20,1,$y_angle,$max_x + ($cell_width * ($i+0.5)/2),$max_y + 300,"100");
		print &text(20,2,$y_angle,$max_x + ($cell_width * ($i+0.6)/2),$max_y + 600,"|");
	    }
	    elsif ($i == $n) { # change $n here to 9 if you want to stop at 20%
		print &text(20,1,$y_angle,$max_x + ($cell_width * ($i+0.5)/2),$max_y + 300,"10");
		print &text(20,2,$y_angle,$max_x + ($cell_width * ($i+0.6)/2),$max_y + 600,"|");
	    }
	    print &box($global_struct->{'black'},$color,$max_x + ($i/2 * $cell_width),$max_y + 600,$cell_width/2,$cell_height,1,1);
	}
	print &text(20,0,$y_angle,$max_x + ($cell_width * ($n+1)/2),$max_y + 300,"% identity");
    }
}

merge_with_legend($global_struct->{'legend_file'}) if ($global_struct->{'merge_with_legend'});

exit 0;

############################################################
#                      end of main
############################################################


sub handle_links {
    for (my $i=0; $i<=$pairs_list->{'count'}; $i++) {    	
	&draw_gene_pairs($pairs_list->{$i}->{'fragmt_id1'},
			 $pairs_list->{$i}->{'fragmt_id2'},
			 $pairs_list->{$i}->{'name1'},
			 $pairs_list->{$i}->{'name2'},
			 $pairs_list->{$i}->{'percent_id'},
			 $pairs_list->{$i}->{'type'},
                         $pairs_list->{$i}->{'color'});
    }
}

sub handle_level {
    my($fragmt_ref) = shift;
    my($level_struct);
    my $level_rend;

    $level_struct = &create_level_struct($fragmt_ref);
    $level_rend = &draw_level($level_struct, $fragmt_ref);
    &draw_ticks($level_struct, $fragmt_ref) if ($global_struct->{'tickmark'});
    &draw_gene_elements($level_struct, $fragmt_ref);
    &draw_chi_square($level_struct, $fragmt_ref);
    $global_struct->{$fragmt_ref->{'level'}} = $level_rend;
}

sub create_level_struct {
    my($fragmt_ref) = shift;
    my(%level_struct);
      
    $level_struct{'level'} = $fragmt_ref->{'level'};
    if ($fragmt_ref->{'lend'} < $fragmt_ref->{'rend'}) {
    	        $level_struct{'x_decrement'} = $fragmt_ref->{'lend'};
		$level_struct{'orientation'} = $global_struct->{'forward'};
		$level_struct{'width'} = $fragmt_ref->{'rend'} - $fragmt_ref->{'lend'};
    }
    else {
		$level_struct{'x_decrement'} = $fragmt_ref->{'rend'};
		$level_struct{'orientation'} = $global_struct->{'reverse'};
		$level_struct{'width'} = $fragmt_ref->{'lend'} - $fragmt_ref->{'rend'};
    }
	
	$level_struct{'lmargin'} = $fragmt_ref->{'lmargin'};
	
    return(\%level_struct);
}

sub merge_with_legend {
    my($f) = shift;
    my($count);

    # very simple minded approach. this function just pulls up the file and 
    #  adds it into the end of the stdout stream. the one thing it does is 
    #  strip out all new color calls, and operates under the ASSUMPTION 
    #  that the color definitions between the legend and the big pic are
    #  the same. the problem is that xfig doesnt like merging two files 
    #  have color definitions in them, and will just create more color
    #  definitions.

    print STDERR "MERGING: $f\n";

    open(LEGEND, $f) || die "hey, file: $f doesnt exist";
    $count = 0;
    while(<LEGEND>) {
	if(/0\s.*\#/) {
	    # dont print these lines
	}
	elsif($count < 9) {
	    # dont print these lines, these are the header. 
	    # a STUPID way to do this, but...
	}
	else {
	    print;
	}
	$count++;
    }
    close(LEGEND);

}


sub process_args {
    my %fragmt_list;
    my %pairs_list;
    my ($opt_A, $opt_F, $opt_I, $opt_D, $opt_x, $opt_P, $opt_C);
    my ($opt_G, $opt_L, $opt_M, $opt_S, $opt_T, $opt_f, $opt_g);
    my ($opt_l, $opt_m, $opt_o, $opt_p, $opt_s, $opt_t, $opt_y);
    my ($opt_nl, $opt_nn, $opt_nt, $opt_oa, $opt_cm, $opt_R);

    $tfobject->TIGR_GetOptions('-A:s' , \$opt_A,
    			       '-F:s' , \$opt_F,
			       '-D'   , \$opt_D,
			       '-G'   , \$opt_G,
			       '-L'   , \$opt_L,
			       '-M:s' , \$opt_M,
			       '-P:s' , \$opt_P,
			       '-C:s' , \$opt_C,
			       '-S:s' , \$opt_S,
			       '-I'   , \$opt_I,
			       '-T'   , \$opt_T,
			       '-f:s' , \$opt_f,
			       '-g:s' , \$opt_g,
			       '-l:s' , \$opt_l,
			       '-m:s' , \$opt_m,
			       '-o:s' , \$opt_o,
			       '-p'   , \$opt_p,
			       '-s:s' , \$opt_s,
			       '-t:s' , \$opt_t,
			       '-x'   , \$opt_x,
			       '-y'   , \$opt_y,
			       '-nl'  , \$opt_nl,
                               '-nn'  , \$opt_nn,
			       '-nt'  , \$opt_nt,
			       '-oa'  , \$opt_oa,
                               '-cm:s', \$opt_cm,
                               '-R:s' , \$opt_R);


    if (length($opt_I) > 0) {
		$global_struct->{'block_toggle'} = 0;
    }

    if (length($opt_D) > 0) {
		$global_struct->{'debug'} = 1;
    }

    if (length($opt_x) > 0) {
		$global_struct->{'xml_toggle'} = 1;
    }

    # User just wants to draw a figure legend.
    if (length($opt_G) > 0) {
		&fig_header();
		make_a_legend();
		exit(1);
    }

    if(length($opt_g) > 0) {
	# this is a way of displaying only a selected set of genes.
		open(F_IN,"$opt_g") || usage("oof. cant open $opt_g");
		while(<F_IN>) {
	    	s/\n//;
	    	$global_struct->{'selected_genes'}->{$_} = 1;
		}
		close(F_IN);
		$global_struct->{'select_genes_toggle'} = 1;
    }

    if ($opt_cm)  {
        $color_match = $opt_cm;
    }
    else {
        $color_match = "p"; # default for drawing gene pairs
    }

    if(length($opt_P) > 0) {
	my ($count) = 0;
	
	open(F_IN,"$opt_P") || usage("oof. cant open $opt_P");
	while(<F_IN>) {
	    chomp;
	    my @l = split /\t+/;
            $pairs_list{$count}->{'type'} = "pairs";
	    $pairs_list{$count}->{'fragmt_id1'} = $l[0];
	    $pairs_list{$count}->{'name1'} = $l[1];
	    $pairs_list{$count}->{'fragmt_id2'} = $l[2];
	    $pairs_list{$count}->{'name2'} = $l[3];
	    
	    # Added by Derrick E. Fouts 7/3/2003 to enable different colored gene pairs depending on BLAST % identity
	    if ($#l > 3)  {
		$pairs_list{$count}->{'percent_id'} = round($l[4]);  # record the percent identity
	    }
	    if ($color_match eq "m") { # to set all gene matches to Herve's USER 55 grey color
		$pairs_list{$count}->{'color'} = $global_struct->{'monogene'};
	    }
	    else {
		$pairs_list{$count}->{'color'} = "";
	    }
	    $count++;
	}
	close(F_IN);
	$pairs_list{'count'} = $count;
        $draw_ident_key = 1;
    }
    
    if(length($opt_R) > 0) { # to handle regional matches from grasta, blast or mummer data [dfouts 04/11/05]
     
        my $frag1 = "";
	my $frag2 = "";
        my $name1 = "";
	my $name2 = "";
        my $count = "";
        my $keycount = keys %pairs_list; # check to see if we also have gene pairs

        if ($keycount > 0)  {
	    $count = $keycount + 1; # start numbering after the last pair
	}
        else {
            $count = 0;
        }
	open(F_IN,"$opt_R") || usage("oof. cant open $opt_R");
	while(<F_IN>) {
	    chomp;
	    my @l = split /\t+/;
            $frag1 = $l[0];
            $name1 = "region\_$frag1\_$count";
	    $pairs_list{$count}->{'fragmt_id1'} = $frag1;
	    $pairs_list{$count}->{'name1'} = $name1;
            $gene{$frag1}{$name1}{'end5'} = $l[1];
	    $gene{$frag1}{$name1}{'type'} = "REGION";
            $gene{$frag1}{$name1}{'end3'} = $l[2];

	    $frag2 = $l[3];
            
            $name2 = "region\_$frag2\_$count";
	    $pairs_list{$count}->{'fragmt_id2'} = $frag2;
	    $pairs_list{$count}->{'name2'} = $name2;
	    $gene{$frag2}{$name2}{'end5'} = $l[4];
            $gene{$frag2}{$name2}{'type'} = "REGION";
            $gene{$frag2}{$name2}{'end3'} = $l[5];
	    
	    if ($#l > 5)  {
		$pairs_list{$count}->{'percent_id'} = round($l[6]);  # record the percent identity
	    }
            if ((length($opt_P) > 0) || ($color_match eq "m")) { # if either gene pairs or monochrome make cornsilk
		$pairs_list{$count}->{'color'} = $global_struct->{'regions'};
	    }
	    else {
		$pairs_list{$count}->{'color'} = "";
	    }
	    $count++;
	}
	close(F_IN);
	
	$pairs_list{'count'} = $count;
	if ($color_match eq "m") {
	    $draw_ident_key = 0;
        }
        else {
            $draw_ident_key = 1;
        }
    }

    if(length($opt_A) > 0) {
	my ($count) = 0;
	my $maxlevel = 0;
	
	open(F_IN,"$opt_A") || usage("oof. cant open $opt_A");
	while(<F_IN>) {
	    chomp;
	    my @l = split /\t+/;
	    $fragmt_list{$count}->{'level'} = $l[0]; 	    
	    $fragmt_list{$count}->{'fragmt_id'} = $l[1];
	    $fragmt_list{$count}->{'lend'} = $l[2];
	    $fragmt_list{$count}->{'rend'} = $l[3];
	    if ($fragmt_list{$count}->{'level'} > $maxlevel)  {  # dfouts
		$maxlevel = $fragmt_list{$count}->{'level'}
	    }
	    if ($l[4]) {
		$fragmt_list{$count}->{'lmargin'} = $l[4];
	    }
	    else {
		$fragmt_list{$count}->{'lmargin'} = 0;
	    } 
	    $count++;
	}
	close(F_IN);
	
	$fragmt_list{'count'} = $count;
	# set structure for the pairs color key
	$fragmt_list{$count+1}->{'level'} = $maxlevel+1; 	    
	$fragmt_list{$count+1}->{'fragmt_id'} = "legend";
	$fragmt_list{$count+1}->{'lend'} = 1;
	$fragmt_list{$count+1}->{'rend'} = 1000;
	$fragmt_list{$count+1}->{'lmargin'} = 0;
    }
    else {
		usage("Hey: must supply -A an assembly file");
    }
    if ($opt_F) {
    	$global_struct->{'gene_file'} = $opt_F;
    }
    else {
    	usage("Hey: must supply -F gene attribute file");
    }
    if ($opt_C) {
    	$global_struct->{'chi_file'} = $opt_C;
    }
    if ($opt_L and $opt_l) {
	usage("cant specify -l (lookup table) and -L (use ident table loci)");
    }
   if ($opt_l) {
       # this is a way of displaying new names for genes
       $global_struct->{'gene_lookup_toggle'} = 1;
       open(F_IN,"$opt_l") || usage("oof. cant open $opt_l");
       while(<F_IN>) {
	   chop;
	   my @l = split(/\s/);
	   $global_struct->{'lookup_table'}->{$l[0]} = $l[1];
       }
       close(F_IN);
    }
    if ($opt_L) {
	# this gets gene names from the ident.locus field
	$global_struct->{'locus_toggle'} = 1;
    }
    if ($opt_y) {
	# this gets gene_syms from ident
	$global_struct->{'gene_sym_toggle'} = 1;
    }
    if ($opt_p) {
	$global_struct->{'page_test_toggle'} = 1;
    }
    if ($opt_t) {
	$global_struct->{'trna_name_toggle'} = 1;
    }
    if ($opt_T) {
	$global_struct->{'truncate_name_toggle'} = 1;
    }
    if ($opt_S) {
	$global_struct->{'short_gene_toggle'} = 1;
	$global_struct->{'short_gene_cutoff'} = $opt_S;
    }
    if ($opt_s) {
	$global_struct->{'line_spacing'} = $opt_s;
    }
    if ($opt_f) {
	$global_struct->{'point_size'} = $opt_f;
	
	$global_struct->{'label_height'} 
	= ($global_struct->{'label_height'} * $global_struct->{'point_size'}) / $global_struct->{'point_size'};

	$global_struct->{'char_width'}  #huh?
	= ($global_struct->{'char_width'} * $global_struct->{'point_size'}) / $global_struct->{'point_size'};

    }
    $global_struct->{'object_spacing'} = $opt_o;
    $global_struct->{'size_supplied_by_user'} = 0;

    if ($opt_m) {
	$global_struct->{'molecule_size'} = $opt_m;
	$global_struct->{'size_supplied_by_user'} = 1;
    }

    $global_struct->{'merge_with_legend'} = 0;

    if ($opt_M) {
	$global_struct->{'legend_file'} = $opt_M;
	$global_struct->{'merge_with_legend'} = 1;
    }

    $global_struct->{'horizontal line'} = 1;

    if ($opt_nl) {
    	$global_struct->{'horizontal line'} = 0;
    }
    
    $global_struct->{'tickmark'} = 1;

    if ($opt_nn) {
	$global_struct->{'no_names'} = 1;
	print STDERR "No names selected\n";
    }
    else {
	$global_struct->{'no_names'} = 0;
    }

    if ($opt_nt) {
    	$global_struct->{'tickmark'} = 0;
    }

    $global_struct->{'arrowstyle'} = 1; # 1 is the new style [DEFAULT]

    if ($opt_oa) {
	$global_struct->{'arrowstyle'} = 0;
    }
    return(\%fragmt_list,\%pairs_list);
}


sub usage {
    my $str = shift;

    print STDERR "$str\n" if (length($str) > 0);

    $tfobject->printHelpInfoAndExit();
    
    exit(1);
}


