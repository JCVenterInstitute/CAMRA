package TIGR::DisplayFunctions;

=pod

=head1 NAME

LinearDisplay.pl - A PERL script that generates publication quality linear genetic (e.g., open reading frame (ORF) maps) in XFig format.

=head1 LICENSE

Copy (C) 2019-2020  The J. Craig Venter Institute (JCVI).  All rights reserved
Written by various people through the years including, but not limited to,
Derrick E. Fouts, Martin Wu, John Quackenbush, Jeremy Peterson, Owen White and members of his team at TIGR.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 CITATION

The following publication used this program to generate linear illustrations of bacteriophage genomes and should be used to site this program:

Duan Y, Llorente C, Brandl K, Chu H, Jian L, Lang S, White RC, Clarke TH, Nguyen K, Torralba M, Sho Y, Liu J, Hernandez-Morales A, Lessor L, Rahman IR, Miyamoto Y, Ly M, Sun W, Kiesel R, Hutmacher F, Lee S, Ventura-Cots M, Bosques-Padilla F, Verna EC, Abraldes JG, Brown Jr RS, Vargas V, Altamirano J, Caballeria J, Shawcross DL, Ho SB, Louvet A, Lucey MR, Mathurin P, Garcia-Tsao G, Bataller R, Tu XM, Eckmann L, van der Donk WA, Young R, Lawley TD, Starkel P, Pride D, Fouts DE, Schnabl B. Bacteriophages against cytolytic Enterococcus faecalis as therapy for alcoholic liver disease. Nature. in press

=head1 SYNOPSIS

  use TIGR::DisplayFunctions;

=head1 DESCRIPTON

This module provides essential subroutines for drawing linear genetic maps.

=head1 CONTACT

    Derick E. Fouts
    dfouts@jcvi.org

=cut

   BEGIN {
      require 5.006_00; # error if using Perl < v5.6.0
   }

use strict;
use Data::Dumper;
use XML::Simple;
require Exporter;

   our @ISA = ('Exporter');
   our @EXPORT = qw(
                   draw_gene_pairs
                   draw_level
                   draw_ticks
                   read_genes
                   read_chi_square
                   draw_gene_elements
                   draw_chi_square
                   position_to_color
                   box
                   text
                   handle_invocation
                   fig_header
                   line2coord
                   x_transform
                   );

    our @EXPORT_OK = qw(
                       draw_gene_pairs
                       draw_level
                       draw_ticks
                       read_genes
                       read_chi_square
                       draw_gene_elements
                       draw_chi_square
                       position_to_color
                       box
                       text
                       handle_invocation
                       fig_header
                       line2coord
                       x_transform
                       );

   our %EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

sub draw_gene_pairs {
  my ($a1, $a2, $name1, $name2, $percent_id, $type, $color_code, $color_match, $color_start, $total_colors, $gene_ref, $global_struct) = @_;
  my $inc = 0;
  my ($x1, $x2, $x3, $y1, $y2, $depth);
  my $max_colors = $total_colors; #default
  return unless (exists $gene_ref->{$a1}{$name1}{'x1'} && exists $gene_ref->{$a2}{$name2}{'x1'});

  # blow thru those structures and fill cells in a 2-dimensional
  #  "dot plot" where each cell that is a match between genes gets
  #  filled with a 1.

  if ((!$color_code) && ($color_match ne "m"))  {
    my $ratio = ($percent_id-1)/100;
    $color_code = &position_to_color(1-$ratio,$total_colors,$color_start); # flip the ratio so that 100% match is 1 and 20% is 0
  }
  if ($type eq "pairs") {
    $depth = $global_struct->{'block_depth'} + $inc++;
    printf("2 3 0 1 %d %d %d 0 20 0.000 0 0 0 0 0 5\n",
	   $color_code,
	   $color_code,
	   $depth);
    if ($gene_ref->{$a2}{$name2}{'y2'} < $gene_ref->{$a1}{$name1}{'y2'}) {
      printf("\t %d %d %d %d %d %d %d %d %d %d\n",
	     $gene_ref->{$a1}{$name1}{'x1'},$gene_ref->{$a1}{$name1}{'y2'},
	     $gene_ref->{$a1}{$name1}{'x2'},$gene_ref->{$a1}{$name1}{'y2'},
	     $gene_ref->{$a2}{$name2}{'x2'},$gene_ref->{$a2}{$name2}{'y3'},
	     $gene_ref->{$a2}{$name2}{'x1'},$gene_ref->{$a2}{$name2}{'y3'},
	     $gene_ref->{$a1}{$name1}{'x1'},$gene_ref->{$a1}{$name1}{'y2'});
    }
    else {
      printf("\t %d %d %d %d %d %d %d %d %d %d\n",
	     $gene_ref->{$a1}{$name1}{'x1'},$gene_ref->{$a1}{$name1}{'y3'},
	     $gene_ref->{$a1}{$name1}{'x2'},$gene_ref->{$a1}{$name1}{'y3'},
	     $gene_ref->{$a2}{$name2}{'x2'},$gene_ref->{$a2}{$name2}{'y2'},
	     $gene_ref->{$a2}{$name2}{'x1'},$gene_ref->{$a2}{$name2}{'y2'},
	     $gene_ref->{$a1}{$name1}{'x1'},$gene_ref->{$a1}{$name1}{'y3'});
    }
  }
  else {
    $depth = 999;
    printf("2 3 0 1 %d %d %d 0 20 0.000 0 0 0 0 0 9\n",
	   $color_code,
	   $color_code,
	   $depth);
    if ($gene_ref->{$a2}{$name2}{'y2'} < $gene_ref->{$a1}{$name1}{'y2'}) {# changed, NOT tested dfouts
      printf("\t %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
	     $gene_ref->{$a1}{$name1}{'x1'},$gene_ref->{$a1}{$name1}{'y1'},
	     $gene_ref->{$a1}{$name1}{'x2'},$gene_ref->{$a1}{$name1}{'y1'},
	     $gene_ref->{$a1}{$name1}{'x2'},$gene_ref->{$a1}{$name1}{'y2'},
	     $gene_ref->{$a2}{$name2}{'x2'},$gene_ref->{$a2}{$name2}{'y3'},
	     $gene_ref->{$a2}{$name2}{'x2'},$gene_ref->{$a2}{$name2}{'y1'},
	     $gene_ref->{$a2}{$name2}{'x1'},$gene_ref->{$a2}{$name2}{'y1'},
	     $gene_ref->{$a2}{$name2}{'x1'},$gene_ref->{$a2}{$name2}{'y3'},
	     $gene_ref->{$a1}{$name1}{'x1'},$gene_ref->{$a1}{$name1}{'y2'},
	     $gene_ref->{$a1}{$name1}{'x1'},$gene_ref->{$a1}{$name1}{'y1'});
    }
    else { # this one works [dfouts]
      printf("\t %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
	     $gene_ref->{$a1}{$name1}{'x1'},$gene_ref->{$a1}{$name1}{'y1'},
	     $gene_ref->{$a1}{$name1}{'x2'},$gene_ref->{$a1}{$name1}{'y1'},
	     $gene_ref->{$a1}{$name1}{'x2'},$gene_ref->{$a1}{$name1}{'y3'},
	     $gene_ref->{$a2}{$name2}{'x2'},$gene_ref->{$a2}{$name2}{'y2'},
	     $gene_ref->{$a2}{$name2}{'x2'},$gene_ref->{$a2}{$name2}{'y1'},
	     $gene_ref->{$a2}{$name2}{'x1'},$gene_ref->{$a2}{$name2}{'y1'},
	     $gene_ref->{$a2}{$name2}{'x1'},$gene_ref->{$a2}{$name2}{'y2'},
	     $gene_ref->{$a1}{$name1}{'x1'},$gene_ref->{$a1}{$name1}{'y3'},
	     $gene_ref->{$a1}{$name1}{'x1'},$gene_ref->{$a1}{$name1}{'y1'});
    }
  }
}

sub fake_genes { # determine the proper level structure and coords for regions (converts end5, end3 to fig units)
  my ($global_struct, $level_ref, $gene_ref) = @_;
  my ($x1, $x2, $y1, $y2, $y3); # draw_gene_pairs only needs x1 and x2, y2 and y3

  $x1 = &x_transform($global_struct, $gene_ref->{'end5'}, $level_ref);
  $x2 = &x_transform($global_struct, $gene_ref->{'end3'}, $level_ref);

  #    $y2------------+
  #      |            |
  #      |            |
  #    $y1            |
  #      |            |
  #      |            |
  #    $y3------------+
  #
  #    $x2           $x1
  #$end5              $end3

  $y1 = &line2coord($global_struct, $level_ref->{'level'});
  $y2 = $y1 - ($global_struct->{'gene_height'} / 2);
  $y3 = $y1 + ($global_struct->{'gene_height'} / 2);
  $gene_ref->{'x1'} = $x1;
  $gene_ref->{'x2'} = $x2;
  $gene_ref->{'y1'} = $y1;
  $gene_ref->{'y2'} = $y2;
  $gene_ref->{'y3'} = $y3;
}

sub draw_level {
  my ($global_struct, $level_ref, $fragmt_ref) = @_;
  my ($x1, $x2, $y);

  $x1 = &x_transform($global_struct, $fragmt_ref->{'lend'}, $level_ref);
  $x2 = &x_transform($global_struct, $fragmt_ref->{'rend'}, $level_ref);

  return $x2 if $global_struct->{'horizontal line'} == 0;

  if ($global_struct->{'tickmark'} == 0) {
    $y = &line2coord($global_struct, $fragmt_ref->{'level'}) + ($global_struct->{'gene_arrow_height'} / 2);
  }
  else {
    $y = &line2coord($global_struct, $fragmt_ref->{'level'});
  }

  # print name of level if the first level of a whole genome or unique level
  &add_label($fragmt_ref->{'fragmt_id'} . " ",$x1,$x2,$y+180,1,2,2,40);

  # print polyline (line for each level)
  printf("2 1 0 1 %d %d %d 0 20 0.000 0 0 0 0 0 2\n",
	 $global_struct->{'black'},
	 $global_struct->{'black'},
	 $global_struct->{'row_depth'});
  #printf("\t %d %d %d %d %d %d %d %d %d %d %d %d\n",
  printf("\t %d %d %d %d\n",
	 $x1,$y,
	 $x2,$y);

  return $x2;
}

sub draw_ticks { # dfouts
  my ($global_struct, $level_ref, $fragmt_ref) = @_;
  my $color = $global_struct->{'black'};
  my $pen = $global_struct->{'black'};
  my $depth = $global_struct->{'gene_depth'}+1; # set depth of tick (larger # = more behind)
  my $filled = 20; # please do fill the box
  my($x1, $y1, $y2);
  my $tick = $global_struct->{'tick'};
  my $major_tick = $global_struct->{'major_tick'};
  my $minor_tick = $global_struct->{'minor_tick'};
  my $position = "";
  my $name = "";
  my $next_major_tick = "";
  my $next_minor_tick = "";
  my $frag_begin = "";
  my $frag_end = "";

  $y1 = &line2coord($global_struct, $level_ref->{'level'})+7;

  #$fragmt_ref->{lend},$fragmt_ref->{rend}

  # F O R  A  T I C K  M A R K
  #
  #    +-----+ $y2
  #    |     |
  #    |     |
  #    |     |
  #    |     |
  #    |     |
  #    |     |
  #    +-----+ $y1
  #      $x1
  #
  #

  $next_major_tick = $position + $major_tick;
  $next_minor_tick = $position + $minor_tick;

  if ($level_ref->{'orientation'} == 1) {
    $frag_begin = $fragmt_ref->{lend};
    $frag_end = $fragmt_ref->{rend};
    if ($fragmt_ref->{lend} == 1)  {
      $position = 0;
    }
    else  {
      $position = $fragmt_ref->{lend};
    }
  }
  else {
    $frag_begin = $fragmt_ref->{rend};
    $frag_end = $fragmt_ref->{lend};
    if ($fragmt_ref->{rend} == 1)  {
      $position = 0;
    }
    else  {
      $position = $fragmt_ref->{rend};
    }
  }

  for ($position; $position <= $frag_end; $position = $position + $tick)  {
    $x1 = &x_transform($global_struct, $position, $level_ref);
    #print STDERR "ticks:  $fragmt_ref->{lend}\t$position\t$tick\t$frag_end\t$fragmt_ref->{rend}\n";
    if ($position == 0)  { # major tick at beginning
      $y2 = $y1 - ($global_struct->{'gene_height'} * 2);
      printf("2 1 0 3 %d %d %d 0 %d 0.000 0 0 0 0 0 2\n",
	     $pen, 
	     $color,
	     $depth,
	     $filled);

      printf("\t %d %d %d %d\n",
	     $x1,$y1,
	     $x1,$y2);
      $name = $position . " bp";
      &add_label($name,$x1-100,$x1-100,$y2-150,1,0,0,30);
      # name; x1; x2; y; position (vert=2, horz = 1, else diagonal); justification (left = 0, center = 1, right = 2); font (Times Roman = 0, Times Italix = 1, Times Bold = 2); font size; color; mid, rotation #
    }
    elsif ($position == $frag_begin)  { # when numbering wraps around
      $y2 = $y1 - ($global_struct->{'gene_height'} * 2);
      printf("2 1 0 3 %d %d %d 0 %d 0.000 0 0 0 0 0 2\n",
	     $pen,
	     $color,
	     $depth,
	     $filled);

      printf("\t %d %d %d %d\n",
	     $x1,$y1,
	     $x1,$y2);
      $name = $position;
      $position = nlowmult($tick, $name);  # get the next position into multiples of $tick
      $next_major_tick = $position + $major_tick;
      $next_minor_tick = $position + $minor_tick;
      &add_label($name,$x1-100,$x1-100,$y2-150,1,0,0,30);
    }
    if (($position == $next_major_tick) && (($frag_end - $position) > $minor_tick)) { # print a major tick
      $y2 = $y1 - ($global_struct->{'gene_height'} * 2);
      printf("2 1 0 3 %d %d %d 0 %d 0.000 0 0 0 0 0 2\n",
	     $pen,
	     $color,
	     $depth,
	     $filled);

      printf("\t %d %d %d %d\n",
	     $x1,$y1,
	     $x1,$y2);
      $name = $position;
      &add_label($name,$x1-100,$x1-100,$y2-150,1,0,0,30);
      $next_major_tick = $position + $major_tick;
    }
    if (($position == $next_minor_tick) && (($frag_end - $position) > $minor_tick)) { # print a minor tick
      $y2 = $y1 - ($global_struct->{'gene_height'});
      printf("2 1 0 3 %d %d %d 0 %d 0.000 0 0 0 0 0 2\n",
	     $pen, 
	     $color,
	     $depth,
	     $filled);

      printf("\t %d %d %d %d\n",
	     $x1,$y1,
	     $x1,$y2);
      $next_minor_tick = $position + $minor_tick;
    }
    else { # print generic small tick
      $y2 = $y1 - ($global_struct->{'gene_height'}/8);
      printf("2 1 0 3 %d %d %d 0 %d 0.000 0 0 0 0 0 2\n",
	     $pen, 
	     $color,
	     $depth,
	     $filled);

      printf("\t %d %d %d %d\n",
	     $x1,$y1,
	     $x1,$y2);
    }
  }
  # now print the last tick
  $position = $frag_end;
  $x1 = &x_transform($global_struct, $position, $level_ref);
  $y2 = $y1 - ($global_struct->{'gene_height'} * 2);
  printf("2 1 0 3 %d %d %d 0 %d 0.000 0 0 0 0 0 2\n",
	 $pen,
	 $color,
	 $depth,
	 $filled);

  printf("\t %d %d %d %d\n",
	 $x1,$y1,
	 $x1,$y2);
  $name = $position . " bp";
  &add_label($name,$x1-100,$x1-100,$y2-150,1,0,0,30);
  # name; x1; x2; y; position (vert=2, horz = 1, else diagonal); justification (left = 0, center = 1, right = 2); font (Times Roman = 0, Times Italix = 1, Times Bold = 2); font size; color; mid, rotation #
}

sub read_genes {
  my ($global_struct, $gene_ref) = @_;
  my $init = 1;
  my $v = "";
  open (IN, $global_struct->{'gene_file'}) || die "can't open $global_struct->{'gene_file'}\n";
  while (<IN>) {
    chomp;
    next if (/^#/); # skip line if commented out (so we can annotate gene_att files
    my $line = $_;
    my ($end5, $end3, $feat_name, $gene_sym, $locus, $type, $role, $frag_id) = split /\t/,$_,8 if ($init);
    if ((($type ne "ORF") && ($type ne "ORFNA") && ($type ne "CLS") && ($type ne "INDEL") && ($type ne "Ti") && ($type ne "Tv") && ($type ne "PMARK") && ($type ne "TERM") && ($type ne "PROM") && ($type ne "ATTP") && ($type ne "TSD") && ($type !~ /RNA/) && ($type ne "RBS")) && ($init))  {
      $v = 1;
    }
    elsif ((($type eq "ORF") || ($type eq "ORFNA") || ($type eq "CLS") || ($type eq "INDEL") || ($type eq "Ti") || ($type eq "Tv") || ($type eq "PMARK") || ($type eq "TERM") || ($type eq "PROM") || ($type eq "ATTP") || ($type eq "TSD") || ($type =~ /RNA/) || ($type eq "RBS")) && ($init)) {
      $v = 2;

    }
    $init = "";
    if ($v == 2) {
      ($end5, $end3, $feat_name, $gene_sym, $locus, $type, $role, $frag_id) = split /\t/,$line,8;
      #print STDERR "--> 1\n";
    }
    else  {
      ($end5, $end3, $feat_name, $gene_sym, $locus, $role, $frag_id, $type) = split /\t/,$line,8;
      #print STDERR "--> 2\n";
    }

    $gene_ref->{$frag_id}{$feat_name}->{'end5'} = $end5;
    $gene_ref->{$frag_id}{$feat_name}->{'end3'} = $end3;
    $gene_ref->{$frag_id}{$feat_name}{'name'} = $feat_name;
    $gene_ref->{$frag_id}{$feat_name}{'main_role'} = $role;
    $gene_ref->{$frag_id}{$feat_name}{'type'} = $type;
    $gene_ref->{$frag_id}{$feat_name}{'final_name'} = &choose_name($global_struct, $type, $feat_name, $gene_sym, $locus);
    #print STDERR "$end5, $end3, $feat_name, $gene_sym, $locus!, $type, $role, $frag_id, $gene_ref->{$frag_id}{$feat_name}{'final_name'}\n";
    #print STDERR "$frag_id, $feat_name, $gene_ref->{$frag_id}{$feat_name}{'end5'}, $gene_ref->{$frag_id}{$feat_name}{'end3'}, $gene_ref->{$frag_id}{$feat_name}{'name'}, $gene_ref->{$frag_id}{$feat_name}{'type'}, $gene_ref->{$frag_id}{$feat_name}{'main_role'}, $gene_ref->{$frag_id}{$feat_name}{'final_name'}\n";
  }
  close(IN);
}

sub read_chi_square {
  my ($global_struct, $chi_ref) = @_;
  my $fragmt_id;
  open (IN, $global_struct->{'chi_file'}) || die "can't open $global_struct->{'chi_file'}\n";
  while (<IN>) {
    chomp;
    if (/^#\s+(\S+.*)/) {
      $fragmt_id = $1;
      $fragmt_id =~ s/\s*$//; # remove any spaces at the end of the line
      print STDERR "FRAGMENT_CHI---> $fragmt_id\n";
    }
    elsif (/(\d+)\s+((\d|\.)+)/) {
      $chi_ref->{$fragmt_id}{$1} = $2;
      #print STDERR "<$fragmt_id>\t<$1>\t<$2>\n";
    }
  }
  close IN;

  for $fragmt_id (keys %{$chi_ref}) {
    my @sort = sort {$a<=>$b} values %{$chi_ref->{$fragmt_id}} ;
    my $range = ($sort[$#sort] - $sort[0]) + 1;
    print STDERR "Range = <$range>\n";

    while (my($key, $value) = each %{$chi_ref->{$fragmt_id}}) {
      $chi_ref->{$fragmt_id}{$key} = ($value-$sort[0])/$range;
    }
  }
}

sub draw_gene_elements {
  my ($level_ref, $fragmt_ref, $global_struct, $gene_ref, $color_lookup) = @_;
  my ($color, $pen, $filled, $end5, $end3);
  my $fragmt_id = $fragmt_ref->{'fragmt_id'};
  my $role;

  foreach my $element (keys %{$gene_ref->{$fragmt_id}}) {
    $end5 = $gene_ref->{$fragmt_id}{$element}{'end5'};
    $end3 = $gene_ref->{$fragmt_id}{$element}{'end3'};

    if (&overlap($fragmt_ref->{lend},$fragmt_ref->{rend},$end5,$end3)) {
      if (!$global_struct->{'select_genes_toggle'} or
	  ($global_struct->{'select_genes_toggle'} and
	   $global_struct->{'selected_genes'}->{$gene_ref->{$fragmt_id}{$element}{'name'}})) {

	$filled = 20;  # filled
	$filled = 41;  # pattern
	$filled = -1;  # no fill

	if (passes_length_test($global_struct, $end5, $end3)) {
	  $role = $gene_ref->{$fragmt_id}{$element}{'main_role'};
	  $color = $color_lookup->{$role};
	  $pen = $global_struct->{'black'};

	  if ($color == $global_struct->{'crosshatch'}) {
	    $color = $global_struct->{'gray'};
	    $filled = 41;  # pattern
	  }
	  else {
	    if ((!defined($color)) && ($gene_ref->{$fragmt_id}{$element}{type} eq "ORF"))  {
	      #unless (($color) | ($gene_ref->{$fragmt_id}{$element}{type} eq "RBS") | ($gene_ref->{$fragmt_id}{$element}{type} eq "TERM") | ($gene_ref->{$fragmt_id}{$element}{type} =~ /RNA/) | ($gene_ref->{$fragmt_id}{$element}{type} eq "PMARK") | ($gene_ref->{$fragmt_id}{$element}{type} eq "PROM")){
	      printf STDERR ("can\'t find color for: %s :: %s drawing in black\n",
			     $role, $gene_ref->{$fragmt_id}{$element}{'name'});
	      $color = $global_struct->{'black'};
	    }
	    $filled = 20;  # pattern
	  }
	}
	else {
	  $pen = $global_struct->{'black'};
	  $color = $global_struct->{'white'};
	  $filled = -1;  # no fill
	}

	if (($gene_ref->{$fragmt_id}{$element}{'type'} eq "ORF") && ($global_struct->{'arrowstyle'})){
	  &draw_genes($global_struct,$level_ref,$gene_ref->{$fragmt_id}{$element},$color,$pen,$filled);
	}
	elsif ($gene_ref->{$fragmt_id}{$element}{'type'} eq "ORF"){
	  &draw_old_genes($global_struct,$level_ref,$gene_ref->{$fragmt_id}{$element},$color,$pen,$filled);
	}
	elsif ($gene_ref->{$fragmt_id}{$element}{'type'} eq "ORFNA") { # no arrowhead desired
	  &draw_gene_block($global_struct,$level_ref,$gene_ref->{$fragmt_id}{$element},$global_struct,$color,$pen,$filled);
	}
	elsif ($gene_ref->{$fragmt_id}{$element}{'type'} eq "CLS"){
	  $filled = 0; # force white fill color independent of role_id
	  &draw_clusters($global_struct,$level_ref,$gene_ref->{$fragmt_id}{$element},$color,$pen,$filled,$role);
	}
	elsif ($gene_ref->{$fragmt_id}{$element}{'type'} =~ /RNA/) {
	  &draw_RNA($global_struct,$level_ref, $gene_ref->{$fragmt_id}{$element}{'type'}, $gene_ref->{$fragmt_id}{$element}{'end5'}, $gene_ref->{$fragmt_id}{$element}{'end3'}, $gene_ref->{$fragmt_id}{$element}{'name'}, $gene_ref->{$fragmt_id}{$element}{'final_name'});
	}elsif ($gene_ref->{$fragmt_id}{$element}{'type'} eq "PROM"){ # if promoter
	  $color = $global_struct->{'green'};
	  $pen = $global_struct->{'green'};
	  &draw_promoters($global_struct,$level_ref,$gene_ref->{$fragmt_id}{$element},$color,$pen,$filled);
	}
	elsif ($gene_ref->{$fragmt_id}{$element}{'type'} eq "PRI"){ # if primer
	  $color = $global_struct->{'red'};
	  $pen = $global_struct->{'red'};
	  &draw_primers($global_struct,$level_ref,$gene_ref->{$fragmt_id}{$element},$color,$pen,$filled);
	}
	elsif ($gene_ref->{$fragmt_id}{$element}{'type'} eq "PRI2"){ # if primer
	  $color = $global_struct->{'magenta'};
	  $pen = $global_struct->{'magenta'};
	  &draw_primers($global_struct,$level_ref,$gene_ref->{$fragmt_id}{$element},$color,$pen,$filled);
	}
	elsif ($gene_ref->{$fragmt_id}{$element}{'type'} eq "TERM"){ # if terminator
	  $color = $global_struct->{'red'};
	  $pen = $global_struct->{'red'};
	  &draw_terminators($global_struct,$level_ref,$gene_ref->{$fragmt_id}{$element},$color,$pen,$filled);
	}
	elsif ($gene_ref->{$fragmt_id}{$element}{'type'} eq "PMARK"){ # if contig break
	  $color = $global_struct->{'red'};
	  $pen = $global_struct->{'red'};
	  &draw_PMARKs($global_struct,$level_ref,$gene_ref->{$fragmt_id}{$element},$color,$pen,$filled);
	}
	elsif ($gene_ref->{$fragmt_id}{$element}{'type'} eq "TSD")  { # if target site duplication
	  $color = $global_struct->{'black'};
	  $pen = $global_struct->{'black'};
	  &draw_TSDs($global_struct,$level_ref,$gene_ref->{$fragmt_id}{$element},$color,$pen,$filled);
	}
	elsif ($gene_ref->{$fragmt_id}{$element}{'type'} eq "ENZ"){ # if restriction enzyme site
	  $color = $global_struct->{'blue'};
	  $pen = $global_struct->{'blue'};
	  &draw_enzymes($global_struct,$level_ref,$gene_ref->{$fragmt_id}{$element},$color,$pen,$filled);
	}
	elsif ($gene_ref->{$fragmt_id}{$element}{'type'} eq "ATTP"){ # if restriction enzyme site
	  $color = $global_struct->{'cyan'};
	  $pen = $global_struct->{'cyan'};
	  &draw_ATTP($global_struct,$level_ref,$gene_ref->{$fragmt_id}{$element},$color,$pen,$filled);
	}
	elsif ($gene_ref->{$fragmt_id}{$element}{'type'} eq "REGION") {
	  &fake_genes($global_struct,$level_ref,$gene_ref->{$fragmt_id}{$element});
	}
	elsif ($gene_ref->{$fragmt_id}{$element}{'type'} eq "INDEL") {
	  $pen = $global_struct->{'black'};
	  $color = $global_struct->{'black'};
	  $filled = 20;
	  &draw_gene_block($global_struct,$level_ref,$gene_ref->{$fragmt_id}{$element},$global_struct,$color,$pen,$filled);
	}
	elsif ($gene_ref->{$fragmt_id}{$element}{'type'} eq "Ti") {
	  $pen = $global_struct->{'red'};
	  $color = $global_struct->{'red'};
	  $filled = 20;
	  &draw_gene_block($global_struct,$level_ref,$gene_ref->{$fragmt_id}{$element},$global_struct,$color,$pen,$filled);
	}
	elsif ($gene_ref->{$fragmt_id}{$element}{'type'} eq "Tv") {
	  $pen = $global_struct->{'blue'};
	  $color = $global_struct->{'blue'};
	  $filled = 20;
	  &draw_gene_block($global_struct,$level_ref,$gene_ref->{$fragmt_id}{$element},$global_struct,$color,$pen,$filled);
	}
      }
    }
  }
}

sub draw_chi_square {
  my ($level_ref, $fragmt_ref, $global_struct, $chi_ref) = @_;
  my $fragmt_id = $fragmt_ref->{'fragmt_id'};
  my ($x, $y1, $y2, $y3);
  for my $position (keys %{$chi_ref->{$fragmt_id}}) {
    if (&overlap($fragmt_ref->{lend},$fragmt_ref->{rend},$position,$position))  {
      $x = &x_transform($global_struct, $position, $level_ref);
      $y1 = &line2coord($global_struct, $level_ref->{'level'})+$global_struct->{'chi_position'};
      $y3 = $y1 + ($global_struct->{'gene_arrow_height'} / 2) + 300 + ($global_struct->{'chi_height'}/4); # needed for correct position of label
      #$y2 = $y3-$chi_ref->{$fragmt_id}{$position}*($global_struct->{'chi_height'}/2); # for upward-facing peaks
      $y2 = $y3+$chi_ref->{$fragmt_id}{$position}*$global_struct->{'chi_height'}; # for downward-facing peaks

      printf("2 1 0 1 %d %d %d 0 20 0.000 0 0 0 0 0 2\n",
	     $global_struct->{'gray'},
	     $global_struct->{'gray'},
	     $global_struct->{'row_depth'});
      printf("\t %d %d %d %d\n",
	     $x,$y3,
	     $x,$y2);
    }
  }
}

sub draw_gene_block {
  my ($level_ref, $gene_ref, $global_struct, $color, $pen, $filled) = @_;
  my ($x1, $x2, $x3, $y1, $y2, $y3, $y4);
  $x1 = &x_transform($global_struct, $gene_ref->{'end5'}, $level_ref);
  $x3 = &x_transform($global_struct, $gene_ref->{'end3'}, $level_ref);

#################################

  #      $y2----------+
  #       |           |
  #      $y4----------+
  #
  #      $x3           $x1

  #$y1 = &line2coord($global_struct, $level_ref->{'level'});
  $y1 = &line2coord($global_struct, $level_ref->{'level'})-7;
  $y4 = $y1;
  $y2 = $y4 - ($global_struct->{'gene_height'})/4;
  #$y3 = $y1 + ($global_struct->{'gene_arrow_height'} / 2); # needed for correct position of label
  #$y2 = $y1 + ($global_struct->{'gene_arrow_height'} / 2) - ($global_struct->{'gene_height'} / 2);
  #$y4 = $y1 + ($global_struct->{'gene_arrow_height'} / 2) + ($global_struct->{'gene_height'} / 2);

  printf("2 3 0 1 %d %d %d 0 %d 0.000 0 0 0 0 0 5\n",
	 $pen,
	 $color,
	 $global_struct->{'gene_depth'},
	 $filled);
  printf("\t %d %d %d %d %d %d %d %d %d %d\n",
	 $x3,$y2,
	 $x1,$y2,
	 $x1,$y4,
	 $x3,$y4,
	 $x3,$y2);

  #&add_label($gene_ref->{'final_name'},$x1,$x3,$y3+800,1,1,0,15); # name, x1, x3, y, 1|0 (below or above), justification (0left|1right|2center), font size in points
  # name; x1; x2; y; position (vert=2, horz = 1, else diagonal); justification (left = 0, center = 1, right = 2); font (Times Roman = 0, Times Italix = 1, Times Bold = 2); font size; color; mid, rotation #

  $gene_ref->{'x1'} = $x1;
  $gene_ref->{'x3'} = $x3;
  $gene_ref->{'y2'} = $y1;
  $gene_ref->{'y4'} = $y4;
}

sub draw_old_genes { # old arrowhead style, pointed box really
  my ($global_struct, $level_ref, $gene_ref, $color, $pen, $filled) = @_;
  my ($x1, $x2, $x3, $y1, $y2, $y3);

  $x1 = &x_transform($global_struct, $gene_ref->{'end5'}, $level_ref);
  $x3 = &x_transform($global_struct, $gene_ref->{'end3'}, $level_ref);
  $x2 = &calc_arrow_x($global_struct,$x1, $x3, $gene_ref->{'end5'}, $gene_ref->{'end3'}, $level_ref->{'orientation'});

  #    $y2------------+
  #    /              |
  #   /               |
  #$y1                |
  #  \                |
  #   \               |
  #    $y3------------+
  #
  #$x3  $x2           $x1
  #$end3              $end5

  $y1 = &line2coord($global_struct, $level_ref->{'level'});
  $y2 = $y1 - ($global_struct->{'gene_height'} / 2);
  $y3 = $y1 + ($global_struct->{'gene_height'} / 2);

  printf("2 3 0 1 %d %d %d 0 %d 0.000 0 0 0 0 0 6\n",
	 $pen, 
	 $color,
	 $global_struct->{'gene_depth'},
	 $filled);

  printf("\t %d %d %d %d %d %d %d %d %d %d %d %d\n",
	 $x3,$y1,
	 $x2,$y2,
	 $x1,$y2,
	 $x1,$y3,
	 $x2,$y3,
	 $x3,$y1);
  if ($global_struct->{'no_names'} == 0) {
    &add_label($gene_ref->{'final_name'},$x1,$x3,$y1-800,1,1,0,30); # name, x1, x2, y, 1|0 (below or above), justification (0left|1right|2center), font, font size in points
    #&add_label($gene_ref->{'final_name'},$x1,$x3,$y1);
  }

  $gene_ref->{'x1'} = $x1;
  $gene_ref->{'x3'} = $x3;
  $gene_ref->{'x2'} = $x2;
  $gene_ref->{'y2'} = $y2;
  $gene_ref->{'y3'} = $y3;
}

sub draw_clusters { #dfouts 06/08/2016
  my ($global_struct, $level_ref, $gene_ref, $color, $pen, $filled, $role) = @_;
  my ($x1, $x2, $x3, $x_radius, $y_radius, $y1, $y2, $y3, $fill_color);

  $x1 = &x_transform($global_struct, $gene_ref->{'end5'}, $level_ref);
  $x3 = &x_transform($global_struct, $gene_ref->{'end3'}, $level_ref);
  $x2 = int(abs(($x3+$x1)/2)); # determine center of circle
  $x_radius = abs($x1-$x2);

  #     ----$y1----
  #   /             \
  #  /               \
  # |       $y2       |
  #  \               /
  #   \             /
  #     ----$y3-----
  #
  #$x3      $x2      $x1

  $y1 = &line2coord($global_struct, $level_ref->{'level'}); #Take -600 off for horizontal line above genes
  $y2 = $y1 + ($global_struct->{'gene_arrow_height'} / 2);
  $y3 = $y1 + $global_struct->{'gene_arrow_height'};
  #$y_radius = abs($y1-$y2); # for more elliptical circles where only width is scaled to gene length
  $y_radius = $x_radius; # for perfect circles scaled to gene length

  if ($role == "1019") {
    $fill_color = $color;
  }
  else {
    $fill_color = $filled;
  }

  &print_circle($x2, $y2, $x_radius, $y_radius, $pen, $fill_color);

  if ($global_struct->{'no_names'} == 0) {
    &add_label($gene_ref->{'final_name'},$x1,$x3,$y2+800,2,2,0,15);
    # name; x1; x2; y; position (vert=2, horz = 1, else diagonal); justification (left = 0, center = 1, right = 2); font (Times Roman = 0, Times Italix = 1, Times Bold = 2); font size; color; mid, rotation #
  }
  # store coords for drawing gene pairs if pairs file is provided
  $gene_ref->{'x1'} = $x1;
  $gene_ref->{'x3'} = $x3;
  $gene_ref->{'x2'} = $x3;
  $gene_ref->{'y2'} = $y1;
  $gene_ref->{'y3'} = $y2;
}

sub draw_genes { #dfouts 07/22/2005
  my ($global_struct, $level_ref, $gene_ref, $color, $pen, $filled) = @_;
  my ($x1, $x2, $x3, $y1, $y2, $y3, $y4, $y5);

  $x1 = &x_transform($global_struct, $gene_ref->{'end5'}, $level_ref);
  $x3 = &x_transform($global_struct, $gene_ref->{'end3'}, $level_ref);
  $x2 = &calc_arrow_x($global_struct,$x1, $x3, $gene_ref->{'end5'}, $gene_ref->{'end3'}, $level_ref->{'orientation'});

  #    $y1
  #    / |
  #   /  $y2----------+
  #$y3                |
  #  \   $y4----------+
  #   \  |
  #    $y5
  #
  #$x3  $x2           $x1

  $y1 = &line2coord($global_struct, $level_ref->{'level'}); #Take -600 off for horizontal line above genes
  $y3 = $y1 + ($global_struct->{'gene_arrow_height'} / 2);
  $y2 = $y1 + ($global_struct->{'gene_arrow_height'} / 2) - ($global_struct->{'gene_height'} / 2);
  $y4 = $y1 + ($global_struct->{'gene_arrow_height'} / 2) + ($global_struct->{'gene_height'} / 2);
  $y5 = $y1 + $global_struct->{'gene_arrow_height'};

  printf("2 3 0 1 %d %d %d 0 %d 0.000 0 0 0 0 0 8\n",
	 $pen,
	 $color,
	 $global_struct->{'gene_depth'},
	 $filled);
  printf("\t %d %d %d %d %d %d %d %d %d %d %d %d\n",
	 $x3,$y3,
	 $x2,$y1,
	 $x2,$y2,
	 $x1,$y2,
	 $x1,$y4,
	 $x2,$y4);
  printf("\t %d %d %d %d\n",
	 $x2,$y5,
	 $x3,$y3);

  if ($global_struct->{'no_names'} == 0) {
    &add_label($gene_ref->{'final_name'},$x1,$x3,$y3+800,2,2,0,15);
    # name; x1; x2; y; position (vert=2, horz = 1, else diagonal); justification (left = 0, center = 1, right = 2); font (Times Roman = 0, Times Italix = 1, Times Bold = 2); font size; color; mid, rotation #
  }

  $gene_ref->{'x1'} = $x1;
  $gene_ref->{'x3'} = $x3;
  $gene_ref->{'x2'} = $x2;
  $gene_ref->{'y2'} = $y1;
  $gene_ref->{'y3'} = $y5;
}

sub draw_thingene {
  my ($global_struct, $level_ref, $gene_ref) = @_;
  my ($x1, $x2, $y);

  $x1 = x_transform($global_struct, $gene_ref->{'end5'}, $level_ref);
  $x2 = x_transform($global_struct, $gene_ref->{'end3'}, $level_ref);

  #
  #   /
  #  /
  #  ---------------- $y
  #  \
  #   \
  #
  #
  # $x2               $x1
  # $end3             $end5

  $y = line2coord($global_struct, $level_ref->{'level'});

  printf("2 1 0 10 %d %d %d 0 -1 0.000 0 0 0 0 0 2\n",
	 $global_struct->{'red'},
	 $global_struct->{'red'},
	 $global_struct->{'gene_depth'}-1);

  printf("\t %d %d %d %d %d %d %d %d %d %d %d %d\n",
	 $x1,$y,
	 $x2,$y);

  if ($global_struct->{'no_names'} == 0) {
    &add_label($gene_ref->{'final_name'},$x1,$x2,$y+800,0,1,1,15);
  }
}

######################################################################
# inserted 03/02/2004 by Derrick E. Fouts from db2fig.dbi
# modified by Derrick Fouts

sub print_squiggle {
  my ($level_ref, $s, $array_ref) = @_;
  #local($s, @l) = @_;
  my ($i, $len, $max_pairs, $str, $max_points);

  $max_pairs = 12;
  #print STDERR ("array = @{$array_ref}\ns = $s\nlevel = $level_ref\n");
  $i = (@{$array_ref} / 2);
  $str = "$s $i";

  $len = @{$array_ref};
  for($i=0;$i<$len-1;$i+=2) {
    $str .= sprintf("\n\t ") if (($i % $max_pairs) == 0);
    $str .= sprintf("%d %d ",
		    ${$array_ref}[$i],
		    ${$array_ref}[$i+1]);
  }
  $str .= sprintf("\n");
  $max_points = 8;

  $str .= sprintf("\t ");
  for($i=1;$i<($len/2) + 1;$i++) {
    if ($i == 1 || $i == ($len/2)) {
      $str .= sprintf("0.000 ");
    }
    else {
      $str .= sprintf("-1.000 ");
    }
    $str .= sprintf("\n\t ") if (($i % $max_points) == 0);
  }
  $str .= sprintf("\n");

  $str =~ s/\s+\n/\n/g;

  printf $str;
}

sub draw_RNA {
  my ($global_struct, $level_ref, $type, $end5, $end3, $name, $final_name) = @_;
  my ($x, $x1, $x2, $y, $width, $height, $label_y);
  my ($i, $squig_params, @array);
  my ($min_y, $max_y, $min_x, $max_x);

  $x1 = &min($end5, $end3);
  $x2 = &max($end5, $end3);

  $width = $x2-$x1;
  $x = &x_transform($global_struct, $x1 + ($width/2), $level_ref);
  $y = &line2coord($global_struct, $level_ref->{'level'});# - ($global_struct->{'gene_height'}/6);
  #$label_y = $y + ($height);
  #if (!$global_struct->{'tickmark'})  { # if no tickmarks
  #$y = &line2coord($global_struct, $level_ref->{'level'});# - ($global_struct->{'gene_height'}*2);
  #$label_y = $y - ($height/4);
  #}

  #$label_y = $y - $global_struct->{'label_offset'};

  #print STDERR ("type = $type\tname = $name\n");

  if (($name =~ /tRNA/) || ($type =~ /tRNA/)) {
    @array = @{$global_struct->{'tRNarray'}};
    $width = $global_struct->{'trna_width'};
    $height = $global_struct->{'trna_height'};
  }
  elsif ($name =~ /23S/) {
    @array = @{$global_struct->{'rRNarray23S'}};
    $height = $global_struct->{'rrna_height'};
  }
  elsif ($name =~ /16S/) {
    @array = @{$global_struct->{'rRNarray16S'}};
    $height = $global_struct->{'rrna_height'};
  }
  elsif ($name =~ /7S|SRP/) {
    @array = @{$global_struct->{'sRNarray7S'}};
    $height = $global_struct->{'rrna_height'};
    print STDERR "\tDRAWING: $name\n";
  }
  elsif ($name =~ /5S/) {
    @array = @{$global_struct->{'rRNarray5S'}};
    $height = $global_struct->{'rrna_height'};
  }
  elsif (($name =~ /tmRNA/) || ($type =~ /tmRNA/)) {
    @array = @{$global_struct->{'tmRNarray'}};
    $width = $global_struct->{'trna_width'};
    $height = $global_struct->{'trna_height'};
    print STDERR "\tDRAWING: $name\n";
  }
  elsif ($name =~ /RNaseP|rnpB/) {
    @array = @{$global_struct->{'sRNarrayP'}};
    $height = $global_struct->{'srna_height'};
    print STDERR "\tDRAWING: $name\n";
  }
  else {
    print STDERR "Could not find structural information for $name!\n";
    return;
  }

  for($i=0; $i<@array-1;$i+=2) {
    $array[$i] = ($array[$i] * $width) + $x;# - ($height/2);
    $array[$i+1] = ($array[$i+1] * $height) + $y - ($height/2);
  }

  $min_y = $y - $height; #- ($global_struct->{'point_size'}*10);
  $max_y = $y;
  $min_x = $x - $width;
  $max_x = int($x + ($width/2) + $global_struct->{'char_width'});

  print "\#$type $name\n";
  print "6 $min_x $min_y $max_x $max_y\n";

  $squig_params = "3 0 0 1 0 7 0 0 -1 0.000 0 0 0";

  &print_squiggle($level_ref, $squig_params, \@array);

  $label_y = $y + ($global_struct->{'gene_arrow_height'} / 2); # same as draw_genes
  #$label_y = $y - ($height+$height/4); #old code
  #if (!$global_struct->{'arrowstyle'})  { # if traditional arrowhead for genes
  #	$label_y = $y - ($height/4);
  #      }
  #print label
  &add_label($final_name, $x, $x, $label_y+800, 2, 2, 0, 15); # modified to print in same area as ORFs

  #&add_label($final_name, $x, $x, $label_y, 1, 1, 0, 15); # original code

  #   $y3 = $y1 + ($global_struct->{'gene_arrow_height'} / 2); from draw_genes sub
  #&add_label($gene_ref->{'final_name'},$x1,$x3,$y3+800,2,2,0,15); #draw_genes code
  # name; x1; x2; y; position (vert=2, horz = 1, else diagonal); justification (left = 0, center = 1, right = 2); font (Times Roman = 0, Times Italix = 1, Times Bold = 2); font size; color; mid, rotation #

  print "-6\n";
  print "\#END OF $type\n";
}

# added 03/04/2004 from spec20functions.pl dfouts

sub draw_promoters { # dfouts
  my ($global_struct, $level_ref, $gene_ref, $color, $pen, $filled) = @_;
  my ($x1, $x2, $x3, $x4, $x5, $y1, $y2, $y3, $y4, $mid);
  my $temp_depth = $global_struct->{'gene_depth'};

  $x1 = &x_transform($global_struct, $gene_ref->{'end5'}, $level_ref);
  $x3 = &x_transform($global_struct, $gene_ref->{'end3'}, $level_ref);
  $x2 = $x3;

  $y1 = &line2coord($global_struct, $level_ref->{'level'});
  $y2 = $y1 - ($global_struct->{'gene_arrow_height'} / 5);
  $y3 = $y1-7;
  if (!$global_struct->{'arrowstyle'})  { # if traditional arrowhead for genes
    $y2 = $y1 - ($global_struct->{'gene_height'} / 2);
    $y3 = $y1 + ($global_struct->{'gene_height'} / 2) + 150;
  }

  if ($x1 > $x2)  {  # if promoter is on minus strand <----
    $mid = $x2 + (($x1 - $x2)/2);
    $x5 = $mid-300;

    # F O R  A N  A R R O W  O N  A  L I N E
    #
    #                +-----+ $y2
    #                |     |
    #                |     |
    #   /            |     |
    #  /|------------+     |
    # / |            |     |
    # \ |            |     |
    #  \|------------+-----+ $y3
    #   \
    #  $x5            $mid
    #
    #          $end3($x2) $end5($x1)
    #
    #
    printf("2 1 0 3 %d %d %d 0 %d 0.000 0 0 0 0 0 2\n",
	   $pen, 
	   $color,
	   $global_struct->{'gene_depth'}-2,
	   $filled);

    if ($global_struct->{'arrowstyle'})  { # if new arrowhead style for genes
      printf("\t %d %d %d %d\n",
	     $mid,$y3,
	     $mid,$y2);
    }
    else {
      printf("\t %d %d %d %d\n",
	     $mid,$y2,
	     $mid,$y3);
    }

    printf("2 1 0 3 %d %d %d 0 %d 0.000 0 0 0 1 0 2\n",
	   $pen, 
	   $color,
	   $global_struct->{'gene_depth'}-2,
	   $filled);

    printf("\t 2 1 1 200 175\n");

    if ($global_struct->{'arrowstyle'})  { # if new arrowhead style for genes
      printf("\t %d %d %d %d\n",
	     $mid+22.5,$y2,
	     $x5,$y2);
    }
    else {
      printf("\t %d %d %d %d\n",
	     $mid+25,$y3,
	     $x5,$y3);
    }
  }
  else  {
    $mid = $x1 + (($x2 - $x1)/2);
    if (!$global_struct->{'arrowstyle'})  { # if traditional arrowhead for genes
      $y2 = $y1 - ($global_struct->{'gene_height'} / 2) - 150;
    }
    $x5 = $mid+300;

    printf("2 1 0 3 %d %d %d 0 %d 0.000 0 0 0 0 0 2\n",
	   $pen, 
	   $color,
	   $global_struct->{'gene_depth'}-2,
	   $filled);

    printf("\t %d %d %d %d\n",
	   $mid,$y3,
	   $mid,$y2);

    printf("2 1 0 3 %d %d %d 0 %d 0.000 0 0 0 1 0 2\n",
	   $pen, 
	   $color,
	   $global_struct->{'gene_depth'}-2,
	   $filled);

    printf("\t 2 1 1 200 175\n");

    printf("\t %d %d %d %d\n",
	   $mid-20,$y2,
	   $x5,$y2);
  }
  #if ($x1 > $x2)  {  # if promoter is on minus strand <----
  #  &add_label($gene_ref->{'final_name'},$x1,$x5,$y3+1100);
  #}
  #else  {
  #  &add_label($gene_ref->{'final_name'},$x1,$x5,$y2+75);
  #}
  $gene_ref->{'x1'} = $mid;
  $gene_ref->{'x3'} = $mid;
  $gene_ref->{'x2'} = $mid;
  $gene_ref->{'y2'} = $y2;
  $gene_ref->{'y3'} = $y3;
}

sub draw_primers { # dfouts
  my ($global_struct, $level_ref, $gene_ref, $color, $pen, $filled) = @_;
  my ($x1, $x2, $x3, $x4, $x5, $y1, $y2, $y3, $y4, $mid);
  my $temp_depth = $global_struct->{'gene_depth'}-3;
  my ($name) = $gene_ref->{'final_name'};

  $x1 = &x_transform($global_struct, $gene_ref->{'end5'}, $level_ref);
  $x3 = &x_transform($global_struct, $gene_ref->{'end3'}, $level_ref);
  $x2 = $x3;

  $y1 = &line2coord($global_struct, $level_ref->{'level'});

  if ($x1 > $x2)  {  # if primer is on minus strand <----
    $mid = $x2 + (($x1 - $x2)/2);
    $y3 = $y1 - ($global_struct->{'gene_height'});
    #$y3 = $y1 - ($global_struct->{'gene_height'}/6);
    if (!$global_struct->{'arrowstyle'}) { # for old style arrowheads
      $y3 = $y1 + ($global_struct->{'gene_height'} / 2) + 150;
    }
    $x5 = $mid-300;

    # F O R  A N  A R R O W
    #
    #   /
    #  /|------------+-----+
    # / |                  |
    # \ |                  |
    #  \|------------+-----+ ($y3, below; $y2, above)
    #   \
    #  $x5            $mid
    #
    #          $end3($x2) $end5($x1)
    #
    #
    printf("2 1 0 3 %d %d %d 0 %d 0.000 0 0 0 1 0 2\n",
	   $pen, 
	   $color,
	   $temp_depth,
	   $filled);

    printf("\t 1 1 1 200 175\n");

    printf("\t %d %d %d %d\n",
	   $mid,$y3,
	   $x5,$y3);
    if ($global_struct->{'arrowstyle'}) { # for old style arrowheads
      &add_label($name,$x5,$x1,$y3-150,1,1,0,10,$color);
    }
    else {
      &add_label($name,$x2,$x1,$y3+550,1,2,0,10,$color);  # name, x1, x2, y, 2|1|0 (vert, horiz or diag), justification (0left|1center|2right), font size in points, color
      # name; x1; x2; y; position (vert=2, horz = 1, else diagonal); justification (left = 0, center = 1, right = 2); font (Times Roman = 0, Times Italix = 1, Times Bold = 2); font size; color; mid, rotation #
    }
  }
  else  {
    $mid = $x1 + (($x2 - $x1)/2);
    $y2 = $y1 - ($global_struct->{'gene_height'});
    #$y2 = $y1 - ($global_struct->{'gene_height'}/6);
    if (!$global_struct->{'arrowstyle'}) { # for old style arrowheads
      $y2 = $y1 - ($global_struct->{'gene_height'} / 2) - 150;
    }
    $x5 = $mid+300;

    printf("2 1 0 3 %d %d %d 0 %d 0.000 0 0 0 1 0 2\n",
	   $pen, 
	   $color,
	   $temp_depth,
	   $filled);

    printf("\t 1 1 1 200 175\n");

    printf("\t %d %d %d %d\n",
	   $mid,$y2,
	   $x5,$y2);
    if ($global_struct->{'arrowstyle'}) { # for old style arrowheads
      &add_label($name,$x5,$x1,$y2-150,1,1,0,10,$color);
    }
    else {
      &add_label($name,$x1+200,$x2+200,$y2-150,2,2,0,10,$color);  # name, x1, x2, y, 2|1|0 (vert, horiz or diag), justification (0left|1center|2right), font size in points, color
      # name; x1; x2; y; position (vert=2, horz = 1, else diagonal); justification (left = 0, center = 1, right = 2); font (Times Roman = 0, Times Italix = 1, Times Bold = 2); font size; color; mid, rotation #
    }
  }
  $gene_ref->{'x1'} = $x1;
  $gene_ref->{'x2'} = $x2;
  $gene_ref->{'y2'} = $y2;
  $gene_ref->{'y3'} = $y3;
}

sub draw_terminators { # dfouts
  my ($global_struct, $level_ref, $gene_ref, $color, $pen, $filled) = @_;
  my ($x1, $x2, $x3, $y1, $y2, $mid);
  my $temp_depth = $global_struct->{'gene_depth'}-4;

  $x1 = &x_transform($global_struct, $gene_ref->{'end5'}, $level_ref);
  $x2 = &x_transform($global_struct, $gene_ref->{'end3'}, $level_ref);
  $x3 = int(abs(($x1+$x2)/2));

  $y1 = &line2coord($global_struct, $level_ref->{'level'})-14;
  $y2 = ($y1 - ($global_struct->{'gene_height'} / 2)) - 75;

  &print_circle($x3, $y2-100, 100, 100, $color, -1);

  #################################s############################
  # Old code for printing terminators on the minus strand
  #
  #if ($x1 > $x2)  {  # if terminator is on minus strand <----
  #  $y3 = $y3+75;
  #  &print_circle($x4, $y3+100, $color, -1);
  #}
  #else  {
  #  $y2 = $y2-75;
  #  &print_circle($x4, $y2-100, $color, -1);
  #}
  #############################################################

  #if ($x1 > $x2)  {  # if terminator is on minus strand <----
  #  &add_label($gene_ref->{'final_name'},$x1,$x3,$y3+1100);
  #}
  #else  {
  #  &add_label($gene_ref->{'final_name'},$x1,$x3,$y2+75);
  #}

  # F O R  T E R M I N A T O R S
    #
    #    +-----+ $y2
    #    |     |
    #    |     |
    #    |     |
    #    |     |
    #    |     |
    #    |     |
    #    +-----+ $y1
    #   $x2   $x1
    #
    #
  printf("2 3 0 1 %d %d %d 0 %d 0.000 0 0 0 0 0 5\n",
         $pen,
	 $color,
         $temp_depth,
         $filled);

  printf("\t %d %d %d %d %d %d %d %d %d %d\n",
	 $x2,$y1,
	 $x2,$y2,
	 $x1,$y2,
	 $x1,$y1,
	 $x2,$y1);

  $gene_ref->{'x1'} = $x1;
  $gene_ref->{'x3'} = $x2;
  $gene_ref->{'x2'} = $x2;
  $gene_ref->{'y1'} = $y1;
  $gene_ref->{'y2'} = $y2;
  $gene_ref->{'y3'} = $y1;
}

sub draw_PMARKs { # dfouts
  my ($global_struct, $level_ref, $gene_ref, $color, $pen, $filled) = @_;
  my ($temp_depth, $x1, $x2, $x3, $x4, $x5, $y1, $y2, $y3, $y4);

  $temp_depth = $global_struct->{'gene_depth'};
  $x1 = &x_transform($global_struct, $gene_ref->{'end5'}, $level_ref);
  $x3 = &x_transform($global_struct, $gene_ref->{'end3'}, $level_ref);
  $x2 = $x3;
  $y1 = &line2coord($global_struct, $level_ref->{'level'})-14;
  $y2 = $y1 - ($global_struct->{'gene_height'} / 2);
  $y3 = $y1;
  if (!$global_struct->{'arrowstyle'}) { # for traditional arrowheads
      $y3 = $y1 + ($global_struct->{'gene_height'} / 2);
  }

  $y4 = abs($y2-$y3);  # determine the height of the "gene"
  $y2 = $y2 - $y4;
  #$y3 = $y3 + $y4;

  # F O R  A  P M A R K
  #
  #    +-----+ $y2
  #    |     |
  #    |     |
  #    |     |
  #    +     + $y1
  #    |     |
  #    |     |
  #    |     |
  #    +-----+ $y3
  #   $x1   $x2,$x3
  #
  #

  printf("2 3 0 1 %d %d %d 0 %d 0.000 0 0 0 0 0 6\n",
         $pen,
	 $color,
         $temp_depth,
         $filled);

  printf("\t %d %d %d %d %d %d %d %d %d %d %d %d\n",
	 $x3,$y1,
	 $x2,$y2,
	 $x1,$y2,
	 $x1,$y3,
	 $x2,$y3,
	 $x3,$y1);

  #&add_label($gene_ref->{'final_name'},$x1,$x3,$y1-150-$y4);
  $gene_ref->{'x1'} = $x1;
  $gene_ref->{'x3'} = $x3;
  $gene_ref->{'x2'} = $x2;
  $gene_ref->{'y2'} = $y2;
  $gene_ref->{'y3'} = $y3;
}

sub draw_enzymes { # dfouts
  my ($global_struct, $level_ref, $gene_ref, $color, $pen, $filled) = @_;
  my ($depth, $x1, $y1, $y2);
  my ($name) = $gene_ref->{'final_name'};

  $x1 = &x_transform($global_struct, $gene_ref->{'end5'}, $level_ref);
  $y1 = &line2coord($global_struct, $level_ref->{'level'});
  $y2 = $y1 - ($global_struct->{'gene_height'});
  $depth = $global_struct->{'gene_depth'}+1; # set depth of tick (larger # = more behind)

    # F O R  A  T I C K  M A R K
    #
    #    +-----+ $y2
    #    |     |
    #    |     |
    #    |     |
    #    |     |
    #    |     |
    #    |     |
    #    +-----+ $y1
    #      $x1
    #
    #

  printf("2 1 0 3 %d %d %d 0 %d 0.000 0 0 0 0 0 2\n",
	 $pen,
	 $color,
	 $depth,
	 $filled);

  printf("\t %d %d %d %d\n",
	 $x1,$y1,
	 $x1,$y2);
  &add_label($name,$x1,$x1,$y2-50,0,0,2,15,1);
  # eg. &add_label($gene_ref->{'final_name'},$x1,$x3,$y3+800,2,2,0,15);
  # name; x1; x2; y; position (vert=2, horz = 1, else diagonal); justification (left = 0, center = 1, right = 2); font (Times Roman = 0, Times Italix = 1, Times Bold = 2); font size; color; mid, rotation # 
  $gene_ref->{'x1'} = $x1;
  $gene_ref->{'y1'} = $y1;
  $gene_ref->{'y2'} = $y2;
}

sub draw_ATTP { # dfouts
  my ($global_struct, $level_ref, $gene_ref, $color, $pen, $filled) = @_;
  my $depth = $global_struct->{'gene_depth'}+1; # set depth of tick (larger # = more behind)
  my ($x1, $x2, $mid, $y1, $y2);
  my ($name) = $gene_ref->{'final_name'};

  $x1 = &x_transform($global_struct, $gene_ref->{'end5'}, $level_ref);
  $x2 = &x_transform($global_struct, $gene_ref->{'end3'}, $level_ref);
  $mid = $x2 + (($x1 - $x2)/2);
  $y1 = &line2coord($global_struct, $level_ref->{'level'});
  $y2 = $y1 + ($global_struct->{'gene_height'}) * 1.5;

  # F O R  A  T I C K  M A R K
  #
  #    +-----+ $y2
  #    |     |
  #    |     |
  #    |     |
  #    |     |
  #    |     |
  #    |     |
  #    +-----+ $y1
  #    $x1   $x2
  #       $mid
  #
  #  printf("2 1 0 3 %d %d %d 0 %d 0.000 0 0 0 1 0 2\n",
  #     $pen,
  #	 $color,
  #        $global_struct->{'gene_depth'}-2,
  #        $filled);

  #      printf("\t 2 1 1 200 175\n");

  #      printf("\t %d %d %d %d\n",
  #         $mid+25,$y3,
  #         $x5,$y3);
  #
  printf("2 1 0 3 %d %d %d 0 %d 0.000 0 0 0 0 1 2\n",
	 $pen,
	 $color,
	 $depth,
	 $filled);

  printf("\t 1 1 1 200 175\n");

  printf("\t %d %d %d %d\n",
	 $mid,$y1,
	 $mid,$y2);
  &add_label($name,$x1,$x2,$y2+150,1,1,0,10,$color);  # name, x1, x2, y, 1|0 (horiz or diag), justification (0left|1center|2right), font size in points, color 
  $gene_ref->{'x1'} = $x1;
  $gene_ref->{'x2'} = $x2;
  $gene_ref->{'y1'} = $y1;
  $gene_ref->{'y2'} = $y2;
}

sub draw_TSDs  {
  my ($global_struct, $level_ref, $gene_ref, $color, $pen, $filled) = @_;
  my ($x1, $x2, $x3, $x4, $y1, $y2, $y3, $y4);

  $x1 = &x_transform($global_struct, $gene_ref->{'end5'}, $level_ref);
  $x3 = &x_transform($global_struct, $gene_ref->{'end3'}, $level_ref);
  $x2 = $x3;
  $x4 = int(abs(($x1+$x2)/2));
  $y1 = &line2coord($global_struct, $level_ref->{'level'});
  $y2 = $y1 - ($global_struct->{'gene_height'} / 2);
  $y3 = $y1 + ($global_struct->{'gene_height'} / 2);

  $y4 = abs($y2-$y3);  # determine the height of the "gene"
  $y2 = $y2 - $y4;

  &print_triangle($x4, $y2, $color, $color, 20);

  $gene_ref->{'x1'} = $x1;
  $gene_ref->{'x3'} = $x3;
  $gene_ref->{'x2'} = $x2;
  $gene_ref->{'y2'} = $y2-75;
  $gene_ref->{'y3'} = $y3;
}

sub print_triangle { #dfouts
  my ($x, $y, $color, $fill_color, $fill_value) = @_;
  my $size = 150;
  #$y = $y-($size*1.5/2);

  # Fields are:
  #    2 - polyline
  #    2 - box
  #    0 - solid
  #    1 - thickness
  #    0 - black (line)
  #    0 - black (fill)
  #  100 - depth
  #    0 - pen style
  #   -1 - area fill
  # 0.000 - style val ??
  #    0 - join style
  #    0 - cap style
  #   -1 - radius
  #    0 - forward arrow
  #    0 - backward arrow
  #    2 - number of points
  print "2 3 0 1 $color $fill_color 100 0 $fill_value 0.000 0 0 -1 0 0 4\n";
  # point coords
  # printf ("\t%d %d %d %d %d %d %d %d\n", $x - $size, $y,
  #	    $x, $y + ($size*1.5), $x + $size, $y, $x - $size, $y);
  printf ("\t%d %d %d %d %d %d %d %d\n", $x - $size, $y - ($size*1.5),
	  $x, $y, $x + $size, $y - ($size*1.5), $x - $size, $y - ($size*1.5));
}

# print_circle (from Mihai Pop's xfig.pm module
#
# draws a circle centered at the given coordinates in the specified color

sub print_circle {
  my ($x, $y, $rx, $ry, $pen, $fill) = @_;

  # Fields are:
  #    1 - ellipse
  #    3 - circle
  #    0 - solid
  #    1 - thickness
  #    0 - black (line)
  #    0 - black (fill)
  #  100 - depth
  #    0 - pen style
  #   -1 - area fill
  # 0.000 - style val ??
  #    1 - direction
  #    0 - angle
  #    cx,cy - center coords
  #    rx,ry - radius
  #    sx,sy - first point
  #    ex,ey - last point
  #    print "1 3 0 2 $color $color 100 0 $fill 0.000 1 0.000 $x $y ",
  #    "100 100 $x $y ", $x + 100, " $y\n";
  print "1 3 0 2 $pen $fill 100 0 $fill 0.000 1 0.000 $x $y ",
    "$rx $ry $x $y ", $x, " $y\n";
}

sub init_color_spectrum {
  my ($total_colors) = shift;
  my (%s);
  my $pi = atan2(1,1) * 4;
  my $r = 0; my $g = 0; my $b = 0;
  #print STDERR "->$total_colors\n";
  for(my $i=0;$i<=$total_colors;$i++) {
    # spectrum algorithm by John Quackenbush.
    # the value of $i must be between 0 and 255
    $r = 255 if ($i <=32);
    $r = abs(int(255*cos(($i-32)*$pi/151))) if ($i>32 && $i <= 107)  ;
    $r = 0 if ($i > 107);

    $g = 0 if ($i < 4);
    $g = abs(int(255*cos(($i-100)*$pi/189))) if ($i > 4 && $i <= 100) ;
    $g = abs(int(255*cos(($i-100)*$pi/294))) if ($i > 100 && $i <= 228) ;
    $g = 0 if ($i > 230) ;

    $b = 0 if ($i < 71);
    $b = abs(int(255*cos(($i-199)*$pi/256))) if ($i >= 71 && $i <= 199) ;
    $b = abs(int(255*cos(($i-199)*$pi/113))) if ($i > 199) ;

    $s{$i} = lc(sprintf("%0.2x%0.2x%0.2x", $r, $g, $b));
  }
  return(\%s);
}

sub position_to_color { # from Jeremy Peterson's script (rewritten by Derrick Fouts to import a ratio)
  my ($ratio, $max_colors, $color_start) = @_;
  my $c = $color_start + (int($ratio * $max_colors));
  return($c);
}

sub box {
  my ($pen_color, $fill_color, $x, $y, $width, $height, $border, $fill) = @_;
  #    (x1,y1)   (x2,y1)
  #           ---
  #           | |
  #           ___
  #    (x1,y2)   (x2,y2)
  #border should be 0 or 1 depending on whether you want a border around the box

  my $x1 = $x;
  my $y1 = $y;
  my $x2 = $x1 + $width;
  my $y2 = $y + $height;

  my $fill_value = 20;
  $fill_value = -1 if($fill == 0);

  my $xfig = sprintf("2 2 0 %d %d %d 50 0 %d 0.000 0 0 -1 0 0 5\n",$border,$pen_color,$fill_color,$fill_value);
  $xfig .= sprintf("         %d %d %d %d %d %d %d %d %d %d\n",$x1,$y1,$x2,$y1,$x2,$y2,$x1,$y2,$x1,$y1);

  return($xfig);
}

sub text {
  my ($global_struct, $text_size, $justification, $angle, $x, $y, $text) = @_;
  # label that tick mark.
  #print STDERR "$text_size, $justification\n";
  #    my $xfig = sprintf("4 %d %d 0 0 16 60 %d 4 765 3225 0", $justification,$black,$angle);
  my $xfig = sprintf("4 %d %d 99 0 0 %d %d 4 765 3225 0", $justification,$global_struct->{'black'},$text_size,$angle);

  $xfig .= sprintf("%d %d %s\\001\n",
		   $x,
		   $y,
		   $text);

  #$xfig .= sprintf("%d %d %s\\001\n",
  #	     $x,
  #	     $y,
  #	     &commas($text));
  return($xfig);
}

######################################################################

sub calc_arrow_x {
  my ($global_struct, $x1, $x3, $end5, $end3, $dir) = @_;
  my ($x);

  if ($dir == $global_struct->{'forward'}) {
    if($end5 < $end3) {
      $x = $x3 - $global_struct->{'gene_arrow_width'};
      $x = max($x,$x1);
    }
    else {
      $x = $x3 + $global_struct->{'gene_arrow_width'};
      $x = min($x,$x1);
    }
  }
  else {
    if($end5 < $end3) {
      $x = $x3 + $global_struct->{'gene_arrow_width'};
      $x = min($x,$x1);
    }
    else {
      $x = $x3 - $global_struct->{'gene_arrow_width'};
      $x = max($x,$x1);
    }
  }
  return($x);
}

sub add_label {
  my ($n) = shift;# ($name,$x1-100,$x1-100,$y2-150,1,0,10)
  my ($x1) = shift;
  my ($x2) = shift;
  my ($y) = shift;
  my ($position) = shift; # vertical = 2, horizontal = 1, else diagonal
  my ($j) = shift; # justification: 0 = Left, 1 = Center, 2 = Right
  my ($f) = shift; # font: 0 = Times Roman, 1 = Times Italic, 2 = Times Bold
  my ($fs) = shift; # font size
  my ($c) = shift; # color (default = 0 or black defined below)
  my ($mid);
  my ($rotation) = 0.7854; # 0.0000 (no rotation)

  ($x1, $x2) = (&min($x1,$x2), &max($x1,$x2));
  #if ($j == "1") {
  $mid = $x1 + (($x2 - $x1)/2);
  #}
  #else  {
  if ($fs == "40") { # printing the fragment_id
    $mid = $x1;
  }
  #}

  if (!defined($c))  { $c = 0;}
  if ($position == 1) {
    # normal
    # code:  object, sub_type(0=left,1=center,2=right), color, depth, pen_style, font, font_size, angle, font_flags, height, length, x, y, string
    #printf("4 $j $c 100 0 $f $fs 0.0000 4 270 675 %d %d %s\\001\n",
    printf("4 $j $c 100 0 $f $fs 0.0000 4 495 3120 %d %d %s\\001\n",
	   $mid,
	   $y,
	   $n);
  }
  elsif ($position == 2) { # vertical
    $rotation = 1.570796327;
    #printf("4 $j $c 100 0 $f $fs %lf 4 270 675 %d %d %s\\001\n",
    printf("4 $j $c 100 0 $f $fs %lf 4 495 3120 %d %d %s\\001\n",
	   $rotation,
	   $mid,
	   $y,
	   $n);
  }
  else {
    #$y = $y - ($global_struct->{'gene_height'} / 2);
    #printf("4 0 $c 100 0 $f $fs %lf 4 270 675 %d %d %s\\001\n",
    printf("4 0 $c 100 0 $f $fs %lf 4 495 3120 %d %d %s\\001\n",
	   $rotation,
	   $mid,
	   $y,
	   $n);
  }
}

sub passes_length_test {
  my ($global_struct, $x, $y) = @_;
  my ($r) = 0;

  if ($global_struct->{'short_gene_toggle'} == 1) { # == 1 means perform the test
    ($x, $y) = ($y, $x) if ($x > $y);
    $r = 1 if (($y - $x + 1) > $global_struct->{'short_gene_cutoff'});
  }
  else {
    $r = 1;
  }
  return($r);
}

sub choose_name {
  my ($global_struct, $type, $orf_name, $sym, $locus) = @_;
  my $name;

  $name = $orf_name;

  if (($global_struct->{'locus_toggle'} == 1) && ($locus)) {
    $name = $locus;
  }
  elsif (($global_struct->{'locus_toggle'} == 1) && (!$locus))  { # if not locus_name provided, do not print anything
    $name = "";
  }

  if ($global_struct->{'truncate_name_toggle'} == 1) {
    if ($type eq "tRNA") {
      $name =~ s/^tRNA-//;
    }
    elsif ($type eq "rRNA") {
      $name =~ s/^[a-zA-Z]*//;
    }
    elsif ($type eq "sRNA") {
      if ($name =~ /tmRNA/) { $name = "tmRNA"; }
      elsif ($name =~ /rnpB|RNaseP/) { $name = "RNaseP"; }
      elsif ($name =~ /7S/) { $name = "7S"; }
    }
    elsif ($type eq "ENZ") {
      $name = $name;
    }
    else  {
      $name =~ s/^[A-Z0]*//;
    }
  }
  if ($global_struct->{'gene_sym_toggle'} == 1 and $sym) {
    $name = $sym;
  }
  return($name);
}

sub handle_invocation {
  my ($go) = @_;
  my ($prog, $invocation, $pwd, $invoke_file);

  $prog = $0;
  $prog =~ s/.*\///;

  $invocation = $prog . " @ARGV";
  $pwd = $ENV{'PWD'};
  $invoke_file = $pwd . "/." . $prog;

  if(-e $invoke_file && $go) {
    open(F_IN,"$invoke_file") || die "can't open $invoke_file";
    while(<F_IN>) {
      print STDERR "last invocation: $_";
    }
    close(F_IN);
    exit(1);
  }

  open(F_OUT,">$invoke_file") || die "coulnt write $invoke_file";
  print F_OUT "$invocation\n";
  close(F_OUT);
}

sub fig_header {
  my ($global_struct, $rgb, $total_colors, $color_match, $color_start) = @_;
  my ($i);

  print "\#FIG 3.2\n";
  print "Portrait\n";	# Portrait, Landscape
  print "Center\n";	# Center, Flush Left
  print "Inches\n";	# Metric, Inches
  print "Letter\n";	# Letter, Legal, Ledger, Tabloid, A, B, C, D, E, A4, A3, A2, A1,
  # A0, B5
  print "100.00\n";	# magnification, percentage
  print "Single\n";	# Single, Multiple (pages)
  print "-2\n";		# -3:background, -2:none, -1:default (transparent color)
  print "1200 2\n";	# units/inch (resolution) 1: lower left (ignored) 2: upper left (coordinate orig)in 
  for($i=$global_struct->{'color_start'};length($rgb->{$i}) != 0;$i++) {
    printf("0 %d \#%s\n",$i,$rgb->{$i});
  }
  if ($color_match eq "s")  {
    print STDERR "----------------------------------------------- ;)\n";
    my $spectrum = &init_color_spectrum($total_colors); # fill hash global %spectrum with rgb colors from 71 up
    for(my $i=0, my $j = $color_start; length($spectrum->{$i}) != 0; $i++, $j++) {
      print STDERR "$spectrum->{$i}\n";
      printf("0 %d \#%s\n",$j,$spectrum->{$i});
    }
  }
}

sub overlap {
  my ($e51, $e31, $e52, $e32) = @_;
  my ($r);

  # find if e51|--------------|e31
  #              e52|--------|e32
  #

  ($e51, $e31) = ($e31, $e51) if ($e31 < $e51);
  ($e52, $e32) = ($e32, $e52) if ($e32 < $e52);

  $r = 0;

  if (($e52 >= $e51) &&
      ($e52 <= $e31) &&
      ($e32 >= $e51) &&
      ($e32 <= $e31)) {
    $r = 1;
  }
  return($r);
}

sub line2coord {
  my ($global_struct, $y) = @_;
  $y = ($y * $global_struct->{'line_spacing'}) + $global_struct->{'y_border'};
  return(int($y));
}

sub x_transform {
  my ($global_struct, $x, $level_ref) = @_;
  my ($width);
  my $level = $level_ref->{'level'};

  if ($level_ref->{'orientation'} == $global_struct->{'forward'}) {
    $x = ($x - $level_ref->{'x_decrement'}) + $global_struct->{'x_border'}+$global_struct->{$level}+$level_ref->{'lmargin'};
  }
  else {
    $x = ($level_ref->{'width'} - ($x - $level_ref->{'x_decrement'})) + $global_struct->{'x_border'}+$global_struct->{$level}+$level_ref->{'lmargin'};
  }
  $x = 1 if ($x < 1);
  return(int($x));
}

sub max {
  my ($x, $y) = @_;
  return ($x >= $y) ? $x :$y;
}

sub min {
  my ($x, $y) = @_;
  return ($x < $y) ? $x :$y;
}
1;
