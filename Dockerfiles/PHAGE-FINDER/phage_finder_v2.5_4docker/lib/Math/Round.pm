package Math::Round;

use strict;
use POSIX;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

require Exporter;

@ISA = qw(Exporter AutoLoader);
@EXPORT = qw(round nearest);
@EXPORT_OK = qw(round nearest round_even round_odd round_rand
   nearest_ceil nearest_floor nearest_rand
   nlowmult nhimult );
$VERSION = '0.05';

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

#--- Determine what value to use for "one-half".  Because of the
#--- perversities of floating-point hardware, we must use a value
#--- slightly larger than 1/2.  We accomplish this by determining
#--- the bit value of 0.5 and increasing it by a small amount in a
#--- lower-order byte.  Since the lowest-order bits are still zero,
#--- the number is mathematically exact.

my $halfhex = unpack('H*', pack('d', 0.5));
if (substr($halfhex,0,2) ne '00' && substr($halfhex, -2) eq '00') {
   #--- It's big-endian.
   substr($halfhex, -4) = '1000';
} else {
   #--- It's little-endian.
   substr($halfhex, 0,4) = '0010';
}

my $half = unpack('d',pack('H*', $halfhex));

sub round {
 my $x;
 my @res = ();
 foreach $x (@_) {
   if ($x >= 0) {
      push @res, POSIX::floor($x + $half);
   } else {
      push @res, POSIX::ceil($x - $half);
   }
 }
 return (wantarray) ? @res : $res[0];
}

sub round_even {
 my $x;
 my @res = ();
 foreach $x (@_) {
   my ($sign, $in, $fr) = _sepnum($x);
   if ($fr == 0.5) {
      push @res, $sign * (($in % 2 == 0) ? $in : $in + 1);
   } else {
      push @res, $sign * POSIX::floor(abs($x) + $half);
   }
 }
 return (wantarray) ? @res : $res[0];
}

sub round_odd {
 my $x;
 my @res = ();
 foreach $x (@_) {
   my ($sign, $in, $fr) = _sepnum($x);
   if ($fr == 0.5) {
      push @res, $sign * (($in % 2 == 1) ? $in : $in + 1);
   } else {
      push @res, $sign * POSIX::floor(abs($x) + $half);
   }
 }
 return (wantarray) ? @res : $res[0];
}

sub round_rand {
 my $x;
 my @res = ();
 foreach $x (@_) {
   my ($sign, $in, $fr) = _sepnum($x);
   if ($fr == 0.5) {
      push @res, $sign * ((rand(4096) < 2048) ? $in : $in + 1);
   } else {
      push @res, $sign * POSIX::floor(abs($x) + $half);
   }
 }
 return (wantarray) ? @res : $res[0];
}

#--- Separate a number into sign, integer, and fractional parts.
#--- Return as a list.
sub _sepnum {
 my $x = shift;
 my ($sign, $i);
 $sign = ($x >= 0) ? 1 : -1;
 $x = abs($x);
 $i = int($x);
 return ($sign, $i, $x - $i);
}

#------ "Nearest" routines (round to a multiple of any number)

sub nearest {
 my ($targ, @inputs) = @_;
 my @res = ();
 my $x;

 $targ = abs($targ) if $targ < 0;
 foreach $x (@inputs) {
   if ($x >= 0) {
      push @res, $targ * int(($x + $half * $targ) / $targ);
   } else {
      push @res, $targ * POSIX::ceil(($x - $half * $targ) / $targ);
   }
 }
 return (wantarray) ? @res : $res[0];
}

# In the next two functions, the code for positive and negative numbers
# turns out to be the same.  For negative numbers, the technique is not
# exactly obvious; instead of floor(x+0.5), we are in effect taking
# ceiling(x-0.5).

sub nearest_ceil {
 my ($targ, @inputs) = @_;
 my @res = ();
 my $x;

 $targ = abs($targ) if $targ < 0;
 foreach $x (@inputs) {
    push @res, $targ * POSIX::floor(($x + $half * $targ) / $targ);
 }
 return (wantarray) ? @res : $res[0];
}

sub nearest_floor {
 my ($targ, @inputs) = @_;
 my @res = ();
 my $x;

 $targ = abs($targ) if $targ < 0;
 foreach $x (@inputs) {
    push @res, $targ * POSIX::ceil(($x - $half * $targ) / $targ);
 }
 return (wantarray) ? @res : $res[0];
}

sub nearest_rand {
 my ($targ, @inputs) = @_;
 my @res = ();
 my $x;

 $targ = abs($targ) if $targ < 0;
 foreach $x (@inputs) {
   my ($sign, $in, $fr) = _sepnear($x, $targ);
   if ($fr == 0.5 * $targ) {
      push @res, $sign * $targ * ((rand(4096) < 2048) ? $in : $in + 1);
   } else {
      push @res, $sign * $targ * int((abs($x) + $half * $targ) / $targ);
   }
 }
 return (wantarray) ? @res : $res[0];
}

#--- Next lower multiple
sub nlowmult {
 my ($targ, @inputs) = @_;
 my @res = ();
 my $x;

 $targ = abs($targ) if $targ < 0;
 foreach $x (@inputs) {
    push @res, $targ * POSIX::floor($x / $targ);
 }
 return (wantarray) ? @res : $res[0];
}

#--- Next higher multiple
sub nhimult {
 my ($targ, @inputs) = @_;
 my @res = ();
 my $x;

 $targ = abs($targ) if $targ < 0;
 foreach $x (@inputs) {
    push @res, $targ * POSIX::ceil($x / $targ);
 }
 return (wantarray) ? @res : $res[0];
}

#--- Separate a number into sign, "integer", and "fractional" parts
#--- for the 'nearest' calculation.  Return as a list.
sub _sepnear {
 my ($x, $targ) = @_;
 my ($sign, $i);
 $sign = ($x >= 0) ? 1 : -1;
 $x = abs($x);
 $i = int($x / $targ);
 return ($sign, $i, $x - $i*$targ);
}

1;

__END__

=head1 NAME

Math::Round - Perl extension for rounding numbers

=head1 SYNOPSIS

  use Math::Round qw(...those desired... or :all);

  $rounded = round($scalar);
  @rounded = round(LIST...);
  $rounded = nearest($target, $scalar);
  @rounded = nearest($target, LIST...);

  # and other functions as described below

=head1 DESCRIPTION

B<Math::Round> supplies functions that will round numbers in different
ways.  The functions B<round> and B<nearest> are exported by
default; others are available as described below.  "use ... qw(:all)"
exports all functions.

=head1 FUNCTIONS

=over 2

=item B<round> LIST

Rounds the number(s) to the nearest integer.  In scalar context,
returns a single value; in list context, returns a list of values.
Numbers that are halfway between two integers are rounded
"to infinity"; i.e., positive values are rounded up (e.g., 2.5
becomes 3) and negative values down (e.g., -2.5 becomes -3).

=item B<round_even> LIST

Rounds the number(s) to the nearest integer.  In scalar context,
returns a single value; in list context, returns a list of values.
Numbers that are halfway between two integers are rounded to the
nearest even number; e.g., 2.5 becomes 2, 3.5 becomes 4, and -2.5
becomes -2.

=item B<round_odd> LIST

Rounds the number(s) to the nearest integer.  In scalar context,
returns a single value; in list context, returns a list of values.
Numbers that are halfway between two integers are rounded to the
nearest odd number; e.g., 3.5 becomes 3, 4.5 becomes 5, and -3.5
becomes -3.

=item B<round_rand> LIST

Rounds the number(s) to the nearest integer.  In scalar context,
returns a single value; in list context, returns a list of values.
Numbers that are halfway between two integers are rounded up or
down in a random fashion.  For example, in a large number of trials,
2.5 will become 2 half the time and 3 half the time.

=item B<nearest> TARGET, LIST

Rounds the number(s) to the nearest multiple of the target value.
TARGET must be positive.
In scalar context, returns a single value; in list context, returns
a list of values.  Numbers that are halfway between two multiples
of the target will be rounded to infinity.  For example:

  nearest(10, 44)    yields  40
  nearest(10, 46)            50
  nearest(10, 45)            50
  nearest(25, 328)          325
  nearest(.1, 4.567)          4.6
  nearest(10, -45)          -50

=item B<nearest_ceil> TARGET, LIST

Rounds the number(s) to the nearest multiple of the target value.
TARGET must be positive.
In scalar context, returns a single value; in list context, returns
a list of values.  Numbers that are halfway between two multiples
of the target will be rounded to the ceiling, i.e. the next
algebraically higher multiple.  For example:

  nearest_ceil(10, 44)    yields  40
  nearest_ceil(10, 45)            50
  nearest_ceil(10, -45)          -40

=item B<nearest_floor> TARGET, LIST

Rounds the number(s) to the nearest multiple of the target value.
TARGET must be positive.
In scalar context, returns a single value; in list context, returns
a list of values.  Numbers that are halfway between two multiples
of the target will be rounded to the floor, i.e. the next
algebraically lower multiple.  For example:

  nearest_floor(10, 44)    yields  40
  nearest_floor(10, 45)            40
  nearest_floor(10, -45)          -50

=item B<nearest_rand> TARGET, LIST

Rounds the number(s) to the nearest multiple of the target value.
TARGET must be positive.
In scalar context, returns a single value; in list context, returns
a list of values.  Numbers that are halfway between two multiples
of the target will be rounded up or down in a random fashion.
For example, in a large number of trials, C<nearest(10, 45)> will
yield 40 half the time and 50 half the time.

=item B<nlowmult> TARGET, LIST

Returns the next lower multiple of the number(s) in LIST.
TARGET must be positive.
In scalar context, returns a single value; in list context, returns
a list of values.  Numbers that are between two multiples of the
target will be adjusted to the nearest multiples of LIST that are
algebraically lower. For example:

  nlowmult(10, 44)    yields  40
  nlowmult(10, 46)            40
  nlowmult(25, 328)          325
  nlowmult(.1, 4.567)          4.5
  nlowmult(10, -41)          -50

=item B<nhimult> TARGET, LIST

Returns the next higher multiple of the number(s) in LIST.
TARGET must be positive.
In scalar context, returns a single value; in list context, returns
a list of values.  Numbers that are between two multiples of the
target will be adjusted to the nearest multiples of LIST that are
algebraically higher. For example:

  nhimult(10, 44)    yields  50
  nhimult(10, 46)            50
  nhimult(25, 328)          350
  nhimult(.1, 4.512)          4.6
  nhimult(10, -49)          -40

=back

=head1 STANDARD FLOATING-POINT DISCLAIMER

Floating-point numbers are, of course, a rational subset of the real
numbers, so calculations with them are not always exact.  In order to
avoid surprises because of this, these routines use a value for
one-half that is very slightly larger than 0.5.  Nevertheless,
if the numbers to be rounded are stored as floating-point, they will
be subject, as usual, to the mercies of your hardware, your C
compiler, etc.  Thus, numbers that are supposed to be halfway between
two others may be stored in a slightly different way and thus behave
surprisingly.

=head1 AUTHOR

Math::Round was written by Geoffrey Rommel E<lt>GROMMEL@cpan.orgE<gt>
in October 2000.

=cut
