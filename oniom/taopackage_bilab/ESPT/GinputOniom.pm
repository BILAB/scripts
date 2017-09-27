package ESPT::GinputOniom;

use base qw(ESPT::ESSfile);
use ESPT::Glib 0.01;
use strict;
use warnings;

our $VERSION = '0.01';

=head1 NAME

ESPT::GinputOniom - Oniom object to treat Gaussian input file.

=head1 SYNOPSIS

   use ESPT::GinputOniom;

   my $ginputoniom = GinputOniom->new();

=head1 DESCRIPTION

This module provides methods to quickly access data and structures contianed in a Gaussian input file
object for ONIOM claculation.  

=cut


### Version History ###
# 0.01	Read in ONIOM Gaussian input file.
#



### To Do List ###

##### NOTE #####
#
################

## the object constructor **


sub new {
	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my $ginputoniom = ESPT::ESSfile->new();

	$ginputoniom->{PROGRAM} = "GAUSSIAN";
	$ginputoniom->{TYPE} = "Input";

	# Link 0 & Route commands
	$ginputoniom->{ROUTE} = undef;
	$ginputoniom->{KEYWORDS} = [];
        $ginputoniom->{OUTPUTFILE} = "";


	# molecular info
	$ginputoniom->{ATOMLIST} = [];
        $ginputoniom->{RESIDUE} = [];
        $ginputoniom->{RESIDUENUM} = [];
        $ginputoniom->{ATOMTYPE} = [];
        $ginputoniom->{LAYERLIST} = [];
        $ginputoniom->{MOVESYMBOL} = [];
        $ginputoniom->{ELEMENTNUM} = [];
        $ginputoniom->{ISLINKATOM} = [];
        $ginputoniom->{LINKATOMTYPE} = [];
        $ginputoniom->{LINKATOMINFO} = [];
        $ginputoniom->{LINKEDATOMNUM} = [];
        $ginputoniom->{ELEMENT} = [];
        $ginputoniom->{PARTIALCHARG} = [];
        $ginputoniom->{CARTCOORD} = [];
        $ginputoniom->{ATOMDEFINE} = [];
        $ginputoniom->{MISSPARMBONDLIST} = [];
        $ginputoniom->{MISSPARMANGLELIST} = [];

	$ginputoniom->{CHARGE} = undef;
	$ginputoniom->{EELEC} = undef;		# SCF electronic energy
	$ginputoniom->{MULTIPLICITY} = undef;
	$ginputoniom->{NATOMS} = undef;
        
        $ginputoniom->{HCHARGE} = 0.0;
        $ginputoniom->{MCHARGE} = 0.0;
        $ginputoniom->{LCHARGE} = 0.0;
        $ginputoniom->{RCHARGE} = 0.0;

        $ginputoniom->{XDIPOLECAL} = 0.0;
        $ginputoniom->{YDIPOLECAL} = 0.0;
        $ginputoniom->{ZDIPOLECAL} = 0.0;
        $ginputoniom->{TOTDIPOLECAL} = 0.0;

        $ginputoniom->{RAW} = 0;
        
	bless($ginputoniom, $class);
	return $ginputoniom;
}


## methods ##

# set filename & spin then digest the file
sub analyze : method {
	my $ginputoniom = shift;
	$ginputoniom->prepare(@_);
	$ginputoniom->digest();
	return;
}


sub digest {

my $ginputoniom = shift;

my $debug = $ginputoniom->{DEBUG} ;
my $raw = $ginputoniom->{RAW};

my $atomflag = 0;
my $atomtot = 0;
my $cartflag = 0;
my $carttot = 0;
my $col = 0;
my $counter = 0;
my $geomcount = 0;
my $row = 0;
my $rparsed = 0;
my $Titleflag = 0;
my $Cflag = 0;
my $weightflag = 0;
my ($curatmtyp, $amberffdum, $curatomname, $curresidname, $curelenum, $curparchg, $curX);
my ($curY,  $curZ, $atomdefined, $line1, $line2, $templine, $isendofatomlist,$structindicator);
my ($atomindex);


$structindicator= 0 ;
$atomindex = 0;

 open(ONIOMINPUTFILE,$ginputoniom->{FILENAME}) || die "Could not read $ginputoniom->{FILENAME}\n$!\n";
 while (<ONIOMINPUTFILE>) {

   if   ($structindicator == 0){
     if (/^\s*-*\d+\s+\d+\s+-*\d+\s+\d+/) {
#     if (/^-*[0-9] [0-9]/) {
       print $_ if $debug >= 1;
       $structindicator = 1;
       next;
     } # if (/^\s*-*\d+\s+\d+\s+-*\d+\s+\d+/) 
     next;
   }  # if   ($structindicator ==0){
   else {
     if (/^\s+$/ or /^$/){ last;}
#     print $_;
     $atomindex++;
#     printf "%d10\n",$atomindex;
#     if ( /^\s+[a-zA-Z\-]+\d+\.\d+\s+-*[01]\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+H/ ) {
#  If the atom is QM region (or allowed to move when -g is turned on), select this atoms.
     if ( /^\s+([a-zA-Z]{1,})-([a-zA-Z0-9\*]{1,})-(-*\d+\.\d+)\s+(-*[01])\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+([HLMhlm])/)
      {
#      print $_;

      if ($ginputoniom->{DEBUG} == 1) {  print "$1\t$2\t$3\t$4\t$5\t$6\t$7\t$8\t$atomindex\n"; }
       $ginputoniom->{ELEMENT}[$atomindex]=$1;
       $ginputoniom->{ATOMTYPE}[$atomindex]=$2;
       $ginputoniom->{PARTIALCHARG}[$atomindex]=$3;
       $ginputoniom->{MOVESYMBOL}[$atomindex]=$4;
       $ginputoniom->{CARTCOORD}[3*$atomindex]=$5;
       $ginputoniom->{CARTCOORD}[3*$atomindex+1]=$6;
       $ginputoniom->{CARTCOORD}[3*$atomindex+2]=$7;
       $ginputoniom->{LAYERLIST}[$atomindex]=$8;
       $ginputoniom->{ISLINKATOM}[$atomindex] = 0;
       $ginputoniom->{XDIPOLECAL} += $5 * $3 * $ginputoniom->{DEBYECONVER};
       $ginputoniom->{YDIPOLECAL} += $6 * $3 * $ginputoniom->{DEBYECONVER};
       $ginputoniom->{ZDIPOLECAL} += $7 * $3 * $ginputoniom->{DEBYECONVER};

# Add up partial charges

       if ( uc($8) eq "L") {$ginputoniom->{LCHARGE} += $3;}
       if ( uc($8) eq "H") {$ginputoniom->{HCHARGE} += $3;}
       if ( uc($8) eq "M") {$ginputoniom->{MCHARGE} += $3;}
       $ginputoniom->{RCHARGE} += $3;

      # If this is a link atom

      if (/^\s+[a-zA-Z]{1,2}-[a-zA-Z0-9]{1,2}--*\d+\.\d+\s+-*[01]\s+-*\d+\.\d+\s+-*\d+\.\d+\s+-*\d+\.\d+\s+[HLMhlm]\s+([a-zA-Z]-\S+)\s+(\d+)/) {
         $ginputoniom->{ISLINKATOM}[$atomindex] = 1;
         $ginputoniom->{LINKATOMINFO}[$atomindex] = $1; 
         $ginputoniom->{LINKEDATOMNUM}[$atomindex] = $1; 
         print " Atom number $atomindex is a link atom info($1).\n" if $debug >= 1;
         
        } # if $line 1

#      push (@qmatomlist, $atomindex);
#      $qmatomnumber++;
      }  # if
   
   } # else
 
 } # while ONIOMINPUTFILE

$ginputoniom->{TOTDIPOLECAL} = sqrt($ginputoniom->{XDIPOLECAL}**2+$ginputoniom->{YDIPOLECAL}**2+$ginputoniom->{ZDIPOLECAL}**2);

} # digest


# Write out total charge

sub writecharges{

 my $ginputoniom = shift;
 my $debug = $ginputoniom->{DEBUG};
 my $raw = $ginputoniom->{RAW};
 
 if ($raw == 0) {
  printf "Total charge of real system is %10.6f.\n",$ginputoniom->{RCHARGE};
  printf "Total charge of high layer is %10.6f.\n",$ginputoniom->{HCHARGE};
  printf "Total charge of medium layer is %10.6f.\n",$ginputoniom->{MCHARGE};
  printf "Total charge of low layer is %10.6f.\n",$ginputoniom->{LCHARGE};
  printf "Total charge of high plus medium layer is %10.6f.\n",$ginputoniom->{HCHARGE}+$ginputoniom->{MCHARGE};
 } else {
  print "Total charge of real system is $ginputoniom->{RCHARGE}\n",
  print "Total charge of high layer is $ginputoniom->{HCHARGE}\n";
  print "Total charge of medium layer is $ginputoniom->{MCHARGE}\n";
  print "Total charge of low layer is $ginputoniom->{LCHARGE}\n";
  print "Total charge of high plus medium layer is ", $ginputoniom->{HCHARGE}+$ginputoniom->{MCHARGE},"\n";
 }

} # sub outputparm{


# Calculate dipolemoment

sub dipolemoment{

 my $ginputoniom = shift;
 my $debug = $ginputoniom->{DEBUG};
 my $raw = $ginputoniom->{RAW};
 
 if ($raw == 0) {
  printf "\nDipole moment (Debye) (X, Y, Z) is \n(%15.4f, %15.4f, %15.4f).\n",$ginputoniom->{XDIPOLECAL},
                                         $ginputoniom->{YDIPOLECAL},$ginputoniom->{ZDIPOLECAL};
  printf "Total Dipole moment (Debye) is %15.4f.\n",$ginputoniom->{TOTDIPOLECAL};
 } else {
  print "\nDipole moment (X, Y, Z) is \n(",$ginputoniom->{XDIPOLECAL},", ", $ginputoniom->{YDIPOLECAL},", ", $ginputoniom->{ZDIPOLECAL},").\n";
  print "Total Dipole moment is $ginputoniom->{TOTDIPOLECAL}\n";
 
 
 }

} # sub outputparm{


1;
__END__

=head1 VERSION

0.01

=head1 AUTHOR

Peng Tao, E<lt>tao.21@osu.edu<gt>

=head1 COPYRIGHT

Copyright (c) 2009~2010 by Peng Tao


=cut

