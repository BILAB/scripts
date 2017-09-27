package ESPT::AMBERFF;

our $VERSION = '0.01';

### Version History ###
# 0.01	Read in AMBER library file and prep file
#

=head1 NAME

ESPT::AMBERFF - AMBER Force Field for Atom type and partial charges lookup

=head1 SYNOPSIS

    use ESPT::AMBERFF;


=head1 DESCRIPTION

This module read in AMBER library file (lib file) and given prep files for 
atom type and partial charges lookup.

=head2 Partial Charges and Atom Type look up

When passing residue name and atom name, this object will be able to return the
atom type and partial charges for this specific atom.

=cut

=head1 ATTRIBUTES

=over 15

=item B<$file-E<gt>new()>

Creates a new AMBERFF object

=cut

use strict;
use warnings;

our %ffatmtyp;
our %ffelenum;
our %ffparchg;
our @ffresiduelist;

our $amberlibfile = "/Users/masa/src/taopackage/ESPT/amberff/all_amino03.lib";
our $prepfolder = "/Users/masa/src/taopackage/ESPT/prepfiles";

my $curresid;
my $endresidue;
my $curatmnam;
my $curatmtyp;
my $curelenum;
my $curparchg;
my $prepfile;
my $prepfilewpath;
my $line1;
my $line2;
my $templine;

my $debug = 4;

 print "Open file $amberlibfile to read in AMBER residue template.\n" if $debug >= 5;

open(AMBERFFFILE,$amberlibfile) || die "Could not read AMBEF library file $amberlibfile\n$!\n";

# Reading AMBER Residue, Atom, Partial charges from AMBER library file $amberlibfile

while (<AMBERFFFILE>) {
  if (/^!entry.(\S\S\S).unit.atoms\s+table/) {

    $curresid = $1; 
    push (@ffresiduelist,$curresid);

    # add current residue as one hash element 
#    print "Current Residue is $1.\n";
    $endresidue=0;
    do {
       $_ = <AMBERFFFILE>;
       if (!/\s+"(\S+)"\s+"(\S+)"\s+-*\d+\s+-*\d+\s+-*\d+\s+-*\d+\s+(\d+)\s+(-*\d+.\d+)/) {$endresidue = 1;}
         else {

#           print "Current atom name is $1, atom type is $2, element number is $3, partial charge is $4.\n";

           $curatmnam = $1;
           $curatmtyp = $2;
           $curelenum = $3,
           $curparchg = $4;
           $ffatmtyp{$curresid}{$curatmnam} =  $curatmtyp ;
           $ffelenum{$curresid}{$curatmnam}  =  $curelenum ;
           $ffparchg{$curresid}{$curatmnam} =  $curparchg ;
           
#           print "Atom with name $curatmnam in residue $curresid has type $ffatmtyp{$curresid}{$curatmnam},";
#           print " element number $ffelenum{$curresid}{$curatmnam}, partial charges $ffparchg{$curresid}{$curatmnam}.\n";

          } # else
    
     } while ($endresidue == 0); # do
    
    
  
   } #  if (/^!entry.\S\S\S.unit.atoms\s+table/) 
   

} # while (<AMBERFFFILE>) 

close (AMBERFFFILE);

# Almost every ONIOM calculation has substrate(s) which is (are) not amino acid. 
# If prepared by AMBER, user will have prep files for these molecules.
# The following code can read these prep files from folder $prepfolder, and 
# take atom name, atom type and partial charges information.
# The code assumes each prep file has this format:
#    0    0    2
#
#This is a remark line
#A5P_dianion.res
#AP2    INT  0
#CORRECT     OMIT DU   BEG
#  0.0000
#   1  DUMM  DU    M    0  -1  -2     0.000      .0        .0      .00000
#   2  DUMM  DU    M    1   0  -1     1.449      .0        .0      .00000
#   3  DUMM  DU    M    2   1   0     1.522   111.1        .0      .00000
#   4  O1    o     M    3   2   1     1.540   111.208   180.000  -0.58161  8
#   5  C1    c     M    4   3   2     1.193    54.434   109.126   0.54347  6
#   6  H1    h4    E    5   4   3     1.093   120.987  -156.388   0.00162  1
#... 
# Please note that original prep file does not have element number. To be safe,
#  user should add element number for each atom at the end of line. In this way
#  it is guaranteed that correct atom symole will be put in the output file

opendir(PREPDIR,$prepfolder) || die "Could not open prep file folder $prepfolder\n$!\n";

foreach $prepfile (readdir(PREPDIR)) {
 if ( $prepfile =~ /^\./ ) {next;}
 print "Reading prep file $prepfile ...\n" if $debug >= 5;
 
 $prepfilewpath = $prepfolder."/".$prepfile;

# print "$prepfilewpath\n";

 open(PREPFILE,$prepfilewpath) || die "Could not read prep file $prepfile in directory $prepfolder\n$!\n";
 
 $line1 = "";
 $line2 = "";
 $templine = "";
 
 while ($templine = <PREPFILE>) {
   $line1 = $line2; 
   $line2 = $templine;
   
   if ( $line2 =~ /^\s*CORR/) {
     
     $line1 =~ /^(\S\S\S)\s+/;
     
    # add current residue as one hash element 
#    print "Current Residue is $1.\n";

     $curresid = $1; 
     push (@ffresiduelist,$curresid);


     # Skip 0.000 and first three DUMM lines after CORR line
     $templine = <PREPFILE>; 
     $templine = <PREPFILE>; 
     $templine = <PREPFILE>; 
     $templine = <PREPFILE>; 
     
     # clear $line1 and $line2;
     $line1 = "";
     $line2 = "";
     
     $endresidue=0;
     do {
        $_ = <PREPFILE>;
        if (!/\s+\d+\s+(\S+)\s+(\S+)\s+\S+\s+-*\d+\s+-*\d+\s+-*\d+\s+-*\d+.\d+\s+-*\d+.\d+\s+-*\d+.\d+\s+(-*\d+.\d+)\s+(\d+)/) {$endresidue = 1;}
          else {

#           print "Current atom name is $1, atom type is $2, element number is $3, partial charge is $4.\n";

            $curatmnam = $1;
            $curatmtyp = uc($2);
            $curparchg = $3;
            $curelenum = $4,
            $ffatmtyp{$curresid}{$curatmnam} =  $curatmtyp ;
            $ffelenum{$curresid}{$curatmnam}  =  $curelenum ;
            $ffparchg{$curresid}{$curatmnam} =  $curparchg ;
           
           print "Atom with name $curatmnam in residue $curresid has type $ffatmtyp{$curresid}{$curatmnam}," if $debug >= 5;
           print " element number $ffelenum{$curresid}{$curatmnam}, partial charges $ffparchg{$curresid}{$curatmnam}.\n" if $debug >= 5;

           } # else
    
      } while ($endresidue == 0); # do
     
   
   
   } # if ( $line2 =~ /^\s+CORR/) {
  
 
  } # while (<PREPFILE>) {

 
 close PREPFILE;
  } # foreach $prepfile (readdir(PREPDIR)) 



sub new {
 our $AMBERFF = {};
 shift;
 $AMBERFF->{DEBUG} = shift;
 
# print "AMBEREE debug is ",$AMBERFF->{DEBUG},"\n";
 
 bless ($AMBERFF);
 return $AMBERFF; 
} # new

# Subroutine atomtype takes atom type for current atom
sub atomtype {
  my $curresidue;
  my $curatomname;
  my $curatomtype;
  
  shift;
  $curresidue = shift;
  $curatomname = shift;
#  print "Atom type inquary for residue $curresidue and atom $curatomname:\n";
  $curatomtype = $ffatmtyp{$curresidue}{$curatomname};
  
  return $curatomtype;

} # atometype


# Subroutine elementnumber takes element number for current atom

sub elementnumber {
  my $curresidue;
  my $curatomname;
  my $curelenum;
  
  shift;
  $curresidue = shift;
  $curatomname = shift;
#  print "Element number inquary for residue $curresidue and atom $curatomname:\n";
  $curelenum = $ffelenum{$curresidue}{$curatomname};
  
  return $curelenum;

} # elementnumber 


# Subroutine partialcharge takes partial charge for current atom

sub partialcharge {
  my $curresidue;
  my $curatomname;
  my $curparchg;
  
  shift;
  $curresidue = shift;
  $curatomname = shift;
#  print "Element number inquary for residue $curresidue and atom $curatomname:\n";
  $curparchg = $ffparchg{$curresidue}{$curatomname};
  
  return $curparchg;

} # partialcharge 


# Subroutine existsresidue checks if the given residue exists or not.

sub existsresidue {
  my $curresidue;
  my $currele;
  my $exists = 0;
  
  shift;
  $curresidue = shift;
  
  foreach $currele (@ffresiduelist) {
     if ( uc($currele) eq uc($curresidue)) {
      $exists = 1;
      last;
     } 
  }

 return $exists;
  
} # existsresidue


1;

__END__

=back

=head1 VERSION

0.01

=head1 SEE ALSO

F<ESPT::PDB>

=head1 AUTHOR

Peng Tao, E<lt>tao.21@osu.eduE<gt>

=head1 COPYRIGHT

Copyright (c) 2009~2010 by Peng Tao


=cut
