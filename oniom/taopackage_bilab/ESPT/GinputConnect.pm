package ESPT::GinputConnect;

use base qw(ESPT::ESSfile);
use ESPT::Glib 0.01;
use strict;
use warnings;

our $VERSION = '0.01';

=head1 NAME

ESPT::GinputConnect - Object to search isolated atoms and atoms with multiple bonds in Gaussian input file.

=head1 SYNOPSIS

   use ESPT::GinputConnect;

   my $ginputconnect = GinputConnect->new();

=head1 DESCRIPTION

This module provides methods to quicklysearch isolated atoms and atoms with multiple bonds in Gaussian input file.

=cut


### Version History ###
# 0.01	Treat Gaussian input file.
#



### To Do List ###

##### NOTE #####
#
################

## the object constructor **


sub new {
	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my $ginputconnect = ESPT::ESSfile->new();

	$ginputconnect->{PROGRAM} = "GAUSSIAN";
	$ginputconnect->{TYPE} = "Input";

	# Link 0 & Route commands
	$ginputconnect->{ROUTE} = undef;
	$ginputconnect->{KEYWORDS} = [];
        $ginputconnect->{OUTPUTFILE} = "";


	# connectivity info
	$ginputconnect->{ATOMLIST} = [];
        $ginputconnect->{RESIDUE} = [];
        $ginputconnect->{RESIDUENUM} = [];
        $ginputconnect->{ATOMTYPE} = [];
        $ginputconnect->{LAYERLIST} = [];
        $ginputconnect->{MOVESYMBOL} = [];
        $ginputconnect->{ELEMENTNUM} = [];
        $ginputconnect->{ISLINKATOM} = [];
        $ginputconnect->{LINKATOMTYPE} = [];
        $ginputconnect->{LINKATOMINFO} = [];
        $ginputconnect->{LINKEDATOMNUM} = [];
        $ginputconnect->{ELEMENT} = [];
        $ginputconnect->{PARTIALCHARG} = [];
        $ginputconnect->{CARTCOORD} = [];
        $ginputconnect->{ATOMDEFINE} = [];
        $ginputconnect->{MISSPARMBONDLIST} = [];
        $ginputconnect->{MISSPARMANGLELIST} = [];
        

	$ginputconnect->{CHARGE} = undef;
	$ginputconnect->{EELEC} = undef;		# SCF electronic energy
	$ginputconnect->{MULTIPLICITY} = undef;
	$ginputconnect->{NATOMS} = undef;
        
        $ginputconnect->{HCHARGE} = 0.0;
        $ginputconnect->{MCHARGE} = 0.0;
        $ginputconnect->{LCHARGE} = 0.0;
        $ginputconnect->{RCHARGE} = 0.0;

        $ginputconnect->{XDIPOLECAL} = 0.0;
        $ginputconnect->{YDIPOLECAL} = 0.0;
        $ginputconnect->{ZDIPOLECAL} = 0.0;
        $ginputconnect->{TOTDIPOLECAL} = 0.0;

        # 2D array to store connectivity table 
        $ginputconnect->{CONNECTTABLE} = [];  
        
        # Store total connections number of each atom
        $ginputconnect->{CONNECTBONDTABLE} = [];
        
        # Store number of connections of each atom from its connectivity line
        $ginputconnect->{TOTALCONNECTIONSTABLE} = [];
        
        # Maximum bonds that we can handle
        $ginputconnect->{MAXBONDS} = 10;
        
        # Threshhold for connections numbers
        $ginputconnect->{Cnumber} = undef;
        
        # Flag to indicate if this atom 
        $ginputconnect->{CONNECTEDFLAG} = []; 
        
        
	bless($ginputconnect, $class);
	return $ginputconnect;
}


## methods ##

# set filename & spin then digest the file
sub analyze : method {
	my $ginputconnect = shift;
	$ginputconnect->prepare(@_);
	$ginputconnect->digest();
	return;
}


sub digest {

my $ginputconnect = shift;

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
my $debug = 0;
my ($curatmtyp, $amberffdum, $curatomname, $curresidname, $curelenum, $curparchg, $curX);
my ($curY,  $curZ, $atomdefined, $line1, $line2, $templine, $isendofatomlist,$structindicator);
my ($atomindex, $ConnetTableIndicator,$LineB4Connect, @currconnetline, $CurrLineElementNumber);
my ($currline, $i,$j,$curratom,$connectedatom, $curconnections);
my ($maxbonds,$atomWmaxbonds,$atomElenWmaxbonds, $curbonds);
my $isONIOM = 1;

$debug = $ginputconnect->{DEBUG} ;

$structindicator= 0 ;
$atomindex = 0;
$ConnetTableIndicator = 0;

 open(ONIOMINPUTFILE,$ginputconnect->{FILENAME}) || die "Could not read $ginputconnect->{FILENAME}\n$!\n";
 while (<ONIOMINPUTFILE>) {

   if   ($structindicator == 0){
#     if (/^\s*-*\d+\s+\d+\s+-*\d+\s+\d+/) {
     if (/^\s*-*\d+\s+\d+/) {
#     if (/^-*[0-9] [0-9]/) {
       print $_ if $debug >= 1;
       $structindicator = 1;
       next;
     } # if (/^\s*-*\d+\s+\d+/) {
     next;
   }  # if   ($structindicator ==0){
   else {
     if (/^\s+$/ or /^$/){ last;}
#     print $_;
     $atomindex++;
#     printf "%d10\n",$atomindex;
#     if ( /^\s*[a-zA-Z\-]+\d+\.\d+\s+-*[01]\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+H/ ) {
#  If the atom is QM region (or allowed to move when -g is turned on), select this atoms.
     if ( /^\s*([a-zA-Z]+)-([a-zA-Z0-9\*]+)-(-*\d+\.\d+)\s+(-*[01])\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+([HLMhlm])/)
      {
       $isONIOM = 1;
       
#      print $_;

      if ($ginputconnect->{DEBUG} == 1) {  print "$1\t$2\t$3\t$4\t$5\t$6\t$7\t$8\t$atomindex\n"; }
       $ginputconnect->{ELEMENT}[$atomindex]=$1;
       $ginputconnect->{ATOMTYPE}[$atomindex]=$2;
       $ginputconnect->{PARTIALCHARG}[$atomindex]=$3;
       $ginputconnect->{MOVESYMBOL}[$atomindex]=$4;
       $ginputconnect->{CARTCOORD}[3*$atomindex]=$5;
       $ginputconnect->{CARTCOORD}[3*$atomindex+1]=$6;
       $ginputconnect->{CARTCOORD}[3*$atomindex+2]=$7;
       $ginputconnect->{LAYERLIST}[$atomindex]=$8;
       $ginputconnect->{ISLINKATOM}[$atomindex] = 0;
       $ginputconnect->{XDIPOLECAL} += $5 * $3 * $ginputconnect->{DEBYECONVER};
       $ginputconnect->{YDIPOLECAL} += $6 * $3 * $ginputconnect->{DEBYECONVER};
       $ginputconnect->{ZDIPOLECAL} += $7 * $3 * $ginputconnect->{DEBYECONVER};

# Add up partial charges

       if ( uc($8) eq "L") {$ginputconnect->{LCHARGE} += $3;}
       if ( uc($8) eq "H") {$ginputconnect->{HCHARGE} += $3;}
       if ( uc($8) eq "M") {$ginputconnect->{MCHARGE} += $3;}
       $ginputconnect->{RCHARGE} += $3;

      # If this is a link atom

      if (/^\s*[a-zA-Z]+-[a-zA-Z0-9\*]+--*\d+\.\d+\s+-*[01]\s+-*\d+\.\d+\s+-*\d+\.\d+\s+-*\d+\.\d+\s+[HLMhlm]\s+([a-zA-Z]-\S+)\s+(\d+)/) {
         $ginputconnect->{ISLINKATOM}[$atomindex] = 1;
         $ginputconnect->{LINKATOMINFO}[$atomindex] = $1; 
         $ginputconnect->{LINKEDATOMNUM}[$atomindex] = $1; 
         print " Atom number $atomindex is a link atom info($1).\n" if $debug >= 1;
         
        } # if $line 1

#      push (@qmatomlist, $atomindex);
#      $qmatomnumber++;
      } elsif (/^\s*[a-zA-Z]+\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)/) {  
 
      $isONIOM = 0;
 
    } # elsif
   
   } # else
 
 } # while ONIOMINPUTFILE
 
 if ($isONIOM == 1) {
     print "This is an ONIOM input file.\n";
   }
  else {
     print "This is a non-ONIOM input file.\n";
  }

$ginputconnect->{TOTDIPOLECAL} = sqrt($ginputconnect->{XDIPOLECAL}**2+$ginputconnect->{YDIPOLECAL}**2+$ginputconnect->{ZDIPOLECAL}**2);

 close(ONIOMINPUTFILE);
 
 print "There are $atomindex atoms.\n";
 
 open(ONIOMINPUTFILE,$ginputconnect->{FILENAME}) || die "Could not read $ginputconnect->{FILENAME}\n$!\n";

 $structindicator= 0;
 $ConnetTableIndicator = 0;
 
 # We got the total number of atoms as $atomindex. Now, initialize connectivity tables
 
 for ($i=1;$i<=$atomindex;$i++){
 
  $ginputconnect->{CONNECTBONDTABLE}[$i] = 0;
  $ginputconnect->{TOTALCONNECTIONSTABLE}[$i] = 0;
  $ginputconnect->{CONNECTEDFLAG}[$i] = 0;
  for ($j=1;$j<=$ginputconnect->{MAXBONDS};$j++) { $ginputconnect->{CONNECTTABLE}[$i][$j] = 0; }
 }


 while (<ONIOMINPUTFILE>) {

   if   ($structindicator == 0){

     if (/^\s*-*\d+\s+\d+/) {

       $structindicator = 1;
       next;
     } # if (/^\s*-*\d+\s+\d+/) {
     next;
   }  # if   ($structindicator ==0){
   elsif ( $ConnetTableIndicator == 0 )  {
   
   if  (/^\s+$/ or /^$/) { $ConnetTableIndicator = 1 }; 
   
   } # elsif ( $ConnetTableIndicator == 0 )  {
   
   elsif (/^\s+$/ or /^$/) { last; }
   
   else {
#    print $_;
    $currline = $_;
    @currconnetline = split;
    $CurrLineElementNumber = @currconnetline;
    print "There are $CurrLineElementNumber elements in current line. They are @currconnetline \n"  if $debug >= 4;
#    print "@currconnetline";

    if ( $CurrLineElementNumber < 1 ) {
     # There is no element in current line. Exit program here, since we don't know how this can happen.
     
     print "There is no element found in currnet line.\n";
     print "Currnet line is ($_)\n";
     die "Encounter strange situation. \nExit.\n";
     
    } # if ( $CurrLineElementNumber > 1 ) {

      elsif ($CurrLineElementNumber < 2) {
      # There is one element in current line. Nothing needs to be done.
      } # elsif ($CurrLineElementNumber < 2) {

      elsif ( $CurrLineElementNumber/2 == int($CurrLineElementNumber / 2) ) {
       # There are even number elements in current line. Exit program here, since we don't know how this can happen.
       print "There are $CurrLineElementNumber (even number) elements found in currnet line.\n";
       print "Currnet line is ($_)\n";
       die "Encounter strange situation. \nExit.\n";
      }
       
       # Read in connectivity information.
      else {
         $curconnections = ($CurrLineElementNumber-1)/2;
         $curratom = $currconnetline[0];
         
         print "Process atom $curratom.\n" if $debug >= 5;

         print "Current atom $curratom connets to $curconnections atoms:" if $debug >= 4;

         $ginputconnect->{CONNECTBONDTABLE}[$curratom] = $ginputconnect->{CONNECTBONDTABLE}[$curratom] + ($CurrLineElementNumber-1)/2;
         $ginputconnect->{TOTALCONNECTIONSTABLE}[$curratom] = ($CurrLineElementNumber-1)/2;
         $ginputconnect->{CONNECTEDFLAG}[$curratom] = 1;
         
         
       # Get 2i elements from @currconnetline (2, 4, 6, $CurrLineElementNumber-1 )
       for ($i=1;$i<= ($CurrLineElementNumber-1)/2;$i++){
         $connectedatom = @currconnetline[2*$i-1];
         $ginputconnect->{CONNECTTABLE}[$curratom][$i] = $connectedatom;
         $ginputconnect->{CONNECTBONDTABLE}[$connectedatom]++;
         $ginputconnect->{CONNECTEDFLAG}[$connectedatom] = 1;
         print "$connectedatom  " if $debug >= 4;
       } # for ($i=1;$i<= ($CurrLineElementNumber-1)/2;$i++){

         print "\n" if $debug >= 4;
      
      } #  else { 
      
   } # else {
    
 
 } # while ONIOMINPUTFILE

 $maxbonds = 0;

 # Check for any isolated atoms ($ginputconnect->{CONNECTEDFLAG} is 0)
 for ($i=1;$i<= $atomindex;$i++) {
   if ($ginputconnect->{CONNECTEDFLAG}[$i] == 0) {
     print "Atom $i is an isolated atom.\n";
     }
   
   if ( $maxbonds < $ginputconnect->{CONNECTBONDTABLE}[$i] ){
     $maxbonds = $ginputconnect->{CONNECTBONDTABLE}[$i];
    }
 } # for ($i=1;$i<= $atomindex;$i++) {
 
 if ($ginputconnect->{Cnumber} >= 0) {$maxbonds = $ginputconnect->{Cnumber};}

 # Found all the atoms with the most connectoins

 for ($i=1;$i<= $atomindex;$i++) {
   if ($maxbonds <= $ginputconnect->{CONNECTBONDTABLE}[$i] ){
     $curbonds = $ginputconnect->{CONNECTBONDTABLE}[$i];
     $atomWmaxbonds = $i;
     $atomElenWmaxbonds = $ginputconnect->{ELEMENT}[$i];
     print "Atom $atomWmaxbonds ($atomElenWmaxbonds) has $curbonds connections.\n";
    
    }
 } # for ($i=1;$i<= $atomindex;$i++) {
 
 
 

} # digest


# Write out total charge

sub writecharges{

 my $ginputconnect = shift;
 my $debug = $ginputconnect->{DEBUG};
 
 printf "Total charge of real system is %10.6f.\n",$ginputconnect->{RCHARGE};
 printf "Total charge of high layer is %10.6f.\n",$ginputconnect->{HCHARGE};
 printf "Total charge of medium layer is %10.6f.\n",$ginputconnect->{MCHARGE};
 printf "Total charge of low layer is %10.6f.\n",$ginputconnect->{LCHARGE};
 printf "Total charge of high plus medium layer is %10.6f.\n",$ginputconnect->{HCHARGE}+$ginputconnect->{MCHARGE};
 

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

