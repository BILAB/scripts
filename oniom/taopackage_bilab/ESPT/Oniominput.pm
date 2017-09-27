package ESPT::Oniominput;

use base qw(ESPT::ESSfile);
use ESPT::Glib 0.01;
use strict;
use warnings;

our $VERSION = '0.02';

=head1 NAME

ESPT::Oniominput - Oniominput object to create Gaussian input file.

=head1 SYNOPSIS

   use ESPT::Oniom;

   my $log = Oniom->new();

=head1 DESCRIPTION

This module provides methods to quickly create Gaussian input file
based on template input file and Gaussian log file. It takes coordinates
from given Guassian log file, and replace the coordinates from the given Gaussian
by these taken coordinates. In this way, the setup of old input file (atom type,
partial charges, link atoms, flag for freeze etc.) can be kept.

=cut

### Version History ###
# 0.01	Take geometry out
# 0.02 Grab Input orientation when there is Z-Matrix orientation 

### To Do ###


## the object constructor **

sub new {
	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my $log = ESPT::ESSfile->new();

	$log->{TYPE} = "log";

	# program info
	$log->{PROGRAM} = "Gaussian";
	$log->{VERSION} = undef;
	$log->{REVISION} = undef;
	$log->{COMPILE} = undef;
	$log->{COMPILEDATE} = undef;
	
	# Link 0 & Route commands
	$log->{LINK0} = [];
	$log->{KEYWORDS} = [];
	$log->{ROUTE} = undef;

	# calc info
	$log->{CPUTIME} = undef;
	$log->{FUNCTIONAL} = undef;
	$log->{NPRIMITIVE} = undef;
	$log->{NCARTESIAN} = undef;
	$log->{OPTIMIZED} = undef;	# flag indicating successful optimization
	$log->{RUNTIME} = undef;
	
	# molecular info
	$log->{BASISLABELS} = [];
	$log->{C} = [];			# coefficient matrix
	$log->{CHARGE} = undef;
	$log->{ETHERM} = undef;		# Thermal corrections to E
	$log->{EELEC} = [];		# SCF electronic energy
	$log->{EIGEN} = [];
	$log->{EINFO} = "E(elec)";	# total energy description
	$log->{ENERGY} = undef; 	# total energy 
	$log->{EZPE} = undef;		# ZPE
	$log->{GSOLV} = undef;		# Delta G of Solvation
	$log->{GTHERM} = undef;		# Thermal corrections to G
	$log->{HOMO} = undef;
	$log->{HTHERM} = undef;		# Thermal corrections to H
	$log->{MOSYMM} = [];
	$log->{MULTIPLICITY} = undef;
	$log->{SADDLEPOINT} = undef;	# Saddle-point Order
	$log->{OCC} =[];		# MO occupation info 
	$log->{PG} = [];
	$log->{S} = [];			# overlap matrix
	$log->{SSQUARED} = [];		# S squared values

	# ONIOM info
	$log->{LOWMODEL} = undef;
	$log->{HIGHMODEL} = undef;
	$log->{LOWREAL} = undef;
	$log->{ONIOMENERGY} = undef;
	$log->{LOWREAL_MODEL} = undef;
	$log->{DIPOLEX} = undef;
	$log->{DIPOLEY} = undef;
	$log->{DIPOLEZ} = undef;
	$log->{DIPOLE} = undef;

        # Operation info
        $log->{GJFTEMPLATE} = undef;
	$log->{OPT} = 0;
	$log->{STRUCTURE} = 0;
        $log->{STRUCTUREFILE} = undef;
        $log->{MOVEFLAG} = 0;
        $log->{VERBOSE} = 0;
        $log->{TEMPLATEINPUT} = 0;

	$log->{GJFINPUT} = undef;
        $log->{GJFOUTPUT} = undef;
        $log->{STRUCTURENUMBER} = 0;
        
	bless($log, $class);
	return $log;
}


## methods ##

=item B<$file-E<gt>analyze(filename [spin])>

Analyze the spin results in file called filename.  Spin defaults to Alpha.

=cut

# set filename & spin then digest the file
sub analyze : method {
	my $log = shift;
	$log->digest();
	return;
}


## subroutines ##



sub digest {
# Files larger than 1Mb should be converted to binary and then
# processed to ensure maximum speed. -D. Ennis, OSC

# For items with multiple occurances, the last value is reported 

my $log = shift;

# flags & counters
my $Cflag = 0;
my $symmflag = 0;
my $rparsed = 0;
my $Sflag = 0;
my $Scount = 0;
my $Ccount = 0;
my $dcount = 0;
my $eigcount = 0;
my $Ecount = 0;
my $ESTATEcount = 0;
my $orbcount = -1;
my $MOcount = 0;
my $PGcount = 0;
my $LowReal_Model = 0;
my $HighModel = 0;
my $ONIOMEnergy = 0;
my $diffLowReal_Modelkcal = 0;
my $diffHighModelkcal = 0;
my $diffONIOMEnergykcal = 0;
my $HFtoKcalmol = 627.509;
my $optsteps = 0;
my @qmatomlist;
my @qmatomsymbollist;
my $atomindex = 0;
my $structindicator = 0;
my $qmatomnumber = 0;
my @currentXYZ;
my $currentX;
my $currentY;
my $currentZ;
my $logatom;
my $logatomsymbol;
my $logatomsymbolUP;
my $gjfatom;
my $gjfatomUP;
my $currentline;
my $i=0;
my $j=0;
my $NAtoms=0;
my $structureindex = 0;
my $totalstructure = 0;





print "Gaussian out file is ",$log->{FILENAME},"\n";

print "Gaussian oniom input file ",$log->{GJFINPUT}," is used as template input file.\n";

print "Gaussian oniom new input file is ", $log->{GJFOUTPUT},".\n";

print "Structure No. ", $log->{STRUCTURENUMBER}," will be taken from log file.\n" if $log->{STRUCTURENUMBER} != 0;

print "The last structure will be taken from log file $log->{FILENAME}.\n" if $log->{STRUCTURENUMBER} == 0;





# test
$totalstructure = 0;
# open filename for reading or display error
open(LOGFILE,$log->{FILENAME}) || die "Could not read $log->{FILENAME}\n$!\n";

# grab Total number of atoms and find out how many structures in LOG file
while (<LOGFILE>){
	# skip blank lines
	next if /^$/;

        if (/^\s+NAtoms=\s*(\d+)/) {
           $NAtoms=$1;
#           print "There are $NAtoms atoms.\n";
          }

      if ((/^\s*Z-Matrix\s+orientation:/) || (/^\s*Standard\s+orientation:/) || (/^\s*Input\s+orientation:/)) {
            $totalstructure++;
           }          
        
} # while (<LOGFILE>){

print "There are $NAtoms atoms in LOG file.\n";

print "There are $totalstructure structures in $log->{FILENAME}.\n";

close (LOGFILE);

if ($log->{STRUCTURENUMBER} == 0) {$log->{STRUCTURENUMBER} = $totalstructure;}

# If the user ask for structure that does not exist, complain.
if ($log->{STRUCTURENUMBER} > $totalstructure) {
  printf "There are %4d of structures in the log file.\n",$totalstructure;
  printf "You asked No. %4d of structures!\n",$log->{STRUCTURENUMBER};
  die "Please double check!";
}  # if ($log->{STRUCTURENUMBER} > $totalstructure) {


$structindicator=0;
open(INPUTFILE,$log->{GJFINPUT}) || die "Could not read $log->{GJFTEMPLATE}\n$!\n";

# Count how many atoms are there in the template guassian input file
while (<INPUTFILE>) {

   if   ($structindicator == 0){
     if (/^\s*-*\d+\s+\d+\s+-*\d+\s+\d+/) {
#     if (/^\s*-*\d+\s+\d+/) {
       print $_;
       $structindicator = 1;
       next;
     } # if (/^[0-9]\s[0-9]) 
     next;
   }  # if   ($structindicator ==0){
   else {
     if (/^\s+$/ or /^$/){ last;}
     $atomindex++;

#     printf "%d10\n",$atomindex;
#     if ( /^\s+[a-zA-Z\-]+\d+\.\d+\s+-*[01]\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+H/ ) {
#  If the atom is QM region (or allowed to move when -g is turned on), select this atoms.
#     if ( (/^\s+([a-zA-Z]{1,2})-.+\d+\.\d+\s+-*[01]\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+H/) 
#          ||  ( (/^\s+([a-zA-Z]{1,2})-.+\d+\.\d+\s+0\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)/)
#                && ($log->{MOVEFLAG}==1))) {
#      print $_;

#      if ($log->{VERBOSE} == 1) {  print "$1\t$2\t$3\t$4\t$atomindex\n"; }
#      push (@qmatomsymbollist,$1);
#      push (@qmatomlist, $atomindex);
#      $qmatomnumber++;
#      }  # if
   
   } # else
 
} # while INPUTFILE

print "There are $atomindex atoms in template input file.\n";

close (INPUTFILE);

if ($atomindex != $NAtoms) {
print "There are $NAtoms atoms in log file, and $atomindex atoms in template input file.\n";
print "They are not equal. Please double check! \n";
die;
} # if ($atomindex != $NAtoms) {
 

if ($log->{VERBOSE} == 1) {
 print "Open gaussian log file  $log->{FILENAME} for processing...\n ";
}

open(LOGFILE,$log->{FILENAME}) || die "Could not read $log->{FILENAME}\n$!\n";

if ($log->{VERBOSE} == 1) {
 print "Open template gaussian input file $log->{GJFINPUT}  for processing...\n ";
}

open(INPUTFILE,$log->{GJFINPUT}) || die "Could not read $log->{GJFTEMPLATE}\n$!\n";

if ($log->{VERBOSE} == 1) {
 print "Open output gaussian input file $log->{GJFOUTPUT} for writing...\n ";
}

open(OUTPUTGJFFILE,">$log->{GJFOUTPUT}") || die "Could not read $log->{GJFOUTPUT}\n$!\n";


while (<INPUTFILE>) {  #NO1

     print OUTPUTGJFFILE $_ ;

     if (/^\s*-*\d+\s+\d+\s+-*\d+\s+\d+/) { last;} 
     
} # while (<INPUTFILE>) { #NO1


$structureindex = 0;

while (<LOGFILE>) {  #No2
     
      if ((/^\s*Z-Matrix\s+orientation:/) || (/^\s*Standard\s+orientation:/) || (/^\s*Input\s+orientation:/)) { #NO3
            $structureindex++;
            if ( $structureindex > $log->{STRUCTURENUMBER} ) {last;}
            if ( $structureindex < $log->{STRUCTURENUMBER} ) { 
              for($i=1;$i<=$NAtoms+4;$i++) {<LOGFILE>;}
              next;
            } else { #NO5
             @currentXYZ = ();
             $_=<LOGFILE>;
             if ( not(/^\s+-+/)) {
              print "Line($_)\n";
              print "Gaussian log format may change. Pay attention to output structure(s).\n";
              }
             $_ = <LOGFILE>;
             if ( not(/^\s+Center/) ) {
              print "Line($_)\n";
              print "Gaussian log format may change. Pay attention to output structure(s).\n";
              }

             $_ = <LOGFILE>;
             if ( not(/^\s+Number\s+Number\s+Type\s+X\s+Y\s+Z/) ) {
              print "Line($_)\n";
              print "Gaussian log format may change. Pay attention to output structure(s).\n";
              }

             $_=<LOGFILE>;
             if ( not(/^\s+-+/)) {
              print "Line($_)\n";
              print "Gaussian log format may change. Pay attention to output structure(s).\n";
              }
              
             for ($i = 1; $i<=$NAtoms;$i++){ #NO6
                $_ = <LOGFILE>;
#                print $_;

                /^\s+.+\s+(\d+)\s+\d+\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)/;
                $logatom  = $1;
                $currentX = $2;
                $currentY = $3;
                $currentZ = $4;
#	        if ($log->{VERBOSE} == 1) {	                
#                  printf "Taking coordinates (%16.10f %16.10f %16.10f) from log file.\n",$2,$3,$4;
#                  }
                $currentline = <INPUTFILE>;
#                print $currentline;
                
#                print "Old input file line:\n",$currentline;
                
#                s/^\s+[a-zA-Z]{1,2}-.+\d+\.\d+\s+-*[01]\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+/$currentX $currentY $currentZ/;
                
                $currentline =~ 
                  s/\s(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+/ $currentX     $currentY     $currentZ    /;
  
#                  print "New input file line:\n",$currentline;
                  print OUTPUTGJFFILE $currentline;
                  
                  $currentline =~  /^\s*([a-zA-Z]+)/;
                  $gjfatom = $1;
#                  print "Current atom from logfile is $logatom, from gaussian input file is $gjfatom.\n";
                  $logatomsymbol=$log->atomconvert($logatom);
#                  print "Current atom from logfile is $logatomsymbol, from gaussian input file is $gjfatom.\n";

#		  change all the symbol to the upper case                    

                  $logatomsymbolUP = $logatomsymbol;
                  $logatomsymbolUP =~ tr/a-z/A-Z/;
                  $gjfatomUP = $gjfatom;
                  $gjfatomUP =~ tr/a-z/A-Z/;

                  if ($logatomsymbolUP ne $gjfatomUP) {
                    print "Atom $i in LOG file is $logatomsymbol, but $gjfatomUP in template Gaussian input file.\n";
                    print "Please double check.\n";
                    die;
                  
                   }
  
                 
                } #NO6
              
              
              
              
            } #NO5
             

           } #NO3          
      



} # while (<LOGFILE>) { #No2

while (<INPUTFILE>) {  #NO9

     print OUTPUTGJFFILE $_ ;

 
} # while (<INPUTFILE>) { #NO9



} #digest

1;
__END__

=head1 VERSION

0.01

=head1 AUTHOR

Peng Tao, E<lt>tao.21@osu.edu<gt>

=head1 COPYRIGHT

Copyright (c) 2009~2010 by Peng Tao

=cut

