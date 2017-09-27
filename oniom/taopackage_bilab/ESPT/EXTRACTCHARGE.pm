package ESPT::EXTRACTCHARGE;

use base qw(ESPT::ESSfile);
use ESPT::Glib 0.01;
use strict;
use warnings;

our $VERSION = '0.03';

=head1 NAME

ESPT::EXTRACTCHARGE - EXTRACTCHARGE object to extract certain type of atomic charges from Gaussian log file.

=head1 SYNOPSIS

   use ESPT::EXTRACTCHARGE;

   my $log = EXTRACTCHARGE->new();

=head1 DESCRIPTION

This module provides methods to extract certain type of atomic charges from Gaussian log file. 
Guassian log files can only be read currently.

=cut

### Version History ###
# 0.01	Extract three type of charges: Mulliken, APT, NPA.


### To Do ###


## the object constructor **

sub new {
	my $invocant          =  shift;
	my $class             =  ref($invocant) || $invocant;
	my $log               =  ESPT::ESSfile->new();

	$log->{TYPE}          =  "log";

	# program info
	$log->{PROGRAM}       =  "Gaussian" ;
	$log->{VERSION}       =  undef      ;
	$log->{REVISION}      =  undef      ;
	$log->{COMPILE}       =  undef      ;
	$log->{COMPILEDATE}   =  undef      ;
	
	# molecular info
	$log->{CHARGE}        =  undef      ;

	# ONIOM info
	$log->{SCFDONE}       =  undef      ;
	$log->{HIGHMODEL}     =  undef      ;
	$log->{LOWREAL}       =  undef      ;
        $log->{MEDMODLE}      =  undef      ;
        $log->{LOWMID}        =  undef      ;
        $log->{MEDMID}        =  undef      ;
	$log->{ONIOMENERGY}   =  undef      ;
	$log->{LOWREAL_MODEL} =  undef      ;
	$log->{DIPOLEX}       =  undef      ;
	$log->{DIPOLEY}       =  undef      ;
	$log->{DIPOLEZ}       =  undef      ;
	$log->{DIPOLE}        =  undef      ;
        $log->{ISTHREELAYER}  =  undef      ;

        # Operation info                    
        $log->{VERBOSE}       =  0          ;         
        $log->{ISNORMAL}      =  0          ;         
        $log->{CHARGETYPE}    =  1          ;              
	$log->{OUTPUTFILE}    =  undef      ;
	$log->{GAUSSFILE}     =  undef      ;
        
                              
	bless($log, $class);  
	return $log;
}


## methods ##

=item B<$file-E<gt>analyze(filename [spin])>

Analyze the spin results in file called filename.  

=cut

# set filename then digest the file
sub analyze : method {
	my $log = shift;
	$log->prepare(@_);
	$log->digest();
	return;
}


## subroutines ##



sub digest {

# For items with multiple occurances, the last value is reported 

my $log = shift;

# flags & counters
my $atomindex = 0;
my $structindicator = 0;
my $qmatomnumber = 0;
my $logatomsymbol;
my $currentline;
my $i=0;
my $j=0;
my $NAtoms=0;
my $structureindex = 0;
my $finalscfenergy;
my $diffscfenergy;
my $verbose;
my $debug;
my $chargetype;
my $outputfile;
my $gfile; 
my $totalcharge;
my @ChargeList;
my @currentatomno;
my @currentAtomSym;


$verbose     =   $log->{VERBOSE}    ;
$debug       =   $log->{DEBUG}      ;
$chargetype  =   $log->{CHARGETYPE} ;
$outputfile  =   $log->{OUTPUTFILE} ;
$gfile       =   $log->{GAUSSFILE}  ;

print "Given Gaussian log file is ",$gfile,"\n" if $debug > 0;

# open filename for reading or display error
open(LOGFILE,$gfile) || die "Could not read $gfile\n$!\n";

# grab everything which may be useful
while (<LOGFILE>){
	# skip blank lines
	next if /^$/;

        if (/^\s*NAtoms=\s*(\d+)/) {
           $NAtoms=$1;
           print "There are $NAtoms atoms.\n" if $debug > 2;
         }
}


 close(LOGFILE);

# If NAtoms information is not available, we need to find it by
# counting in standard orientation.

   if ( $NAtoms == 0 ) {
      print "NAtoms info is not presented in $gfile. \n\nTrying to find number of atoms in another way...\n\n" if $debug >= 0;
      
      open(LOGFILE,$gfile) || die "Could not read $gfile\n$!\n";

      while (<LOGFILE>){
 
           if ((/^\s*Z-Matrix\s+orientation:/i) || (/^\s*Standard\s+orientation:/i) || (/^\s+Input\s+orientation:/i) ) {

               $_ = <LOGFILE>;
               if ( not(/^\s+-+/)) {
                  print "Line($_)\n" if $debug >= 3;
                  print "Gaussian log format may change. Pay attention to output structure(s).\n" if  $debug >= 0;
                  }
               $_ = <LOGFILE>;
               if ( not(/^\s*Center/i) ) {
                  print "Line($_)\n" if $debug >= 3;
                  print "Gaussian log format may change. Pay attention to output structure(s).\n"  if $debug >= 0;
                  }

               $_ = <LOGFILE>;
               if ( not(/^\s*Number\s+Number\s+Type\s+X\s+Y\s+Z/i) ) {
                  print "Line($_)\n" if $debug >= 3;
                  print "Gaussian log format may change. Pay attention to output structure(s).\n"  if $debug >= 0;
                  }

               $_ = <LOGFILE>;

               if ( not(/^\s+-+/)) {
                  print "Line($_)\n" if $debug >= 3;
                  print "Gaussian log format may change. Pay attention to output structure(s).\n" if  $debug >= 0;
                  }
               
               $_ = <LOGFILE>;
               
               while (not(/^\s+-+/)) {
                  $_ = <LOGFILE>;
                  $NAtoms++;
                } # while (not(/^\s+-+/) {
              }

          } # while LOGFILE

       } # if  ( $NAtoms == 0 ) {
   
 close(LOGFILE);

 print "There are $NAtoms atoms in file $gfile.\n\n" if $debug >= 0;

 
# Extract Mulliken charges. 
# Mulliken Charges example:
#   Mulliken atomic charges:
#                1
#       1  C   -0.351292
#       2  H    0.125741
#       3  H    0.373141
#       4  H    0.135956
#       5  F   -0.908440
#       6  F   -0.375106
#   Sum of Mulliken atomic charges =  -1.00000
#

if ($chargetype == 1) {

 close(LOGFILE);

 open(LOGFILE,$gfile) || die "Could not read $gfile\n$!\n";

 while (<LOGFILE>){ # tag 1
   
   if (/^\s+Mulliken\s+atomic\s+charges:\s*$/) { # tag 2
    
    # Skip first line 
    $_=<LOGFILE>;
    
    # Empty several arraies
    @ChargeList     =  ()  ;
    @currentatomno  =  ()  ;
    @currentAtomSym =  ()  ;
    
    # Read in all Mulliken charges
    for ($i = 0; $i<$NAtoms;$i++){
    
       $_ = <LOGFILE>;

           /^\s+(\d+)\s+(\S+)\s+(-*\d+\.\d+)/;
           push (@currentatomno,  $1);
           push (@currentAtomSym, $2);
           push (@ChargeList,     $3);
           
           print "Current line is \n$_"                           if $debug >= 2;
           print "$i\t\t",$currentatomno[$i],"\t"                 if $debug >= 1;
           print $currentAtomSym[$i],"\t",$ChargeList[$i],"\n"    if $debug >= 1;
           print "\t\t$1\t$2\t$3\n"                               if $debug >= 3;

      } # for
    
    $_ = <LOGFILE>;
    
    /^\s+Sum\s+of\s+Mulliken\s+atomic\s+charges\s*=\s*(-*\d+\.\d+)/;
    $totalcharge = $1;
    print "Total Mulliken charges are $totalcharge\n"                    if $debug >= 1;
    
    
    } # if (/^\s+Mulliken\s+atomic\s+charges/)  tag 2
    
   
  
  } #  while (<LOGFILE>){ # tag 1
 

 } # if ($chargetype == 1) {Extract Mulliken charges 
 
 


# Extract APT charges. 
# APT Charges example:
#   APT atomic charges:           
#                1                
#       1  C    0.484463          
#       2  H   -0.100277          
#       3  H    0.292772          
#       4  H   -0.084655          
#       5  F   -1.012599          
#       6  F   -0.579704          
#   Sum of APT charges=  -1.00000 
#

if ($chargetype == 2) { # tag 5

 close(LOGFILE);

 open(LOGFILE,$gfile) || die "Could not read $gfile\n$!\n";

 while (<LOGFILE>){ # tag 4
   
   if (/^\s+APT\s+atomic\s+charges:\s*$/) { # tag 3
    
    # Skip first line 
    $_=<LOGFILE>;
    
    # Empty several arraies
    @ChargeList     =  ()  ;
    @currentatomno  =  ()  ;
    @currentAtomSym =  ()  ;
    
    # Read in all APT charges
    for ($i = 0; $i<$NAtoms;$i++){
    
       $_ = <LOGFILE>;

           /^\s+(\d+)\s+(\S+)\s+(-*\d+\.\d+)/;
           push (@currentatomno,  $1);
           push (@currentAtomSym, $2);
           push (@ChargeList,     $3);
           
           print "Current line is \n $_"                           if $debug >= 2;
           print "$i\t\t",$currentatomno[$i],"\t"                  if $debug >= 1;
           print $currentAtomSym[$i],"\t",$ChargeList[$i],"\n"     if $debug >= 1;
           print "\t\t$1\t$2\t$3\n"                                if $debug >= 3;

      } # for
    
    $_ = <LOGFILE>;
    
    /^\s+Sum\s+of\s+APT\s+charges\s*=\s*(-*\d+\.\d+)/;
    $totalcharge = $1;
    print "Total APT charges are $totalcharge\n"                if $debug >= 1;
    
    
    } # if (/^\s+APT\s+atomic\s+charges:\s*$/)  tag 3
    
  } #  while (<LOGFILE>){ # tag 4
 

 } # if ($chargetype == 2) {Extract APT charges tag 5
 
 


# Extract NPA charges. 
# NPA Charges example:
#  
#                                         Natural Population 
#                  Natural  -----------------------------------------------
#      Atom  No    Charge         Core      Valence    Rydberg      Total
#   -----------------------------------------------------------------------
#        C    1   -0.23316      1.99930     4.21495    0.01891     6.23316
#        H    2    0.15052      0.00000     0.84819    0.00129     0.84948
#        H    3    0.30852      0.00000     0.68986    0.00162     0.69148
#        H    4    0.15881      0.00000     0.83998    0.00121     0.84119
#        F    5   -0.94532      2.00000     7.94442    0.00090     9.94532
#        F    6   -0.43937      1.99996     7.43027    0.00915     9.43937
#   =======================================================================
#     * Total *   -1.00000      5.99926    21.96766    0.03308    28.00000
#

if ($chargetype == 3) { # tag 6

 close(LOGFILE);

 open(LOGFILE,$gfile) || die "Could not read $gfile\n$!\n";

 while (<LOGFILE>){ # tag 7
    
   if (/^\s+Summary\s+of\s+Natural\s+Population\s+Analysis:/) { # tag 8

    $_=<LOGFILE>;
    if ( not(/^\s*$/) ) {
       print "Line($_)\n";
       print "Natural Population Analysis (NPA) output format may change. Pay attention to output charges.\n";
       }

    $_=<LOGFILE>;
    if ( not(/^\s+Natural\s+Population\s+/) ) {
       print "Line($_)\n";
       print "Natural Population Analysis (NPA) output format may change. Pay attention to output charges.\n";
    }
    
    $_=<LOGFILE>;
    if ( not(/^\s+Natural\s+-+/) ) {
       print "Line($_)\n";
       print "Natural Population Analysis (NPA) output format may change. Pay attention to output charges.\n";
       }
       
    $_=<LOGFILE>;
    if ( not(/^\s+Atom\s+No\s+Charge\s+Core\s+Valence\s+Rydberg\s+Total/) ) {
       print "Line($_)\n";
       print "Natural Population Analysis (NPA) output format may change. Pay attention to output charges.\n";
       }

    $_=<LOGFILE>;
    if ( not(/^\s*-+/)) {
       print "Line($_)\n";
       print "Gaussian log format may change. Pay attention to output structure(s).\n";
       }
    
    # Empty several arraies
    @ChargeList     =  ()  ;
    @currentatomno  =  ()  ;
    @currentAtomSym =  ()  ;
    
    # Read in all APT charges
    for ($i = 0; $i<$NAtoms;$i++){
    
       $_ = <LOGFILE>;

           /^\s+(\S+)\s+(\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)/;
           push (@currentAtomSym, $1);
           push (@currentatomno,  $2);
           push (@ChargeList,     $3);
           
           print "Current line is \n$_"                          if $debug >= 2;
           print "$i\t\t",$currentatomno[$i],"\t"                if $debug >= 1;
           print $currentAtomSym[$i],"\t",$ChargeList[$i],"\n"   if $debug >= 1;
           print "\t\t$1\t$2\t$3\n"                              if $debug >= 3;

      } # for
    
    $_ = <LOGFILE>;  # Skip line ========================
    print "Line after APT charge list:\n" if $debug >= 3;
    print $_                              if $debug >= 3;
    
    $_ = <LOGFILE>;
    
    /(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)/;
    $totalcharge = $1;
    print "Total NPA charges are $totalcharge\n"                if $debug >= 1;
    
    
    } # if (/^\s+Natural\s+Population\s*$/)  tag 8
    
  } #  while (<LOGFILE>){ # tag 7
 

 } # if ($chargetype == 3) {Extract NPA charges tag 6
 
 
# Write out charges
    
                 
   open(OUTPUTFILE,">$outputfile") || die "Could not write to $outputfile\n$!\n";
    
    $j = 0;
    
    for ($i = 0; $i<$NAtoms;$i++){
      printf OUTPUTFILE " %9.6f", $ChargeList[$i];
      
      $j++;
      if ($j == 8) {
          print OUTPUTFILE "\n";
          $j = 0;
        }
     }
                     
    print OUTPUTFILE "\n";

    print "\nSuccessfully wrote charge file $outputfile.\n\n" if $debug >= 0;

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

