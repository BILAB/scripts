#!/usr/bin/perl 

use lib "/Users/YoshitakaM/apps/taopackage/";
use ESPT::Glog 0.07;
use ESPT::Gfchk 0.01;
use ESPT::EXTRACTCHARGE 0.03;
use strict;
use warnings;

=head1 NAME

EXTRACTCHARGE - Extracts extracts atomic charges from a Gaussian log file and save 
them to a file in RESP format.

=head1 SYNOPSIS

B<extractcharge>  [ B<-c?> (?:1~3) ] [ B<-d?> (?:0~3) ] [ B<-g> F<filename> ] [ B<-h> ] [ B<-o> F<filename> ] [ B<-q> ]

=head1 DESCRIPTION

This program is an extraction tool for Gaussian output files. Its main purpose
is to extract the charges of specified atoms from a Gaussian log file. 

=begin comment


### Version History ###
# 1.0 	Extract Mulliken, APT, NPA charges.
#
### 2 Do List ###


=end comment

=cut

### Main Program ###
our $version = 1.0;

# check for arguments
&usage if ( $#ARGV < 0 );
&help if $ARGV[0] eq "-h";
&help if $ARGV[0] eq "--help";

our( $debug, $verbose);
our( $chargetype, $HaveOutputFile, $OutputFilename, $HaveGFile, $gfile, $IsGaussianLogFile, $IsComplete);

# parse arguments
for (my $i=0; $i<=$#ARGV; $i++) {
	$debug          = -1           if $ARGV[$i] eq "-q"   ;
	$debug          =  1           if $ARGV[$i] eq "-d"   ;
	$debug          =  0           if $ARGV[$i] eq "-d0"  ;
        $debug          =  1           if $ARGV[$i] eq "-d1"  ;
        $debug          =  2           if $ARGV[$i] eq "-d2"  ;
        $debug          =  3           if $ARGV[$i] eq "-d3"  ;
        $debug          =  4           if $ARGV[$i] eq "-d4"  ;
        $debug          =  5           if $ARGV[$i] eq "-d5"  ;
        $debug          =  6           if $ARGV[$i] eq "-d6"  ;
	                                                      
        $chargetype     =  1           if $ARGV[$i] eq "-c"   ;
	$chargetype     =  1           if $ARGV[$i] eq "-c1"  ;
	$chargetype     =  2           if $ARGV[$i] eq "-c2"  ;
	$chargetype     =  3           if $ARGV[$i] eq "-c3"  ;
        
        $HaveOutputFile =  1           if $ARGV[$i] eq "-o"   ;
        $OutputFilename =  $ARGV[$i+1] if $ARGV[$i] eq "-o"   ;

        $HaveGFile      =  1           if $ARGV[$i] eq "-g"   ;
        $gfile          =  $ARGV[$i+1] if $ARGV[$i] eq "-g"   ;
                                                              
        $verbose        =  1           if $ARGV[$i] eq "-v"   ;
                                                              
}                                                             
                                                              
# set defaults                                                
$debug           ||=   0;                                     
$verbose         ||=   0;                                     
$chargetype      ||=   1;                                     
$HaveOutputFile  ||=   0;
$OutputFilename  ||=   "Chargefile.dat";
$HaveGFile       ||=   0;
our $calc = ESPT::EXTRACTCHARGE->new();

print "\nEXTRACTCHARGE $version : Atomic charges extraction tool for Gaussian log files.\n\n" if $debug >= 0;



=head1 OPTIONS

Command line option specifications are processed from left to right and may 
be specified more than once. If conflicting options are specified, later  
specifications override earlier ones. 

=over 16

=item B<-c?> (?:1~3)

Type of charges that will be extracted from the given Gaussian log file.
Currently available charge types:
1 for Mulliken charges,
2 for atomic polar tensor (APT) charges,
3 for natural population analysis (NPA) charges. 
Default is Mulliken charges.


=item B<-d?> (?:0~3)

Turn on debug printing. The printing level can be controlled by a given number.
The larger the number, the more information will be printed when program is running.


=item B<-g> F<filename>

Guassian log file from which the charges will be extracted.


=item B<-h>

=item B<--help>

Print usage information. Cannot be used with other options.


=item B<-o> F<filename>

Output file which will contain extracted charges in the following format:

  0.494560 -0.595262  0.442097  0.071576  0.220557 -0.114917 -0.521308  0.235333
  0.084777 -0.048257  

Default file name is Chargefile.dat

=item B<-q>

Run in quiet mode and do not print progress messages.


=back

=cut

# determine generating ESS and create log object


if ( $HaveGFile == 0 ) {
   print "WARNING: there is no GAUSSIAN log file provided. No action can be taken without log file.\n" if $debug >= 0;
   die;
  }


if ($gfile eq $OutputFilename) {
   print "WARNING: you accidentally assigned Gaussian log file $gfile as charge output file.\n Program stopped.\n" if $debug >= 0;
   die;
  }

open(LOGFILE,$gfile) || die "Could not read $gfile\n$!\n";

$IsGaussianLogFile = 0  ;
$IsComplete        = 0  ;

while (<LOGFILE>){
	# skip blank lines
        
       
	next if /^$/;
	
	# Gaussian 
	if ( /^\s+Entering\s+Gaussian\s+System/ ){
                $calc->debug($debug);
                $calc->{VERBOSE}        =   $verbose         ; 
                $calc->{CHARGETYPE}     =   $chargetype      ; 
                $calc->{OUTPUTFILE}     =   $OutputFilename  ;
                $calc->{GAUSSFILE}      =   $gfile           ;
                $IsGaussianLogFile      =   1                ;
	}

	if ( /^\s+Normal\s+termination/ ){
             $IsComplete   = 1  ;
          }
        
}

$calc->{ISNORMAL}   =  $IsComplete ;

  
if ( $IsComplete == 0 ) {
   print "WARNING: Given Gaussian log file does end noramally. This program will keep running to extract\n" if $debug >= 0;
   print "WARNING: requested charges. But please be cautious about these charges.\n\n" if $debug >= 0;
  }
  

$calc->analyze;

print "\nEXTRACTCHARGE ends.\n\n" if $debug >= 0;

print "\n" if $debug >= 0;

exit;



## Subroutines ##

# diplay program info
sub proginfo {
        print "EXTRACTCHARGE version : $version\tInput file: $gfile\n\n";
}

# display help on usage
sub help {
	system("perldoc /Users/YoshitakaM/apps/taopackage/scripts/extractcharge.pl");
	exit;
}

# display help on usage
sub usage {
        print "\nEXTRACTCHARGE $version : Atomic charges extraction tool for Gaussian log files.\n";
        print "\nUsage: extractcharge [options] \n";
        print "\t-c? \t\tTye of charges to extract (?:1~3). Default is 1.\n";
        print "\t\t\t\tCharge type: 1:Mulliken; 2:APT; 3:NPA\n";
        print "\t-d? \t\t\tDebug print (?:0~3)\n";
        print "\t-g Filename\tESS output file name\n";
	print "\t-h \t\tPrint usage information\n";
        print "\t-o Filename\tOutput file for extracted charges. Default is Chargefile.dat\n";
        print "\t-q \t\tQuiet mode for less running info.\n\n";
        exit;
}

1;

=head1 EXAMPLES

=over

=item extractcharge

Called with no parameters at all, EXTRACTCHARGE will display usage information. 
If B<-h> or B<--help> is passed, then the full EXTRACTCHARGE documentation is 
displayed via perldoc.

=item extractcharge -c3 -g foo.log -o foo.qout

EXTRACTCHARGE extracts NPA charges from foo.log file and writes them to foo.qout.

=back


=head1 VERSION

1.0


=head1 AUTHOR

Peng Tao, E<lt>tao.21@osu.eduE<gt>


=head1 COPYRIGHT

Copyright (c) 2009~2010 by Peng Tao


=cut
