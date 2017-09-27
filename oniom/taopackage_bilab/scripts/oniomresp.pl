#!/usr/bin/perl

use lib "/Users/YoshitakaM/apps/taopackage/";
use File::Basename;
use ESPT::ONIOMRESP 0.1;
use strict;

=head1 NAME

ONIOMRESP - Utility to assist in the RESP charge fitting process of ONIOM calculations
   
=head1 SYNOPSIS

B<oniomresp>  [ B<-d?> (?:0~3) ] [ B<-m?> 1~3 ] [ B<-q> ] [ B<-h> ] [ B<-mid> ] 
[ B<-g> F<Gaussian_file_name> ][ B<-o> F<output_file_name> ] [ B<-#> number ] 
[ B<-s> step_number ] [ B<-qin> F<charge_qout_file> ] [ B<-p> F<map_file> ] 
[ B<-c> F<charge_compare_file> ] [ B<-list> F<atom_list_file> ]

=head1 DESCRIPTION

This program generates a Gaussian input file for an ESP calculation or a RESP fitting 
input file from a Gaussian log file or takes partial charges and adds them to an 
ONIOM input file.

There are three different tasks that ONIOMRESP can do, therefore, there are three modes (controlled by flag -m?) to run
ONIOMRESP.

mode -m1
  Generate a Gaussian input file for an ESP calculation from a Gaussian log file. 
  
mode -m2
 Generate a RESP input file from a Gaussian log file

mode -m3
 Take the partial charges from a given charge file and add them to a Gaussian ONIOM input file.


=cut

### Version History ###
# 1.0  Generates RESP fitting input file for given log file,and takes RESP fitted 
#      charges to assigned oniom input file.
#
# 2008 Peng Tao

### To Do List ###
# 1    Generate Gaussian input file for RESP fitting from ONIOM log file
#

### Main Program ###
our $version = "1.0";

# check for arguments
usage() if ( $#ARGV < 0 );

help() if $ARGV[0] eq "-h";
help() if $ARGV[0] eq "--help";

our ($debug, $extension, $mode, $base, $dir, $ext);
our ($gfile, $capHnum, $outputfile, $qinfile,$mapfile, $chgcompfile, $stepnum, $AtomListFile);
our ($isJobfile, $isLogfile);
our ($setcapH, $setoutput,$setmode);
our ($haveInput, $haveqinfile, $havemapfile, $havecompchgfile, $setmid,$HaveListFile);


$setoutput = 0;
$isJobfile = 0;
$isLogfile = 0;
$haveInput = 0;
$setcapH   = 0;
$haveqinfile = 0;
$havemapfile = 0;
$mode = 1;
$setmode = 0;
$capHnum = 0;
$havecompchgfile = 0;
$stepnum = 0;
$setmid = 0;
$HaveListFile = 0;

# parse arguments
for (my $i=0; $i<=$#ARGV; $i++) {
	$debug      = -1 if $ARGV[$i] eq "-q";
	$debug      =  1 if $ARGV[$i] eq "-d";
	$debug      =  0 if $ARGV[$i] eq "-d0";
        $debug      =  1 if $ARGV[$i] eq "-d1";
        $debug      =  2 if $ARGV[$i] eq "-d2";
        $debug      =  3 if $ARGV[$i] eq "-d3";
        $debug      =  4 if $ARGV[$i] eq "-d4";
        $debug      =  5 if $ARGV[$i] eq "-d5";
        $debug      =  6 if $ARGV[$i] eq "-d6";
        
        $setmid     = 1  if $ARGV[$i] eq "-mid";

        $stepnum    =  $ARGV[$i + 1] if $ARGV[$i] eq "-s";
        
        if ($ARGV[$i] eq "-m")  {
           $mode       =  1;
           $setmode    =  1;
          }    
        
        if ($ARGV[$i] eq "-m1")  {
           $mode       =  1;
           $setmode    =  1;
          }    

        if ($ARGV[$i] eq "-m2")  {
           $mode       =  2;
           $setmode    =  1;
          }    

        if ($ARGV[$i] eq "-m3")  {
           $mode       =  3;
           $setmode    =  1;
          }    
        
	if ($ARGV[$i] eq "-g") {
             $gfile = $ARGV[$i + 1] ;
             $haveInput = 1;
          }
        
	if ($ARGV[$i] eq "-#") {
             $capHnum = $ARGV[$i + 1] ;
             $setcapH = 1;
          }
          
        if ($ARGV[$i] eq "-o") {
          $outputfile  = $ARGV[$i + 1];
          $setoutput  = 1;
          }

        if ($ARGV[$i] eq "-qin") {
          $qinfile    = $ARGV[$i + 1];
          $haveqinfile  = 1;
          }
          
        if ($ARGV[$i] eq "-p") {
          $mapfile      = $ARGV[$i + 1];
          $havemapfile  = 1;
          }
          
        if ($ARGV[$i] eq "-c") {
          $chgcompfile      = $ARGV[$i + 1];
          $havecompchgfile  = 1;
          }

        if ($ARGV[$i] eq "-list") {
          $AtomListFile     = $ARGV[$i + 1];
          $HaveListFile     = 1;
          }
        
} 

# set defaults

$debug ||= 0;

print "\nONIOMRESP $version : Utility that helps RESP charge fitting process of ONIOM calculation\n" if $debug >= 0;

if ( $stepnum < 0 ) {
  print "\n Negative step number ($stepnum) obtained, the last geometry will be used.\n" if $debug >= 0;
  $stepnum = 0;
 }

print "\n" if $debug >= 0;

# system('date') if $debug >= 0;

print "\n" if $debug >= 0;

if ($haveInput == 0) {
  print "Gaussian file is missing.\n\n" if $debug >= 0;
  die "Exit.\n\n$!";
 }

if ( $mode == 1 ) {
  if ( $setmode ==  0) {
      print "No -m flag set in command line. Will use default.\n\n" if $debug >= 0; 
     }    
     
  print "Mode one is selected.\n\n" if $debug >= 0;
  print "Will generate Gaussian input file for ESP calculation from given Gaussian log file.\n\n" if $debug >= 0;
 }

if ( $mode == 2 ) {
  print "Mode two is selected.\n\n" if $debug >= 0;
  print "Will generating RESP setup .in file for given log file.\n\n" if $debug >= 0;
 }

if ( $mode == 3 ) {
  print "Mode three is selected.\n\n" if $debug >= 0;
  print "Will take RESP charges, and put it to ONIOM input file.\n\n" if $debug >= 0;
 }
 
 
if ( ($haveqinfile == 0) && ($mode == 3) ) {
  print "You seem want to assign RESP charge to your ONIOM file.\n" if $debug >= 0;
  print "But you didn't tell me which one is your charge file by -qin flag.\n" if $debug >= 0;
  print "Please try again using -qin to tell me which one is your charge file. Thanks.\n\n" if $debug >= 0;
  die "Exit.\n$!";
 }


=head1 OPTIONS

Command line option specifications are processed from left to right and may 
be specified more than once. If conflicting options are specified, later  
specifications override earlier ones.

=over 16

=item B<-d?> (?:0~3)

Turn on debug printing. The printing level can be controlled by a given number.
The larger the number, the more information will be printed when the program is running.


=item B<-m?> (?:1~3)

Running mode setup. There are three modes: -m1, -m2, -m3. 
 Default (either -m or no flag) is -m1.

mode -m1
  Generate a Gaussian input file for an ESP calculation from a given Gaussian log file. 
   The reason for using a the Gaussian log file is that scaling factors 
  (listed in the log file on lines "ONIOM: Cut beween") for positioning 
  the hydrogen capping atoms are available in the log file.
  usage:
  oniomresp -m1 -g gaussianlogfile.log -o gaussianinput_4ESP.gjf

  Default output file name is defined as following.
  If no -o flag is given, then _4ESP.gjf will be added to the base of 
   the given gaussian log file.
   e.g. gaussianlogfile_4ESP.gjf will be generated 
        from gaussianlogfile.log
  
mode -m2
 Generate a RESP setup.in file from a Gaussian log file
 usage:
 oniomresp -m2 -g gaussianlogfile.log -o respsetupfile.in
 
 Default output file name is defined as following.
 If no -o flag is given, then _RESP.in will be added to the 
   base of the given gausian log file.
  e.g. gaussianlogfile_4ESP_RESP.in will be generated from 
   gaussianlogfile_4ESP.log
 
 This mode will take the log file and generate a RESP setup.in file.

mode -m3
 Take the partial charges from given charge file and add them to a Gaussian ONIOM input file.
 usage:
 oniomresp -m3 -g oniomfilename.gjf -qin RESPchargefile.qout -o newoniomfile.gjf -c chargecomparefile.txt
 
 Default output file name is defined as following.
 If no -o flag is given, then _NEWchg.gjf will be added to the base of 
   the given gaussian ONIOM input file.
  e.g. gaussianinput_NEWchg.gjf will be generated 
    from gaussianinput.gjf

 chargecomparefile.txt has both old and new partial charges 
    for comparison purpose
 Default output file name is defined as following.
 If no -c flag is given, then _CompCHG.txt will be added to the base of 
   the given gaussian ONIOM input file.
  e.g. gaussianinput_CompCHG.txt will be generated from gaussianinput.gjf

Please note that depending on which mode is specified, either a Gaussian ONIOM input file or
Gaussian ONIOM log file is needed.
 
 
=item B<-h>

=item B<--help>

Print full ONIOMRESP documentation via perldoc. Can not be used with other options.

=item B<-q>

Run in quiet mode and do not print progress messages.

=item B<-g> F<Gaussian_file_name>

Gaussian file. This file can be either an ONIOM input file or a general log file depending 
on the -m flag. For -m1, an ONIOM log file will be used to generate an ESP single point 
calculation job input file. For -m2, a Gaussian log file will be used to generate a RESP 
input file. For -m3, an ONIOM input file will be used to generate a new input file with 
new charges on the QM atoms.


=item B<-o> F<Output_file_name>

Output file. This file can be either a Gaussian or RESP input file.


=item B<-#> number

number represents capping hydrogen atoms in the given Gaussian log file (from ESP calculation). 
With this option, additional input will be added in a RESP input file to add zero charge 
constraint on capping hydrogen atoms. The capping hydrogens are assumed to be the last
atoms listed in the coordinates.

=item B<-mid>

Used in mode 1.
When set, the middle layer (instead of the high layer) of a three-layer ONIOM job
will be extracted with link atoms for further (RESP) calculations. By default the high
layer (usually QM region) will be extracted with link atoms.


=item B<-s> step_number

When an optimization ONIOM log file is given, user can choose the structure along 
the optimization path that will be used for the RESP charge fitting.
Default is the last geometry.

=item B<-c> F<charge_compare_file>

Charge comparison file. It will display old and new charges together
for comparison purpose.

=item B<-qin> F<charge_qout_file>

Charge file that contains all the partial charges.


=item B<-p> F<map_file>

map_file builds the mapping relationship between the atoms in the ONIOM input file and the atoms in the
the partial charge file.

When this file is provided, ONIOMRESP will take the charge of each atom in a given charge file 
F<charge_qout_file> and assign it to its corresponding ONIOM atom.

Example:

3450 E<10>
...E<10>
3459 E<10>
3582 E<10>
3583 E<10>
3584 E<10>
3585 

Partial charge 1 in F<charge_qout_file> will be assigned to atom 3450 in 
the ONIOM input file "..." means partial charges of atoms 1 through 10 in 
the charge file will be assigned to atoms 3450 through 3459 in ONIOM input file. 
This option saves typing time.

If this file is not provided, ONIOMRESP will just take the partial charge of 
each atom from the partial charge file and assign it to the QM atom in 
the ONIOM input file in sequence.

=item B<-list> F<atom_list_file>

File contains the list of atoms that need to be extracted from the given ONIOM log 
file for the ESP charge calculation.

When this file is provided, ONIOMRESP will extract corresponding atoms from the given ONIOM log file
and write them to a Gaussian input file for the ESP calculation.

Format:

AtomNumber

or

AtomNumber  LinkHostNumber

or 

AtomNumber LinkHostNumber ScalingFactor

or 

AtomNumber LinkHostNumber ScalingFactor Element

For "AtomNumber  LinkHostNumber" input, a hydrogen atom will be created as a link atom using 
a default scaling factor (0.723886, which is for a C-H bond scaled from a C-C bond) based 
on coordinates of AtomNumber and LinkHostNumber.

For "AtomNumber  LinkHostNumber ScalingFactor" input, a hydrogen atom will be created as 
link using the given sacling factor based on coordinates of AtomNumber and LinkHostNumber.


For "AtomNumber  LinkHostNumber ScalingFactor Element" input, a given Element atom 
will be created as a link atom using the given scaling factor based on the coordinates of 
AtomNumber and LinkHostNumber.


Example:

3450 E<10>
...E<10>
3459 E<10>
3582 3590 E<10>
3583 E<10>
3584 3587  0.761575 HE<10>
3585 E<10>
4380 E<10>
... E<10>
4390 E<10>
4391 E<10>
4392 4399 0.761575 E<10>
4393 E<10>
4394 E<10>


=back

=cut



# Generate new ONIOMtoPDB object
our $oniomresp = ESPT::ONIOMRESP->new(); 

# check for input or log file
open(FILEIN,$gfile) || die "Could not read $gfile\n$!\n";

 

# determine file type and set extension

while (<FILEIN>){
        # skip blank lines
        next if /^$/;
	
	# Gaussian 
	if ( /^\s+Entering\s+Gaussian\s+System/ ){
                $isLogfile = 1;
		last;
	}
    }
    
close(FILEIN);

# If input file is not log file, it should be job input file, since we do not
#  use fchk file yet.
if ($isLogfile == 1) {
   $isJobfile = 0;
   print "Input file $gfile is Gaussian log file.\n\n" if $debug >= 0 ;
  } else {
   print "Input file $gfile is Gaussian job file.\n\n" if $debug >= 0 ;
   $isJobfile = 1;
  }

# Check if we have desired input files for the job
# When mode 3 is set

if ( ($isLogfile == 1) && ( $mode == 3 ) )
 {
   print "When mode 3 is set, an ONIOM input file is desired.\n" if $debug >= 0;
   print "But file $gfile is a Gaussian Log file.\n" if $debug >= 0;
   print "Please rethink what you want to do.\n\n" if $debug >= 0;
   die "Exit.\n\n";
 
 }

# When either mode 1 or 2 is set

if ( ($isJobfile == 1) && ( $mode != 3 ) )
 {
   print "When either mode 1 or 2 is set, a Gaussian log file is desired.\n" if $debug >= 0;
   print "But file $gfile is probably a Gaussian input file.\n" if $debug >= 0;
   print "Please rethink what you want to do.\n\n" if $debug >= 0;
   die "Exit.\n\n$!";
 
 }

# Get information ready passing to ONIOMRESP package for different jobs.

# When mode 1 is set, one output file name is needed

if ( ( $mode == 1 ) && ( $setoutput  == 0 )) {
   
  ($base, $dir, $ext) = fileparse($gfile, qr/\.[log|out]*/);
  $outputfile = $base."_4ESP.gjf";
 }

# When mode 2 is set, one output file name is needed

if ( ($mode == 2) &&  ($setoutput  == 0) ) {
    ($base, $dir, $ext) = fileparse($gfile, qr/\.[log|out]*/);
    $outputfile = $base."_RESP.in";
 }

# When mode 3 is set, one output file name is needed

if ( ($mode == 3) ) {
  if ( $setoutput  == 0 ) {
     ($base, $dir, $ext) = fileparse($gfile, qr/\.[com|gjf]*/);
     $outputfile = $base."_NEWchg.gjf";
    }
  
  if ( $havecompchgfile == 0 ) {
     ($base, $dir, $ext) = fileparse($gfile, qr/\.[com|gjf]*/);
     $chgcompfile = $base."_CompCHG.txt";
     print "Charge comparison information will be written to file $chgcompfile.\n\n" if $debug >= 0 ; 
    }
  
  }

print "File $outputfile will be used as output file.\n\n" if $debug >= 0 ; 


 
# passing parameters to package
$oniomresp->debug($debug);
$oniomresp->{G_FILE_NAME}      = $gfile;
$oniomresp->{OUTPUT_FILE}      = $outputfile;
$oniomresp->{MODE}             = $mode;
$oniomresp->{CapH_Num}         = $capHnum;

$oniomresp->{Step_Num}         = $stepnum;
$oniomresp->{Map_File}         = $mapfile;
$oniomresp->{Have_Map_File}    = $havemapfile;
$oniomresp->{Charge_Comp_File} = $chgcompfile;
$oniomresp->{Qin_File}         = $qinfile;


$oniomresp->{IS_LOG_FILE}      = $isLogfile;
$oniomresp->{IS_JOB_FILE}      = $isJobfile;

$oniomresp->{Set_Mid}          = $setmid;

$oniomresp->{HaveListFile}     = $HaveListFile;

$oniomresp->{AtomListFile}     = $AtomListFile if $HaveListFile == 1;


$oniomresp->analyze;

#$oniomresp->writeoutputpdb;

print "\nSuccessfully wrote output file.\n\n" if $debug >= 0;

#system('date') if $debug >= 0;

print "\nONIOMRESP ends.\n\n" if $debug >= 0;

print "\n" if $debug >= 0;

## Subroutines ##

# display help on usage
sub help {
	system("perldoc /Users/YoshitakaM/apps/taopackage/scripts/oniomresp.pl");
	exit;
}

sub usage {
        print "\nONIOMRESP $version : Utility that helps RESP charge fitting process of ONIOM calculation\n";
        print "\nUsage: oniomresp [options] \n";
        print "\t-d? \t\t\tDebug print (?:0~3) \n";
        print "\t-m? \t\t\tChoose mode (?:1~3) \n";
        
	print "\t-h\t\t\tprint full documentation\n";
	print "\t-q\t\t\tQuiet mode\n";

        print "\t-mid\t\t\tExtract middle and high layers from a three-layer ONIOM log file (in mode 1)\n";
	print "\t-g    File Name\t\tGaussian file.\n";
	print "\t-#    Interger Number \tNumber of capping hydrogen atoms\n";
	print "\t-s    Interger Number \tStep Number in ONIOM opt log file for RESP charge fitting\n";
	print "\t-qin  File Name\t\tCharge qout file\n";
	print "\t-p    File Name\t\tMapping atoms in ONIOM input file and RESP fitting file.\n";
        print "\t-list File Name\t\tAtom list that need to be extracted from given ONIOM log file for ESP calculation\n";
        print "\t-c    File Name\t\tChargecomparefile: Old and New RESP charge comparison.\n";
	print "\t-o    Output File Name\n\n";
        
        print "Example for mode 1:\n";
        print "  oniomresp -m1 -g gaussianlogfilename.log -o gaussianinput4ESP.gjf\n\n";
        
        print "Example for mode 2:\n";
        print "  oniomresp -m2 -# CappingHNum -g gaussianlogfilename.log -o respsetupfile.in\n\n";
        
        print "Example for mode 3:\n";
        print "  oniomresp -m3 -g oniominput.gjf  -qin RESPchargefile.qout -o newoniomfile.gjf -c chargecomparefile.txt\n\n";
        
        
	print "\n";
        exit;
}


1;

=head1 EXAMPLES

=over

=item oniomresp

Called without any parameters, ONIOMRESP will display usage information. 
If B<-h> or B<--help> is passed, then the full ONIOMRESP documentation is displayed via perldoc.


=item oniomresp -m1 -g foo.log -o foo_4ESP.in

ONIOMRESP reads a Gaussian ONIOM log file foo.log, then creates a Gaussian QM single point calculation
file named foo_4ESP.in for the electrostatic potential calculation.


=item oniomresp -m2 -g foo.log -# 3 -o fooRESP.in

ONIOMRESP reads foo.log, then creates a RESP input file named fooRESP.in.
The last three atoms will be treated as capping atoms with zero partial charge constraints.

=item oniomresp -m3 -g foo.gjf -qin fooRESP.qout -o fooNEWchg.gjf -c fooOldNewchg.txt

ONIOMRESP reads foo.gjf and fooRESP.qout, then takes charges from fooRESP.qout, assigns
them to the QM atoms in foo.gjf, and saves fooNEWchg.gjf with new RESP charges on QM atoms.
The old and new charges are listed together in fooOldNewchg.txt for the user's reference.


=back

=head1 NOTES

ONIOMRESP

=head1 VERSION

1.0

=head1 AUTHOR

Peng Tao, E<lt>tao.21@osu.eduE<gt>


=head1 COPYRIGHT

Copyright (c) 2009~2010 by Peng Tao


=cut


