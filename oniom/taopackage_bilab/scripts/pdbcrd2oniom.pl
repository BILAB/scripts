#!/usr/bin/perl

use lib "/Users/YoshitakaM/apps/taopackage";
use File::Basename;
use ESPT::PDBCRD2ONIOM 0.1;
use strict;

=head1 NAME

PDBCRD2ONIOM - Takes the coordinates from a PDB file, and puts them into a Gaussian ONIOM file
   
=head1 SYNOPSIS

B<pdbcrd2oniom> [ B<-d?> (?:0~3) ] [ B<-g> F<ONIOM_file_name> ] [ B<-h> ]
[ B<-o> F<output_ONIOM_file_name> ] [ B<-pdb> F<PDB_file_name> ] [ B<-q> ]

=head1 DESCRIPTION

This program generates a new Gaussian ONIOM input file from a Gaussian ONIOM input file 
with coordinates from a given PDB file.
 
=cut

### Version History ###
# 1.0   Take ONIOM template files and PDB file and generate ONIOM file with 
#        new coordinates from given PDB.
#
# 2008 Peng Tao

### To Do List ###
# 

### Main Program ###
our $version = "1.0";

# check for arguments
usage() if ( $#ARGV < 0 );

help() if $ARGV[0] eq "-h";
help() if $ARGV[0] eq "--help";

our ($debug, $extension, $file, $isPDB);
our ($name,  $listresid, $near);
our ($gfile, $pdbfile, $outputfile,$setoutput);
our ($isJobfile, $isLogfile, $isPDB,$wholepdbfile);
our ($haveInput, $havePDB);


$setoutput = 0;
$isJobfile = 0;
$isLogfile = 0;
$haveInput = 0;
$havePDB   = 0;
$wholepdbfile = 0;


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
	if ($ARGV[$i] eq "-g") {
             $gfile = $ARGV[$i + 1] ;
             $haveInput = 1;
          }
        
	if ($ARGV[$i] eq "-pdb") {
             $pdbfile = $ARGV[$i + 1] ;
             $havePDB = 1;
          }
          
        if ($ARGV[$i] eq "-o") {
          $outputfile  = $ARGV[$i + 1];
          $setoutput   = 1;
          }
} 

# set defaults

$debug ||= 0;

print "\nPDBCRD2ONIOM $version : PDB Coordinates --> Gaussian ONIOM file\n\n" if $debug >= 0;

if ($haveInput == 0) {
  print "Gaussian job is missing.\n" if $debug >= 0;
  die "Exit.\n$!";
 }

if ($havePDB == 0) {
  print "Template PDB file is missing.\n" if $debug >= 0;
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


=item B<-g> F<ONIOM_file_name>

Gaussian ONIOM input file.


=item B<-h>

=item B<--help>

Print full PDBCRD2ONIOM documentation via perldoc. Cannot be used with other options.


=item B<-o> F<output_ONIOM_file_name>

Output ONIOM file. Default name will be F<output_ONIOM_file_name>_NEW.gjf.


=item B<-pdb> F<PDB_file_name>

PDB file. Coordinates from this PDB file will be used to generate a
new ONIOM input file using the given Gaussian ONIOM input file as a template.


=item B<-q>

Run in quiet mode and do not print progress messages.



=back

=cut

# Generate new ONIOMtoPDB object
our $input = ESPT::PDBCRD2ONIOM->new(); 

# check for input or log file
open(FILEIN,$gfile) || die "Could not read $gfile\n$!\n";

# If no output file set, use INPUTFILENAME.pdb as output file
 if ( $setoutput == 0) {
  $outputfile = $gfile."_NEW.gjf";
  print "\nNo output PDB file name obtained, use $outputfile.\n" if $debug >= 0;
 }
 

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
   print "\nInput file $gfile is Gaussian Log file.\n" ;
   print "\nCurrently, we cannot use ONIOM log file to create ONIOM input file.\n";
   die "Exit!\n\n";
  } else {
   print "\nInput file $gfile is a Gaussian job file.\n" if $debug >= 0 ;
   $isJobfile = 1;
  }

$isPDB = 0;

# Check template PDB file
open(FILEPDB,$pdbfile) || die "Could not read $pdbfile\n$!\n";

while (<FILEPDB>){
        # skip blank lines
        next if /^$/;
        # PDB file
        if ( /^ATOM\s+\d+/ or /^HETATM\s+\d+/ ) {
                $isPDB = 1;
         	last;
        }
}
close(FILEPDB);

if ( $isPDB == 0 ) { 
  print "\nGiven PDB file $pdbfile is NOT a PDB file. Please double check.\n";
  die "PDBCRD2ONIOM ended. $!\n";
}


print "\nGenerate ONIOM input file $outputfile using $gfile as template and cooridnate from $pdbfile.\n" 
      if $debug >= 0;
 
# read file contents
$input->debug($debug);
$input->{G_FILE_NAME}      = $gfile;
$input->{PDB_FILE_NAME}    = $pdbfile;
$input->{OUTPUT_FILE_NAME} = $outputfile;
$input->{IS_LOG_FILE}      = $isLogfile;
$input->{IS_JOB_FILE}      = $isJobfile;

$input->analyze;

$input->writeoutput;

print "\nPDBCRD2ONIOM exits.\n\n" if $debug >= 0;



## Subroutines ##

# display help on usage
sub help {
	system("perldoc /Users/YoshitakaM/apps/taopackage/scripts/pdbcrd2oniom.pl");
	exit;
}

sub usage {
        print "\nPDBCRD2ONIOM $version : PDB Coordinates --> Gaussian ONIOM file\n";
        print "\nUsage: pdbcrd2oniom [options] \n";
        print "\t-d? \t\t\tDebug print (?:0~3) \n";
	print "\t-g File Name\t\tGaussian ONIOM file name\n";
	print "\t-h\t\t\tprint full documentation\n";
        print "\t-o Output File Name\tNew ONIOM input file name\n";
	print "\t-pdb File Name\t\tPDB file name\n";
	print "\t-q\t\t\tQuiet mode\n";
	print "\n";
        exit;
}


1;

=head1 EXAMPLES

=over

=item pdbcrd2oniom

Called without any parameters, PDBCRD2ONIOM will display usage information. 
If B<-h> or B<--help> is passed then the full PDBCRD2ONIOM documentation is displayed via perldoc.

=item pdbcrd2oniom -g foo.gjf -pdb new.pdb -o foo_New.gjf

PDBCRD2ONIOM reads a Gaussian ONIOM input file foo.gjf, then creates a new Gaussian 
ONIOM input file named foo_New.gjf using the coordinates from the PDB file new.pdb.

=back

=head1 NOTES

PDBCRD2ONIOM is very different from PDB2ONIOM. PDB2ONIOM creates a new ONIOM input file using
a PDB file. PDBCRD2ONIOM needs a PDB file and an existing Gaussian ONIOM input file as a template. 
These two files should have exactly the same number and types of atoms and in the exact same order.
PDBCRD2ONIOM then takes ONLY coordinates from the given PDB file and replaces ONLY the coordinates
in the given Gaussian ONIOM input file, nothing else.

The reason for developing this tool is that some users may find it easier to manipulate 
biomolecule structures when using PDB files. Using both PDBCRD2ONIOM and ONIOM2PDB, users
can easily create a PDB file with an optimized geometry from a Gaussian ONIOM file using
ONIOM2PDB, make some changes to the structure in the PDB file using their favorite tools (to 
create transition state or product, etc), and create a new ONIOM input file with modified coordinates
using PDBCRD2ONIOM.

Please keep in mind that a PDB file only keeps 3 decimal places for the coordinate value. This is 
less than in a Gaussian ONIOM job file (6 decmial places). 


=head1 VERSION

1.0

=head1 AUTHOR

Peng Tao, E<lt>tao.21@osu.eduE<gt>

=head1 COPYRIGHT

Copyright (c) 2009~2010 by Peng Tao


=cut


