#!/usr/bin/perl

use lib "/Users/YoshitakaM/apps/taopackage/";
use File::Basename;
use ESPT::PDB 0.03;
use strict;

=head1 NAME

PDB2ONIOM - Generates a Gaussian ONIOM input file from a PDB file
   
=head1 SYNOPSIS

B<pdb2oniom> [ B<-d?> (?:0~3) ] [ B<-h> ] [ B<-i> F<filename> ] [ B<-near> number ] 
[ B<-o> F<filename> ] [ B<-q> ] [ B<-resid> F<listfile> ] 

=head1 DESCRIPTION

This program generates Gaussian ONIOM input files from PDB files.   
 

=cut

### Version History ###
# 1.0   Take PDB file and generate ONIOM input file.
# 1.1	Set freeze flag according to user's request
# 1.2	Need flag -i for input file
#
# 2009 Peng Tao

### To Do List ###
#

### Main Program ###
our $version = "1.2";

# check for arguments
usage() if ( $#ARGV < 0 );

#help() if $ARGV[0] eq "-h";
#help() if $ARGV[0] eq "--help";

our ($debug, $file, $isPDB);
our ($outputname,  $listresid, $near);
our ($topexist, $topfile, $isTOP, $topcheck); # M.Yamada
our ($help); # M.Yamada

# parse arguments
for (my $i=0; $i<=$#ARGV; $i++) {
	$debug =  1 if $ARGV[$i] eq "-d";
	$debug = -1 if $ARGV[$i] eq "-q";
	$debug =  0 if $ARGV[$i] eq "-d0";
        $debug =  1 if $ARGV[$i] eq "-d1";
        $debug =  2 if $ARGV[$i] eq "-d2";
        $debug =  3 if $ARGV[$i] eq "-d3";
        $debug =  4 if $ARGV[$i] eq "-d4";
        $debug =  5 if $ARGV[$i] eq "-d5";
        $debug =  6 if $ARGV[$i] eq "-d6";
	$outputname = $ARGV[$i + 1] if ( $ARGV[$i] eq "-o");
        $listresid=$ARGV[$i + 1] if $ARGV[$i] eq "-resid";
        $near=$ARGV[$i + 1] if $ARGV[$i] eq "-near";
        $file = $ARGV[$i + 1] if $ARGV[$i] eq "-i";
        $topfile = $ARGV[$i + 1] if $ARGV[$i] eq "-t"; # M.Yamada
	$topexist = 1 if $ARGV[$i] eq "-t"; # M.Yamada
	$help = 1 if $ARGV[$i] eq "-h"; # M.Yamada
	$help = 1 if $ARGV[$i] eq "--help"; # M.Yamada
} 

# set defaults
$help ||= 0; # M.Yamada
$debug ||= 0;
$near ||= 6.0;
$listresid ||= "THERE_IS_NO_CORE_RESIDUE_SET_HERE";
$outputname ||= "newGaussianfile.gjf";

help() if $help eq 1; # M.Yamada

print "\nPDB2ONIOM $version PDB -> Gaussian ONIOM input file\n\n" if $debug >= 0;

if (not (defined $file)) {
 print "No PDB file was given. Please use flag -i to provide a PDB file.\n";
 die "Exit.\n$!";
}
# M.Yamada
if (not (defined $topfile)) {
 $topexist = 0;
 print "No .top file\n" if $debug >= 0;
}
# M.Yamada END

=head1 OPTIONS

Command line option specifications are processed from left to right and may 
be specified more than once. If conflicting options are specified, later  
specifications override earlier ones.

=over 16

=item B<-d?> (?:0~3)

Turn on debug printing. The printing level can be controlled by a given number.
The larger the number, the more information will be printed when the program is running.


=item B<-h>

=item B<--help>

Print full PDB2ONIOM documentation via perldoc. Cannot be used with other options.

=item B<-i> F<filename>

PDB file to be used to generate the Gaussian ONIOM input file.

=item B<-near> number

number is in Angstroms. Any residue with at least one atom within this distance 
of any atom in the core residues will be free to move during the optimization process 
(with flag 0 in the Gaussian ONIOM job file). Default is 6.0 A.

=item B<-o> F<filename>

New Gaussian ONIOM input file. The default is newGaussianfile.gjf.

=item B<-q>

Run in quiet mode and do not print progress messages.

=item B<-resid> F<listfile>

File containing the list of core residues. The distance between atoms and these core residues
will be used to set up the optimization flags. The format of file is 
[Residue Name] "Residue ID", eg:

[GLU] "345"

[ASP] "233"

[GLY] "344"

Information for each core residue needs to be on a separate line.

Please note that Residue ID is a string, which can contain both the residue number and chain name, 
for example "A 345" (residue 345 in chain A) may also be used as Residue ID.


=back

=cut

# check for prmtop or fchk file
open(FILEIN,$file) || die "Could not read $file\n$!\n";

# get path and base filename

our ($input);

$isPDB = 0;

# determine file type and set extension
while (<FILEIN>){
        # skip blank lines
        next if /^$/;
	
        # PDB file
        if ( /^ATOM\s+\d+/ or /^HETATM\s+\d+/ ) {
                $input = ESPT::PDB->new();   
                $input->{DEBUG} = $debug;
                $input->{RESIDFILE} = $listresid;
                $input->{NEAR} = $near;
                print "It is a PDB file \n" if $debug>0;
                $isPDB = 1;
         	last;
        }
}
close(FILEIN);

if ( $isPDB == 0 ) { 
  print "Input file $file is NOT a PDB file. Please double check. \n";
  die "PDB2ONIOM ended.\n";
}

#print "Processing a ", $input->get("PROGRAM"), " ", $input->get("TYPE"), " file.\n" if $debug >= 0;

# M.Yamada
$input->{TOPEXIST} = $topexist;
if ( $input->{TOPEXIST} == 1 ) {
 open(TOPIN,$topfile) || die "Could not read $topfile\n$!\n";
 $isTOP = 0;
 $topcheck = 0;
 while (<TOPIN>){
         # skip blank lines
         next if /^$/;
         next if /^\s*;/;
 	
         # TOP file
         if ( /[\s*moleculetype\s*]/ and $topcheck == 0 ) {
                 $topcheck = 1;
               	 next;
         }
         if ( /[\s*atoms\s*]/ and $topcheck == 1 ) {
                 print "It is a TOP file \n" if $debug>0;
                 $topcheck = 2;
                 $isTOP = 1;
               	 next;
         }
         if ( /[\s*bonds\s*]/ and $topcheck == 2 ) {
                 $topcheck = 3;
               	 last;
         }
 }
 close(TOPIN);

 if ( $isTOP == 0 ) { 
   print "Input file $topfile is NOT a topology file. Please double check. \n";
   die "PDB2ONIOM ended.\n";
 }
 if ( $topcheck == 2 ) { 
   print "Input file $topfile has no bonds. Please double check. \n";
   die "PDB2ONIOM ended.\n";
 }
}
# M.Yamada END

# read file contents
$input->debug($debug);
$input->analyze($file, $topfile); # M.Yamada : $topfile

$input->writeoniom($outputname);

#print $input->{ORIGINALLINE}[0];

print "Successfully wrote ", $outputname, " file.\n" if $debug >= 0;

print "\nPDB2ONIOM ends.\n\n" if $debug >= 0;


## Subroutines ##

# display help on usage
sub help {
	system("perldoc /Users/YoshitakaM/apps/taopackage/scripts/pdb2oniom.pl");
	exit;
}

sub usage {
        print "\nPDB2ONIOM $version : PDB -> Gaussian ONIOM input file\n";
        print "\nUsage: pdb2oniom [options]\n";
        print "\t-d? \t\tDebug print (?:0~3) \n";
	print "\t-h\t\tPrint full documentation\n";
	print "\t-o name\t\tOutput file name, defaults to newGaussianfile.gjf\n";
	print "\t-q\t\tQuiet mode\n";
        print "\t-resid listfile\tFile that contains the list of core residues\n";
        print "\t-near number\tDistance in angstrom to define moving region within core residues (Default 6.0)\n";
	print "\t-i filename\tA valid pdb file.\n";
	print "\n";
        exit;
}


1;

=head1 EXAMPLES

=over

=item pdb2oniom

Called without any parameters, PDB2ONIOM will display usage information. 
If B<-h> or B<--help> is passed, then the full PDB2ONIOM documentation is displayed via perldoc.

=item pdb2oniom -i foo.pdb -o foo.gjf

PDB2ONIOM reads foo.pdb, then creates a Gaussian ONIOM input file named foo.gjf.
Since no core residue list file was provided, all atoms in foo.gjf are free to move during
geometry optimization.

=item pdb2oniom -resid residfile -near 6.5 foo.pdb -o foo.gjf

PDB2ONIOM reads foo.pdb, then creates a Gaussian ONIOM input file named foo.gjf.
All atoms in any non-core residue within 6.5 angstrom of any core 
residue listed in residfile are are free to move during geometry optimization run by foo.gjf.


=back

=head1 NOTES

PDB2ONIOM also generates a .onb file. This file contains the Gaussian ONIOM job setup and
residue information. Do not delete this file. This file will be needed for programs SETMVFLG 
and TRANSGEOM.

=head1 SEE ALSO

SETMVFLG 
TRANSGEOM

=head1 VERSION

1.2

=head1 AUTHOR

Peng Tao, E<lt>tao.21@osu.eduE<gt>

=head1 COPYRIGHT

Copyright (c) 2009~2010 by Peng Tao


=cut


