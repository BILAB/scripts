#!/usr/bin/perl

use lib "/Users/YoshitakaM/apps/taopackage/";
use File::Basename;
use ESPT::PDB 0.03;
use strict;

=head1 NAME

PDBCORE - Extracts part of a PDB structure surrounding core residues
   
=head1 SYNOPSIS

B<pdbcore> [ B<-d?> (?:0~3) ] [ B<-h> ] [ B<-i> F<filename> ] [ B<-near> number ] 
[ B<-o> F<filename> ] [ B<-q> ] [ B<-resid> F<listfile> ]

=head1 DESCRIPTION

This program extracts part of a structure from a given PDB file based on 
a specified distance between the core region and the surrounding atoms 
of the given PDB file.
 
=cut

### Version History ###
# 1.0   Take PDB file and generate PDB file.
#
# 2007 Peng Tao

### To Do List ###
#

### Main Program ###
our $version = "1.0";

# check for arguments
usage() if ( $#ARGV < 0 );

help() if $ARGV[0] eq "-h";
help() if $ARGV[0] eq "--help";

our ($debug, $file, $isPDB);
our ($outputname,  $listresid, $near);

# parse arguments
for (my $i=0; $i<=$#ARGV; $i++) {
	$debug = 1 if $ARGV[$i] eq "-d";
	$debug = -1 if $ARGV[$i] eq "-q";
	$debug =  0 if $ARGV[$i] eq "-d0";
        $debug =  1 if $ARGV[$i] eq "-d1";
        $debug =  2 if $ARGV[$i] eq "-d2";
        $debug =  3 if $ARGV[$i] eq "-d3";
        $debug =  4 if $ARGV[$i] eq "-d4";
        $debug =  5 if $ARGV[$i] eq "-d5";
        $debug =  6 if $ARGV[$i] eq "-d6";
	$outputname = $ARGV[$i + 1] if  $ARGV[$i] eq "-o";
        $listresid=$ARGV[$i + 1] if $ARGV[$i] eq "-resid";
        $near=$ARGV[$i + 1] if $ARGV[$i] eq "-near";
        $file = $ARGV[$i + 1] if $ARGV[$i] eq "-i";
} 

# set defaults
$debug ||= 0;
$near ||= 12.0;
$listresid ||= "THERE_IS_NO_CORE_RESIDUE_SET_HERE";
$outputname ||= "newONIOMFile.gjf";

print "\nPDBCORE $version PDB -> smaller PDB\n\n" if $debug >= 0;

if (not (defined $file)) {
 print "No PDB file was given. Please use flag -i to provide a PDB file.\n";
 die "Exit.\n$!";
}

if ($listresid eq "THERE_IS_NO_CORE_RESIDUE_SET_HERE") {
 print "No core residues list file was given. Please use flag -resid to provide such file.\n";
 print "Refer to manual page for more information: pdbcore -h\n";
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


=item B<-h>

=item B<--help>

Print full PDBCORE documentation via perldoc. Cannot be used with other options.


=item B<-i> F<filename>

Input PDB file. 


=item B<-near> number

number is in Angstroms. Atoms within number distance of core residue atoms will
be extracted. Default is 12.0.


=item B<-o> F<filename>

Output PDB file. The default is newONIOMFile.gjf.


=item B<-q>

Run in quiet mode and do not print progress messages.



=item B<-resid> F<listfile>

File containing the list of core residues. The distance between surrounding 
atoms and these core residues will be used to set up the optimization flags. 
The file format is 
[Residue Name] "Residue ID", eg:

[GLU] "345"

[ASP] "233"

[GLY] "344"

Information for each core residue needs to be on a separate line.

Please note that Residue ID is a string, which can contain both the residue number and chain name, 
for example "A 345" (residue 345 in chain A) may also be used as Residue ID.


=back

=cut

# check file
open(FILEIN,$file) || die "Could not read $file\n$!\n";


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
                print "It is a PDB file \n" if $debug>1;
                $isPDB = 1;
         	last;
        }
}
close(FILEIN);

if ( $isPDB == 0 ) { 
  print "Input file $file is NOT a PDB file. Please double check. \n";
  die "PDBCORE ended.\n";
}

#print "Processing a ", $input->get("PROGRAM"), " ", $input->get("TYPE"), " file.\n" if $debug >= 0;

# read file contents
$input->debug($debug);
$input->analyze($file);


$input->writepdb($outputname);


print "Successfully wrote ", $outputname, " file.\n" if $debug >= 0;

print "\nPDBCORE ends.\n\n" if $debug >= 0;


## Subroutines ##

# display help on usage
sub help {
	system("perldoc /Users/YoshitakaM/apps/taopackage/scripts/pdbcore.pl");
	exit;
}

sub usage {
        print "\nPDBCORE $version : PDB -> smaller PDB\n";
        print "\nUsage: pdbcore [options]\n";
        print "\t-d? \t\tDebug print (?:0~3) \n";
	print "\t-h\t\tPrint full documentation\n";
	print "\t-i filename\tA valid pdb file.\n";
        print "\t-near number\tDistance in angstrom to define move region within core residues\n";
	print "\t-o filename\tOutput file name, defaults to newONIOMFile.gjf\n";
	print "\t-q\t\tQuiet mode\n";
        print "\t-resid listfile\tFile that contains the list of core residue numbers\n";
	print "\n";
        exit;
}


1;

=head1 EXAMPLES

=over

=item pdbcore

Called without any parameters, PDBCORE will display usage information. 
If B<-h> or B<--help> is passed then the full PDBCORE documentation is displayed via perldoc.

=item pdbcore -resid listfile -near 14 -i foo.pdb -o foosmall.pdb

PDBCORE takes all the residues which have at least one atom within 14 Angstroms of any atom 
in the core region and saves to foosmall.pdb.

=back

=head1 NOTES

When running PDBCORE, there may be some warning signs that some atoms or residues cannot
be recognized.  Therefore, atom types or partial charges cannot be assigned. These warnings
usually do not affect the output PDB file. Please check the output PDB file.

=head1 VERSION

1.0

=head1 AUTHOR

Peng Tao, E<lt>tao.21@osu.eduE<gt>

=head1 COPYRIGHT

Copyright (c) 2009~2010 by Peng Tao


=cut


