#!/usr/bin/perl 

use lib "/Users/YoshitakaM/apps/taopackage/";
use ESPT::GTable 0.01;
use strict;
use warnings;

=head1 NAME

GAUSSIANTABLE - Generates a comparison table for given Gaussian log files.

=head1 SYNOPSIS

B<gaussiantable> [ B<-d?> (?:1~3) ] [ B<-h> ] [ B<-i> F<input_file> ] [ B<-o> F<output_file> ] [ B<-q> ] [ B<-r> ]

=head1 DESCRIPTION

This program generates thermodynamic data based on the given Gaussian output file list.

=begin comment

### Structure ###


### Version History ###
 	
 1.0 	Treat ONIOM files

### 2 Do List ###


=end comment

=cut

### Main Program ###
our $version = "1.0";

# check for arguments
&usage if ( $#ARGV < 0 );
&help if $ARGV[0] eq "-h";
&help if $ARGV[0] eq "--help";

our( $debug, $spin, $energy, $freq, $firstMO, $lastMO, $mo, $units, $xyz );
our( $file, $calc, $proginfo, $oniom, $gjftemplate, $opt, $structure, $raw);
our( $structurefilename, $moveflag, $verbose, $templateinput);
our( $oniominput, $gjfinput, $gjfoutput, $structurenumber,$layerlevel);
our($inputfile, $outputfile);



$opt = 0;
$structure = 0;


# parse arguments
for (my $i=0; $i<=$#ARGV; $i++) {
	$debug = -1 if $ARGV[$i] eq "-q";
	$debug = 1 if $ARGV[$i] eq "-d";
	$debug = 0 if $ARGV[$i] eq "-d0";
        $debug = 1 if $ARGV[$i] eq "-d1";
        $debug = 2 if $ARGV[$i] eq "-d2";
        $debug = 3 if $ARGV[$i] eq "-d3";
        $debug = 4 if $ARGV[$i] eq "-d4";
        $debug = 5 if $ARGV[$i] eq "-d5";
        $debug = 6 if $ARGV[$i] eq "-d6";
        $raw   = 1 if $ARGV[$i] eq "-r";

        $inputfile = $ARGV[$i+1] if $ARGV[$i] eq "-i";
        $outputfile = $ARGV[$i+1] if $ARGV[$i] eq "-o";
}

# set defaults
$debug ||= 0;
$inputfile ||= "gtable.in";
$outputfile ||= "gtable.out";
$raw   ||= 0;

print "\nGAUSSIANTABLE $version : Generating comparing table for given gaussian files.\n\n" if $debug >= 0;

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

Print usage information. Cannot be used with other options.


=item B<-i> F<input_file>

Input file. If not given, GAUSSIANTABLE will try to read from gtable.in (Default).


=item B<-o> F<output_file>

Output file, Default is gtable.out.


=item B<-q> 

Quiet mode with least output during running.


=item B<-r>

When used, all raw values from the Gaussian file are displayed.
By default, energies are displayed with up to six decimal places.
Dipole moments are displayed with up to four decimal places.
Partial charges are displayed with up to six decimal places.


=back

=head1 Input File Format

GAUSSIANTABLE can take one input file in text form. Here are rules for the input file:
all blank lines are ignored, all lines beginning with # are ignored, and spaces before
# are overlooked.

Gaussian log files are organized by blocks.

The format for each block looks like:

<Label>Sulfoxide Thiirane 6-31G(d) ME </Label>E<13>E<10> <ONIOM> E<13>E<10> File1 E<13>E<10> File2 (full path) E<13>E<10> File3 (full path) E<13>E<10> </ONIOM>


<Label> and </Label> need to be in the same line.

GAUSSIANTABLE extracts B<Sulfoxide Thiirane 6-31G(d) ME> as the label for this block.
B<<ONIOM>> tells GAUSSIANTABLE that this is a block of ONIOM jobs.
B<</ONIOM>> tells GAUSSIANTABLE that the ONIOM block ends here.

GAUSSIANTABLE extracts ONIOM energies for each job, uses the result from File1 as
the reference, and outputs absolute energies (in I<hartree>) and relative energies
(in I<kcal/mol>) for all the files (including F<File1>).

There is no limitation on the number of files listed in each block. GAUSSIANTABLE will 
output absolute and relative energies (with respect to the first file in this block) 
for each file.
The files in each block should be either all 2-layer or all 3-layer ONIOM jobs. Mixture of
2- and 3-layer ONIOM jobs in a block will cause GAUSSIANTABLE to exit.

The file name with full path can be used, so multiple files in different
folders can be compared. Path information is not necessary for files located in the folder
in which GAUSSIANTABLE runs.

There is no limitation on the number of blocks listed in each input file. 
GAUSSIANTABLE will treat each block in the order that they appear in the input file 
and write output for each block in the same order.

The current version 1.0 of GAUSSIANTABLE can only analyze Gaussian ONIOM job files. 
In the future, it will be updated to treat regular Gaussian job (non-ONIOM) files. 

=cut

# Generate Gaussian Table object
 our $gtable = ESPT::GTable->new();
 
 $gtable->{InputFilename} = $inputfile;
 
 $gtable->{OutputFilename} = $outputfile;
 
 $gtable->{DEBUG} = $debug;
 
 $gtable->{RAW} = $raw;



# analyze file and extract standard data

print "Input file $inputfile is used to generate Gaussian Table. \n" if $debug > 0;

$gtable->analyze($file, $spin, $debug);


print "\nGAUSSIANTABLE ends.\n\n" if $debug >= 0;

print "\n" if $debug >= 0;

exit;



## Subroutines ##


# diplay program info
sub proginfo {
        print "GAUSSIANTABLE version: $version\tInput file: $file\n\n";
}

# Job completion check
sub chkcomplete {
	print "\nThis calculation is not complete!!\n" if $calc->get("COMPLETE") == 0;
}

# display help on usage
sub usage {
        print "\nGAUSSIANTABLE $version : Generating comparing table for given gaussian files.\n";
        print "\nUsage: gaussiantable [options] \n";
        print "\t-d? \t\tDebug print (?:0~3) \n";
	print "\t-h \t\tPrint usage information\n";
        print "\t-i Filename\tInput file.\n";
        print "\t-o Filename\tOutput file.\n";
        print "\t-q \t\tQuiet model\n";
	print "\t-r\t\tDisplay raw data. (Keep as many decimal places as available)\n\n";
        exit;
}


# display help using perldoc option
sub help {
	system("perldoc /Users/YoshitakaM/apps/taopackage/scripts/gaussiantable.pl");
	exit;
}


1;

=head1 EXAMPLES

=over

=item gaussiantable

Called without any parameters, GAUSSIANTABLE will display usage information. 
If B<-h> or B<--help> is passed, then the full GAUSSIANTABLE documentation is displayed via perldoc.

=item gaussiantable -i foo.in -o foo.out

GAUSSIANTABLE reads foo.in, then writes output to foo.out. 


=back


=head1 VERSION

1.0


=head1 AUTHOR

Peng Tao, E<lt>tao.21@osu.eduE<gt>

=head1 COPYRIGHT

Copyright (c) 2009~2010 by Peng Tao

=cut
