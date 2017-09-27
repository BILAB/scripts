#!/usr/bin/perl

use lib "/Users/YoshitakaM/apps/taopackage/";
use File::Basename;
use ESPT::GoutOniom 0.01;
use ESPT::GinputOniom 0.01;
use strict;

=head1 NAME

CHARGESUM - Provides a summary of the partial charge information for an ONIOM output or input file.

=head1 SYNOPSIS

B<chargesum> [ B<-d?> (?:0~3) ] [ B<-g> F<Gaussian_ONIOM_file> ] [ B<-h> ] [ B<-q> ] [ B<-r> ] 

=head1 DESCRIPTION

This program provides a summary of the partial charge information for an ONIOM output or input file.
 
=cut

### Version History ###
# 1.0   Add up partial charges of all atoms from Gaussian log file.
# 1.1   Add option to treat Gaussian input file
#
# 2007 Peng Tao

### Main Program ###
our $version = "1.1";

# check for arguments
usage() if ( $#ARGV < 0 );

help() if $ARGV[0] eq "-h";
help() if $ARGV[0] eq "--help";

our ($gfile, $outparmfile, $input, $isgout, $debug, $raw);


# parse arguments
for (my $i=0; $i<=$#ARGV; $i++) {
	$gfile = $ARGV[$i + 1] if $ARGV[$i] eq "-g";
	$debug = -1 if $ARGV[$i] eq "-q";
        $debug = 1 if $ARGV[$i] eq "-d";
	$debug =  0 if $ARGV[$i] eq "-d0";
        $debug =  1 if $ARGV[$i] eq "-d1";
        $debug =  2 if $ARGV[$i] eq "-d2";
        $debug =  3 if $ARGV[$i] eq "-d3";
        $debug =  4 if $ARGV[$i] eq "-d4";
        $debug =  5 if $ARGV[$i] eq "-d5";
        $debug =  6 if $ARGV[$i] eq "-d6";
        $raw   =  1 if $ARGV[$i] eq "-r";
        
} 

# set defaults
$debug ||= 0;
$raw   ||= 0;

print "\nCHARGESUM $version : Provides a summary of the partial charge information for an ONIOM output or input file.\n\n"  if $debug >= 0;

if (not (defined $gfile)) {
 print "No Gaussian ONIOM file was given. Please use flag -g to provide such a file.\n";
 die "Exit.\n$!";
}


=head1 OPTIONS

Command line option specifications are processed from left to right and may 
be specified more than once. If conflicting options are specified, later  
specifications override earlier ones. 

=over 16

=item B<-d?> (?:0~3)

Turn on debug printing. The printing level can be controlled by a given number.
The larger the number, the more information will be printed when program is running.

=item B<-g> F<Gaussian_ONIOM_file>

Gaussian ONIOM file, either input or output file.


=item B<-h>

=item B<--help>

Print full CHARGESUM documentation via perldoc. Cannot be used with other options.


=item B<-q>

Run in quiet mode and do not print progress messages.


=item B<-r>

When used, all raw values from the Gaussian file are displayed.
By default, energies are displayed with up to six decimal places.
Dipole moments are displayed with up to four decimal places.
Partial charges are displayed with up to six decimal places.

=back

=cut

# get path and base filename
#(my $base, my $dir, my $ext) = fileparse($gfile, qr/\.[log]*/);

# check for Gaussian output file
open(FILEIN,$gfile) || die "Could not read $gfile\n$!\n";

print "Opening $gfile for process\n\n" if $debug >= 0;

$isgout = 0;
# determine file type and set extension
while (<FILEIN>){
        # skip blank lines
        next if /^$/;

	# Gaussian 
	if ( /^\s+Entering\s+Gaussian\s+System/ ){
                print "Given file $gfile is Gaussian log file.\n"  if $debug >= 1;
		$input = ESPT::GoutOniom->new();
                $isgout = 1 ;
		last;
	}
	
}
close(FILEIN);

if ( $isgout == 0 )  { 
   
   print "Given file $gfile is not Gaussian log file.\nTreat it as Gaussian input file.\n"  if $debug >= 0;
   $input = ESPT::GinputOniom->new();

} # if ( $isgout == 0 )  

$input->{RAW} = $raw;

#print "Processing a ", $input->get("PROGRAM"), " ", $input->get("TYPE"), " file.\n" if $debug >= 0;

$input->debug($debug);

$input->analyze($gfile);

$input->writecharges;

$input->dipolemoment;

# $input->writeoutput($goutfile);

print "\nCHARGESUM ends.\n\n" if $debug >= 0;

## Subroutines ##

# display help on usage
sub help {
	system("perldoc /Users/YoshitakaM/apps/taopackage/scripts/chargesum.pl");
	exit;
}

sub usage {
        print "\nCHARGESUM $version : Provides a summary of the partial charge information for an ONIOM output or input file.\n";
        print "\nUsage: chargesum [options] \n";
        print "\t-d? \t\tDebug print (?:0~3) \n";
	print "\t-g Gaussian_ONIOM_file_name \tGaussian ONIOM input or output file.\n";
	print "\t-h\t\tPrint full documentation\n";
	print "\t-q\t\tQuiet mode\n";
	print "\t-r\t\tDisplay raw data. (Keep as many decimal places as available)\n";
	print "\n";
        exit;
}


1;

=head1 EXAMPLES

=over

=item chargesum

Called without any parameters, CHARGESUM will display usage information. 
If B<-h> or B<--help> is passed, then the full CHARGESUM documentation is 
displayed via perldoc.

=item chargesum -g foo.gjf

CHARGESUM reads foo.gif, then sums the total partial charges of each layer of 
the ONIOM calculation.


=back

=head1 NOTES

CHARGESUM can treat both two-layer and three-layer ONIOM calculations. Please interpret results based on
the setup of your system.

=head1 VERSION

1.1

=head1 AUTHOR

Peng Tao, E<lt>tao.21@osu.eduE<gt>

=head1 COPYRIGHT

Copyright (c) 2009~2010 by Peng Tao


=cut




