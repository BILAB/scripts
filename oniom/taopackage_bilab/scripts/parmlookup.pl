#!/usr/bin/perl

use lib "/Users/YoshitakaM/apps/taopackage/";
use File::Basename;
use ESPT::GoutOniom 0.01;
use strict;

=head1 NAME

PARMLOOKUP - Looks up missing force field parameters for an ONIOM calculation.

=head1 SYNOPSIS

B<parmlookup> [ B<-d?> (?:0~3) ] [ B<-g> F<Gaussian_ONIOM_log_file> ] [ B<-h> ] [ B<-q> ]  [ B<-o> F<parmfilename> ]

=head1 DESCRIPTION

This program looks up missing AMBER force field parameters for ONIOM calculations 
based on a given Gaussian output and AMBER parameters files. It goes through the AMBER
force field parameter file for amino acids first, then the GAFF (General AMBER Force Field).
 
=cut

### Version History ###
# 1.0   Look up missing paramters reported by Gaussian.
#	
#
# 2009 Peng Tao

### Main Program ###
our $version = "1.0";

# check for arguments
usage() if ( $#ARGV < 0 );

help() if $ARGV[0] eq "-h";
help() if $ARGV[0] eq "--help";

our ($goutfile, $outparmfile, $input, $isgout, $debug);
our $isall = 0; # M.Yamada

# parse arguments
for (my $i=0; $i<=$#ARGV; $i++) {
	$goutfile = $ARGV[$i + 1] if $ARGV[$i] eq "-g";
	$isall = 1 if $ARGV[$i] eq "-all"; # M.Yamada
	$debug = -1 if $ARGV[$i] eq "-q";
        $debug = 1 if $ARGV[$i] eq "-d";
	$debug = 0 if $ARGV[$i] eq "-d0";
        $debug = 1 if $ARGV[$i] eq "-d1";
        $debug = 2 if $ARGV[$i] eq "-d2";
        $debug = 3 if $ARGV[$i] eq "-d3";
        $debug = 4 if $ARGV[$i] eq "-d4";
        $debug = 5 if $ARGV[$i] eq "-d5";
        $debug = 6 if $ARGV[$i] eq "-d6";
        $outparmfile = $ARGV[$i+1] if $ARGV[$i] eq "-o";
} 

# set defaults
$debug ||= 0;
$outparmfile ||= "NewParm.txt";
        
print "\nPARMLOOKUP $version Look up missing parameters for ONIOM calculation.\n\n" if $debug >= 0;

if (not (defined $goutfile)) {
 print "No Gaussian ONIOM log file was given. Please use flag -g to provide such a file.\n";
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


=item B<-g> F<Gaussian_ONIOM_log_file> 

Gaussian log file. The missing parameter messages should be located 
at the end of this file. Other information obtained from this file includes all the atoms and
their types.


=item B<-h>

=item B<--help>

Print full PARMLOOKUP documentation via perldoc. Cannot be used with other options.


=item B<-q>

Run in quiet mode and do not print progress messages.


=item B<-o> F<parmfilename>

Output file containing missing parameters. Default is NewParm.txt.


=back

=cut

# get path and base filename
(my $base, my $dir, my $ext) = fileparse($goutfile, qr/\.[log]*/);

# check for gaussian output file
open(FILEIN,$goutfile) || die "Could not read $goutfile\n$!\n";

print "Opening $goutfile for process\n" if $debug >= 0;

$isgout = 0;
# determine file type and set extension
while (<FILEIN>){
        # skip blank lines
        next if /^$/;

	# Gaussian 
	if ( /^\s+Entering\s+Gaussian\s+System/ ){
		$input = ESPT::GoutOniom->new();
                $isgout = 1 ;
                $input->{PARMLOOKUP} = 1;
		last;
	}
	
}
close(FILEIN);

if ( $isgout == 0 and $isall == 0 )  { die "Given file $goutfile is NOT gaussian output file.\n Please check.\n";}


#print "Processing a ", $input->get("PROGRAM"), " ", $input->get("TYPE"), " file.\n" if $debug >= 0;

$input->debug($debug);
$input->{isall} = $isall; # M.Yamada

$input->analyze($goutfile);

$input->lookupparms;

$input->parmreduction;

$input->outputparm($outparmfile);

# $input->writeoutput($goutfile);


print "\nSuccessfully wrote ", $outparmfile, " file.\n" if $debug >= 0;

print "\nPARMLOOKUP ends.\n\n" if $debug >= 0;

## Subroutines ##

# display help on usage
sub help {
	system("perldoc /Users/YoshitakaM/apps/taopackage/scripts/parmlookup.pl");
	exit;
}

sub usage {
        print "\nPARMLOOKUP $version : Look up missing parameters for ONIOM calculation.\n";
        print "\nUsage: parmlookup [options] \n";
        print "\t-d? \t\tDebug print (?:0~3) \n";
	print "\t-g log file \tGaussian ONIOM output file.\n";
	print "\t-h\t\tPrint full documentation\n";
	print "\t-q\t\tQuiet mode\n";
	print "\t-o parmfilename\tOutput file contains found parameters. Default is NewParm.txt\n";
	print "\n";
        exit;
}



1;

=head1 EXAMPLES

=over

=item parmlookup

Called with no parameters at all, PARMLOOKUP will display usage information. 
If B<-h> or B<--help> is passed, then the full PARMLOOKUP documentation is displayed via perldoc.

=item parmlookup -g foo.log -o fooparm.txt

PARMLOOKUP reads foo.log, then looks up missing parameters in the AMBER Force Field files and 
saves hem to fooparm.txt.


=back

=head1 NOTES

If PARMLOOKUP cannot find a parameter in the AMBER force field files, a zero value is
displayed. A user can estimate values for parameters that cannot be found by looking
up corresponding parameters for atoms with similar atom types. 

=head1 VERSION

1.0

=head1 AUTHOR

Peng Tao, E<lt>tao.21@osu.eduE<gt>

=head1 COPYRIGHT

Copyright (c) 2009~2010 by Peng Tao


=cut




