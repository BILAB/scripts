#!/usr/bin/perl

use lib "/Users/YoshitakaM/apps/taopackage/";
use File::Basename;
use ESPT::GinputConnect 0.01;
use strict;

=head1 NAME

CHECKCONNECT - Searches for isolated atoms and atoms with multiple connections in an Gaussian ONIOM input file.

=head1 SYNOPSIS

B<checkconnect> [ B<-c> number ] [ B<-d?> (?:0~3) ] [ B<-g> F<Gaussian_job_file_name> ] [ B<-h> ] [ B<-q> ]  
                 
=head1 DESCRIPTION

This program looks through the connectivity table of the given Gaussian input file
and searches for isolated atoms and atoms with multiple connections as sanity 
check for Gaussian ONIOM jobs.
 
=cut

### Version History ###
# 1.0   Read into Gaussian input file
#
# 2009 Peng Tao

### Main Program ###
our $version = "1.0";

# check for arguments
usage() if ( $#ARGV < 0 );

help() if $ARGV[0] eq "-h";
help() if $ARGV[0] eq "--help";

our ($gfile, $outparmfile, $input, $isgout, $debug, $cnumber);


# parse arguments
for (my $i=0; $i<=$#ARGV; $i++) {
	$gfile = $ARGV[$i + 1] if $ARGV[$i] eq "-g";
        $cnumber = $ARGV[$i + 1] if $ARGV[$i] eq "-c";
	$debug = -1 if $ARGV[$i] eq "-q";
	$debug = 1 if $ARGV[$i] eq "-d";
	$debug = 0 if $ARGV[$i] eq "-d0";
        $debug = 1 if $ARGV[$i] eq "-d1";
        $debug = 2 if $ARGV[$i] eq "-d2";
        $debug = 3 if $ARGV[$i] eq "-d3";
        $debug = 4 if $ARGV[$i] eq "-d4";
        $debug = 5 if $ARGV[$i] eq "-d5";
        $debug = 6 if $ARGV[$i] eq "-d6";
} 

# set defaults
$debug ||= 0;
$cnumber ||= -1;

print "\nCHECKCONNECT $version : Gaussian ONIOM input file sanity check for connections.\n\n" if $debug >= 0;

if (not (defined $gfile)) {
 print "No Gaussian ONIOM input file was given. Please use flag -g to provide such a file.\n";
 die "Exit.\n$!";
}



=head1 OPTIONS

Command line option specifications are processed from left to right and may 
be specified more than once. If conflicting options are specified, later  
specifications override earlier ones.

=over 16

=item B<-c> number

Any atoms with the same or more connections are printed. If this option is 
not used, then all the atoms with the most connections are printed.

=item B<-d?> (?:0~3)

Turn on debug printing. The printing level can be controlled by a given number.
The larger the number, the more information will be printed when the program is running.


=item B<-g> Gaussian_job_file_name  

Gaussian ONIOM input file.


=item B<-h>

=item B<--help>

Print full CHECKCONNECT documentation via perldoc. Cannot be used with other options.


=item B<-q>

Run in quiet mode and do not print progress messages.

=back

=cut

# get path and base filename
# (my $base, my $dir, my $ext) = fileparse($gfile, qr/\.[log]*/);

# check for Gaussian output file
open(FILEIN,$gfile) || die "Could not read $gfile\n$!\n";

print "Opening $gfile for processing...\n" if $debug >= 0;

$isgout = 0;

# determine file type and set extension
while (<FILEIN>){
        # skip blank lines
        next if /^$/;

	# Gaussian 
	if ( /^\s+Entering\s+Gaussian\s+System/ ){
                print "Given file $gfile is Gaussian log file.\n";
                print "Current version of CheckConnect $version only deal with Gaussian input file. \n";
                print "Use Gaussian input file for connectivity check.\n";
                die "Exit\n";
	}
	
}
close(FILEIN);

if ( $isgout == 0 )  { 
   
   print "Treat the file as a Gaussian input file.\n";
   $input = ESPT::GinputConnect->new();
   $input->{Cnumber} = $cnumber;

} # if ( $isgout == 0 )  


#print "Processing a ", $input->get("PROGRAM"), " ", $input->get("TYPE"), " file.\n" if $debug >= 0;

$input->debug($debug);

$input->analyze($gfile);

print "\nCHECKCONNECT ends.\n\n" if $debug >= 0;


## Subroutines ##

# display help on usage
sub help {
	system("perldoc /Users/YoshitakaM/apps/taopackage/scripts/checkconnect.pl");
	exit;
}

sub usage {
        print "\nCHECKCONNECT $version : Look up isolated (not connected to any other atom) atoms in Gaussian input file.\n";
        print "\nUsage: chargeaddup [options]\n";
        print "\t-c\t\tconnection number\n";
        print "\t-d? \t\tDebug print (?:0~3) \n";
	print "\t-g Gaussian_file_name \tGaussian input file.\n";
	print "\t-h\t\tPrint full documentation\n";
	print "\t-q\t\tQuiet mode\n";
	print "\n";
        exit;
}


1;

=head1 EXAMPLES

=over

=item checkconnect

Called without any parameters, CHECKCONNECT will display usage information. 
If B<-h> or B<--help> is passed, then the full CHECKCONNECT documentation 
is displayed via perldoc.

=item checkconnect -g foo.gjf

CHECKCONNECT reads foo.gif, then finds any isolated atoms and atoms with
the most connections.

=item checkconnect -g foo.gjf -c 5

CHECKCONNECT reads foo.gif, then finds any isolated atoms and atoms with
connections equal to or larger than 5.

=back

=head1 NOTES

CHECKCONNECT

=head1 VERSION

1.0

=head1 AUTHOR

Peng Tao, E<lt>tao.21@osu.eduE<gt>

=head1 COPYRIGHT

Copyright (c) 2009~2010 by Peng Tao


=cut



