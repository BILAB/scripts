#!/usr/bin/perl

use lib "/Users/YoshitakaM/apps/taopackage/";
use File::Basename;
use ESPT::ONIOMGEOM 0.02;
use strict;

=head1 NAME

TRANSGEOM - Transfers optimized geometry from a model system to a production run input
   
=head1 SYNOPSIS

B<Transgeom> [ B<-d?> (?:0~3) ] [ B<-h> ] [ B<-q> ] [ B<-atmp> ] [ B<-atmtypp> ] [ B<-chgp> ] [ B<-coordp> ]
[ B<-layerp> ] [ B<-movep> ] [ B<-t> tolerance ] F<filename>

=head1 DESCRIPTION

This program extracts the atomic setup information from calculations of a model system 
and prepares a production run input
 
=cut

### Version History ###
# 1.0   Open setup file and read 5 file names:
#     Model system ONIOM input with optimized geometry.
#     ONIOM .onb file for model system.
#     Production ONIOM input file.
#     ONIOM .onb file for production system.
#     ONIOM input file name for production system with optimized geometry
#
# 1.1 User can choose take either coordinate, atom, atom type, partial charge,
#     moving flag, layer information, or any combination of them from model
#     system to production system.
#
# #
# 2009 Peng Tao

### To Do List ###
#

### Main Program ###
our $version = "1.1";

# check for arguments
usage() if ( $#ARGV < 0 );

help() if $ARGV[0] eq "-h";
help() if $ARGV[0] eq "--help";


our ($debug, $extension, $file, $tol,$atomprod,$atmtypprod,$chgprod,$moveprod,$layerprod,$coordprod,$name);

# parse arguments
for (my $i=0; $i<=$#ARGV; $i++) {
	$debug = 1 if $ARGV[$i] eq "-d";
	$debug = 0 if $ARGV[$i] eq "-d0";
        $debug = 1 if $ARGV[$i] eq "-d1";
        $debug = 2 if $ARGV[$i] eq "-d2";
        $debug = 3 if $ARGV[$i] eq "-d3";
        $debug = 4 if $ARGV[$i] eq "-d4";
        $debug = 5 if $ARGV[$i] eq "-d5";
        $debug = 6 if $ARGV[$i] eq "-d6";
	$debug = -1 if $ARGV[$i] eq "-q";
        $tol = $ARGV[$i+1] if  $ARGV[$i] eq "-t";
        $atomprod = 1 if $ARGV[$i] eq "-atmp";   # if set, use element from production ONIOM input file
        $atmtypprod = 1 if $ARGV[$i] eq "-atmtypp"; # if set, use atom type from production ONIOM input file
        $chgprod = 1 if $ARGV[$i] eq "-chgp"; # if set, use partial charges from production ONIOM input file
        $moveprod = 1 if $ARGV[$i] eq "-movep"; # if set, use moving flag from production ONIOM input file
        $layerprod = 1 if $ARGV[$i] eq "-layerp"; # if set, use layer setup from production ONIOM input file
        $coordprod = 1 if $ARGV[$i] eq "-coordp"; # if set, use coordinates from production ONIOM input file
        $file = $ARGV[$i] if $i == $#ARGV;
} 

# set defaults
$debug ||= 0;
$tol ||= 0.00001;
$atomprod ||= 0;
$atmtypprod ||= 0;
$chgprod ||= 0; 
$moveprod ||= 0;
$layerprod ||= 0;
$coordprod ||= 0;

print "\nTRANSGEOM $version : Put optimized geometry from model system to production run input\n" if  $debug >= 0;
 
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

Print full TRANSGEOM documentation via perldoc. Cannot be used with other options.


=item B<-q>

Run in quiet mode and do not print progress messages.


=item B<-atmp>

If set, use the element information from the production ONIOM input file.
Default uses model system input file.


=item B<-atmtypp>

If set, use the atom type information from the production ONIOM input file.
Default uses model system input file.


=item B<-chgp>

If set, use the partial charges information from the production ONIOM input file.
Default uses model system input file.


=item B<-coordp> 

If set, use the coordinate information from the production ONIOM input file.
Default uses model system input file.


=item B<-movep>

If set, use the move flag information from the production ONIOM input file.
Default uses model system input file.


=item B<-layerp>

If set, use the layer setup information from the production ONIOM input file.
Default uses model system input file.


=item B<-t> tolerance

Tolerance to decide if an atom has been moved during the optimization.
Default is 0.00001;


=item F<filename>

  Format:
   Modelgjf         ONIOM input file name for model system
   Modelonb         ONIOM ONB file for model system
   Productiongjf    ONIOM input file name for production run system
   Productiononb    ONIOM ONB file for production run system
   Productioninput  New ONIOM input file for production run system (Default is Prod_ONIOM_NEW.gjf)

These files can be specified in any order.

The parameter filename must be the last command line argument when running transgeom.

=back

=cut

open(FILEIN,$file) || die "Could not read $file\n$!\n";

# get path and base filename

our ($input);

$input = ESPT::ONIOMGEOM->new();   
$input->{DEBUG} = $debug;

$input->{tolerance}=$tol;
$input->{atomprod} = $atomprod;
$input->{atmtypprod} = $atmtypprod;
$input->{chgprod} = $chgprod;
$input->{moveprod} = $moveprod;
$input->{layerprod} = $layerprod;
$input->{coordprod} = $coordprod;

my $defaultprodname = $input->{DefaultProdNewInput};

print "\nPut optimized geometry from model system to production run ONIOM input file.\n\n" if $debug >= 0;

print "\nTolerance to decide if one atom moved or not during optimization is $tol.\n" if $debug >= 0;

# Read model and production oniom input and ONB files
while(<FILEIN>) {
       next if /^$/;
       
       if (/^\s*modelgjf/i) {
         /^\s*modelgjf\s+(\S+)/i;
         $input->{ModelGjfFile}=$1;
         print "Model system ONIOM input file name is \n\"$1\"\n\n" if $debug >= 1;
         next;
        }

       if (/^\s*modelonb/i) {
         /^\s*modelonb\s+(\S+)/i;
         $input->{ModelONBFile}=$1;
         print "Model system ONB template file name is \n\"$1\"\n\n" if $debug >= 1;
         next;
        }

       if (/^\s*productiongjf/i) {
         /^\s*productiongjf\s+(\S+)/i;
         $input->{ProdGjfFile}=$1;
         print "Production system ONIOM input file name is \n\"$1\"\n\n" if $debug >= 1;
         next;
        }
        
       if (/^\s*productiononb/i) {
         /^\s*productiononb\s+(\S+)/i;
         $input->{ProdONBFile}=$1;
         print "Production system ONB template file name is \n\"$1\"\n\n" if $debug >= 1;
         next;
        }

       if (/^\s*productioninput/i) {
         /^\s*productioninput\s+(\S+)/i;
         $input->{ProdNewInput}=$1;
         print "Production system ONIOM input file name with optimized geometry from model system is \n\"$1\"\n\n" if $debug >= 1;
         next;
        }

} # while() {

close(FILEIN);

# Check if all four necessary file names are obtained.

if ( not defined($input->{ModelGjfFile}) ) 
  {print "Model system ONIOM input file name not found. \nExit!\n";die;}

if ( not defined($input->{ModelONBFile})) 
  {print "Model system ONB template file name not found. \nExit!\n";die;}

if ( not defined($input->{ProdGjfFile})) 
  {print "Production system ONIOM input file name not found. \nExit!\n";die;}

if ( not defined($input->{ProdONBFile}) ) 
  {print "Production system ONB template file name not found. \nExit!\n";die;}
  
if ( not defined($input->{ProdNewInput}) ) 
  {print "New Production system ONIOM input file name not found. Use default ($defaultprodname)\n";
   $input->{ProdNewInput}=$defaultprodname;}
  

# read file contents
$input->debug($debug);
$input->analyze();


$input->writeprodoniom();


print "\nSuccessfully wrote ", $input->{ProdNewInput}, " file.\n" if $debug >= 0;


print "\nTRANSGEOM ends.\n\n" if $debug >= 0;

## Subroutines ##

# display help on usage
sub help {
	system("perldoc /Users/YoshitakaM/apps/taopackage/scripts/transgeom.pl");
	exit;
}

sub usage {
        print "\nTRANSGEOM $version : Put optimized geometry from model system to production run input\n";
        print "\nUsage: Transgeom [options] filename \n";
        print "\t-d? \t\tDebug print (?:0~3) \n";
	print "\t-h\t\tprint full documentation\n";
	print "\t-q\t\tQuiet mode\n";
        print "\t-atmp\t\tuse element from production ONIOM input file\n";
        print "\t-atmtypp\tuse atom type from production ONIOM input file\n";
        print "\t-chgp\t\tuse partial charges from production ONIOM input file\n";
        print "\t-coordp\t\tuse coordinates from production ONIOM input file\n";
        print "\t-layerp\t\tuse layer setup from production ONIOM input file\n";
        print "\t-movep\t\tuse moving flag from production ONIOM input file\n";
        print "\t-t\t\ttolerance\n";
      
	print "\tfilename\tSetup file.\n";

	print "\t\t\tSetup file format:\n";
	print "\t\t\t   Modelgjf         ONIOM input file name for model system\n";
	print "\t\t\t   Modelonb         ONIOM ONB file for model system\n";
	print "\t\t\t   Productiongjf    ONIOM input file name for production run system\n";
	print "\t\t\t   Productiononb    ONIOM ONB file for production run system\n";
	print "\t\t\t   Productioninput  New ONIOM input file for production run system\n";
	print "\n";
        exit;
}


1;

=head1 EXAMPLES

=over

=item transgeom

Called without any parameters, TRANSGEOM will display usage information. 
If B<-h> or B<--help> is passed, then the full TRANSGEOM documentation is displayed via perldoc.

=item transgeom foo.in

TRANSGEOM reads foo.in, then creates a production run input. 
Please note that there is no flag needed before setup file name in command line.
Sample foo.in
   Modelgjf         ONIOMmodelinput.gjf
   Modelonb         ONIOMmodelinput.onb
   Productiongjf    ONIOMproductioninput.gjf
   Productiononb    ONIOMproductioninput.onb
   Productioninput  ONIOMproductioninputwithmodelgeom.gjf


=back

=head1 NOTES

The model system is not necessary smaller than the production system. A full protein model 
can be used as a model system, and a reduced size model for the same protein can be used 
as the production system. TRANSGEOM will extract corresponding portion of the geometry 
from the full size model and build a partial model.


=head1 VERSION

1.1

=head1 AUTHOR

Peng Tao, E<lt>tao.21@osu.eduE<gt>


=head1 COPYRIGHT

Copyright (c) 2009~2010 by Peng Tao



=cut



