#!/usr/bin/perl 

use lib "/Users/YoshitakaM/apps/taopackage";
use ESPT::Glog 0.07;
use ESPT::Oniom 0.04;
use ESPT::Oniominput 0.02;
use strict;
use warnings;

=head1 NAME

ONIOMLOG - Gaussian ONIOM job monitoring tool.

=head1 SYNOPSIS

B<oniomlog> [ B<-d?> (?:0~3) ] [ B<-fn> number ] [ B<-fo> F<filename> ] [ B<-g> ] [ B<-h> ] [ B<--help> ]  [ B<-ha> ]
[ B<-i> F<filename> ] [ B<-l> [123] ] [ B<-o> ] [ B<-oi> ] [ B<-q> ] [ B<-r> ] [ B<-s> F<filename> ] [ B<-t> F<filename> ] 

=head1 DESCRIPTION

This program is a summary and extraction tool for Gaussian ONIOM jobs. Its main purpose
is to monitor and analyze running and finished Gaussian ONIOM jobs. 


=begin comment

### Version History ###
 1.0    Deal with ONIOM input and output file	

### 2 Do List ###

=end comment

=cut

### Main Program ###
our $version = "1.0";

# check for arguments
&usage if ( $#ARGV < 0 );
&help  if $ARGV[0] eq "-h";
&help if $ARGV[0] eq "--help";

our( $debug, $energy, $units, $xyz, $InHatree, $savestruct, $haveInput);
our( $file, $calc, $gjftemplate, $opt, $structure, $raw);
our( $structurefilename, $moveflag, $verbose, $templateinput);
our( $oniominput, $gjfoutput, $structurenumber,$layerlevel);



$opt = 0;
$structure = 0;
$InHatree  = 0;
$oniominput = 0 ;
$savestruct = 0;
$haveInput  = 0;

# parse arguments
for (my $i=0; $i<=$#ARGV; $i++) {
	$debug = 1 if $ARGV[$i] eq "-d";
        $debug = -1 if $ARGV[$i] eq "-q";
	$debug = 0 if $ARGV[$i] eq "-d0";
        $debug = 1 if $ARGV[$i] eq "-d1";
        $debug = 2 if $ARGV[$i] eq "-d2";
        $debug = 3 if $ARGV[$i] eq "-d3";
        $debug = 4 if $ARGV[$i] eq "-d4";
        $debug = 5 if $ARGV[$i] eq "-d5";
        $debug = 6 if $ARGV[$i] eq "-d6";
        $raw   =  1 if $ARGV[$i] eq "-r";
#       $xyz = 1  if $ARGV[$i] eq "-xyz";
        $layerlevel = $ARGV[$i+1] if $ARGV[$i] eq "-l";
        $opt = 1 if $ARGV[$i] eq "-o";
        $InHatree = 1 if $ARGV[$i] eq "-ha";
        $structure = 1 if $ARGV[$i] eq "-s";
        $structurefilename = $ARGV[$i+1] if $ARGV[$i] eq "-s";
        $savestruct = 1 if $ARGV[$i] eq "-s";
        $gjftemplate = $ARGV[$i+1] if $ARGV[$i] eq "-t";
        $templateinput = 1 if $ARGV[$i] eq "-t";
        $moveflag = 1 if $ARGV[$i] eq "-g";
        $verbose = 1 if $ARGV[$i] eq "-v";
        $oniominput = 1 if $ARGV[$i] eq "-oi";
        $gjfoutput = $ARGV[$i+1] if $ARGV[$i] eq "-fo";
        $structurenumber = $ARGV[$i+1] if $ARGV[$i] eq "-fn";
        
        $file = $ARGV[$i+1]  if $ARGV[$i] eq "-i";
        $haveInput = 1   if $ARGV[$i] eq "-i";
}

# set defaults
$debug ||= 0;
$energy ||= 0;
$moveflag ||= 0;
$verbose ||= 0;
$templateinput ||= 0;
$gjfoutput ||="newONIOMinput.gjf";
$structurenumber ||= 0;
$layerlevel ||= 0;
$raw   ||= 0;

print "\nONIOMLOG $version : Summary and monitor tool for Gaussian ONIOM log files.\n\n"  if $debug >= 0;

if ($haveInput == 0) { die "Please use flag -i to provide a Gaussian ONIOM log file.\n$!\n";}

# If none of -l or -g are used, choose QM region (layer 1) 
if (($moveflag == 0) && ($layerlevel==0)) {$layerlevel=1};

if (($structure == 1) && ($templateinput == 0)) {
  print "Structure needs to be extracted. But no template file is provided.\n"  if $debug >= 0;
  print "Try to extract ONIOM job info from given Gaussian Log file $file.\n"  if $debug >= 0;
}

if (($oniominput == 1) && ($templateinput == 0)) {
  print "A Gaussian ONIOM job file is needed when generating a new input file.\n"  if $debug >= 0;
  print "Use flag -t to provide a template file.\n"  if $debug >= 0;
  die "$!\n";
}

=head1 OPTIONS

Command line option specifications are processed from left to right and may
be specified more than once. If conflicting options are specified, later
specifications override earlier ones. The input filename passed to ONIOMLOG must be
a Gaussian ONIOM log file (partial or complete). 

=over 16

=item B<-d?> (?:0~3)

Turn on debug printing. The printing level can be controlled by a given number.
The larger the number, the more information will be printed when the program is running.


=item B<-fn> number

Specify the structrue from geometry optimization job log file. Default is the last output structure (0).


=item B<-fo> F<filename>

Write a Gaussian ONIOM input file with the new coordinates. Default is newONIOMinput.gjf.


=item B<-g>

Extract atoms which are allowed to move during the optimization. 


=item B<-h>

=item B<--help>

Print usage information. Cannot be used with other options.


=item B<-ha>

Using absolute (not relative) energies in hatree for optimization path output


=item B<-i> F<filename>

Gaussian ONIOM log file.


=item B<-l> [123]

Choose layers for structure extraction. 1: H layer, 2: M layer, 3: L layer.
It can be any combination of 123 in any order:1, 2, 3, 12, 213, 13, 23, etc.
Default is 0 when -g is used. If -l or -g are not used, level will be set to 1 (QM region).
In this way, users have more flexibility to choose any part of the system. 


=item B<-o>

Energy of each optimized step will be extracted and printed out. 
If both -s and -o are used, structures along the optimization path will also be 
extracted and saved to the file (file name given by -s flag).
The structure file can be loaded by a common visualization program (VMD for example) 
and played as a movie.


=item B<-oi>

Write a new Gaussian ONIOM input file. This provides a method to quickly 
create a Gaussian input file based on a template input file and Gaussian 
log file. It takes coordinates from the Guassian log file and replaces 
the coordinates in the template file with those from the log file. 
In this way, the setup of the old input file (atom type,
partial charges, link atoms, flag for freeze etc.) can be kept.


=item B<-q>

Run in quiet mode and do not print progress messages.


=item B<-r>

When used, all raw figures from Gaussian log file are displayed.
By default, energies are displayed up to six decimal places.
Dipole moments are displayed up to four decimal places.
Partial charges are displayed up to six decimal places.


=item B<-s> F<filename>

Specifies the portion of the structure to be extracted and output (as XYZ file) 
to a file with the given file name. When this flag is not used, no structure 
will be extracted. When flag -s is used, a filename must be provided. 

=item B<-t> F<filename>

Gaussian ONIOM input file to be used as a template.


=back


=cut

       
if ($oniominput == 1){
    $calc = ESPT::Oniominput->new();
    $calc->{FILENAME}=$file;
    $calc->{GJFINPUT}=$gjftemplate;
    $calc->{GJFOUTPUT}=$gjfoutput;
    $calc->{STRUCTURENUMBER}=$structurenumber;
    $calc->{VERBOSE}=$verbose;
    $calc->{DEBUG}=$debug;
} else {
    $calc = ESPT::Oniom->new();
    $calc->{FILENAME}=$file;
    $calc->{GJFTEMPLATE}=$gjftemplate;
    $calc->{OPT}=$opt;
    $calc->{InHatree} = $InHatree;
    $calc->{STRUCTURE}=$structure;
    $calc->{STRUCTUREFILE}=$structurefilename;
    $calc->{MOVEFLAG}=$moveflag;
    $calc->{VERBOSE}=$verbose;
    $calc->{TEMPLATEINPUT}=$templateinput;
    $calc->{LAYERINFO}=$layerlevel;
    $calc->{DEBUG}=$debug;
    $calc->{RAW}=$raw;
    }


#print $calc->get("PROGRAM"), " ", $calc->get("TYPE"), "\n" if $debug >= 1;

#print $gjftemplate;

# analyze file and extract standard data

$calc->analyze($file);



print "\nONIOMLOG ends.\n\n" if $debug >= 0;

exit;



## Subroutines ##



# Job completion check
sub chkcomplete {
	print "\nThis calculation is not complete!!\n" if $calc->get("COMPLETE") == 0;
}

# display help on usage
sub usage {
        print "\nONIOMLOG $version : Summary and monitor tool for Gaussian ONIOM log files.\n";
        print "\nUsage: oniomlog [options] \n";
        print "\t-d? \t\tDebug print (?:0~3) \n";
        print "\t-fn number\tNumber of structrue used from log file. Default is the last one (0).\n";
        print "\t-fo Filename\tGaussian ONIOM input file with new coordidates. Default(newONIOMinput.gjf)\n";
        print "\t-g \t\tExtract atoms which are allowed to move during optimization.\n";
	print "\t-h \t\tPrint usage information\n";
        print "\t-ha \t\tUsing absolute (not relative) energies in hatree in optimization path output\n";
        print "\t-i  filename \tGaussian ONIOM job log file\n\n";
        print "\t-l [123]\tChoose layers for structure extraction: 1:H, 2:M, 3:L\n";
        print "\t-o \t\tTreat Gaussian ONIOM calculation with multiple optimization steps\n";
        print "\t-oi \t\tTreat Gaussian ONIOM input file.\n";
        print "\t-q \t\tTurn on quiet mode.\n";
	print "\t-r\t\tDisplay raw data. (Keep as many decimal places as available)\n";
        print "\t-s Filename\tOutput extracted structure from ONIOM calculations\n";
        print "\t-t Filename \tGaussian ONIOM input file as template\n\n";
        exit;
}



# display help using perldoc option
sub help {
	system("perldoc /Users/YoshitakaM/apps/taopackage/scripts/oniomlog.pl");
        exit;
}

1;

=head1 EXAMPLES

=over 

=item oniomlog

Called without any parameters, ONIOMLOG will display usage information. 
If B<-h> or B<--help> is passed then the full ONIOMLOG documentation is displayed via perldoc.

=item oniomlog -i long.log

Check the energy of a Gaussian ONIOM job


=item oniomlog -o -i long.log

Check the energies of a Gaussian ONIOM job with optimization information


=item oniomlog -s oniomqm.xyz -o -i oniom.log

Extract QM region geometries for all of the structures along an optimization from 
a Gaussian ONIOM job (oniom.log). and save to oniom.xyz. The ONIOM setup 
information is taken from the Gaussian ONIOM log file (oniom.log).


=item oniomlog -s oniom.xyz -t oniominput.gjf -g -i oniom.log

Extract all the moving parts from the last structure in a Gaussian ONIOM job (oniom.log)
and save to oniom.xyz. The ONIOM setup information is taken from the Gaussian ONIOM 
input file (oniominput.gjf). oniominput.gjf will be used as a template for
the Gaussian ONIOM job setup (e.g. the moving flag). If the flag -t is not used to
provide a template file, the Gaussian ONIOM job setup information will be obtained
from oniom.log.


=item oniomlog -s oniom.xyz -l 2 -t oniominput.gjf -g -o -i oniom.log

Extract only layer 2 for all of the structures along the optimization from a 
Gaussian ONIOM job (oniom.log) and save to oniom.xyz. 


=item oniomlog -oi -t oniominput.gjf -fo new.gjf -i oniom.log

Create a new Gaussian ONIOM input file (new.gjf) from a Gaussian job  (oniom.log). 
A Gaussian ONIOM input file (oniominput.gjf) is needed as a template.


=item oniomlog -oi -fn 3 -t oniominput.gjf -fo new3.gjf -i oniom.log

Create a new Gaussian ONIOM input file (new.gjf) using the third structure along
the optimization from a Gaussian job  (oniom.log). A Gaussian ONIOM input file 
(oniominput.gjf) is needed as a template.


=back


=head1 VERSION

1.0


=head1 AUTHOR

Peng Tao, E<lt>Tao.21@osu.eduE<gt>


=head1 COPYRIGHT

Copyright (c) 2009~2010 by Peng Tao


=cut
