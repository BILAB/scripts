package ESPT::ONIOMLog;

use base qw(ESPT::ESSfile);
use ESPT::Glib 0.01;
use strict;
use warnings;

our $VERSION = '0.02';

=head1 NAME

ESPT::ONIOMLog - ONIOMLog object to treat Gaussian log file.

=head1 SYNOPSIS

   use ESPT::ONIOMLog;

   my $log = ONIOMLog->new();

=head1 DESCRIPTION

This module provides methods to quickly access data and structures contianed in a Gaussian log file
object for ONIOM claculation.  Guassian log files can only be read currently.

=cut

### Version History ###
# 0.01	Export ONIOM energy information
# 0.02  Treat three layer ONIOM calculation.


### To Do ###


## the object constructor **

sub new {
	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my $log = ESPT::ESSfile->new();

	$log->{TYPE} = "log";

	# program info
	$log->{PROGRAM} = "Gaussian";
	$log->{VERSION} = undef;
	$log->{REVISION} = undef;
	$log->{COMPILE} = undef;
	$log->{COMPILEDATE} = undef;
	
	# Link 0 & Route commands
	$log->{LINK0} = [];
	$log->{KEYWORDS} = [];
	$log->{ROUTE} = undef;

	# calc info
	$log->{CPUTIME} = undef;
	$log->{FUNCTIONAL} = undef;
	$log->{NPRIMITIVE} = undef;
	$log->{NCARTESIAN} = undef;
	$log->{OPTIMIZED} = undef;	# flag indicating successful optimization
	$log->{RUNTIME} = undef;

	
	# ONIOM info
	$log->{LOWMODEL} = undef;
	$log->{HIGHMODEL} = undef;
	$log->{LOWREAL} = undef;
        $log->{MEDMODLE} = undef;
        $log->{LOWMID} = undef;
        $log->{MEDMID} = undef;
	$log->{ONIOMENERGY} = undef;
	$log->{LOWREAL_MODEL} = undef;

        # Operation info
        $log->{GJFTEMPLATE} = undef;
	$log->{OPT} = 0;
	$log->{STRUCTURE} = 0;
        $log->{STRUCTUREFILE} = undef;
        $log->{MOVEFLAG} = 0;
        $log->{VERBOSE} = 0;
        $log->{TEMPLATEINPUT} = 0;
        $log->{ISNORMAL} = 0;
	bless($log, $class);
	return $log;
}


## methods ##

=item B<$file-E<gt>analyze(filename [spin])>

Analyze the spin results in file called filename.  Spin defaults to Alpha.

=cut

# set filename & spin then digest the file
sub analyze : method {
	my $log = shift;
	$log->{FILENAME} = shift;
        $log->{DEBUG} = shift;
	$log->digest();
	return;
}


## subroutines ##



sub digest {
# Files larger than 1Mb should be converted to binary and then
# processed to ensure maximum speed. -D. Ennis, OSC

# For items with multiple occurances, the last value is reported 

my $log = shift;

# flags & counters
my $Cflag = 0;
my $symmflag = 0;
my $rparsed = 0;
my $Sflag = 0;
my $Scount = 0;
my $Ccount = 0;
my $dcount = 0;
my $eigcount = 0;
my $Ecount = 0;
my $ESTATEcount = 0;
my $orbcount = -1;
my $MOcount = 0;
my $PGcount = 0;
my $LowReal_Model = 0;
my $HighModel = 0;
my $ONIOMEnergy = 0;
my $diffLowReal_Modelkcal = 0;
my $diffHighModelkcal = 0;
my $diffONIOMEnergykcal = 0;
my $HFtoKcalmol = 627.509;
my $optsteps = 0;
my @qmatomlist;
my @qmatomsymbollist;
my $atomindex = 0;
my $structindicator = 0;
my $qmatomnumber = 0;
my @currentXYZ;
my $currentline;
my $i=0;
my $j=0;
my $NAtoms=0;
my $structureindex = 0;
my $debug;

my $infocontrol = 1;


$debug = $log->{DEBUG};

print "Gaussian out file is ",$log->{FILENAME},"\n" if $debug >= $infocontrol ;

$log->{ISNORMAL} = 0;

# open filename for reading or display error
open(LOGFILE,$log->{FILENAME}) || die "Could not read $log->{FILENAME}\n$!\n";

# grab everything which may be useful
while (<LOGFILE>){
	# skip blank lines
	next if /^$/;

        if (/^\s+NAtoms=\s*(\d+)/) {
           $NAtoms=$1;
#           print "There are $NAtoms atoms.\n";
         }


        # ONIOM gridpoint  1 method:  low   system:  model energy
        if ( /^\s+ONIOM:\s+gridpoint\s+1\s+method:\s+low\s+system:\s+model\s+energy:\s+(.+)/ ) {
		$log->{LOWMODEL} = $1;
#                print "Current Low Model energy is ", $log->{LOWMOEDL},"\n";
		next;
        }
        
        # ONIOM gridpoint  2 method:  high   system:  model energy
        if ( /^\s+ONIOM:\s+gridpoint\s+2\s+method:\s+high\s+system:\s+model\s+energy:\s+(.+)/ ) {
		$log->{HIGHMODEL} = $1;
#                print "Current high Model energy is ", $log->{HIGHMOEDL},"\n";
#                print $_;
		next;
        }

        # ONIOM gridpoint  3 method:  low   system:  real energy
        if ( /^\s+ONIOM:\s+gridpoint\s+3\s+method:\s+low\s+system:\s+real\s+energy:\s+(.+)/ ) {
		$log->{LOWREAL} = $1;
#                print "Current low Real energy is ", $log->{LOWREAL},"\n";
#                print $_;
                $log->{LOWREAL_MODEL} = $log->{LOWREAL} - $log->{LOWMODEL};
#                print "Current low Real energy - low Model Energy is ",$LowReal_Model ,"\n";
                $log->{ISTHREELAYER} = 0;
		next;
        }

        # ONIOM gridpoint  2 method:  med   system:  model energy  for three layer ONIOM
        if ( /^\s+ONIOM:\s+gridpoint\s+2\s+method:\s+med\s+system:\s+model\s+energy:\s+(.+)/ ) {
		$log->{MEDMODEL} = $1;
		next;
        }
        
        # ONIOM gridpoint  3 method:  low   system:  mid energy  for three layer ONIOM
        if ( /^\s+ONIOM:\s+gridpoint\s+3\s+method:\s+low\s+system:\s+mid\s+energy:\s+(.+)/ ) {
		$log->{LOWMID} = $1;
		next;
        }

        # ONIOM gridpoint  4 method:  high  system:  model energy  for three layer ONIOM
        if ( /^\s+ONIOM:\s+gridpoint\s+4\s+method:\s+high\s+system:\s+model\s+energy:\s+(.+)/ ) {
		$log->{HIGHMODEL} = $1;
		next;
        }

        # ONIOM gridpoint  5 method:  med  system:  mid energy  for three layer ONIOM
        if ( /^\s+ONIOM:\s+gridpoint\s+5\s+method:\s+med\s+system:\s+mid\s+energy:\s+(.+)/ ) {
		$log->{MEDMID} = $1;
		next;
        }

        # ONIOM gridpoint  6 method:  low  system:  real energy  for three layer ONIOM
        if ( /^\s+ONIOM:\s+gridpoint\s+6\s+method:\s+low\s+system:\s+real\s+energy:\s+(.+)/ ) {
		$log->{LOWREAL} = $1;
                $log->{ISTHREELAYER} = 1;
		next;
        }
        
        # ONIOM  energy
        if ( /^\s+ONIOM:\s+extrapolated\s+energy\s+=(.+)/ ) {
		$log->{ONIOMENERGY} = $1;
		next;
        }

        # Dipole moment 
        if ( /^\s+ONIOM:\s+Dipole\s+moment\s+\(Debye\):/ ) {
		$_ = <LOGFILE>;
#                print $_;
                /^\s+X=\s+(.+)\s+Y=\s+(.+)\s+Z=\s+(.+)\s+Tot=\s+(.+)/;
                $log->{DIPOLEX} = $1;
                $log->{DIPOLEY} = $2;
                $log->{DIPOLEZ} = $3;
                $log->{DIPOLE}  = $4;
#                print "Current Dipole moments are ", $log->{DIPOLEX},$log->{DIPOLEY},$log->{DIPOLEZ},$log->{DIPOLE},"\n";
        	next;
        
        }

        # Check if is normal termination 
        if ( /^\s+Normal\s+termination/ ) {
                $log->{ISNORMAL}  = 1;
        	next;
        
        }
        
        
}

if ($log->{ISTHREELAYER} == 0) {
  print "Final ONIOM energy Report for two layer ONIOM calculation(Hartree):\n\n" if $debug >= $infocontrol ;
  print "High Level Model:\t",$log->{HIGHMODEL},"\n" if $debug >=$infocontrol ;
  print "Low  Level Model:\t",$log->{LOWMODEL},"\n" if $debug >=$infocontrol ;
  print "Low  Level Real :\t",$log->{LOWREAL},"\n" if $debug >=$infocontrol ;
  print "Low  Level Real-Model :\t",$log->{LOWREAL_MODEL},"\n\n" if $debug >=$infocontrol ;
  print "ONIOM Energy :\t",$log->{ONIOMENERGY},"\n\n" if $debug >=$infocontrol ;
  print "Dipole moment on X :\t",$log->{DIPOLEX},"\n" if $debug >=$infocontrol ;
  print "Dipole moment on Y :\t",$log->{DIPOLEY},"\n" if $debug >=$infocontrol ;
  print "Dipole moment on Z :\t",$log->{DIPOLEZ},"\n" if $debug >=$infocontrol ;
  print "Dipole moment (Total):\t",$log->{DIPOLE},"\n\n\n" if $debug >=$infocontrol;
 } else {
  print "Final ONIOM energy Report for three layer ONIOM calculation(Hartree):\n\n" if $debug >=$infocontrol;
  print "Low  Level Model:\t",$log->{LOWMODEL},"\n" if $debug >=$infocontrol;
  print "Med  Level Model:\t",$log->{MEDMODEL},"\n" if $debug >=$infocontrol;
  print "Low  Level Mid:\t\t",$log->{LOWMID},"\n" if $debug >=$infocontrol;
  print "High Level Model:\t",$log->{HIGHMODEL},"\n" if $debug >=$infocontrol;
  print "Med  Level Mid:\t\t",$log->{MEDMID},"\n" if $debug >=$infocontrol;
  print "Low  Level Real :\t",$log->{LOWREAL},"\n" if $debug >=$infocontrol;
  print "ONIOM Energy :\t",$log->{ONIOMENERGY},"\n\n" if $debug >=$infocontrol;
  print "Dipole moment on X :\t",$log->{DIPOLEX},"\n" if $debug >=$infocontrol;
  print "Dipole moment on Y :\t",$log->{DIPOLEY},"\n" if $debug >=$infocontrol;
  print "Dipole moment on Z :\t",$log->{DIPOLEZ},"\n" if $debug >=$infocontrol;
  print "Dipole moment (Total):\t",$log->{DIPOLE},"\n\n\n" if $debug >=$infocontrol;
 }


if ( $log->{ISNORMAL}  == 0 ) {
 print "Attention: this Gaussian calculation $log->{FILENAME} did not terminate normally!\n" if $debug >= 0 ;
}


 close (LOGFILE);

} #digest

1;
__END__

=head1 VERSION

0.01

=head1 AUTHOR

Peng Tao, E<lt>tao.21@osu.edu<gt>

=head1 COPYRIGHT

Copyright (c) 2009~2010 by Peng Tao

=cut

