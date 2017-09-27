package ESPT::Oniom;

use base qw(ESPT::ESSfile);
use ESPT::Glib 0.01;
use strict;
use warnings;

our $VERSION = '0.04';

=head1 NAME

ESPT::Oniom - Oniom object to treat Gaussian log file.

=head1 SYNOPSIS

   use ESPT::Oniom;

   my $log = Oniom->new();

=head1 DESCRIPTION

This module provides methods to quickly access data and structures contianed in a Gaussian log file
object for ONIOM claculation.  Guassian log files can only be read currently.

=cut

### Version History ###
# 0.01	Export ONIOM energy information
# 0.02  Treat three layer ONIOM calculation.
# 0.03  Grab Input orientation when there is Z-Matrix orientation 
# 0.04  Display energy in hatree using flag -ha
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
	
	# molecular info
	$log->{BASISLABELS} = [];
	$log->{C} = [];			# coefficient matrix
	$log->{CHARGE} = undef;
	$log->{ETHERM} = undef;		# Thermal corrections to E
	$log->{EELEC} = [];		# SCF electronic energy
	$log->{EIGEN} = [];
	$log->{EINFO} = "E(elec)";	# total energy description
	$log->{ENERGY} = undef; 	# total energy 
	$log->{EZPE} = undef;		# ZPE
	$log->{GSOLV} = undef;		# Delta G of Solvation
	$log->{GTHERM} = undef;		# Thermal corrections to G
	$log->{HOMO} = undef;
	$log->{HTHERM} = undef;		# Thermal corrections to H
	$log->{MOSYMM} = [];
	$log->{MULTIPLICITY} = undef;
	$log->{SADDLEPOINT} = undef;	# Saddle-point Order
	$log->{OCC} =[];		# MO occupation info 
	$log->{PG} = [];
	$log->{S} = [];			# overlap matrix
	$log->{SSQUARED} = [];		# S squared values

	# ONIOM info
	$log->{LOWMODEL} = undef;
	$log->{HIGHMODEL} = undef;
	$log->{LOWREAL} = undef;
        $log->{MEDMODLE} = undef;
        $log->{LOWMID} = undef;
        $log->{MEDMID} = undef;
	$log->{ONIOMENERGY} = undef;
	$log->{LOWREAL_MODEL} = undef;
	$log->{DIPOLEX} = undef;
	$log->{DIPOLEY} = undef;
	$log->{DIPOLEZ} = undef;
	$log->{DIPOLE} = undef;
        $log->{ISTHREELAYER} = undef;

        # Operation info
        $log->{GJFTEMPLATE} = undef;
	$log->{OPT} = 0;
        $log->{InHatree} = 0;
	$log->{STRUCTURE} = 0;
        $log->{STRUCTUREFILE} = undef;
        $log->{MOVEFLAG} = 0;
        $log->{VERBOSE} = 0;
        $log->{TEMPLATEINPUT} = 0;
        $log->{ISNORMAL} = 0;
        $log->{LAYERINFO} = 0;
        $log->{RAW} = 0;
        
	bless($log, $class);
	return $log;
}


## methods ##

=item B<$file-E<gt>analyze(filename [spin])>


=cut

# set filename & spin then digest the file
sub analyze : method {
	my $log = shift;
	$log->prepare(@_);
	$log->digest();
	return;
}


## subroutines ##



sub digest {

# For items with multiple occurances, the last value is reported 

my $log = shift;

# flags & counters
my $layerinfo = $log->{LAYERINFO};
my $debug     =  $log->{DEBUG};
my $raw       =  $log->{RAW};
my $Hlayer = 0;
my $Mlayer = 0;
my $Llayer = 0;

my $LowReal_Model = 0;
my $HighModel = 0;
my $ONIOMEnergy = 0;
my $diffLowReal_Modelkcal = 0;
my $diffLowReal_Modelhatree = 0;
my $diffHighModelkcal = 0;
my $diffHighModelhatree = 0;
my $diffONIOMEnergykcal = 0;
my $diffONIOMEnergyhatree = 0;
my $HFtoKcalmol = 627.509;
my $optsteps = 0;
my @outatomlist;
my @outatomsymbollist;
my $atomindex = 0;
my $structindicator = 0;
my $outatomnumber = 0;
my @currentXYZ;
my $currentline;
my $i=0;
my $j=0;
my $NAtoms=0;
my $structureindex = 0;

my ($counter, $line1, $line2, $templine,$isendofatomlist);


print "\nGaussian out file is ",$log->{FILENAME},"\n" if $debug >= 0;

#print "\nDebug level is $debug\n" if $debug >= 0;

if ($debug >= 0) {

print "\nGaussian oniom input file ",$log->{GJFTEMPLATE}," is used as template.\n" 
      if $log->{TEMPLATEINPUT} ==1;

print "\nOptimization -o was used\n" if $log->{OPT} ==1;

print "\nOutput xyz structure file is ", $log->{STRUCTUREFILE},".\n" if $log->{STRUCTURE} == 1;

print "\nMoveflag -g is turned on\n" if $log->{MOVEFLAG} ==1;

} #if ($log->{VERBOSE} == 1) 


# test

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
  print "Final ONIOM energy report for a two-layer ONIOM calculation (Hartree):\n";
  print "  This corresponds to the last complete step of an unfinished geometry optimization job.\n\n";
  if ($raw == 1) {
     print "High Level Model:\t",$log->{HIGHMODEL},"\n";
     print "Low  Level Model:\t",$log->{LOWMODEL},"\n";
     print "Low  Level Real :\t",$log->{LOWREAL},"\n";
     print "Low  Level Real-Model :\t",$log->{LOWREAL_MODEL},"\n\n";
     print "ONIOM Energy :\t",$log->{ONIOMENERGY},"\n\n";
     print "\nDipole moment (Debye) (X, Y, Z) is \n(",$log->{DIPOLEX},", ", $log->{DIPOLEY},", ", $log->{DIPOLEZ},").\n";
     print "Total Dipole moment (Debye) is $log->{DIPOLE}\n\n";
  } else {
     printf "High Level Model:\t%16.6f\n",$log->{HIGHMODEL};
     printf "Low  Level Model:\t%16.6f\n",$log->{LOWMODEL};
     printf "Low  Level Real :\t%16.6f\n",$log->{LOWREAL};
     printf "Low  Level Real-Model :\t%16.6f\n\n",$log->{LOWREAL_MODEL};
     printf "ONIOM Energy :\t\t%16.6f\n\n",$log->{ONIOMENERGY};
     printf "\nDipole moment (Debye) (X, Y, Z) is \n(%15.4f, %15.4f, %15.4f).\n",
            $log->{DIPOLEX},$log->{DIPOLEY},$log->{DIPOLEZ};
     printf "Total Dipole moment (Debye) is %15.4f.\n\n",$log->{DIPOLE};
  }
 } else {
  print "Final ONIOM energy report for a three-layer ONIOM calculation(Hartree):\n";
  print "  This corresonds to the last complete step of an unfinished geometry optimization job.\n\n";
  if ($raw == 1) {
    print "Low  Level Model:\t",$log->{LOWMODEL},"\n";
    print "Med  Level Model:\t",$log->{MEDMODEL},"\n";
    print "Low  Level Mid:\t\t",$log->{LOWMID},"\n";
    print "High Level Model:\t",$log->{HIGHMODEL},"\n";
    print "Med  Level Mid:\t\t",$log->{MEDMID},"\n";
    print "Low  Level Real :\t",$log->{LOWREAL},"\n";
    print "ONIOM Energy :\t",$log->{ONIOMENERGY},"\n\n";
    print "\nDipole moment (Debye) (X, Y, Z) is \n(",$log->{DIPOLEX},", ", $log->{DIPOLEY},", ", $log->{DIPOLEZ},").\n";
    print "Total Dipole moment (Debye) is $log->{DIPOLE}\n\n";
  } else {
    printf "Low  Level Model:\t%16.6f\n",$log->{LOWMODEL};
    printf "Med  Level Model:\t%16.6f\n",$log->{MEDMODEL};
    printf "Low  Level Mid:\t\t%16.6f\n",$log->{LOWMID};
    printf "High Level Model:\t%16.6f\n",$log->{HIGHMODEL};
    printf "Med  Level Mid:\t\t%16.6f\n",$log->{MEDMID};
    printf "Low  Level Real :\t%16.6f\n",$log->{LOWREAL};
    printf "ONIOM Energy :\t\t%16.6f\n\n",$log->{ONIOMENERGY};
    printf "\nDipole moment (Debye) (X, Y, Z) is \n(%15.4f, %15.4f, %15.4f).\n",
           $log->{DIPOLEX},$log->{DIPOLEY},$log->{DIPOLEZ};
    printf "Total Dipole moment (Debye) is %15.4f.\n\n",$log->{DIPOLE};
  }
 }




if ( $log->{ISNORMAL}  == 0 ) {
 print "Attention: this Gaussian calculation did not terminate normally!\n\n"  if $debug >= 0;
}

$LowReal_Model = $log->{LOWREAL_MODEL} ;
$HighModel     = $log->{HIGHMODEL}     ;
$ONIOMEnergy   = $log->{ONIOMENERGY}   ;

if ( $log->{OPT} == 1 ) {
 $optsteps = 0;
 print "OPT flag -o was turned on. Energy of each step will be printed.\n"  if $debug >= 0;
 
 if ( $log->{InHatree} == 0 ) {
    print "Energies along the optimization path (kcal/mol):\n";
   } else{
    print "Energies along the optimization path (hatree):\n";
   } 
 
 if ($log->{ISTHREELAYER} == 0) {
    print "Step Number\t     ONIOM \t\tHigh Model\t    Low (Real-Model)\n";
   } else {
    print "Step Number\t   ONIOM \t\tHigh Model\n";
   
   }

  
 close (LOGFILE);
 open(LOGFILE,$log->{FILENAME}) || die "Could not read $log->{FILENAME}\n$!\n";
 
 while (<LOGFILE>){
	# skip blank lines
	next if /^$/;

        # ONIOM gridpoint  1 method:  low   system:  model energy
        if ( /^\s+ONIOM:\s+gridpoint\s+1\s+method:\s+low\s+system:\s+model\s+energy:\s+(.+)/ ) {
		$log->{LOWMODEL} = $1;
		next;
        }
        
        # ONIOM gridpoint  2 method:  high   system:  model energy
        if ( /^\s+ONIOM:\s+gridpoint\s+2\s+method:\s+high\s+system:\s+model\s+energy:\s+(.+)/ ) {
		$log->{HIGHMODEL} = $1;
                $diffHighModelkcal = ($log->{HIGHMODEL} - $HighModel) * $HFtoKcalmol;
                $diffHighModelhatree = $log->{HIGHMODEL} - $HighModel;
		next;
        }

        # ONIOM gridpoint  3 method:  low   system:  real energy
        if ( /^\s+ONIOM:\s+gridpoint\s+3\s+method:\s+low\s+system:\s+real\s+energy:\s+(.+)/ ) {
		$log->{LOWREAL} = $1;
                $log->{LOWREAL_MODEL} = $log->{LOWREAL} - $log->{LOWMODEL};
                $diffLowReal_Modelkcal = ($log->{LOWREAL_MODEL} - $LowReal_Model) * $HFtoKcalmol;
                $diffLowReal_Modelhatree = $log->{LOWREAL_MODEL} - $LowReal_Model;
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
                $diffHighModelkcal = ($log->{HIGHMODEL} - $HighModel) * $HFtoKcalmol;
                $diffHighModelhatree = $log->{HIGHMODEL} - $HighModel;
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
                $diffONIOMEnergykcal = ($log->{ONIOMENERGY} - $ONIOMEnergy) * $HFtoKcalmol;
                $diffONIOMEnergyhatree = $log->{ONIOMENERGY} - $ONIOMEnergy;
                
                 $optsteps++;
                 if  ( $log->{ISTHREELAYER} == 0 ) {
                    if ( $log->{InHatree} == 0 ) {
                        if ($raw == 0){
                         printf " %3d     %20.6f %20.6f %20.6f\n",$optsteps,$diffONIOMEnergykcal,$diffHighModelkcal,$diffLowReal_Modelkcal;
                        } else {
                         print $optsteps,"                  ",$diffONIOMEnergykcal,"     ",$diffHighModelkcal,"     ",$diffLowReal_Modelkcal,"\n";
                        }
                       }  else {
                        if ($raw == 0){
                          printf " %3d     %20.6f %20.6f %20.6f\n",$optsteps,$log->{ONIOMENERGY},$log->{HIGHMODEL},$log->{LOWREAL_MODEL};
                         } else {
                          print $optsteps,"             ",$log->{ONIOMENERGY},"     ",$log->{HIGHMODEL},"      ",$log->{LOWREAL_MODEL},"\n";
                         }
                        }
                    } else {
                     if ( $log->{InHatree} == 0 ) {
                      if ($raw == 0) {
                         printf " %3d     %20.6f %20.6f\n",$optsteps,$diffONIOMEnergykcal,$diffHighModelkcal;
                       } else {
                         print $optsteps,"                  ",$diffONIOMEnergykcal,"     ",$diffHighModelkcal,"\n";
                       }
                      } else {
                        if ($raw == 0) {
                         printf " %3d     %20.6f %20.6f\n",$optsteps,$log->{ONIOMENERGY},$log->{HIGHMODEL};
                        } else {
                         print $optsteps,"               ",$log->{ONIOMENERGY},"   ",$log->{HIGHMODEL},"\n";
                        }
                       }
                    }

#                 print "$log->{ONIOMENERGY}\t$log->{HIGHMODEL}\t$log->{LOWREAL_MODEL}\n";
                 
                } 
		next;
        } # while LOGFILE
 print "There are $optsteps steps of optimization in this job.\n" if $debug >= 0;
 print "\n\n" if $debug >= 0;
 } # opt if



if ( $log->{STRUCTURE} == 1 ) {

 if ($log->{DEBUG} >= 1) {

  print "STRUCTURE flag -s was used.\n";
  print "Certain part of structure will be extracted and saved to file $log->{STRUCTUREFILE}.\n";
  } # if ($log->{VERBOSE} == 1) {
 
 close (LOGFILE);


# Findout which layers need to be extracted.

if ( ($layerinfo == 1)   || ($layerinfo == 12)  || ($layerinfo == 13)  || 
     ($layerinfo == 123) || ($layerinfo == 21)  || ($layerinfo == 31)  || 
     ($layerinfo == 132) || ($layerinfo == 213) || ($layerinfo == 312) || 
     ($layerinfo == 321) || ($layerinfo == 231) )
   {$Hlayer = 1;}

if ( ($layerinfo == 2)   || ($layerinfo == 12)  || ($layerinfo == 23)  || 
     ($layerinfo == 123) || ($layerinfo == 21)  || ($layerinfo == 32)  || 
     ($layerinfo == 132) || ($layerinfo == 213) || ($layerinfo == 312) || 
     ($layerinfo == 321) || ($layerinfo == 231) )
   {$Mlayer = 1;}

if ( ($layerinfo == 3)   || ($layerinfo == 13)  || ($layerinfo == 23)  || 
     ($layerinfo == 123) || ($layerinfo == 31)  || ($layerinfo == 32)  || 
     ($layerinfo == 132) || ($layerinfo == 213) || ($layerinfo == 312) || 
     ($layerinfo == 321) || ($layerinfo == 231) )
   {$Llayer = 1;}


# If user need level 2 in a two-layer ONIOM job, it should be L layer, not M layer
if  ( $log->{ISTHREELAYER} == 0 ) {
   if ($Mlayer == 1) {
   $Llayer = 1; 
   }
 }

# Construct list of atoms that need to be extracted and output
 if ($log->{TEMPLATEINPUT} ==1) {
   open(TEMPLATEFILE,$log->{GJFTEMPLATE}) || die "Could not read $log->{GJFTEMPLATE}\n$!\n";
   print "Open $log->{GJFTEMPLATE} for process.\n" if $debug >= 0;
   while (<TEMPLATEFILE>) {

     if   ($structindicator == 0){
       if (/^\s*-*\d+\s+\d+\s+-*\d+\s+\d+/) {
#       if (/^\s*-*[0-9]\s+[0-9]/) {

         print $_ if $debug >= 0;
         $structindicator = 1;
         next;
       } # if (/^[0-9]\s[0-9]) 
       next;
     }  # if   ($structindicator ==0){
     else {
       if (/^\s+$/ or /^$/){ last;}
#       print $_;
       $atomindex++;
#       printf "%d10\n",$atomindex;
#       if ( /^\s+[a-zA-Z\-]+\d+\.\d+\s+-*[01]\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+H/ ) {
#    Atom selection 
       if ( ((/^\s*([a-zA-Z]+)-[a-zA-Z0-9\*]+--*\d*.\d*\s+-*\d\s+(-*\d*\.\d*)\s+(-*\d*\.\d*)\s+(-*\d*\.\d*)\s+([Hh])/)
             && ($Hlayer==1) )
            || ((/^\s*([a-zA-Z]+)-[a-zA-Z0-9\*]+--*\d*.\d*\s+-*\d\s+(-*\d*\.\d*)\s+(-*\d*\.\d*)\s+(-*\d*\.\d*)\s+([Mm])/)
             && ($Mlayer==1) )       
            || ((/^\s*([a-zA-Z]+)-[a-zA-Z0-9\*]+--*\d*.\d*\s+-*\d\s+(-*\d*\.\d*)\s+(-*\d*\.\d*)\s+(-*\d*\.\d*)\s+([Ll])/)
             && ($Llayer==1) ) 
            ||  ((/^\s*([a-zA-Z]+)-[a-zA-Z0-9\*]+--*\d*.\d*\s+0\s+(-*\d*\.\d*)\s+(-*\d*\.\d*)\s+(-*\d*\.\d*)/)
                  && ($log->{MOVEFLAG}==1))) {
#        print $_;

        if ($log->{DEBUG} >= 1) {  print "$1\t$2\t$3\t$4\t$atomindex\n"; }
        push (@outatomsymbollist,$1);
        push (@outatomlist, $atomindex);
        $outatomnumber++;
        }  # if
     
     } # else
   
   } # while TEMPLATEFILE

   print "There are $outatomnumber atoms selected for output structures.\n" if $debug >= 0;
   } else {

      # No template file is provided. Try to extract ONIOM job info from given log file.   
   
      open(GOUTFILE,$log->{FILENAME}) || die "Could not read $log->{FILENAME}\n$!\n";

      # set atom counter as 0
      $atomindex = 0; 
      
      # Process ONIOM Gaussian output file 
      
       $line1 = "";
       $line2 = "";
       $templine = "";
      
      
      while ($templine = <GOUTFILE>) {
        
         $line1 = $line2; 
         $line2 = $templine;
         
         # looking for the starting of atom list in given log file.
         if (( $line1 =~ /^\s*Charge\s+=\s+-*\d+\s+Multiplicity/) and (not($line2 =~ /^\s*Charge\s+=\s+-*\d+\s+Multiplicity/))) {
          
          $isendofatomlist  = 0; 
          #  $line1 =~ /^(\S\S\S)\s+/;
          print "Start reading atom list from Gaussian output file $log->{FILENAME} \nat line:\n$line2" if $debug >= 0;
                
          $atomindex = 1;
      
          
          while ($isendofatomlist == 0) {
          
            $line1 = $line2; 
            $line2 = <GOUTFILE>;
            
            $_ = $line1;
            
            if ( ((/^\s*([a-zA-Z]+)-[a-zA-Z0-9\*]+--*\d*.\d*\s+-*\d\s+(-*\d*\.\d*)\s+(-*\d*\.\d*)\s+(-*\d*\.\d*)\s+([Hh])/)
                  && ($Hlayer==1) )
                  || ((/^\s*([a-zA-Z]+)-[a-zA-Z0-9\*]+--*\d*.\d*\s+-*\d\s+(-*\d*\.\d*)\s+(-*\d*\.\d*)\s+(-*\d*\.\d*)\s+([Mm])/)
                   && ($Mlayer==1) )       
                  || ((/^\s*([a-zA-Z]+)-[a-zA-Z0-9\*]+--*\d*.\d*\s+-*\d\s+(-*\d*\.\d*)\s+(-*\d*\.\d*)\s+(-*\d*\.\d*)\s+([Ll])/)
                   && ($Llayer==1) ) 
                  ||  ((/^\s*([a-zA-Z]+)-[a-zA-Z0-9\*]+--*\d*.\d*\s+0\s+(-*\d*\.\d*)\s+(-*\d*\.\d*)\s+(-*\d*\.\d*)/)
                   && ($log->{MOVEFLAG}==1))) {
              print "$1\t$2\t$3\t$4\t$atomindex\n" if $debug >= 2; 
              push (@outatomsymbollist,$1);
              push (@outatomlist, $atomindex);
              $outatomnumber++;
              }  # if
             $atomindex++;
            
             if ( $line2 =~ /^\s*$/) {$isendofatomlist = 1;}
      
           } # while ($isendofatomlist != 0) 
    
          
         $log->{NATOMS} = $atomindex-1;
         print "There are $log->{NATOMS} atoms in gaussian output file $log->{FILENAME}.\n" if $debug >= 1;
         last;
         
         } # if (( $line1 =~ 
       
       } # while ($templine = <GOUTFILE>)
      close (GOUTFILE);
      print "There are $outatomnumber atoms selected for output structures.\n" if $debug >= 0;
         
 }
       
 open(LOGFILE,$log->{FILENAME}) || die "Could not read $log->{FILENAME}\n$!\n";
 open(XYZFILE,">$log->{STRUCTUREFILE}") || die "Could not read $log->{STRUCTUREFILE}\n$!\n";
 
 
 $structindicator = 0;
 
 while (<LOGFILE>){
 
#     if ( $structindicator==0) {
      if ((/^\s*Z-Matrix\s+orientation:/) || (/^\s*Standard\s+orientation:/) || (/^\s*Input\s+orientation:/) ) {
         @currentXYZ = ();
#         print $_;
         $_=<LOGFILE>;
         if ( not(/^\s+-+/)) {
            print "Line($_)\n" if $debug >= 0;
            print "Gaussian log format may change. Pay attention to output structure(s).\n" if $debug >= 0;
            }
         $_ = <LOGFILE>;
         if ( not(/^\s+Center/) ) {
            print "Line($_)\n" if $debug >= 0;
            print "Gaussian log format may change. Pay attention to output structure(s).\n" if $debug >= 0;
            }

         $_ = <LOGFILE>;
         if ( not(/^\s+Number\s+Number\s+Type\s+X\s+Y\s+Z/) ) {
            print "Line($_)\n" if $debug >= 0;
            print "Gaussian log format may change. Pay attention to output structure(s).\n" if $debug >= 0;
            }

         $_ = <LOGFILE>;
         $j = 0;
         for ($i = 1; $i<=$NAtoms;$i++){
     
            $_ = <LOGFILE>;

#            if (/^\s+.+\s+.+\s+.+\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)/){
#               print $_;
#               print "is $i and j is $j\n";
#               print "Current atom for pick up is ",$outatomlist[$j],"\n";
                
               if ($i == $outatomlist[$j]) {
                /^\s+.+\s+.+\s+.+\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)/;
                push (@currentXYZ,$1);
                push (@currentXYZ,$2);
                push (@currentXYZ,$3);
#                print "$i\t\t$j\n";
#                print $currentXYZ[3*($j)],"\t",$currentXYZ[3*($j)+1],"\t",$currentXYZ[3*($j)+2],"\n";
#                print "$1\t$2\t$3\n";
                $j++;
                if ($j==$outatomnumber) {last;}
               } #if ($i == @outatomlist[$j]) {

#              } # if (/^\s+.+\s+.+\s+.+\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)/){

         
          } # for
          
# Write out every structure from optimization  if -o is turned on

         if ($log->{OPT}==1) {

           $structureindex++;
           print XYZFILE "$outatomnumber\n";
           print XYZFILE "ONIOM structure from Gaussian File $log->{FILENAME} structure Number $structureindex.\n";
           for ($j = 0; $j<$outatomnumber;$j++) {
           printf XYZFILE  "%10s %20.6f%20.6f%20.6f\n",$outatomsymbollist[$j],$currentXYZ[3*$j],$currentXYZ[3*$j+1],$currentXYZ[3*$j+2];
      
        } # for ($j = 1; $j<=$outatomnumber;$j++) {
         
         
         
         } # if ($log->{OPT}==1) {
           
         
 
        
        } # if ((/^\s+Z-Matrix\s+orientation:/) || (/^Standard\s+orientation:/)|| (/^\s*Input\s+orientation:/) ) 
        
      
#      }  # if ( $structindicator==0) 
       
       

    } # while LOGFILE

# Write the last structure from current file,  when -o is NOT turned on
# the format is in XYZ. It can be transformed to other format by BABEL.
    
    if ($log->{OPT}!=1) { 
      print XYZFILE "$outatomnumber\n";
      print XYZFILE "ONIOM structure from Gaussian File $log->{FILENAME}\n";
      for ($j = 0; $j<$outatomnumber;$j++) {
         printf XYZFILE  "%10s %20.6f%20.6f%20.6f\n",$outatomsymbollist[$j],$currentXYZ[3*$j],$currentXYZ[3*$j+1],$currentXYZ[3*$j+2];
      
        } # for ($j = 1; $j<=$outatomnumber;$j++) {
     
     }  # if ($log->{OPT}!=1) { 
    
 } # if ( $log->{STRUCTURE} == 1 ) {

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

