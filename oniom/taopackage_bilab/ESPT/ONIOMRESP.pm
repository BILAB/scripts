package ESPT::ONIOMRESP;

use base qw(ESPT::ESSfile);
use ESPT::Glib 0.01;
use strict;
use warnings;

our $VERSION = '0.1';

=head1 NAME

ESPT::ONIOMRESP - Oniom RESP object to facilitate RESP charge fitting.

=head1 SYNOPSIS

   use ESPT::ONIOMRESP;

   my $oniomresp = ONIOMRESP->new();

=head1 DESCRIPTION

=over 15

This module provides methods to facilitate RESP charge fitting for ONIOM calculation

=cut

### Version History ###
# 0.1	Mode two and three work now (refer to oniomresp.pl help)


### To Do ###
# Mode one

## the object constructor **

sub new {
	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my $oniomresp = ESPT::ESSfile->new();

	$oniomresp->{TYPE} = "oniom";

	# program info
	$oniomresp->{PROGRAM}  = "Gaussian";
	$oniomresp->{VERSION}  = undef;
	$oniomresp->{REVISION} = undef;
	$oniomresp->{COMPILE}  = undef;
	$oniomresp->{COMPILEDATE} = undef;
        
        # Operation info
        $oniomresp->{G_FILE_NAME} = undef;
        $oniomresp->{OUTPUT_FILE} = undef;
        $oniomresp->{MODE} = undef;
        $oniomresp->{CapH_Num} = undef;
        $oniomresp->{Map_File} = undef;
        $oniomresp->{Have_Map_File} = undef;
        $oniomresp->{Charge_Comp_File} = undef;
        $oniomresp->{Qin_File} = undef;
        $oniomresp->{Step_Num} = 0;
        $oniomresp->{IS_LOG_FILE} = undef;
	$oniomresp->{IS_JOB_FILE} = 0;
        $oniomresp->{NAtoms} = 0;

        # Genereal ONIOM atoms info
	$oniomresp->{ATOMSGfile} = [];
	$oniomresp->{ATOMTYPEGfile} = [];
	$oniomresp->{PARTIALCHARGGfile} = [];
        $oniomresp->{ATOMFREEZEFLAGGfile} = [];
	$oniomresp->{CARTCOORDGfile} = [];
	$oniomresp->{AtomLayerGfile} = [];
	$oniomresp->{MovedGfile} = [];
	$oniomresp->{LineNumGfile} = [];
        $oniomresp->{Original_Line} = [];
        $oniomresp->{IsLinkAtom} = [];
        $oniomresp->{LinkAtomEle} = [];
        $oniomresp->{LinkAtomType} = [];
        $oniomresp->{LinkedAtomNum} = [];


        # Variables used in Mode One
        $oniomresp->{Set_Mid} = undef;
        $oniomresp->{ISTWOLAYER} = undef;
        $oniomresp->{ISTHREELAYER} = undef;
        $oniomresp->{Layers} = undef;
        $oniomresp->{HaveListFile} = undef;
        $oniomresp->{AtomListFile} = undef;
        $oniomresp->{DefaultScale} = 0.723886;
        $oniomresp->{DefaultLinkAtom} = "H";

        # Variables used in Mode Two
        $oniomresp->{Element_NUM} = [];
        $oniomresp->{CARTCOORD} = [];
        $oniomresp->{Charge} = 0;
        $oniomresp->{Multiplicity} = 1;

        # Variables used in Mode Three with ONIOM input file (Gfile)
        $oniomresp->{RESP_charges} = [];
        $oniomresp->{ATOMNUMGfile} = undef;
        $oniomresp->{NeedNewRESPGfile} = [];
        $oniomresp->{RESPChargeIndexGfile} = [];


	# System info
	$oniomresp->{ATOMSModelGjf} = [];    # Element symbol from input file
        $oniomresp->{RESIDUEModelGjf} = [];
        $oniomresp->{RESIDUENUMModelGjf} = [];
        $oniomresp->{ATOMTYPEModelGjf} = [];
        $oniomresp->{ELEMENTNUMModelGjf} = [];
        $oniomresp->{ELEMENTModelGjf} = [];
        $oniomresp->{PARTIALCHARGModelGjf} = [];
        $oniomresp->{CARTCOORDModelGjf} = [];
        $oniomresp->{ATOMDEFINEModelGjf} = [];
        $oniomresp->{ATOMFREEZEFLAGModelGjf} = [];
        $oniomresp->{AtomLayerModelGjf} = [];
        $oniomresp->{ATOMINRESIDUEModelGjf} = [];
        $oniomresp->{ORIGINALLINEModelGjf} = [];
        $oniomresp->{LineNumModelGjf} = [];
        $oniomresp->{ATOMNUMModelGjf} = 0;
        $oniomresp->{MovedModelGjf} = [];  # flag to label if this atom moved from its original place
        $oniomresp->{ProdIndexModelGjf} = []; # This atom's corresponding atom index in production system

        
        $oniomresp->{ProdNewInput} = undef;
        $oniomresp->{DefaultProdNewInput} = "Prod_ONIOM_NEW.gjf";

	$oniomresp->{GJFINPUT} = undef;
        $oniomresp->{GJFOUTPUT} = undef;
        $oniomresp->{STRUCTURENUMBER} = 0;
        
        
	bless($oniomresp, $class);
	return $oniomresp;
}

# set filename & spin then digest the file
sub analyze : method {
	my $oniomresp = shift;
	$oniomresp->digest();
	return;
}

## subroutines ##

sub digest {

my $oniomresp = shift;
my $debug     = $oniomresp->{DEBUG};
my $gfile     = $oniomresp->{G_FILE_NAME};
my $outputfile = $oniomresp->{OUTPUT_FILE};
my $mode = $oniomresp->{MODE};
my $capHnum = $oniomresp->{CapH_Num} ;
my $stepnum = $oniomresp->{Step_Num} ;
my $mapfile = $oniomresp->{Map_File} ;
my $havemapfile = $oniomresp->{Have_Map_File} ;
my $chgcompfile = $oniomresp->{Charge_Comp_File};
my $qinfile = $oniomresp->{Qin_File} ;
my $isLogfile = $oniomresp->{IS_LOG_FILE};
my $isJobfile = $oniomresp->{IS_JOB_FILE};
my $atomindex;
my $structindicator;
my $lineindex;
my $i;
my $j;
my $distance;
my $mismatch;

my ($prodgjf,$prodonb,$modelonb,$modelgjf,$tol);


# Running for mode 1
if ($mode == 1) {
   $oniomresp->modeone;
 }

# Running for mode 2
if ($mode == 2) {
   $oniomresp->modetwo;
 }


# Running for mode 3
if ($mode == 3) {
   $oniomresp->modethree;
 }


}  # sub digest {


# Subroutine for running mode one
# Generate gaussian input file for ESP calculation from given gaussian Log file.

sub modeone {
   my $oniomresp    = shift;
   my $debug        = $oniomresp->{DEBUG};
   my $gfile        = $oniomresp->{G_FILE_NAME};
   my $outputfile   = $oniomresp->{OUTPUT_FILE};
   my $mode         = $oniomresp->{MODE};
   my $capHnum      = $oniomresp->{CapH_Num} ;
   my $stepnum      = $oniomresp->{Step_Num} ;
   my $mapfile      = $oniomresp->{Map_File} ;
   my $havemapfile  = $oniomresp->{Have_Map_File} ;
   my $chgcompfile  = $oniomresp->{Charge_Comp_File};
   my $qinfile      = $oniomresp->{Qin_File} ;
   my $isLogfile    = $oniomresp->{IS_LOG_FILE};
   my $isJobfile    = $oniomresp->{IS_JOB_FILE};
   my $setmid       = $oniomresp->{Set_Mid};
   my $HaveListFile = $oniomresp->{HaveListFile};
   my $AtomListFile = $oniomresp->{AtomListFile};
   my $DefaultScale =$oniomresp->{DefaultScale};
   my  $DefaultLinkAtom = $oniomresp->{DefaultLinkAtom};
   my $NAtoms = 0;
   my $totalSteps = 0;
   my $isnormal = 0;
   my $undefined = undef;
   my ($i,$j,$atomindex);
   my ($is_twolayer,$is_threelayer,$linkmatch,$skipindicator);
   
   my ($logatom,$currentX,$currentY,$currentZ,$currentline, $oribondlen, $linkbondlen);
   
   my ($currLinkHost,$currHost);
   
   my ($HighMultiplicity,$HighCharge,$MidMultiplicity,$MidCharge);

   my (@LinkHostEle,@LinkAtomEle,@LinkIndex,@HostEle,@HostIndex,@Scale1,@Scale2,@Scale3);
   # 
   # About @LinkHostEle,@LinkAtomEle,@LinkIndex,@HostEle,@HostIndex,@Scale
   #
   # ONIOM: Cut between C /H  4300 and C  4303 factor= 0.723886 0.723886
   #
   # ONIOM: Cut between C            /H          4300 and   C        4303 factor= 0.723886 0.723886
   #
   #                    C              H             4300       C        4303         0.723886 0.723886
   #                    @LinkHostEle  @LinkAtomEle   @LinkIndex @HostEle @HostIndex   @Scale1  @Scale2   @Scale3
   #
   #
   
   # Same arrays again, these info obtained from Input section in log file. Obtained for sanity check
   my (@LinkHostEleIn,@LinkAtomEleIn,@LinkIndexIn,@HostEleIn,@HostIndexIn,@LinkHostLayer);
   
   my (@ListAtom,@ListHost,@ListScale,@ListLinkElement,@HasLink,$List_Index,$skipfrom, $listLinkAtomNum);
   
   my ($tempAtom,$tempHost,$tempScale,$tempLinkElement,$tempHasLink);
   
   
   my ($link_index,$link_indexIn,$line_index,$stepcounter, $chosenLAYER, $chosenlayer);
   
   # Variables for output geometry
   my (@modelatom,@modelX,@modelY,@modelZ,$modelindex);
   
   @modelatom = ();
   @modelX = ();
   @modelY = ();
   @modelZ = ();
      
   $oniomresp->{ISTWOLAYER} = 0;
   $oniomresp->{ISTHREELAYER} = 0;
   $oniomresp->{Layers} = 1;
   
# open Gaussian Log file for reading or display error
  
  print "Open Gaussian Log file $gfile ......\n\n"  if $debug >= 1;
  open(LOGFILE,$gfile) || die "Could not read $gfile\n$!\n";

# Process gaussian log file for first time to obtain necessary information
  
  while (<LOGFILE>){
    
 	# skip blank lines
	next if /^$/;

        if (/^\s*NAtoms=\s*(\d+)/i) {
           $NAtoms=$1;
        
        # print "There are $NAtoms atoms.\n";
         }

        # Check if is normal termination 
        if ( /^\s*Normal\s+termination/i ) {
                $isnormal  = 1;
        	next;
         }

        # Find out charge and multiplicity
        # Try to match line like:
        #  Charge =  1 Multiplicity = 1 for high  level calculation on model system.
        
        if ( /^\s*Charge\s*=\s*(-*\d+)\s+Multiplicity\s*=\s*(\d+)\s+for\s+high\s+level\s+calculation\s+on\s+model\s+system/i ) {
           print "Charge is $1, multiplicity is $2 for high level calculation on model system.\n\n" if $debug >= 1;
           $HighCharge = $1;
           $HighMultiplicity = $2;
          }

        # Try to match line like:
        #  Charge =  0 Multiplicity = 1 for med   level calculation on mid   system.

        if ( /^\s*Charge\s*=\s*(-*\d+)\s+Multiplicity\s*=\s*(\d+)\s+for\s+med\s+level\s+calculation\s+on\s+mid\s+system/i ) {
           print "Charge is $1, multiplicity is $2 for med   level calculation on mid   system.\n\n" if $debug >= 1;
           $MidCharge = $1;
           $MidMultiplicity = $2;
          }


        # Check if this is a two-layer ONIOM job
        if ( /^\s*ONIOM:\s+gridpoint\s+1\s+method:\s+low\s+system:\s+model\s+energy:\s+(.+)/i ) {
		$oniomresp->{ISTWOLAYER} = 1;
                $oniomresp->{Layers} = 2;
		next;
        }


        # Check if this is a three-layer ONIOM job
        if ( /^\s*ONIOM:\s+gridpoint\s+6\s+method:\s+low\s+system:\s+real\s+energy:\s+(.+)/i ) {
                $oniomresp->{ISTHREELAYER} = 1;
                $oniomresp->{ISTWOLAYER} = 0;
                $oniomresp->{Layers} = 3;
		next;
        }
        
        # Count total step number in given log file
        if ((/^\s*Z-Matrix\s+orientation:/i) || (/^\s*Standard\s+orientation:/i)) {
               $totalSteps++;
              }
        
        

          
   } # while (<LOGFILE>){

 close (LOGFILE);
 
 # If given log file is not a ONIOM job
 if ( ($oniomresp->{ISTWOLAYER} + $oniomresp->{ISTHREELAYER}) == 0 ) {
   print "\n Attention: given file $gfile may not be an ONIOM job.\n" if $debug >= 0;
   print "Still try to proceed, but please double check log file and output files.\n" if $debug >= 0;
   
  }
  
  $is_twolayer = $oniomresp->{ISTWOLAYER};
  $is_threelayer = $oniomresp->{ISTHREELAYER}; 
  
# If not terminated normally
   if ( $isnormal == 0 )  {
     print "Warning: It occurs that $gfile does not terminate normally.\n\n" if $debug >= 0;
     print "Please find out why before you proceed on RESP charges fitting.\n\n" if $debug >= 0;
   
    }
  
  if ( $totalSteps > 1 ) {
    print "\n There are ",$totalSteps-1," steps optimization in $gfile.\n"  if $debug >= 0;
    
    }

# If NAtoms information is not available, we need to find it by
# counting in standard orientation.

   if ( $NAtoms == 0 ) {
      print "NAtoms info is not presented in $gfile. \n\nTrying to find number of atoms in another way...\n\n" if $debug >= 0;
      
      open(LOGFILE,$gfile) || die "Could not read $gfile\n$!\n";

      while (<LOGFILE>){
 
           if ((/^\s*Z-Matrix\s+orientation:/i) || (/^\s*Standard\s+orientation:/i)) {

               $_ = <LOGFILE>;
               if ( not(/^\s+-+/)) {
                  print "Line($_)\n" if $debug >= 3;
                  print "Gaussian log format may change. Pay attention to output structure(s).\n" if  $debug >= 0;
                  }
               $_ = <LOGFILE>;
               if ( not(/^\s*Center/i) ) {
                  print "Line($_)\n" if $debug >= 3;
                  print "Gaussian log format may change. Pay attention to output structure(s).\n"  if $debug >= 0;
                  }

               $_ = <LOGFILE>;
               if ( not(/^\s*Number\s+Number\s+Type\s+X\s+Y\s+Z/i) ) {
                  print "Line($_)\n" if $debug >= 3;
                  print "Gaussian log format may change. Pay attention to output structure(s).\n"  if $debug >= 0;
                  }

               $_ = <LOGFILE>;

               if ( not(/^\s+-+/)) {
                  print "Line($_)\n" if $debug >= 3;
                  print "Gaussian log format may change. Pay attention to output structure(s).\n" if  $debug >= 0;
                  }
               
               $_ = <LOGFILE>;
               
               while (not(/^\s+-+/)) {
                  $_ = <LOGFILE>;
                  $NAtoms++;
                } # while (not(/^\s+-+/) {
              }

          } # while LOGFILE

       } # if  ( $NAtoms == 0 ) {
   
   print "There are $NAtoms atoms in file $gfile.\n\n" if $debug >= 0;
   
   $oniomresp->{NAtoms} = $NAtoms;
   
   close (LOGFILE);
   
   if ( $stepnum == 0 ) 
    {
     $stepnum = $totalSteps ;
    }  elsif ( $stepnum > $totalSteps - 1 )  {
     print "\n Given step number $stepnum is greater than found step number ",$totalSteps-1,". Use found step number intead.\n"  if $debug >= 0;
     $stepnum = $totalSteps;
    } else {
     # Since in optimization job, number of "Z-Matrix orientation" lines is one larger than real optimization steps
     # The obtained step number needs to increase one to match the desired step number.
     $stepnum++;
    }
   
   if ( $stepnum > 1 ) {
     print "\n The actual step number ",$stepnum-1,", geometry after number $stepnum Z-Matrix orientation or Standard orientation is used.\n" if $debug >= 2;
    } else {
     print "\nGiven file $gfile has only one geometry, and may not be an optimization file.\n\n" if $debug >= 2;
    }
     
# Process gaussian log file again for distance scaling factor for link atoms
#  ONIOM set up (atom type, layer, move, etc) and coordinatesfor each atom.
 
  open(LOGFILE,$gfile) || die "Could not read $gfile\n$!\n";
  
  $link_index = 0;
  
  $line_index = 0;

  $link_indexIn = 0;  
  
  $stepcounter = 0;
  
  while (<LOGFILE>){  #M1NO1
    
    $line_index++;
    
   
    if ( $is_twolayer == 1 ) {
       #
       # Read line like:
       #   ONIOM: Cut between C /H  4300 and C  4303 factor= 0.723886 0.723886
       #
       if ( /^\s*ONIOM:\s+Cut between\s+([a-zA-Z]+)\s+\/([a-zA-Z]+)\s+(\d+)\s+and\s+([a-zA-Z]+)\s+(\d+)\s+factor=\s+(\d*.\d*)\s+(\d*.\d*)/i ) 
          {
          
          print $_  if $debug >= 4;
          
          $LinkHostEle[$link_index] = $1;
          $LinkAtomEle[$link_index] = $2;
          $LinkIndex[$link_index]   = $3;
          $HostEle[$link_index]     = $4;
          $HostIndex[$link_index]   = $5;
          $Scale1[$link_index]      = $6;
          $Scale2[$link_index]      = $7;

          $link_index++;
          
          if ( $link_index == 1 ) {
              print "Link atoms information obtained from \"ONIOM: Cut\" section.\n\n" if $debug >= 3;
             }


          print " $link_index: Link atom $1 change to $2 ($3); host $4 ($5). Scale Factor: $6 $7.\n" if $debug >= 3;
          
        }
      
      } else {
      
       if ( /^\s*ONIOM:\s+Cut between\s+([a-zA-Z]+)\s+\/([a-zA-Z]+)\s+(\d+)\s+and\s+([a-zA-Z]+)\s+(\d+)\s+factor=\s+(\d*.\d*)\s+(\d*.\d*)\s+(\d*.\d*)/i ) 
          {
          #
          # Read line like:
          #   ONIOM: Cut between C /H   241 and C   244 factor= 0.723886 0.723886 0.723886 
          #
          
          print $_  if $debug >= 4;
          
          $LinkHostEle[$link_index] = $1;
          $LinkAtomEle[$link_index] = $2;
          $LinkIndex[$link_index]   = $3;
          $HostEle[$link_index]     = $4;
          $HostIndex[$link_index]   = $5;
          $Scale1[$link_index]      = $6;
          $Scale2[$link_index]      = $7;
          $Scale3[$link_index]      = $8;
          $link_index++;

          if ( $link_index == 1 ) {
              print "Link atoms information obtained from \"ONIOM: Cut\" section.\n\n" if $debug >= 3;
             }

          print " $link_index: Link atom $1 change to $2 ($3); host $4 ($5). Scale Factor: $6 $7 $8.\n" if $debug >= 3;
          
         }
      }
   
   if ( /\s*Symbolic\s+Z-matrix:/i ) { #M1NO2
      
     # When encounter line
     #  Symbolic Z-matrix: 
     # it is time to find ONIOM input part

     
     if ($is_twolayer == 1) {
        
        # Skip three lines:
        #  Charge =-11 Multiplicity = 1 for low   level calculation on real  system.
        #  Charge =  1 Multiplicity = 1 for high  level calculation on model system.
        #  Charge =  1 Multiplicity = 1 for low   level calculation on model system.
        #
        $_ = <LOGFILE>;
        $line_index++;
        if (not/^\s*Charge\s*=\s*-*\d+\s+Multiplicity\s*=\s*\d+\s+for\s+low\s+level\s+calculation\s+on\s+real\s+system/i) {
            print $_ if $debug >= 0;
            print "Gaussian log format may change for two-layer ONIOM job. Pay attention to output file.\n" if $debug >= 0;
          }
        
        $line_index++;
        $_ = <LOGFILE>;
        if (not/^\s*Charge\s*=\s*-*\d+\s+Multiplicity\s*=\s*\d+\s+for\s+high\s+level\s+calculation\s+on\s+model\s+system/i) {
            print $_ if $debug >= 0;
            print "Gaussian log format may change for two-layer ONIOM job. Pay attention to output file.\n" if $debug >= 0;
          }
        
        $line_index++;
        $_ = <LOGFILE>;
        if (not/^\s*Charge\s*=\s*-*\d+\s+Multiplicity\s*=\s*\d+\s+for\s+low\s+level\s+calculation\s+on\s+model\s+system/i) {
            print $_ if $debug >= 0;
            print "Gaussian log format may change for two-layer ONIOM job. Pay attention to output file.\n" if $debug >= 0;
          }
        
     
      } else { #M1NO3

        # Skip six lines:
        # Charge = -3 Multiplicity = 1 for low   level calculation on real  system.
        # Charge =  0 Multiplicity = 1 for med   level calculation on mid   system.
        # Charge =  0 Multiplicity = 1 for low   level calculation on mid   system.
        # Charge = -3 Multiplicity = 1 for high  level calculation on model system.
        # Charge = -3 Multiplicity = 1 for med   level calculation on model system.
        # Charge = -3 Multiplicity = 1 for low   level calculation on model system.
        #
      
        $_ = <LOGFILE>;
        $line_index++;
        if (not/^\s*Charge\s*=\s*-*\d+\s+Multiplicity\s*=\s*\d+\s+for\s+low\s+level\s+calculation\s+on\s+real\s+system/i) {
            print $_ if $debug >= 0;
            print "Gaussian log format may change for three-layer ONIOM job. Pay attention to output file.\n" if $debug >= 0;
          }

        $_ = <LOGFILE>;
        $line_index++;
        if (not/^\s*Charge\s*=\s*-*\d+\s+Multiplicity\s*=\s*\d+\s+for\s+med\s+level\s+calculation\s+on\s+mid\s+system/i) {
            print $_ if $debug >= 0;
            print "Gaussian log format may change for three-layer ONIOM job. Pay attention to output file.\n" if $debug >= 0;
          }

        $_ = <LOGFILE>;
        $line_index++;
        if (not/^\s*Charge\s*=\s*-*\d+\s+Multiplicity\s*=\s*\d+\s+for\s+low\s+level\s+calculation\s+on\s+mid\s+system/i) {
            print $_ if $debug >= 0;
            print "Gaussian log format may change for three-layer ONIOM job. Pay attention to output file.\n" if $debug >= 0;
          }
        
        $_ = <LOGFILE>;
        $line_index++;
        if (not/^\s*Charge\s*=\s*-*\d+\s+Multiplicity\s*=\s*\d+\s+for\s+high\s+level\s+calculation\s+on\s+model\s+system/i) {
            print $_ if $debug >= 0;
            print "Gaussian log format may change for three-layer ONIOM job. Pay attention to output file.\n" if $debug >= 0;
          }

        $_ =  <LOGFILE>;
        $line_index++;
        if (not/^\s*Charge\s*=\s*-*\d+\s+Multiplicity\s*=\s*\d+\s+for\s+med\s+level\s+calculation\s+on\s+model\s+system/i) {
            print $_ if $debug >= 0;
            print "Gaussian log format may change for three-layer ONIOM job. Pay attention to output file.\n" if $debug >= 0;
          }
        
        $_ = <LOGFILE>;
        $line_index++;
        if (not/^\s*Charge\s*=\s*-*\d+\s+Multiplicity\s*=\s*\d+\s+for\s+low\s+level\s+calculation\s+on\s+model\s+system/i) {
            print $_ if $debug >= 0;
            print "Gaussian log format may change for three-layer ONIOM job. Pay attention to output file.\n" if $debug >= 0;
          }
      } #M1NO3
      
     
     # Now, it's time to read in ONIOM setup information.
     # Read in $NAtoms lines like:
     # H-H-0.216500         -1    40.56     22.841    52.249   L 
     
     print "Link atoms information obtained from Z-Matrix input section.\n\n" if $debug >= 3;
  
     for ($i = 0; $i < $NAtoms; $i++){ #M1NO4
     
            $_ = <LOGFILE>;
            $line_index++;
            $oniomresp->{Original_Line}[$i] = $_ ;
            print $_ if $debug >= 6;                        
            
             if (/\s*([a-zA-Z]+)-([a-zA-Z0-9\*]+)-(-*\d*.\d*)\s+(-*\d)\s+-*\d*\.\d*\s+-*\d*\.\d*\s+-*\d*\.\d*\s+([HMLhml])\s*$/) 
              { #M1NO5
              
               # If no link atom info in current atom
               print "$1\t$2\t$3\t$4\t$5\t$i\n" if $debug >= 5; 
               $oniomresp->{ATOMSGfile}[$i]          = $1;
               $oniomresp->{ATOMTYPEGfile}[$i]       = $2;
               $oniomresp->{PARTIALCHARGGfile}[$i]   = $3;
               $oniomresp->{ATOMFREEZEFLAGGfile}[$i] = $4;
               $oniomresp->{AtomLayerGfile}[$i]      = $5;
               $oniomresp->{LineNumGfile}[$i]        = $line_index;
               $oniomresp->{IsLinkAtom}[$i]          = 0;
               $oniomresp->{LinkAtomEle}[$i]         = undef;
               $oniomresp->{LinkAtomType}[$i]        = undef;
               $oniomresp->{LinkedAtomNum}[$i]       = undef;
              
              }  elsif  #M1NO5
               (/\s*([a-zA-Z]+)-([a-zA-Z0-9\*]+)-(-*\d*.\d*)\s+(-*\d)\s+-*\d*\.\d*\s+-*\d*\.\d*\s+-*\d*\.\d*\s+([HMLhml])\s*([a-zA-Z]+)-([a-zA-Z0-9\*]+)\s+(\d+)/)
                # With link atom info
                # e.g. C-CT--0.122638       0     50.42455   35.15176   58.6587  L    H-H4     4447  0.        0. 
               { #M1NO6
                print "$1\t$2\t$3\t$4\t$5\t$6\t$7\t$8\t$i\n" if $debug >= 5; 
                $oniomresp->{ATOMSGfile}[$i]          = $1;
                $oniomresp->{ATOMTYPEGfile}[$i]       = $2;
                $oniomresp->{PARTIALCHARGGfile}[$i]   = $3;
                $oniomresp->{ATOMFREEZEFLAGGfile}[$i] = $4;
                $oniomresp->{AtomLayerGfile}[$i]      = $5;
                $oniomresp->{LineNumGfile}[$i]        = $line_index;
                $oniomresp->{IsLinkAtom}[$i]          = 1;
                $oniomresp->{LinkAtomEle}[$i]         = $6;
                $oniomresp->{LinkAtomType}[$i]        = $7;
                $oniomresp->{LinkedAtomNum}[$i]       = $8;
                
                $LinkHostEleIn[$link_indexIn] = $1;
                $LinkAtomEleIn[$link_indexIn] = $6;
                $LinkIndexIn[$link_indexIn]   = $i+1;
                $HostIndexIn[$link_indexIn]   = $8;
                $LinkHostLayer[$link_indexIn] = $5;
                
                $link_indexIn++;

                print " From Z-Matrix input: $link_indexIn. link atom $1 change to $6"," (", $i+1,"), host atom $8.\n" if $debug >= 3;
                
                
               }  #M1NO6
              
              if ( $debug >= 3 ) {
                 print " Atom ",          $i+1," (line:",$line_index,") undefined.\n" unless defined $oniomresp->{ATOMSGfile}[$i]; 
                 print " Atom Type ",     $i+1," (line:",$line_index,") undefined.\n" unless defined $oniomresp->{ATOMTYPEGfile}[$i]; 
                 print " Partial Charge ",$i+1," (line:",$line_index,") undefined.\n" unless defined $oniomresp->{PARTIALCHARGGfile}[$i]; 
                 print " Move flag ",     $i+1," (line:",$line_index,") undefined.\n" unless defined $oniomresp->{ATOMFREEZEFLAGGfile}[$i]; 
                 print " Layer ",         $i+1," (line:",$line_index,") undefined.\n" unless defined $oniomresp->{AtomLayerGfile}[$i]; 
                }


         
          } #M1NO4
          
         print "\n" if $debug >= 3;
   
    }  #M1NO2
      
   # Locate the desired geometries   
   if ( (/^\s*Z-Matrix\s+orientation:/i) || (/^\s*Standard\s+orientation:/i) ) { #M1NO7
      $stepcounter++;
      
      if ( $stepcounter != $stepnum ) {
        # Not desired geometry, need to skip it.
        for ($i = 0; $i < $NAtoms+4; $i++) { $_ = <LOGFILE>; }
       } else { #M1NO8
        # Desired geometry, need to obtain coordinates
        # Skip next four lines:
        # ---------------------------------------------------------------------
        # Center     Atomic     Atomic              Coordinates (Angstroms)
        # Number     Number      Type              X           Y           Z
        # ---------------------------------------------------------------------

         $_ = <LOGFILE>;
         if ( not(/^\s+-+/)) {
            print $_ if $debug >= 1;
            print "Gaussian log format may change. Pay attention to output structure(s).\n" if $debug >= 1;
            }
            
         $_ = <LOGFILE>;
         if ( not(/^\s*Center/i) ) {
            print $_ if $debug >= 1;
            print "Gaussian log format may change. Pay attention to output structure(s).\n" if $debug >= 1;
            }

         $_ = <LOGFILE>;
         if ( not(/^\s*Number\s+Number\s+Type\s+X\s+Y\s+Z/i) ) {
            print $_ if $debug >= 1;
            print "Gaussian log format may change. Pay attention to output structure(s).\n" if $debug >= 1;
            }

         $_ = <LOGFILE>;
         
         # Read in coordinate for each atom
         for ($i = 0; $i < $NAtoms;$i++){
     
                $_ = <LOGFILE>;
                
                # Match line to obtain:             X           Y           Z 
                #  1          7      20001020       41.507000   23.202000   52.275000
                #
                
                /^\s+.+\s+.+\s+.+\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)/;
                $currentX = $1;
                $currentY = $2;
                $currentZ = $3;
                print $_  if $debug >= 6;
                print "XYZ for atom ",$i+1,": $1\t$2\t$3\n" if $debug >= 5;
          
                $oniomresp->{CARTCOORDGfile}[3*$i]   = $currentX;
                $oniomresp->{CARTCOORDGfile}[3*$i+1] = $currentY;
                $oniomresp->{CARTCOORDGfile}[3*$i+2] = $currentZ;
              }
       
       } #M1NO8
   
     } #M1NO7
  
  
  } #M1NO1
  
  close (LOGFILE);


 if ( $HaveListFile == 0 ) { #M1N15
 
  # Got everything, now it's time to generate structure with capping atoms
  
  # Grab all atoms that desired for output structure
  # without -mid flag: H layer only for both 2- and 3-layer ONIOM jobs
  # with    -mid flag: H layer for 2-layer ONIOM job
  #                    H and M layers for 3-layer ONIOM job
  #
  
  $modelindex = 0;
  
  for ($i = 0; $i < $NAtoms; $i++) { #M1NO9
      # Grab H layer for sure
      if ( ($oniomresp->{AtomLayerGfile}[$i] eq "H") ||
           ($oniomresp->{AtomLayerGfile}[$i] eq "h")  ) {
           
             $modelatom[$modelindex] = $oniomresp->{ATOMSGfile}[$i] ;
             $modelX[$modelindex]    = $oniomresp->{CARTCOORDGfile}[3*$i];
             $modelY[$modelindex]    = $oniomresp->{CARTCOORDGfile}[3*$i+1];
             $modelZ[$modelindex]    = $oniomresp->{CARTCOORDGfile}[3*$i+2];
             $modelindex++;
             print "\n Grab atom ",$oniomresp->{ATOMSGfile}[$i]," (",$i+1,") in H layer as atom $modelindex in output structure."  if $debug >= 3;
           }
           
       # Grab M layer if necessary
       if ( ( $setmid == 1 ) && 
            ( ($oniomresp->{AtomLayerGfile}[$i] eq "M") ||
              ($oniomresp->{AtomLayerGfile}[$i] eq "m")    ) ) {
             $modelatom[$modelindex] = $oniomresp->{ATOMSGfile}[$i] ;
             $modelX[$modelindex]    = $oniomresp->{CARTCOORDGfile}[3*$i];
             $modelY[$modelindex]    = $oniomresp->{CARTCOORDGfile}[3*$i+1];
             $modelZ[$modelindex]    = $oniomresp->{CARTCOORDGfile}[3*$i+2];
             $modelindex++;
             print "\n Grab atom ",$oniomresp->{ATOMSGfile}[$i]," (",$i+1,") in M layer as atom $modelindex in output structure."  if $debug >= 3;
           }
   } #M1NO9
   
   print "\n\n"  if $debug >= 3;
  
  # Start writing to output file for Gaussian calculation.
  
  print "Open output file $outputfile to generate Gaussian input for ESP calculation......\n\n"  if $debug >= 1;
  open(OUTPUTFILE,">$outputfile") || die "Could not write $outputfile\n$!\n";
  
  print  OUTPUTFILE "%chk=$outputfile.chk\n";
  print  OUTPUTFILE "%mem=2000MB\n";
  print  OUTPUTFILE "%nprocshared=2\n";
  print  OUTPUTFILE "#p hf/6-31g(d) iop(6/33=2) pop=mk sp test scf=(xqc,tight)\n\n";
  print  OUTPUTFILE "Gaussian Calc for RESP fitting based on $gfile\n";
  
  if ( $is_twolayer == 1 ) {
     print  OUTPUTFILE "two-layer oniom calculation. H(igh) layer is extracted.\n";
   } elsif ($setmid == 1) {
      print  OUTPUTFILE "three-layer oniom calculation. H(igh) and M(iddle) region are extracted.\n";
     } else {
      print  OUTPUTFILE "three-layer oniom calculation. H(igh) region is extracted.\n";
     }
  
  
   
  # Treat capping atom. Calculate coordinates using scaling factors obtained form log file
  # Find each link atom
  
  # We use $chosenLAYER  and $chosenlayer to control link host atoms in which layers are used to generate link atoms.
  #  In two layer ONIOM ($is_twolayer == 1), we use atoms in L layer anyway
  #  In 3-layer ONIOM with -mid flag ($setmid == 1), we also use atoms in L layer 
  #  In 3-layer ONIOM without -mid flag ($setmid == 0), we also use atoms in M layer 
  
  $capHnum = 0;
  
  if ( ($is_twolayer == 1) || ( $setmid == 1) ) {
      $chosenLAYER = "L";
      $chosenlayer = "l";
    } else {
      $chosenLAYER = "M";
      $chosenlayer = "m";
    }
  
    # Locate All Link host atoms
    # $i runs from 0 to $link_indexIn
    
    for ($i = 0;$i < $link_indexIn; $i++) { #M1N11
        
       # Sanity check for these link related atoms
       
       if ( $debug >= 1) { #M1N12
       
         # Variable for sanity check
         $linkmatch = 0;
       
         # Check if info indexed by $link_index matches with info indexed by $link_indexIn
         for ($j = 0; $j < $link_indexIn; $j++) { #M1N13
           
           if ( $LinkIndex[$i] == $LinkIndexIn[$j] ) { #M1N14
               
               $linkmatch = 1;
               
               if (  $HostIndex[$i] != $HostIndexIn[$j]  ) {
                  $linkmatch = 0;
                  print "\n Link atoms sanity check:\n";
                  print " Atom ",$LinkIndex[$i]," in MM layer is connected to atom ",$HostIndex[$i]," according to ONIOM: Cut between section.\n";
                  print "                          but to atom ",$HostIndexIn[$j]," according to input section.\n";
                } elsif ( $HostEle[$i] ne $oniomresp->{ATOMSGfile}[ $HostIndexIn[$j]-1 ]) {
                  $linkmatch = 0;
                  print "\n Link atoms sanity check:\n";
                  print " Host for atom ",$LinkIndex[$i]," is a ",$HostEle[$i]," according to ONIOM: Cut between section.\n";
                  print "                          but is a ",$oniomresp->{ATOMSGfile}[ $HostIndexIn[$j]-1 ]," according to input section.\n";
                }
                
               if (  $LinkHostEle[$i] ne $LinkHostEleIn[$j]  ) {
                  $linkmatch = 0;
                  print "\n Link atoms sanity check:\n";
                  print " Atom ",$LinkIndex[$i]," in MM layer is element ",$LinkHostEle[$i]," according to ONIOM: Cut between section.\n";
                  print "                          but is element ",$LinkHostEleIn[$j]," according to input section.\n";
                }
             
                
               if (  $LinkAtomEle[$i] ne $LinkAtomEleIn[$j]  ) {
                  $linkmatch = 0;
                  print "\n Link atoms sanity check:\n";
                  print " Atom ",$LinkIndex[$i]," in MM layer will be replaced by ",$LinkAtomEle[$i]," according to ONIOM: Cut between section.\n";
                  print "                          but by ",$LinkAtomEleIn[$j]," according to input section.\n";
                }
                
             
             } #M1N14
        
          } #M1N13
          
        } #M1N12
        
        # Output warning if there is mismatch
        if ( ($debug >= 1) && ( $linkmatch == 0 ) ) {
           print "\n There are mismatch information about link atoms ",$LinkIndex[$i],".\n";
         }
       

       $currLinkHost = $LinkIndex[$i] - 1;
       
       $currHost     = $HostIndex[$i] - 1;

       # Check if we need to use this atom to generate link atom
       
       if ( ($oniomresp->{AtomLayerGfile}[$currLinkHost] ne $chosenLAYER) && 
            ($oniomresp->{AtomLayerGfile}[$currLinkHost] ne $chosenlayer)    ) {next;}
       
       
       $capHnum++;
       
       # Calculate new coordinate of link atom
       if ($is_twolayer == 1) {
         if ( $Scale1[$i] != $Scale2[$i] )  {
           print "\n Atom ",$LinkIndex[$i]," has two different scaling facotr for bond length: ",$Scale1[$i]," and ",$Scale2[$i],"\n";
           print "\n  Use ",$Scale1[$i]," to calculate coordinates of link atom.\n";
           }
         }
       
       if ($is_threelayer == 1) {
         if ( ($Scale1[$i] != $Scale2[$i]) || ($Scale1[$i] != $Scale3[$i]) )  {
           print "\n Atom ",$LinkIndex[$i]," has two different scaling facotr for bond length: ",
                 $Scale1[$i],", ",$Scale2[$i]," and ",$Scale3[$i],"\n";
                 
           print "\n  Use ",$Scale1[$i]," to calculate coordinates of link atom.\n";
           }
         }
       
       $currentX     = ( $oniomresp->{CARTCOORDGfile}[3*$currLinkHost] - $oniomresp->{CARTCOORDGfile}[3*$currHost] ) * $Scale1[$i] 
                         + $oniomresp->{CARTCOORDGfile}[3*$currHost];
                         
       $currentY     = ( $oniomresp->{CARTCOORDGfile}[3*$currLinkHost+1] - $oniomresp->{CARTCOORDGfile}[3*$currHost+1] ) * $Scale1[$i] 
                         + $oniomresp->{CARTCOORDGfile}[3*$currHost+1];
                         
       $currentZ    = ( $oniomresp->{CARTCOORDGfile}[3*$currLinkHost+2] - $oniomresp->{CARTCOORDGfile}[3*$currHost+2] ) * $Scale1[$i] 
                         + $oniomresp->{CARTCOORDGfile}[3*$currHost+2];
       
       # Original bond distance
       $oribondlen = sqrt( ($oniomresp->{CARTCOORDGfile}[3*$currLinkHost] - $oniomresp->{CARTCOORDGfile}[3*$currHost]) ** 2 +
                           ( $oniomresp->{CARTCOORDGfile}[3*$currLinkHost+1] - $oniomresp->{CARTCOORDGfile}[3*$currHost+1] ) ** 2 +
                           ( $oniomresp->{CARTCOORDGfile}[3*$currLinkHost+2] - $oniomresp->{CARTCOORDGfile}[3*$currHost+2] ) ** 2 );

       # Bond distance after scaling

       $linkbondlen = sqrt( ($currentX - $oniomresp->{CARTCOORDGfile}[3*$currHost]) ** 2 +
                            ($currentY - $oniomresp->{CARTCOORDGfile}[3*$currHost+1]) ** 2 +
                            ($currentZ - $oniomresp->{CARTCOORDGfile}[3*$currHost+2]) ** 2 );
       
                         
       $modelatom[$modelindex] = $LinkAtomEle[$i];
       $modelX[$modelindex]    = $currentX;
       $modelY[$modelindex]    = $currentY;
       $modelZ[$modelindex]    = $currentZ;
       $modelindex++;
       
       if ($debug >= 3) {
            
            printf "\n Coordinate of cap atom %3s at X: %13.6f  Y: %13.6f  Z: %13.6f\n",$LinkAtomEle[$i],$currentX,$currentY,$currentZ;
            
            printf " replaces atom %3s (%6d) at X: %13.6f  Y: %13.6f  Z: %13.6f, in layer %2s\n",$oniomresp->{ATOMSGfile}[ $currLinkHost ],
                                                                                   $LinkIndex[$i],
                                                                                   $oniomresp->{CARTCOORDGfile}[3*$currLinkHost],
                                                                                   $oniomresp->{CARTCOORDGfile}[3*$currLinkHost+1],
                                                                                   $oniomresp->{CARTCOORDGfile}[3*$currLinkHost+2],
                                                                                   $oniomresp->{AtomLayerGfile}[ $currLinkHost ];
                                                                                        
            printf " connects host %3s (%6d) at X: %13.6f  Y: %13.6f  Z: %13.6f, in layer %2s\n",$oniomresp->{ATOMSGfile}[ $HostIndex[$i]-1 ],
                                                                                   $HostIndex[$i],
                                                                                   $oniomresp->{CARTCOORDGfile}[3*$currHost],
                                                                                   $oniomresp->{CARTCOORDGfile}[3*$currHost+1],
                                                                                   $oniomresp->{CARTCOORDGfile}[3*$currHost+2],
                                                                                   $oniomresp->{AtomLayerGfile}[ $currHost ];
                                                                                   
            printf " Original bond length: %13.6f, scaled bond length: %13.6f, scale factor: %13.6f.\n", $oribondlen, $linkbondlen,$Scale1[$i];
         }
       
       # Put link atom into $modelatom array
      
      
      } #M1N11
   
   
   print  OUTPUTFILE "There are $capHnum capping atoms.\n\n";
   
   if ( $capHnum > 1 ) {
       print  "\n There are $capHnum capping atoms added.\n\n" if $debug >= 0;
      } elsif (  $capHnum > 0 ) {
       print  "\n There is one capping atom added.\n\n" if $debug >= 0;
      } else {
       print  "\n No capping atom added.\n\n" if $debug >= 0;
      }
   
   
   if ( $setmid == 1 ) { 
       print OUTPUTFILE "$MidCharge  $MidMultiplicity\n";
     } else {
      print OUTPUTFILE "$HighCharge  $HighMultiplicity\n";
     }
 
  for ($i = 0; $i< $modelindex; $i++ )   {
      printf OUTPUTFILE "%4s %15.7f %15.7f %15.7f\n",$modelatom[$i],$modelX[$i],$modelY[$i],$modelZ[$i];
   
   }
  print OUTPUTFILE "\n\n";


 } #M1N15

 # When there is a atom list file provided, we use it to generate output file


 if ( $HaveListFile == 1 ) { #M1N16

     # Start writing to output file for Gaussian calculation.
     
     print "Open output file $outputfile to generate Gaussian input for ESP calculation......\n\n"  if $debug >= 1;
     open(OUTPUTFILE,">$outputfile") || die "Could not write $outputfile\n$!\n";
     
     print  OUTPUTFILE "%chk=$outputfile.chk\n";
     print  OUTPUTFILE "%mem=2000MB\n";
     print  OUTPUTFILE "%nprocshared=2\n";
     print  OUTPUTFILE "#p hf/6-31g(d) iop(6/33=2) pop=mk sp test scf=(xqc,tight)\n\n";
     print  OUTPUTFILE "Gaussian Calc for RESP fitting based on $gfile\n";
     print  OUTPUTFILE "Atoms are extracted using atom list file $AtomListFile.\n";
     print  OUTPUTFILE "No charge and multiplicity info available. Using 0 1, please correct it.\n";
     

   
   # Read in Atom List file 
   #
   #
   
   
  # open given Atom List File
  #  Read in atom list file like this
  # 4303  4300 
  # ...
  # 4310
  # 4320 4317 0.9000
  # 4321
  # ...
  # 4325
  # 4362 4359 0.723886 H
  # ...
  # 4369
  # 4447
  # ...
  # 4454
  # 5010
  # 5014
  # ...
  # 5023  5026  0.823886 F
  # 5024
  # 5025
  # 5054
  # 5055
  # 5056
  
  print "Open atom list file $AtomListFile ......\n\n"  if $debug >= 1;
  open(ATOMLISTFILE,$AtomListFile) || die "Could not read $AtomListFile\n$!\n";
  
  $List_Index = 0;
  $skipindicator = 0;
  @ListAtom = ();
  @ListHost = ();
  @ListScale= ();
  @ListLinkElement = ();
  @HasLink = ();

  while (<ATOMLISTFILE>) { #M3N17
  
   next if /^\s*$/;
    
     if (/^\.\.\./) {
       # Need to treat skip mode: e.g. 4447 ... 4454
       #
       print " Set \$skipindicator as 1\n"  if $debug >= 5;
       $skipindicator = 1;
       } elsif ($skipindicator == 0) { #M3N18
        # Not in skip mode, just read in atom info
        # Atom id only
        if (/\s*(\d+)\s*$/) {
           $ListAtom[$List_Index] = $1;
           $HasLink[$List_Index] = 0;
           $ListHost[$List_Index] = $undefined;
           $ListScale[$List_Index] = $undefined;
           $ListLinkElement[$List_Index] = $undefined;
           
         }
        
        # Atom ID   Host Atom ID

        if (/\s*(\d+)\s+(\d+)\s*$/) {
           $ListAtom[$List_Index] = $1;
           $ListHost[$List_Index] = $2;
           $ListScale[$List_Index] = $DefaultScale;
           $ListLinkElement[$List_Index] = $DefaultLinkAtom;
           $HasLink[$List_Index] = 1;
         }
         
        # Atom ID   Host Atom ID    Scaling Factor

        if (/\s*(\d+)\s+(\d+)\s+(\d*.\d+)\s*$/) {
           $ListAtom[$List_Index] = $1;
           $ListHost[$List_Index] = $2;
           $ListScale[$List_Index] = $3;
           $ListLinkElement[$List_Index] = $DefaultLinkAtom;
           $HasLink[$List_Index] = 1;
         }
         
        # Atom ID   Host Atom ID    Scaling Factor   Lint Element

        if (/\s*(\d+)\s+(\d+)\s+(\d*.\d+)\s+([a-zA-Z]+)\s*$/) {
           $ListAtom[$List_Index] = $1;
           $ListHost[$List_Index] = $2;
           $ListScale[$List_Index] = $3;
           $ListLinkElement[$List_Index] = $4;
           $HasLink[$List_Index] = 1;
         }
        
        if ($List_Index == 0) {
           print "\nRead in atom list file $AtomListFile.\n\n" if $debug >= 3;
         }
        
        if ( $HasLink[$List_Index] == 1 ) {
           printf " %4d: atom %6d connects to atom %6d, using scaling factor %10.7f for link atom %3s.\n",
                  $List_Index+1,$ListAtom[$List_Index],$ListHost[$List_Index],
                  $ListScale[$List_Index],$ListLinkElement[$List_Index] if $debug >= 3;
          } else {
            printf " %4d: atom %6d\n",$List_Index+1,$ListAtom[$List_Index] if $debug >= 3;
          
          }
         
        
        $List_Index++;
       
       }   #M3N18
         else {  #M3N19
            
           # In skip mode, need to read this line, add missing atom list, then put current atom in
           # 
           #
           $skipfrom = $ListAtom[$List_Index-1];
           $skipindicator = 0;
           
             # Atom id only
             if (/\s*(\d+)\s*$/) {
                 $tempAtom = $1;
                 $tempHost = $undefined;
                 $tempScale = $undefined;
                 $tempLinkElement = $undefined;
                 $tempHasLink = 0;
                }
             
             # Atom ID   Host Atom ID
     
             if (/\s*(\d+)\s+(\d+)\s*$/) {
                $tempAtom = $1;
                $tempHost = $2;
                $tempScale = $DefaultScale;
                $tempLinkElement = $DefaultLinkAtom;
                $tempHasLink = 1;
              }
              
             # Atom ID   Host Atom ID    Scaling Factor
     
             if (/\s*(\d+)\s+(\d+)\s+(\d*.\d+)\s*$/) {
                 $tempAtom = $1;          
                 $tempHost = $2;          
                 $tempScale = $3;
                 $tempLinkElement = $DefaultLinkAtom;
                 $tempHasLink = 1;               
              }
              
             # Atom ID   Host Atom ID    Scaling Factor   Lint Element
     
             if (/\s*(\d+)\s+(\d+)\s+(\d*.\d+)\s+([a-zA-Z]+)\s*$/) {
                 $tempAtom = $1;                 
                 $tempHost = $2;                 
                 $tempScale = $3;                
                 $tempLinkElement = $4;
                 $tempHasLink = 1;      
              }
             
                          
         if ( $tempAtom < $skipfrom ) { #M3N20
       
             print "\n Attention: atom index $skipfrom larger than $tempAtom, will run from $skipfrom down to $tempAtom.\n\n"  if $debug >= 0;

             for ($i = $skipfrom-1; $i > $tempAtom; $i--) {
                  printf " %4d: atom %6d\n",$List_Index+1, $i if $debug >= 3;
                  $ListAtom[$List_Index] = $i;
                  $ListHost[$List_Index] = $undefined;
                  $ListScale[$List_Index] = $undefined;
                  $ListLinkElement[$List_Index] = $undefined;
                  $HasLink[$List_Index] = 0;
                  $List_Index++;
                }
         
          
          } #M3N20
            elsif ( $tempAtom > $skipfrom) { #M3N21

             for ($i = $skipfrom+1; $i < $tempAtom; $i++) {
                  printf " %4d: atom %6d\n",$List_Index+1, $i if $debug >= 3;
                  $ListAtom[$List_Index] = $i;
                  $ListHost[$List_Index] = $undefined;
                  $ListScale[$List_Index] = $undefined;
                  $ListLinkElement[$List_Index] = $undefined;
                  $HasLink[$List_Index] = 0;
                  $List_Index++;
                }
            
             } #M3N21
               else {
                  print "\n Warning: two numbers before and after ... are the same as $skipfrom.\n" if $debug >= 0;
                  print "\n  This will create duplicate atoms in output file. Please double check your atom list file $AtomListFile.\n" if $debug >= 0;
            }
           
           
            # After adding skipped atoms, we need to add the last one into the list
            if ( $tempHasLink == 1 ) {
               printf " %4d: atom %6d connects to atom %6d, using scaling factor %10.7f for link atom %3s.\n",
                      $List_Index+1,$tempAtom,$tempHost,
                      $tempScale,$tempLinkElement if $debug >= 3;
              } else {
                printf " %4d: atom %6d\n",$List_Index+1,$tempAtom if $debug >= 3;
              }
                
            $ListAtom[$List_Index] =      $tempAtom ;
            $ListHost[$List_Index] =        $tempHost ;
            $ListScale[$List_Index] =       $tempScale ;
            $ListLinkElement[$List_Index] =   $tempLinkElement;
            $HasLink[$List_Index] =    $tempHasLink ;
            $List_Index++;
         
         }  #M3N19

  
   } #M3N17
   
   print "\n" if $debug >= 3;
   
  # Add atoms listed in atom list file to final model atom list first
  #  if there is any, the capping atom will be added later.
  
  $modelindex = 0;
  $listLinkAtomNum = 0;

  printf "\n Structure extracted using atom list file $AtomListFile.\n\n" if  $debug >= 4;
  
  for ($i = 0; $i < $List_Index; $i++) { #M1N22
     $modelatom[$modelindex] = $oniomresp->{ATOMSGfile}[ $ListAtom[$i]-1 ];
     $modelX[$modelindex]    = $oniomresp->{CARTCOORDGfile}[ 3*($ListAtom[$i]-1) ];
     $modelY[$modelindex]    = $oniomresp->{CARTCOORDGfile}[ 3*($ListAtom[$i]-1) + 1 ];
     $modelZ[$modelindex]    = $oniomresp->{CARTCOORDGfile}[ 3*($ListAtom[$i]-1) + 2 ];
     $modelindex++;
     printf "  %5d: %4s  %15.7f  %15.7f  %15.7f\n",$modelindex, $modelatom[$modelindex-1],
            $modelX[$modelindex-1],$modelY[$modelindex-1],$modelZ[$modelindex-1] if  $debug >= 4;
            
     printf "  From line %11d in $gfile as\n",$oniomresp->{LineNumGfile}[$ListAtom[$i]-1] if  $debug >= 5;
     
     print $oniomresp->{Original_Line}[$ListAtom[$i]-1] if  $debug >= 5;
     
     $listLinkAtomNum += $HasLink[$i];
     
   } #M1N22
  
  # We need to create link atoms, if there are any.
  if ( $listLinkAtomNum > 0 ) { #M1N23
    printf "\n There are $listLinkAtomNum link atoms need to be constructed according to $AtomListFile.\n\n" if  $debug >= 1;
    
    # Go through list and treat link atoms 
    for ($i = 0; $i < $List_Index; $i++) { #M1N24
       
       if ( $HasLink[$i] == 0 ) { next;}

       $currHost = $ListAtom[$i] - 1;

       $currLinkHost     = $ListHost[$i] - 1;


       $currentX     = ( $oniomresp->{CARTCOORDGfile}[3*$currLinkHost] - $oniomresp->{CARTCOORDGfile}[3*$currHost] ) * $ListScale[$i] 
                         + $oniomresp->{CARTCOORDGfile}[3*$currHost];
                         
       $currentY     = ( $oniomresp->{CARTCOORDGfile}[3*$currLinkHost+1] - $oniomresp->{CARTCOORDGfile}[3*$currHost+1] ) * $ListScale[$i] 
                         + $oniomresp->{CARTCOORDGfile}[3*$currHost+1];
                         
       $currentZ    = ( $oniomresp->{CARTCOORDGfile}[3*$currLinkHost+2] - $oniomresp->{CARTCOORDGfile}[3*$currHost+2] ) * $ListScale[$i] 
                         + $oniomresp->{CARTCOORDGfile}[3*$currHost+2];


       # Original bond distance
       $oribondlen = sqrt( ($oniomresp->{CARTCOORDGfile}[3*$currLinkHost] - $oniomresp->{CARTCOORDGfile}[3*$currHost]) ** 2 +
                           ( $oniomresp->{CARTCOORDGfile}[3*$currLinkHost+1] - $oniomresp->{CARTCOORDGfile}[3*$currHost+1] ) ** 2 +
                           ( $oniomresp->{CARTCOORDGfile}[3*$currLinkHost+2] - $oniomresp->{CARTCOORDGfile}[3*$currHost+2] ) ** 2 );

       # Bond distance after scaling

       $linkbondlen = sqrt( ($currentX - $oniomresp->{CARTCOORDGfile}[3*$currHost]) ** 2 +
                            ($currentY - $oniomresp->{CARTCOORDGfile}[3*$currHost+1]) ** 2 +
                            ($currentZ - $oniomresp->{CARTCOORDGfile}[3*$currHost+2]) ** 2 );

       
       $modelatom[$modelindex] = $ListLinkElement[ $i ];
       $modelX[$modelindex]    =$currentX ;
       $modelY[$modelindex]    =$currentY ;
       $modelZ[$modelindex]    =$currentZ ;
       $modelindex++;
       
       if ($debug >= 4) {
            
            printf "\n Coordinate of cap atom %3s at X: %13.6f  Y: %13.6f  Z: %13.6f\n",$modelatom[$modelindex-1],$currentX,$currentY,$currentZ;
            
            printf " uses atom     %3s (%6d) at X: %13.6f  Y: %13.6f  Z: %13.6f, in layer %2s\n",$oniomresp->{ATOMSGfile}[ $currLinkHost ],
                                                                                   $currLinkHost+1,
                                                                                   $oniomresp->{CARTCOORDGfile}[3*$currLinkHost],
                                                                                   $oniomresp->{CARTCOORDGfile}[3*$currLinkHost+1],
                                                                                   $oniomresp->{CARTCOORDGfile}[3*$currLinkHost+2],
                                                                                   $oniomresp->{AtomLayerGfile}[ $currLinkHost ];
                                                                                        
            printf " connects host %3s (%6d) at X: %13.6f  Y: %13.6f  Z: %13.6f, in layer %2s\n",$oniomresp->{ATOMSGfile}[ $currHost ],
                                                                                   $currHost+1,
                                                                                   $oniomresp->{CARTCOORDGfile}[3*$currHost],
                                                                                   $oniomresp->{CARTCOORDGfile}[3*$currHost+1],
                                                                                   $oniomresp->{CARTCOORDGfile}[3*$currHost+2],
                                                                                   $oniomresp->{AtomLayerGfile}[ $currHost ];
                                                                                   
            printf " Original bond length: %13.6f, scaled bond length: %13.6f, scale factor: %13.6f.\n", $oribondlen, $linkbondlen,$ListScale[$i];
         }
       
       
       printf "  %5d: %4s  %15.7f  %15.7f  %15.7f\n",$modelindex, $modelatom[$modelindex-1],
              $modelX[$modelindex-1],$modelY[$modelindex-1],$modelZ[$modelindex-1] if  $debug >= 4;
              
       printf "  From line %11d in $gfile as\n",$oniomresp->{LineNumGfile}[$ListAtom[$i]-1] if  $debug >= 5;
       
       print $oniomresp->{Original_Line}[$ListAtom[$i]-1] if  $debug >= 5;
       
       
     } #M1N24
     
   
   } #M1N23
   
   print  OUTPUTFILE "There are $listLinkAtomNum capping adtoms added.\n\n0 1\n";
 
   for ($i = 0; $i< $modelindex; $i++ )   {
      printf OUTPUTFILE "%4s %15.7f %15.7f %15.7f\n",$modelatom[$i],$modelX[$i],$modelY[$i],$modelZ[$i];
   
    }
   print OUTPUTFILE "\n\n";


 
  } #M1N16


 } # sub modeone {


# Subroutine for running mode two
# Generating RESP setup .in file for given log file

sub modetwo {
   my $oniomresp = shift;
   my $debug     = $oniomresp->{DEBUG};
   my $gfile     = $oniomresp->{G_FILE_NAME};
   my $outputfile = $oniomresp->{OUTPUT_FILE};
   my $mode = $oniomresp->{MODE};
   my $capHnum = $oniomresp->{CapH_Num} ;
   my $stepnum = $oniomresp->{Step_Num} ;
   my $mapfile = $oniomresp->{Map_File} ;
   my $havemapfile = $oniomresp->{Have_Map_File} ;
   my $chgcompfile = $oniomresp->{Charge_Comp_File};
   my $qinfile = $oniomresp->{Qin_File} ;
   my $isLogfile = $oniomresp->{IS_LOG_FILE};
   my $isJobfile = $oniomresp->{IS_JOB_FILE};
   my $NAtoms = 0;
   my $isnormal = 0;
   my ($i,$atomindex);
   my ($logatom,$currentX,$currentY,$currentZ,$currentline);

   
# open Gaussian Log file for reading or display error
  
  print "Open Gaussian Log file $gfile ......\n\n"  if $debug >= 1;
  open(LOGFILE,$gfile) || die "Could not read $gfile\n$!\n";

# grab everything which may be useful
  while (<LOGFILE>){
    
 	# skip blank lines
	next if /^$/;

        if (/^\s*NAtoms=\s*(\d+)/i) {
           $NAtoms=$1;
#           print "There are $NAtoms atoms.\n";
         }

        # Check if is normal termination 
        if ( /^\s*Normal\s+termination/i ) {
                $isnormal  = 1;
        	next;
         }

# Find out charge and multiplicity
        
        if ( /^\s*Charge\s*=\s*(-*\d+)\s+Multiplicity\s*=\s*(\d+)/i ) {
           print "In $gfile, charge is $1, multiplicity is $2.\n\n" if $debug >= 1;
           $oniomresp->{Charge} = $1;
           $oniomresp->{Multiplicity} = $2;
          }
          
   } # while (<LOGFILE>){

 close (LOGFILE);
   
# If not terminated normally
   if ( $isnormal == 0 )  {
     print "Warning: It occurs that $gfile does not terminate normally.\n\n" if $debug >= 0;
     print "Please find out why before you proceed on RESP charges fitting.\n\n" if $debug >= 0;
   
    }

# If NAtoms information is not available, we need to find it by
# counting in standard orientation.

   if ( $NAtoms == 0 ) {
      print "NAtoms info is not presented in $gfile. \n\nTrying to find number of atoms in another way...\n\n" if $debug >= 0;
      
      open(LOGFILE,$gfile) || die "Could not read $gfile\n$!\n";

      while (<LOGFILE>){
 
           if ((/^\s*Z-Matrix\s+orientation:/i) || (/^\s*Standard\s+orientation:/i)) {

               $_ = <LOGFILE>;
               if ( not(/^\s+-+/)) {
                  print "Line($_)\n" if $debug >= 3;
                  print "Gaussian log format may change. Pay attention to output structure(s).\n" if  $debug >= 0;
                  }
               $_ = <LOGFILE>;
               if ( not(/^\s*Center/i) ) {
                  print "Line($_)\n" if $debug >= 3;
                  print "Gaussian log format may change. Pay attention to output structure(s).\n"  if $debug >= 0;
                  }

               $_ = <LOGFILE>;
               if ( not(/^\s*Number\s+Number\s+Type\s+X\s+Y\s+Z/i) ) {
                  print "Line($_)\n" if $debug >= 3;
                  print "Gaussian log format may change. Pay attention to output structure(s).\n"  if $debug >= 0;
                  }

               $_ = <LOGFILE>;

               if ( not(/^\s+-+/)) {
                  print "Line($_)\n" if $debug >= 3;
                  print "Gaussian log format may change. Pay attention to output structure(s).\n" if  $debug >= 0;
                  }
               
               $_ = <LOGFILE>;
               
               while (not(/^\s+-+/)) {
                  $_ = <LOGFILE>;
                  $NAtoms++;
                } # while (not(/^\s+-+/) {
              }

          } # while LOGFILE

       } # if  ( $NAtoms == 0 ) {
   
   print "There are $NAtoms atoms in file $gfile.\n\n" if $debug >= 0;
   
   $oniomresp->{NAtoms} = $NAtoms;
   
   close (LOGFILE);
   
   
# Read in atom info from log file

      open(LOGFILE,$gfile) || die "Could not read $gfile\n$!\n";

      while (<LOGFILE>){
 
           if ((/^\s*Z-Matrix\s+orientation:/i) || (/^\s*Standard\s+orientation:/i)) { #NO5

               $_ = <LOGFILE>;
               if ( not(/^\s+-+/)) {
                  print "Line($_)\n" if $debug >= 3;
                  print "Gaussian log format may change. Pay attention to output structure(s).\n" if  $debug >= 0;
                  }
               $_ = <LOGFILE>;
               if ( not(/^\s*Center/i) ) {
                  print "Line($_)\n" if $debug >= 3;
                  print "Gaussian log format may change. Pay attention to output structure(s).\n"  if $debug >= 0;
                  }

               $_ = <LOGFILE>;
               if ( not(/^\s*Number\s+Number\s+Type\s+X\s+Y\s+Z/i) ) {
                  print "Line($_)\n" if $debug >= 3;
                  print "Gaussian log format may change. Pay attention to output structure(s).\n"  if $debug >= 0;
                  }

               $_ = <LOGFILE>;

               if ( not(/^\s+-+/)) {
                  print "Line($_)\n" if $debug >= 3;
                  print "Gaussian log format may change. Pay attention to output structure(s).\n" if  $debug >= 0;
                  }

# Read into atom info
               
             for ($i = 0; $i<$NAtoms;$i++){ #NO6
                $_ = <LOGFILE>;
                print $_ if  $debug >= 6;

                /^\s+.+\s+(\d+)\s+\d+\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)/;

                if ( $debug >= 5) {printf "Taking atom %3d coordinates (%16.10f %16.10f %16.10f) from log file.\n",$1,$2,$3,$4;}
                $oniomresp->{Element_NUM}[$i]   = $1;
                $oniomresp->{CARTCOORD}[3*$i]   = $2;
                $oniomresp->{CARTCOORD}[3*$i+1] = $3;
                $oniomresp->{CARTCOORD}[3*$i+2] = $4;
                } #NO6

               last;

              } #NO5

          } # while LOGFILE

# Now create RESP setup file for this Gaussian Log file.

  print "Open RESP setup file $outputfile for writing ......\n\n"  if $debug >= 1;
  open(RESPFILE,">$outputfile") || die "Could not read $outputfile\n$!\n";
  
  print RESPFILE "Resp charges fitting for $gfile with zero charge on last $capHnum atoms\n";
  print RESPFILE "\n";
  print RESPFILE " &cntrl\n";
  print RESPFILE "\n";
  print RESPFILE " nmol = 1,\n";
  print RESPFILE " ihfree = 1,\n";
  print RESPFILE " ioutopt = 1,\n";
  print RESPFILE " qwt = 0.0005,\n";
  print RESPFILE "\n";
  print RESPFILE " &end\n";
  print RESPFILE "    1.0\n";
  print RESPFILE "Resp charges fitting setup\n";
  printf RESPFILE "%5d%5d\n",$oniomresp->{Charge}, $NAtoms;
  
# print out each atom  
  for ($i = 0; $i<$NAtoms;$i++){ #NO7
    printf RESPFILE "%5d%5d\n",$oniomresp->{Element_NUM}[$i], 0;
   } #NO7
  
# print zero constraint for capping atom  
  if ($capHnum > 0 ) {
    for ($i = $NAtoms-$capHnum; $i<$NAtoms;$i++) {
        print RESPFILE "    1    0.0\n";
        printf RESPFILE "%5d%5d\n",1,$i+1;
      }
   }
  print RESPFILE "\n";
  print RESPFILE "\n";
   
  close(RESPFILE);
  
  print "Finish writing $outputfile!\n\n "  if $debug >= 0;
  
} # sub modetwo {


# Subroutine for running mode three
# Taking RESP charge, and put it to ONIOM input file.

sub modethree {

   my $oniomresp    = shift;
   my $debug        = $oniomresp->{DEBUG};
   my $gfile        = $oniomresp->{G_FILE_NAME};
   my $outputfile   = $oniomresp->{OUTPUT_FILE};
   my $mode         = $oniomresp->{MODE};
   my $capHnum      = $oniomresp->{CapH_Num} ;
   my $stepnum      = $oniomresp->{Step_Num} ;
   my $mapfile      = $oniomresp->{Map_File} ;
   my $havemapfile  = $oniomresp->{Have_Map_File} ;
   my $chgcompfile  = $oniomresp->{Charge_Comp_File};
   my $qinfile      = $oniomresp->{Qin_File} ;
   my $isLogfile    = $oniomresp->{IS_LOG_FILE};
   my $isJobfile    = $oniomresp->{IS_JOB_FILE};
   my $NAtoms       = 0;
   my $isnormal     = 0;
   my ($i,$j);
   my ($atomindex,$lineindex,$structindicator,$structurestart,$curlayer);
   my ($logatom,$currentX,$currentY,$currentZ,$currentline);
   my %layers = ();
   my %layersatm = ();
   my (@currentchg,@RESPchg,$totalchgnum,$RESPchgindex,$currnewchg,$currnewline);
   my (@map,$mapindex,$skipindicator,$currmapline,$skipfrom);
   

   $layers{'H'} = 0;
   $layers{'M'} = 0;
   $layers{'L'} = 0;
   
   $layersatm{'H'} = 0;
   $layersatm{'M'} = 0;
   $layersatm{'L'} = 0;

# open Gaussian ONIOM Input file for reading or display error

   print "Open Gaussian ONIOM input file $gfile ......\n\n"  if $debug >= 1;
   open(ONIOMINFILE,$gfile) || die "Could not read $gfile\n$!\n";

   $atomindex       = 0;
   $lineindex       = 0;
   $structindicator = 0;

# Read in all atom info from Gaussian ONIOM Input file

   while (<ONIOMINFILE>) {  #M3NO1
   
      $lineindex++;

# Locate start of geometry info      
      if   ($structindicator == 0){
        if (/^\s*-*\d+\s+\d+\s+(-*\d+\s+\d+)/) {
   
          print $_ if $debug >= 3;
          $structindicator = 1;
          $structurestart = $lineindex;
          next;
        } # if (/^\s*-*\d+\s+\d+\s+-*\d+\s+\d+/) 
     
        next;
      }  # if   ($structindicator ==0){
      else {
        if (/^\s+$/ or /^$/){ last;}
        print $_ if $debug >= 6;
        
         $oniomresp->{Original_Line}[$atomindex] = $_;

         #  Read in each atom information
         /^\s*([a-zA-Z]+)-([a-zA-Z0-9\*]+)-(-*\d*\.\d+)\s+(-*[01])\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+([HML])/;
            
         print "$1\t$2\t$3\t$4\t$5\t$6\t$7\t$8\t$atomindex\n" if $debug >= 5; 
         $oniomresp->{ATOMSGfile}[$atomindex]          = $1;
         $oniomresp->{ATOMTYPEGfile}[$atomindex]       = $2;
         $oniomresp->{PARTIALCHARGGfile}[$atomindex]   = $3;
         $oniomresp->{ATOMFREEZEFLAGGfile}[$atomindex] = $4;
         $oniomresp->{CARTCOORDGfile}[3*$atomindex]    = $5;
         $oniomresp->{CARTCOORDGfile}[3*$atomindex+1]  = $6;
         $oniomresp->{CARTCOORDGfile}[3*$atomindex+2]  = $7;
         $oniomresp->{AtomLayerGfile}[$atomindex]      = $8;
         $oniomresp->{MovedGfile}[$atomindex]          = 0;
         $oniomresp->{LineNumGfile}[$atomindex]        = $lineindex;
         $oniomresp->{NeedNewRESPGfile}[$atomindex]    = 0;
         $oniomresp->{RESPChargeIndexGfile}[$atomindex]= -1;     
         # Count each layer
         $curlayer = $8;
         $atomindex++;
         $layersatm{$curlayer}++;
      
      } # else
   
    } # while (<ONIOMFILE>) { #M3NO1
   
   close(ONIOMINFILE);
   
   $oniomresp->{ATOMNUMGfile} = $atomindex;
   
   print "There are $atomindex atoms in file $gfile. \n"       if $debug >= 1; 
   print "There are ",$layersatm{'H'},"\t atoms in H layer.\n" if $debug >= 1;
   print "There are ",$layersatm{'M'},"\t atoms in M layer.\n" if $debug >= 1;
   print "There are ",$layersatm{'L'},"\t atoms in L layer.\n" if $debug >= 1;

   # Read in RESP charges 

   print "\nOpen RESP charge file $qinfile ......\n\n"  if $debug >= 0;
   open(CHARGEFILE,$qinfile) || die "Could not read $qinfile\n$!\n";
   
   $oniomresp->{RESP_charges} = ();
   
 # Read in all charges from RESP charge file

   while (<CHARGEFILE>) {  #M3NO2
     if (/^\s*$/){ last;}

     print "Read in charge file:\n" if $debug >= 5;
     print $_ if $debug >= 5; 
     
     @currentchg = split;
     
     push(@RESPchg,@currentchg);
     
     print "Split into individual charge\n" if $debug >= 5;
     
     for ($i=0;$i < @currentchg;$i++) {
         
         print "\"",$currentchg[$i]."\"  " if $debug >= 5;
       
        }

     print "\n" if $debug >= 5;

   } # end of #M3NO2
   
   $totalchgnum = @RESPchg;
   print "\nRead in $totalchgnum RESP charges.\n"  if $debug >= 4;
   
  
 # If no mapping file provided, RESP charges will be assigned to QM atoms in sequence.
 
 if ($havemapfile == 0)  {  #M3NO3
   
   print "No map file provided for atoms need new charges.\n"  if $debug >= 0;
   print "Charges will be taken from $qinfile, and assigned to QM atoms in $gfile in sequence. \n\n"  if $debug >= 0;
   
# open Gaussian ONIOM Input file for generating new oniom input file

   print "Open Gaussian ONIOM input file $gfile again......\n\n"  if $debug >= 1;
   open(ONIOMINFILE,$gfile) || die "Could not read $gfile\n$!\n";
   
   print "Open new Gaussian ONIOM input file $outputfile for writing......\n\n"  if $debug >= 1;
   open(ONIOMOUTFILE,">$outputfile") || die "Could not write to $outputfile\n$!\n";
   
   print "Open file $chgcompfile for writing......\n\n"  if $debug >= 1;
   open(CHGCOMPFILE,">$chgcompfile") || die "Could not write to $chgcompfile\n$!\n";
   
   print CHGCOMPFILE "Compare old and new charges for following atoms\n";
   print CHGCOMPFILE "Old charges from $gfile. New charges from $qinfile\n";
#   printf CHGCOMPFILE "%12s%12s%12s%12s%12s%12s%12s%12s\n","Atom","Type","OldCharge","NewCharge","Layer","X","Y","Z";
   printf CHGCOMPFILE "%5s%5s%15s%15s%10s%10s%11s%7s%15s%15s%15s",
                      "Atom","Type","OldCharge","NewCharge","SerialNum","LineNum","ChargeNum","Layer","X","Y","Z";
   print CHGCOMPFILE "\n";
   
   for ($i=0;$i < $structurestart; $i++) {  
     $_= <ONIOMINFILE>;
     print ONIOMOUTFILE;
     
     }  
   $RESPchgindex=0;
   for ($i=0;$i < $atomindex; $i++) {  #M3NO4
      
      $_ = <ONIOMINFILE>;

      if ( ($oniomresp->{AtomLayerGfile}[$i] eq 'H') ||
           ($oniomresp->{AtomLayerGfile}[$i] eq 'h')  ) { #M3NO5
           
           print "-"x66,"\n"  if $debug >= 4;

           print $oniomresp->{Original_Line}[$i] if $debug >= 4;

           $currnewchg=$RESPchg[$RESPchgindex];
           print " New charge is $currnewchg\n" if $debug >= 5;
           $RESPchgindex++;
           
           $currnewline = $oniomresp->{Original_Line}[$i];
           
           if ( $currnewchg == 0.000000 )  {
               print "\nWarning: zero charge 0.000000 will be used for atom ",$i+1,".\n" if $debug >= 0;
               print " Original line:\n ".$currnewline if $debug >= 0;
               print " located at line ",$structurestart+$i+1 if $debug >= 0;
               print " in $gfile.\n\n"  if $debug >= 0;
             }

           $currnewline =~ s/-(-*\d*\.\d+)/-$currnewchg/;
           print $currnewline if $debug >= 4;
           print ONIOMOUTFILE $currnewline;
           printf CHGCOMPFILE  "%5s%5s%15.7f%15.7f%10d%10d%11d%7s%15.7f%15.7f%15.7f\n",
                 $oniomresp->{ATOMSGfile}[$i],
                 $oniomresp->{ATOMTYPEGfile}[$i],
                 $oniomresp->{PARTIALCHARGGfile}[$i],
                 $currnewchg,
                 $i+1,
                 $oniomresp->{LineNumGfile}[$i],
                 $RESPchgindex,
                 $oniomresp->{AtomLayerGfile}[$i],
                 $oniomresp->{CARTCOORDGfile}[3*$i],
                 $oniomresp->{CARTCOORDGfile}[3*$i+1],
                 $oniomresp->{CARTCOORDGfile}[3*$i+2];
      
          } else {
           print ONIOMOUTFILE $oniomresp->{Original_Line}[$i];
          } #M3NO5
   
   
     }  #M3NO4 
    
    print "-"x66,"\n"  if $debug >= 4;

    while (<ONIOMINFILE>) { print ONIOMOUTFILE;}
   
 
 
 
   }  #M3NO3
  
 if ($havemapfile == 1)  {  #M3NO6

   print "Map file $mapfile was provided for atoms need new charges.\n"  if $debug >= 0;
   print "Charges will be taken from $qinfile, and assigned to atoms in $gfile using this file. \n\n"  if $debug >= 0;
   
   print "\nOpen map file $mapfile ......\n\n"  if $debug >= 0;
   open(MAPFILE,$mapfile) || die "Could not read $mapfile\n$!\n";
   
   # Read in mapping relationships from map file
   
   $mapindex = 0;
   $skipindicator = 0;
   @map = ();
   
   while (<MAPFILE>) { #M3NO7
#     print $_;

     next if /^\s*$/;
     
     chop;
     $currmapline = $_;
     
     # Treat index represented separatedly or connects by ...
     #______________
     #  3425
     #  3427
     #  3459
     #______________
     # will be read separatedly
     #______________
     #  3426
     #  ...
     #  3440
     #______________
     #  3426 through 3440 will be read
     #  
     #  or
     #______________
     #  3440
     #  ...
     #  3426
     #______________
     #  3440 through 3426 will be read
     #

     if (/^\.\.\./) {
       print " Set \$skipindicator as 1\n"  if $debug >= 4;
       $skipindicator = 1;
       } elsif ($skipindicator == 0) {
         print " Set $mapindex element in map array as $currmapline\n"  if $debug >= 4;
         $map[$mapindex] = $currmapline;
         $oniomresp->{NeedNewRESPGfile}[$currmapline-1] = 1;
         $oniomresp->{RESPChargeIndexGfile}[$currmapline-1] = $mapindex;
         $mapindex++;   
       } else {
         $skipindicator = 0;
         $skipfrom = $map[$mapindex-1];
         
         if ( $currmapline < $skipfrom ) {
             print "\n Attention: map index $skipfrom larger than $currmapline, will run from $skipfrom down to $currmapline.\n\n"  if $debug >= 0;

             for ($i = $skipfrom-1; $i >= $currmapline; $i--) {
                  print " Set $mapindex element in map array as $i\n" if $debug >= 4;
                  $map[$mapindex] = $i;
                  $oniomresp->{NeedNewRESPGfile}[$i-1] = 1;
                  $oniomresp->{RESPChargeIndexGfile}[$i-1] = $mapindex;
                  $mapindex++;   
                }
                
            } elsif ( $currmapline > $skipfrom) {
                     for ($i = $skipfrom+1; $i <= $currmapline; $i++) {
                       print " Set $mapindex element in map array as $i\n" if $debug >= 4;
                       $map[$mapindex] = $i;
                       $oniomresp->{NeedNewRESPGfile}[$i-1] = 1;
                       $oniomresp->{RESPChargeIndexGfile}[$i-1] = $mapindex;
                       $mapindex++;   
                      }
            } else {
              print "\n Warning: two numbers before and after ... are the same as $skipfrom.\n" if $debug >= 0;
            
            }
       }
    } #M3NO7
   
   if ($debug >= 3) {
     print " Map index read from $mapfile:\n";
     print " Mapping atom index in charge file $qinfile and ONIOM file $gfile.\n";
     for ($i = 0; $i < $mapindex; $i++) {
        printf "%5d%10d\n", $i+1,$map[$i];
      }
     print "\n There are $mapindex relationship read.\n\n";
    }
   
   # Sanity check: compare total charge number from RESP charge file and total atom numbers in map file
   if ( $totalchgnum < $mapindex ) {
      print " Warning: there are $mapindex atoms in map file $mapfile, but only $totalchgnum charges in charge file $qinfile.\n";
      print " Don't know where to find missing ",$mapindex - $totalchgnum," charges.\n";
      die " Exit.\n\n$!";
    
    }
   
   # After obtaining map relationship, now ready to extact new charges.
   # open Gaussian ONIOM Input file for generating new oniom input file

   print "Open Gaussian ONIOM input file $gfile again......\n\n"  if $debug >= 1;
   open(ONIOMINFILE,$gfile) || die "Could not read $gfile\n$!\n";
   
   print "Open new Gaussian ONIOM input file $outputfile for writing......\n\n"  if $debug >= 1;
   open(ONIOMOUTFILE,">$outputfile") || die "Could not write to $outputfile\n$!\n";
   
   print "Open file $chgcompfile for writing......\n\n"  if $debug >= 1;
   open(CHGCOMPFILE,">$chgcompfile") || die "Could not write to $chgcompfile\n$!\n";
   
   print CHGCOMPFILE "Compare old and new charges for following atoms\n";
   print CHGCOMPFILE "Old charges from $gfile. New charges from $qinfile\n";
   printf CHGCOMPFILE "%5s%5s%15s%15s%10s%10s%11s%7s%15s%15s%15s",
                      "Atom","Type","OldCharge","NewCharge","SerialNum","LineNum","ChargeNum","Layer","X","Y","Z";
   print CHGCOMPFILE "\n";
   
   for ($i=0;$i < $structurestart; $i++) {  
     $_= <ONIOMINFILE>;
     print ONIOMOUTFILE;
     
     }  
     
   
   # Run through all atoms
   
   for ($i = 0; $i < $atomindex; $i++) {  #M3NO8
      
      $_ = <ONIOMINFILE>;
       
      if ( $oniomresp->{NeedNewRESPGfile}[$i] == 0 ) {

          print ONIOMOUTFILE $oniomresp->{Original_Line}[$i];

       } else { #M3NO9
           print "-"x66,"\n"  if $debug >= 4;

           print $oniomresp->{Original_Line}[$i] if $debug >= 4;

           $currnewchg = $RESPchg[ $oniomresp->{RESPChargeIndexGfile}[$i] ];
           print " New charge is $currnewchg\n" if $debug >= 5;
           
           $currnewline = $oniomresp->{Original_Line}[$i];
           
           if ( $currnewchg == 0.000000 )  {
               print "\nWarning: zero charge 0.000000 will be used for atom ",$i+1,".\n" if $debug >= 0;
               print " Original line:\n ".$currnewline if $debug >= 0;
               print " located at line ",$structurestart+$i+1 if $debug >= 0;
               print " in $gfile.\n\n"  if $debug >= 0;
             }

           $currnewline =~ s/-(-*\d*\.\d+)/-$currnewchg/;
           print $currnewline if $debug >= 4;
           print ONIOMOUTFILE $currnewline;
           printf CHGCOMPFILE  "%5s%5s%15.7f%15.7f%10d%10d%11d%7s%15.7f%15.7f%15.7f\n",
                 $oniomresp->{ATOMSGfile}[$i],
                 $oniomresp->{ATOMTYPEGfile}[$i],
                 $oniomresp->{PARTIALCHARGGfile}[$i],
                 $currnewchg,
                 $i+1,
                 $oniomresp->{LineNumGfile}[$i],
                 $oniomresp->{RESPChargeIndexGfile}[$i]+1,
                 $oniomresp->{AtomLayerGfile}[$i],
                 $oniomresp->{CARTCOORDGfile}[3*$i],
                 $oniomresp->{CARTCOORDGfile}[3*$i+1],
                 $oniomresp->{CARTCOORDGfile}[3*$i+2];
       
       } #M3NO9
      
   
     }  #M3NO8
    
    print "-"x66,"\n"  if $debug >= 4;
    
    # Finishing printing rest part of ONIOM input file.
    while (<ONIOMINFILE>) { print ONIOMOUTFILE;}
   
 
   
   }  #M3NO6
   

   
 } # sub modethree {


1;

__END__

=back

=head1 VERSION

0.1

=head1 SEE ALSO

F<ESPT::ONIOMinput>

=head1 AUTHOR

Peng Tao, E<lt>tao.21@osu.eduE<gt>

=head1 COPYRIGHT

Copyright (c) 2009~2010 by Peng Tao


=cut
