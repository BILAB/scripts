package ESPT::SETMVFLG;

use base qw(ESPT::ESSfile);
use ESPT::Glib 0.01;
use strict;
use warnings;

our $VERSION = '0.01';

=head1 NAME

ESPT::SETMVFLG - Oniom files object to reset moving flag in a Gaussian ONIOM input file.

=head1 SYNOPSIS

   use ESPT::SETMVFLG;

   my $oniom = SETMVFLG->new();

=head1 DESCRIPTION

This module provides a method to reset moving flag in a Gaussian ONIOM input

=cut

### Version History ###
# 0.01	Reset moving flag


### To Do ###


## the object constructor **

sub new {
	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my $oniom = ESPT::ESSfile->new();

	$oniom->{TYPE} = "oniom";

	# program info
	$oniom->{PROGRAM} = "Gaussian";
	$oniom->{VERSION} = undef;
	$oniom->{REVISION} = undef;
	$oniom->{COMPILE} = undef;
	$oniom->{COMPILEDATE} = undef;
	
	# ONIOM info
	$oniom->{LOWMODEL} = undef;
	$oniom->{HIGHMODEL} = undef;
	$oniom->{LOWREAL} = undef;
	$oniom->{LOWREAL_MODEL} = undef;

        # Operation info
        $oniom->{GJFTEMPLATE} = undef;
	$oniom->{OPT} = 0;
	$oniom->{STRUCTURE} = 0;
        $oniom->{STRUCTUREFILE} = undef;
        $oniom->{MOVEFLAG} = 0;
        $oniom->{VERBOSE} = 0;
        $oniom->{TEMPLATEINPUT} = 0;
        
        
        $oniom->{UseONBCoord} = undef;
        
        $oniom->{ModelGjfFile}  = undef;
        $oniom->{ModelONBFile}  = undef;
        $oniom->{ModelNewInput} = undef;


	# Model system input file info
        $oniom->{ModelGjfFile} = undef;
	$oniom->{ATOMSModelGjf} = [];    # Element symbol from input file
        $oniom->{RESIDUEModelGjf} = [];
        $oniom->{RESIDUENUMModelGjf} = [];
        $oniom->{ATOMTYPEModelGjf} = [];
        $oniom->{ELEMENTNUMModelGjf} = [];
        $oniom->{ELEMENTModelGjf} = [];
        $oniom->{PARTIALCHARGModelGjf} = [];
        $oniom->{CARTCOORDModelGjf} = [];
        $oniom->{ATOMDEFINEModelGjf} = [];
        $oniom->{ATOMFREEZEFLAGModelGjf} = [];
        $oniom->{AtomLayerModelGjf} = [];
        $oniom->{ATOMINRESIDUEModelGjf} = [];
        $oniom->{ORIGINALLINEModelGjf} = [];
        $oniom->{LineNumModelGjf} = [];
        $oniom->{ATOMNUMModelGjf} = 0;

        $oniom->{ATOMNewFREEZEFLAG} = [];

	# Model system ONB file info
        $oniom->{ModelONBFile} = undef;
	$oniom->{ATOMSModelONB} = []; # Element symbol from input file
        $oniom->{RESIDUEModelONB} = [];
        $oniom->{RESIDUENUMModelONB} = [];
        $oniom->{ATOMTYPEModelONB} = [];
        $oniom->{ELEMENTNUMModelONB} = [];
        $oniom->{ELEMENTModelONB} = [];
        $oniom->{PARTIALCHARGModelONB} = [];
        $oniom->{CARTCOORDModelONB} = [];
        $oniom->{ATOMDEFINEModelONB} = [];
        $oniom->{ATOMFREEZEFLAGModelONB} = [];
        $oniom->{ATOMINRESIDUEModelONB} = [];
        $oniom->{ORIGINALLINEModelONB} = [];
        $oniom->{AtomLayerModelONB} = [];
        $oniom->{AtomNameModelONB} = [];
        
        $oniom->{ATOMINRESIDUE} = [];
        $oniom->{RESIDUEWITHINCORE} = [];
        
   	# These variables are related to core residues set up
        $oniom->{CORERESIDUE} = [];
        $oniom->{CORERESIDUENAME} = [];
        $oniom->{CORERESIDUENO} = undef;
        $oniom->{RESIDFILE} = undef;
        $oniom->{NEAR} = 6.0;
        $oniom->{RESIDUENEARCORE} = undef;
     
        # These variables are list of core atoms
        $oniom->{COREATOMSID} = [];
        $oniom->{COREATOMSTOT} = undef;
        
        # Variables for labeling the residues to assist freeze flag setup
        $oniom->{RESIDUEIDLIST} = [];
        $oniom->{RESIDUENAMELIST} = [];
        $oniom->{RESIDUEWITHINCORE} = [];
        $oniom->{TOTALRESIDUENUM} = undef;
        
        $oniom->{CARTCOORD4SetFlag} = [];

        # Summary for flag setup
        $oniom->{TotalMovingAtomsBefore} = 0;
        $oniom->{TotalMovingResiduesBefore} = 0;
        $oniom->{TotalMovingAtomsAfter} = 0;
        $oniom->{TotalMovingResiduesAfter} = 0;

        
        $oniom->{ModelNewInput} = undef;
        $oniom->{DefaultModelNewInput} = "ONIOM_NEW_Moving_Flag.gjf";

	$oniom->{GJFINPUT} = undef;
        $oniom->{GJFOUTPUT} = undef;
        
	bless($oniom, $class);
	return $oniom;
}

# set filename & spin then digest the file
sub analyze : method {
	my $oniom = shift;
	$oniom->digest();
	return;
}

## subroutines ##

sub digest {

my $oniom = shift;
my $debug = $oniom->{DEBUG};
my $modelgjf = $oniom->{ModelGjfFile};
my $modelonb = $oniom->{ModelONBFile};
my $modelnewgjf = $oniom->{ModelNewInput};
my $useONBcoord = $oniom->{UseONBCoord};
my $near = $oniom->{NEAR};
my $residfile = $oniom->{RESIDFILE};
my ($atomindex, $resindex);
my $residref = -999999;
my $resnameref = "";
my ($curresidname, $curresidnum);
my $structindicator;
my $lineindex;
my $i;
my $j;
my $distance;
my $mismatch;
my ($coreresindex, $coreatomindex, $curcoreresid, $curcoreresname);

# Read in Model GJF input file

print "\nReading ONIOM input file $modelgjf...\n" if $debug >= 0;

open(MODELGJFFILE,$modelgjf) || die "Could not read $modelgjf\n$!\n";

$structindicator = 0;
$atomindex = 0;
$lineindex = 0;

while (<MODELGJFFILE>) {  #NO1
   $lineindex++;
   if   ($structindicator == 0){
     if (/^\s*-*\d+\s+\d+\s+-*\d+\s+\d+/) {

       print $_ if $debug >= 1;
       $structindicator = 1;
       next;
     } # if (/^\s*-*\d+\s+\d+\s+-*\d+\s+\d+/) 
     next;
   }  # if   ($structindicator ==0){
   else {
     if (/^\s+$/ or /^$/){ last;}
#     print $_;
     
     $oniom->{ORIGINALLINEModelGjf}[$atomindex] = $_;
#  Read in each atom information
     /^\s*([a-zA-Z]+)-([a-zA-Z0-9\*]+)-(-*\d*\.\d+)\s+(-*[01])\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+([HMLhml])/;
      
      print "$1\t$2\t$3\t$4\t$5\t$6\t$7\t$8\t$atomindex\n" if $debug >= 3; 
      $oniom->{ATOMSModelGjf}[$atomindex] = $1;
      $oniom->{ATOMTYPEModelGjf}[$atomindex] = $2;
      $oniom->{PARTIALCHARGModelGjf}[$atomindex] = $3;
      $oniom->{ATOMFREEZEFLAGModelGjf}[$atomindex] = $4;
      $oniom->{CARTCOORDModelGjf}[3*$atomindex] = $5;
      $oniom->{CARTCOORDModelGjf}[3*$atomindex+1] = $6;
      $oniom->{CARTCOORDModelGjf}[3*$atomindex+2] = $7;
      $oniom->{AtomLayerModelGjf}[$atomindex] = $8;
      $oniom->{MovedModelGjf}[$atomindex] = 0;
      $oniom->{LineNumModelGjf}[$atomindex] = $lineindex;
      
      $oniom->{ATOMNewFREEZEFLAG}[$atomindex] = -1;
      
      # If -b is not set, use coordinates from ONIOM job file to set optimization flag
      if ($useONBcoord == 0) {
        $oniom->{CARTCOORD4SetFlag}[3*$atomindex]   = $5;
        $oniom->{CARTCOORD4SetFlag}[3*$atomindex+1] = $6;
        $oniom->{CARTCOORD4SetFlag}[3*$atomindex+2] = $7;
       }
      
      $atomindex++;
   
   } # else

 } # while (<MODELGJFFILE>) { #NO1

$oniom->{ATOMNUMModelGjf} = $atomindex;

print "There are $atomindex atoms in file $modelgjf. \n" if $debug >= 0; 
print "Done!\n\n" if $debug >= 0;

close (MODELGJFFILE);



# Read in Model ONB input file

print "Reading ONIOM ONB file $modelonb...\n" if $debug >= 0;

open(MODELONBFILE,$modelonb) || die "Could not read $modelonb\n$!\n";

# set atom and residues counters as 0
$atomindex = 0;
$resindex = -1;

while (<MODELONBFILE>) { 

     if (/^\s+$/ or /^$/){ last;}
#     print $_;
     
     $oniom->{ORIGINALLINEModelONB}[$atomindex] = $_;
#  Read in each atom information
     /^\s*([a-zA-Z]+)-([a-zA-Z0-9\*]+)-(-*\d*\.\d+)\s+(-*[01])\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+([HML])\s+(\S+)\s+(\S+)\s+(\d+)/;
      
      print "$1\t$2\t$3\t$4\t$5\t$6\t$7\t$8\t$9\t$10\t$11\t$atomindex\n" if $debug >= 3; 
      $oniom->{ATOMSModelONB}[$atomindex] = $1;
      $oniom->{ATOMTYPEModelONB}[$atomindex] = $2;
      $oniom->{PARTIALCHARGModelONB}[$atomindex] = $3;
      $oniom->{ATOMFREEZEFLAGModelONB}[$atomindex] = $4;
      $oniom->{CARTCOORDModelONB}[3*$atomindex] = $5;
      $oniom->{CARTCOORDModelONB}[3*$atomindex+1] = $6;
      $oniom->{CARTCOORDModelONB}[3*$atomindex+2] = $7;
      $oniom->{AtomLayerModelONB}[$atomindex] = $8;
      $oniom->{AtomNameModelONB}[$atomindex] = $9;
      $oniom->{RESIDUEModelONB}[$atomindex] = $10;
      $oniom->{RESIDUENUMModelONB}[$atomindex] = $11;
      $oniom->{ATOMINRESIDUE}[$atomindex] = $resindex;
      
      $curresidname = $oniom->{RESIDUEModelONB}[$atomindex];
      $curresidnum = $oniom->{RESIDUENUMModelONB}[$atomindex];


      # If -b is set, use coordinates from ONB file to set optimization flag
      if ($useONBcoord == 1) {
        $oniom->{CARTCOORD4SetFlag}[3*$atomindex]   = $5;
        $oniom->{CARTCOORD4SetFlag}[3*$atomindex+1] = $6;
        $oniom->{CARTCOORD4SetFlag}[3*$atomindex+2] = $7;
       }

      $atomindex++;
     #Check if encounter a new residue, if yes, store residue number and name. 
     
     if ( ($residref eq $curresidnum) and ( $resnameref eq $curresidname) )
      { next; }
     
     $resindex++;
     $oniom->{RESIDUEIDLIST}[$resindex] = $curresidnum;
     $oniom->{RESIDUENAMELIST}[$resindex] = $curresidname;
     $oniom->{RESIDUEWITHINCORE}[$resindex] = 0;
     $residref = $curresidnum;
     $resnameref = $curresidname;
     $oniom->{ATOMINRESIDUE}[$atomindex-1] = $resindex;
     print "Encounter new residue $resindex, name $curresidname, index $curresidnum.\n" if $debug >= 2;

 } # while (<MODELGJFFILE>) {

if ($oniom->{ATOMNUMModelGjf} != $atomindex)  {
   print "There are $atomindex atoms in $modelonb, but ",$oniom->{ATOMNUMModelGjf}," atoms in $modelgjf!\n";
   print "Please check!\n";
   die;
  }


print "There are $atomindex atoms in file $modelonb. \n" if $debug >= 0; 

$resindex++;
print "There are $resindex residues in file $modelonb.\n" if $debug >= 0;
$oniom->{TOTALRESIDUENUM} =$resindex; 

close (MODELONBFILE);

print "Done!\n\n" if $debug >= 0;


# Sanity check for Model system

print "Sanity check for Model system...\n" if $debug >= 1;

for ($i=0; $i < $oniom->{ATOMNUMModelGjf}; $i++) {
  
  if ($oniom->{ATOMSModelGjf}[$i] ne $oniom->{ATOMSModelONB}[$i]) {
    print "Atom ",$i+1," line ",$oniom->{LineNumModelGjf}[$i],
        " in file $modelgjf is different \n  between $modelgjf and $modelonb.\n" if $debug >= 1;
    print " Original line in $modelgjf is\n",$oniom->{ORIGINALLINEModelGjf}[$i] if $debug >= 1;
    print " Original line in $modelonb is\n",$oniom->{ORIGINALLINEModelONB}[$i] if $debug >= 1;
    print "\n\n" if $debug >= 1;
  
   }
   

 } # for ($i=0; $i < $oniom->{ATOMNUMModelGjf}; $i++)

#print "Done!\n\n" if $debug >= 0;



# Read in core residues information from core residues list file

 print "Open $residfile file for process.\n" if $debug >= 1;
 open(CORERESFILE,$residfile) || die "Could not read $residfile\n$!\n";

 $coreresindex = 0;

 while (<CORERESFILE>) {
 
  next if /^$/;
  /\[(.+)\]\s+"(.+)"/;
  $oniom->{CORERESIDUENAME}[$coreresindex] = $1;
  $oniom->{CORERESIDUE}[$coreresindex] = $2;
  
  $curcoreresid = $oniom->{CORERESIDUE}[$coreresindex];
  $curcoreresname = $oniom->{CORERESIDUENAME}[$coreresindex];
  $coreresindex++;

  print "Read in residue number $coreresindex with residue name as ($curcoreresname) residue number as ($curcoreresid) \n" if $debug >= 1;
    
  } # while (<CORERESFILE>) {


  $oniom->{CORERESIDUENO} = $coreresindex;
  
  if ( $coreresindex == 0) {
   print "There is no residue number found in file $residfile, please double check.\n";
   die;
   }
 
  print "There are $coreresindex residues in the core region. \n" if $debug >= 1;

 
 # Find list of atoms which belong to core region
 $coreatomindex = 0;
 for ($i=0; $i < $oniom->{ATOMNUMModelGjf};$i++) {
   
   if ( $oniom->atombelong2core($i) == 0) {next;}
    
   $oniom->{COREATOMSID}[$coreatomindex] = $i;
   $coreatomindex++;
   print "Atom $i (in residue $oniom->{RESIDUEModelONB}[$i] $oniom->{RESIDUENUMModelONB}[$i]) belongs to core region.\n" if $debug >= 2;
 
  } # for (my $i=0; $i < $pdb->{NATOMS};$i++)
 
 $oniom->{COREATOMSTOT} = $coreatomindex;
 print "There are total $coreatomindex atoms in the core region.\n" if $debug >= 1;

 $oniom->decideflag();


} # sub digest 



#Subroutine that decides the freeze flag of each atom according to its distance
# to core residues.

sub decideflag {
 my $oniom = shift;
 my $curdistance;
 my $curresidue;
 my $debug = $oniom->{DEBUG};
 my $totresiduewithin = 0;
 my $frzbackboneatoms = 0;
 
 
 # go through all the atoms to check if each residue are within certain distance of core or not.
 for (my $i=0; $i < $oniom->{ATOMNUMModelGjf};$i++) {
   
#   print "Process (",$i,") (",$oniom->{ATOMINRESIDUE}[$i],") (",$oniom->{RESIDUEWITHINCORE}[$oniom->{ATOMINRESIDUE}[$i]],")\n";
   if ($oniom->{RESIDUEWITHINCORE}[$oniom->{ATOMINRESIDUE}[$i]] == 1) {next;}
   
   if ( $oniom->atomwithincore($i) == 1) {
     $oniom->{RESIDUEWITHINCORE}[$oniom->{ATOMINRESIDUE}[$i]] = 1;
     $curresidue = $oniom->{ATOMINRESIDUE}[$i];
     $totresiduewithin++;
     print "Residue ($curresidue) is close to core region.\n" if $debug >= 2;
     
    }
   
   
  } # for (my $i=0; $i < $oniom->{ATOMINRESIDUE};$i++)

 print "There are $totresiduewithin residues close to core region (including those residues in the core region). \n" if $debug >= 1;

 $oniom->{RESIDUENEARCORE} = $totresiduewithin;

 # go through all the atoms to setup freeze flag.
 for (my $i=0; $i < $oniom->{ATOMNUMModelGjf};$i++) {
  
#   if ($oniom->{RESIDUEWITHINCORE}[$oniom->{ATOMINRESIDUE}[$i]] == 0) {next;}
   
   # If this atom belongs to core region, it is free to move.
   if ( $oniom->atombelong2core($i) == 1 ) {
    $oniom->{ATOMNewFREEZEFLAG}[$i] = 0;
    next;
    }
   
   # If this atom belongs to a residue that is within core, it is free to move.
   
   if ($oniom->{RESIDUEWITHINCORE}[$oniom->{ATOMINRESIDUE}[$i]] == 1){
      $oniom->{ATOMNewFREEZEFLAG}[$i] = 0;
      next;
    }
   
  }# for (my $i=0; $i < $oniom->{ATOMNUMModelGjf};$i++)
  

} # sub decideflag 


# Subroutine atomwithincore: Check if one atom is within certain distance of core or not.

sub atomwithincore{

 my $oniom = shift;
 my $atomid = shift;
 my $curcoreatom;
 my $distance;
 my $debug = $oniom->{DEBUG};
 
 for (my $i=0; $i < $oniom->{COREATOMSTOT};$i++) {
   $curcoreatom = $oniom->{COREATOMSID}[$i];
   $distance = $oniom->atomdistance($atomid,$curcoreatom);
   if ( $distance < $oniom->{NEAR}) {
   
      print "Atom $atomid is at least $distance A to core region. \n" if  $debug >= 2;
      return 1;
     } 
   
 
  } # for (my $i=0; $i < $oniom->{COREATOMSTOT};$i++)

 print "Atom $atomid is NOT close to core region. \n" if  $debug >= 2;

 return 0;
 
} # sub atomwithincore


# Subroutine atombelong2core: check if given atom belong to core residues or not
sub atombelong2core{

 my $oniom = shift;
 my $atomid = shift;
 my $debug = $oniom->{DEBUG};
 
 my $oriname = $oniom->{RESIDUEModelONB}[$atomid];
 my $oriid = $oniom->{RESIDUENUMModelONB}[$atomid];
 my $curname;
 my $curid;
 
 
 for (my $j=0; $j < $oniom->{CORERESIDUENO}; $j++) {
   $curname  = $oniom->{CORERESIDUENAME}[$j];
   $curid    = $oniom->{CORERESIDUE}[$j];
   
   print "Compare (".$oriname.") with (".$curname.") and (".
                    $oriid.") with (".$curid.")\n" if $debug >= 2;

   if ( (($oriname =~ /$curname/) or ($curname =~ /$oriname/)) and 
        (($oriid =~ /$curid/) and ($curid =~ /$oriid/)) ){ 
          print "Match (".$oriname.") with (".$curname.") and (".
                    $oriid.") with (".$curid.")\n" if $debug >= 2;
          return 1;
        }
  }
 
 return 0;

} # sub atombelong2core



#Subroutine AtomDistance: calculate distance between two atoms
sub atomdistance{

 my $oniom = shift;
 my $atomid1 = shift;
 my $atomid2 = shift;
 my $distance;
 
 $distance = sqrt( ( $oniom->{CARTCOORD4SetFlag}[3*$atomid1] - $oniom->{CARTCOORD4SetFlag}[3*$atomid2] )**2 + 
                   ( $oniom->{CARTCOORD4SetFlag}[3*$atomid1+1] - $oniom->{CARTCOORD4SetFlag}[3*$atomid2+1] )**2 + 
                   ( $oniom->{CARTCOORD4SetFlag}[3*$atomid1+2] - $oniom->{CARTCOORD4SetFlag}[3*$atomid2+2] )**2  );
 return $distance;


} # sub atomdistance


#Subroutine that summarizes optimization flag set before and after.

sub flagsetsummary {
 my $oniom = shift;
 my $debug = $oniom->{DEBUG};
 my $totMvAtmB4   = 0;
 my $totMvResB4   = 0;
 my $totMvAtmAft  = 0;
 my $totMvResAft  = 0;
 my $residref = -999999;
 my $resnameref = "";
 my $resindex = 0;
 my ($curresidname, $curresidnum);
 my $i;
 
  for ($i=0; $i < $oniom->{ATOMNUMModelGjf};$i++) {

   $curresidname = $oniom->{RESIDUEModelONB}[$i];
   $curresidnum = $oniom->{RESIDUENUMModelONB}[$i];
   
     if ($oniom->{ATOMFREEZEFLAGModelGjf}[$i] == 0 ) {$totMvAtmB4++;}
     if ($oniom->{ATOMNewFREEZEFLAG}[$i] == 0 ) {$totMvAtmAft++;}
     
     if ( ($residref eq $curresidnum) and ( $resnameref eq $curresidname) )
      { next; }
     
     $resindex++;
     $residref = $curresidnum;
     $resnameref = $curresidname;
     if ($oniom->{ATOMFREEZEFLAGModelGjf}[$i] == 0 ) {$totMvResB4++;}
     if ($oniom->{ATOMNewFREEZEFLAG}[$i] == 0 ) {$totMvResAft++;}
     
   
  } # for (my $i=0; $i < $oniom->{ATOMINRESIDUE};$i++)

  $oniom->{TotalMovingAtomsBefore} = $totMvAtmB4;
  $oniom->{TotalMovingResiduesBefore} = $totMvResB4;
  $oniom->{TotalMovingAtomsAfter} =$totMvAtmAft;
  $oniom->{TotalMovingResiduesAfter} = $totMvResAft;

print "\nIn file $oniom->{ModelGjfFile}, there are $totMvAtmB4 atoms, $totMvResB4 residues subject to geometry optimization.\n"  if $debug >= 0;

print "\nAfter reset, there are $totMvAtmAft atoms, $totMvResAft residues subject to geometry optimization.\n\n"  if $debug >= 0;


} # sub flagsetsummary {


#Subroutine WriteModelONIOM: write new model ONIOM input file with new moving flag

sub writemodeloniom{

 my $oniom = shift;
 my $debug = $oniom->{DEBUG};

 my $modelnewgjf      = $oniom->{ModelNewInput}; 
 my $i;
 my $currline;
 my ($atom,$atomtyp,$partialcharge,$moveflag,$layer,$x,$y,$z);
 my $modelindex;

 
 open(INPUTFILE,$oniom->{ModelGjfFile}) || die "Could not read $oniom->{ModelGjfFile}\n$!\n";
 
 print "Open $modelnewgjf to write new ONIOM input file...\n" if $debug >= 0;

 open(OUTPUTGJFFILE,">$modelnewgjf") || die "Could not open $modelnewgjf\n$!\n";
 
 while (<INPUTFILE>) {  
     
     # Keep everything before atom list

     print OUTPUTGJFFILE $_ ;

     if (/^\s*-*\d+\s+\d+\s+-*\d+\s+\d+/) { last;} 
     
  } # while (<INPUTFILE>) { 

 for ($i=0; $i< $oniom->{ATOMNUMModelGjf}; $i++) {
   
   $currline = <INPUTFILE>;
   
   # Use original file if do not need to be replaced;
   
   if ( $oniom->{ATOMFREEZEFLAGModelGjf}[$i] == $oniom->{ATOMNewFREEZEFLAG}[$i]){
     print OUTPUTGJFFILE $currline ;
     next;
    } 
   
   $moveflag = $oniom->{ATOMNewFREEZEFLAG}[$i];
   
   print "\nAtom ", $i+1," ",$oniom->{ATOMSModelGjf}[$i], " in residue ", $oniom->{RESIDUEModelONB}[$i]," change flag from ", $oniom->{ATOMFREEZEFLAGModelGjf}[$i],
         " to $moveflag.\n" if $debug >=1;
   
   print "Before the line is \n",$currline if $debug >=2;
   
   if ($moveflag == 0) {
      $currline =~ s/\s(-*[01])\s/  $moveflag /;
     } else {
      $currline =~ s/\s(-*[01])\s/ $moveflag /;
     }

   print "After substitution the line is \n",$currline if $debug >=2;
#   print $atom,$atomtyp, $partialcharge, $moveflag, $x, $y, $z, $layer,"\n" if $debug >=1;
 
   print OUTPUTGJFFILE $currline ;
 
 
  } # for ($i=0; $i< $oniom->{ATOMNUMProdGjf}; $i++)
 
 while (<INPUTFILE>) {  
 
     print OUTPUTGJFFILE $_ ;

 
  } # while (<INPUTFILE>) { 
 

} #sub writeprodoniom{




1;

__END__

=back

=head1 VERSION

0.01

=head1 SEE ALSO

F<ESPT::PDB>

=head1 AUTHOR

Peng Tao, E<lt>tao.21@osu.eduE<gt>

=head1 COPYRIGHT

Copyright (c) 2009~2010 by Peng Tao

=cut
