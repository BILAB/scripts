package ESPT::ONIOMGEOM;

use base qw(ESPT::ESSfile);
use ESPT::Glib 0.01;
use strict;
use warnings;

our $VERSION = '0.02';

=head1 NAME

ESPT::OniomGeom - Oniom input and template files object to create new ONIOM gaussian input file.

=head1 SYNOPSIS

   use ESPT::OniomGeom;

   my $oniom = OniomGeom->new();

=head1 DESCRIPTION

This module provides methods to put optimized geometry from model system to production run ONIOM input

=cut

### Version History ###
# 0.01	Use four template files to create new files


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
	$oniom->{ONIOMENERGY} = undef;
	$oniom->{LOWREAL_MODEL} = undef;
	$oniom->{DIPOLEX} = undef;
	$oniom->{DIPOLEY} = undef;
	$oniom->{DIPOLEZ} = undef;
	$oniom->{DIPOLE} = undef;

        # Operation info
        $oniom->{GJFTEMPLATE} = undef;
	$oniom->{OPT} = 0;
	$oniom->{STRUCTURE} = 0;
        $oniom->{STRUCTUREFILE} = undef;
        $oniom->{MOVEFLAG} = 0;
        $oniom->{VERBOSE} = 0;
        $oniom->{TEMPLATEINPUT} = 0;

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
        $oniom->{MovedModelGjf} = [];  # flag to label if this atom moved from its original place
        $oniom->{ProdIndexModelGjf} = []; # This atom's corresponding atom index in production system

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

	# Production system input file info
	$oniom->{ProdGjfFile} = undef;
	$oniom->{ATOMSProdGjf} = []; # Element symbol from input file
        $oniom->{RESIDUEProdGjf} = [];
        $oniom->{RESIDUENUMProdGjf} = [];
        $oniom->{ATOMTYPEProdGjf} = [];
        $oniom->{ELEMENTNUMProdGjf} = [];
        $oniom->{ELEMENTProdGjf} = [];
        $oniom->{PARTIALCHARGProdGjf} = [];
        $oniom->{CARTCOORDProdGjf} = [];
        $oniom->{ATOMDEFINEProdGjf} = [];
        $oniom->{ATOMFREEZEFLAGProdGjf} = [];
        $oniom->{ATOMINRESIDUEProdGjf} = [];
        $oniom->{ORIGINALLINEProdGjf} = [];
        $oniom->{AtomLayerProdGjf} = [];
        $oniom->{LineNumProdGjf} = [];
        $oniom->{ATOMNUMProdGjf} = 0;
        $oniom->{InModelProdGjf} = [];   # Flag to show if this atom is in model system or not
        $oniom->{NeedReplaceProdGjf} = []; # Flag if this atom needs to new coordinates or not
        $oniom->{ModelAtmIndexProdGjf} = []; # This atom's corresponding atom index in model system
        

	# Production system ONB file info
        $oniom->{ProdONBFile} = undef;
	$oniom->{ATOMSProdONB} = []; # Element symbol from input file
        $oniom->{RESIDUEProdONB} = [];
        $oniom->{RESIDUENUMProdONB} = [];
        $oniom->{ATOMTYPEProdONB} = [];
        $oniom->{ELEMENTNUMProdONB} = [];
        $oniom->{ELEMENTProdONB} = [];
        $oniom->{PARTIALCHARGProdONB} = [];
        $oniom->{CARTCOORDProdONB} = [];
        $oniom->{ATOMDEFINEProdONB} = [];
        $oniom->{ATOMFREEZEFLAGProdONB} = [];
        $oniom->{ATOMINRESIDUEProdONB} = [];
        $oniom->{ORIGINALLINEProdONB} = [];
        $oniom->{AtomLayerProdONB} = [];
        $oniom->{AtomNameProdONB} = [];
        
        $oniom->{ProdNewInput} = undef;
        $oniom->{DefaultProdNewInput} = "Prod_ONIOM_NEW.gjf";

	$oniom->{GJFINPUT} = undef;
        $oniom->{GJFOUTPUT} = undef;
        $oniom->{STRUCTURENUMBER} = 0;
        
        $oniom->{tolerance} = 0.0001;
	$oniom->{atomprod} = undef;
        $oniom->{atmtypprod}= undef;
        $oniom->{chgprod}  = undef;
        $oniom->{moveprod} = undef;
        $oniom->{layerprod}   = undef;  
        $oniom->{coordprod}   = undef;  
        
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
my $prodgjf = $oniom->{ProdGjfFile};
my $prodonb = $oniom->{ProdONBFile};
my $prodnewgjf = $oniom->{ProdNewInput};
my $atomindex;
my $structindicator;
my $lineindex;
my $tol;
my $i;
my $j;
my $distance;
my $mismatch;

# Read in Model GJF input file

$tol = $oniom->{tolerance};

print "Reading Model system ONIOM input file $modelgjf...\n" if $debug >= 0;

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
      $oniom->{ProdIndexModelGjf}[$atomindex] = 0;
      $oniom->{LineNumModelGjf}[$atomindex] = $lineindex;
      $atomindex++;
   
   } # else

 } # while (<MODELGJFFILE>) { #NO1

$oniom->{ATOMNUMModelGjf} = $atomindex;

print "There are $atomindex atoms in file $modelgjf. \n" if $debug >= 0; 
print "Done!\n" if $debug >= 0;

close (MODELGJFFILE);


# Read in Production GJF input file

print "Reading Production system ONIOM input file $prodgjf...\n" if $debug >= 0;

open(PRODGJFFILE,$prodgjf) || die "Could not read $prodgjf\n$!\n";

$structindicator = 0;
$atomindex = 0;
$lineindex = 0;

while (<PRODGJFFILE>) {  #NO1
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
     
     $oniom->{ORIGINALLINEProdGjf}[$atomindex] = $_;
#  Read in each atom information
     /^\s*([a-zA-Z]+)-([a-zA-Z0-9\*]+)-(-*\d*\.\d+)\s+(-*[01])\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+([HMLhml])/;
      
      print "$1\t$2\t$3\t$4\t$5\t$6\t$7\t$8\t$atomindex\n" if $debug >= 3; 
      $oniom->{ATOMSProdGjf}[$atomindex] = $1;
      $oniom->{ATOMTYPEProdGjf}[$atomindex] = $2;
      $oniom->{PARTIALCHARGProdGjf}[$atomindex] = $3;
      $oniom->{ATOMFREEZEFLAGProdGjf}[$atomindex] = $4;
      $oniom->{CARTCOORDProdGjf}[3*$atomindex] = $5;
      $oniom->{CARTCOORDProdGjf}[3*$atomindex+1] = $6;
      $oniom->{CARTCOORDProdGjf}[3*$atomindex+2] = $7;
      $oniom->{AtomLayerProdGjf}[$atomindex] = $8;
      $oniom->{LineNumProdGjf}[$atomindex] = $lineindex;
      $oniom->{InModelProdGjf}[$atomindex] = 0;
      $oniom->{NeedReplaceProdGjf}[$atomindex] = 0;
      $oniom->{ModelAtmIndexProdGjf}[$atomindex] = 0;
      $atomindex++;
   
   } # else

 } # while (<PRODGJFFILE>) { #NO1

$oniom->{ATOMNUMProdGjf} = $atomindex;

print "There are $atomindex atoms in file $prodgjf. \n" if $debug >= 0; 
print "Done!\n" if $debug >= 0;

close (PRODGJFFILE);



# Read in Model ONB input file

print "Reading Model system ONIOM ONB file $modelonb...\n" if $debug >= 0;

open(MODELONBFILE,$modelonb) || die "Could not read $modelonb\n$!\n";

$atomindex = 0;

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
      
      $atomindex++;
   

 } # while (<MODELGJFFILE>) {

if ($oniom->{ATOMNUMModelGjf} != $atomindex)  {
   print "There are $atomindex atoms in $modelonb, but ",$oniom->{ATOMNUMModelGjf}," atoms in $modelgjf!\n";
   print "Please check!\n";
   die;
  }


print "There are $atomindex atoms in file $modelonb. \n" if $debug >= 0; 
print "Done!\n" if $debug >= 0;

close (MODELONBFILE);


# Read in production ONB input file

print "Reading production system ONIOM ONB file $prodonb...\n" if $debug >= 0;

open(PRODONBFILE,$prodonb) || die "Could not read $prodonb\n$!\n";

$atomindex = 0;

while (<PRODONBFILE>) { 

     if (/^\s+$/ or /^$/){ last;}
#     print $_;
     
     $oniom->{ORIGINALLINEProdONB}[$atomindex] = $_;
#  Read in each atom information
     /^\s*([a-zA-Z]+)-([a-zA-Z0-9\*]+)-(-*\d*\.\d+)\s+(-*[01])\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+([HML])\s+(\S+)\s+(\S+)\s+(\d+)/;
      
      print "$1\t$2\t$3\t$4\t$5\t$6\t$7\t$8\t$9\t$10\t$11\t$atomindex\n" if $debug >= 3; 
      $oniom->{ATOMSProdONB}[$atomindex] = $1;
      $oniom->{ATOMTYPEProdONB}[$atomindex] = $2;
      $oniom->{PARTIALCHARGProdONB}[$atomindex] = $3;
      $oniom->{ATOMFREEZEFLAGProdONB}[$atomindex] = $4;
      $oniom->{CARTCOORDProdONB}[3*$atomindex] = $5;
      $oniom->{CARTCOORDProdONB}[3*$atomindex+1] = $6;
      $oniom->{CARTCOORDProdONB}[3*$atomindex+2] = $7;
      $oniom->{AtomLayerProdONB}[$atomindex] = $8;
      $oniom->{AtomNameProdONB}[$atomindex] = $9;
      $oniom->{RESIDUEProdONB}[$atomindex] = $10;
      $oniom->{RESIDUENUMProdONB}[$atomindex] = $11;
      
      $atomindex++;
   

 } # while (<PRODGJFFILE>) {

if ($oniom->{ATOMNUMProdGjf} != $atomindex)  {
   print "There are $atomindex atoms in $prodonb, but ",$oniom->{ATOMNUMProdGjf}," atoms in $prodgjf!\n";
   print "Please check!\n";
   die;
  }


print "There are $atomindex atoms in file $prodonb. \n" if $debug >= 0; 
print "Done!\n" if $debug >= 0;

close (PRODONBFILE);


# Sanity check for Model system

print "Sanity check for Model system...\n" if $debug >= 0;

for ($i=0; $i < $oniom->{ATOMNUMModelGjf}; $i++) {
  
  if ($oniom->{ATOMSModelGjf}[$i] ne $oniom->{ATOMSModelONB}[$i]) {
    print "Atom ",$i+1," line ",$oniom->{LineNumModelGjf}[$i],
        " in file $modelgjf is different \n  between $modelgjf and $modelonb.\n" if $debug >= 0;
    print " Original line in $modelgjf is\n",$oniom->{ORIGINALLINEModelGjf}[$i] if $debug >= 0;
    print " Original line in $modelonb is\n",$oniom->{ORIGINALLINEModelONB}[$i] if $debug >= 0;
    print "\n\n" if $debug >= 0;
  
   }
   

 } # for ($i=0; $i < $oniom->{ATOMNUMModelGjf}; $i++)

print "Done!\n" if $debug >= 0;

# Sanity check for Production system

print "Sanity check for Production system...\n" if $debug >= 0;

for ($i=0; $i < $oniom->{ATOMNUMProdGjf}; $i++) {
  
  if ($oniom->{ATOMSProdGjf}[$i] ne $oniom->{ATOMSProdONB}[$i]) {
    print "Atom ",$i+1," line ",$oniom->{LineNumProdGjf}[$i],
        " in file $prodgjf is different \n  between $prodgjf and $prodonb.\n" if $debug >= 0;
    print " Original line in $prodgjf is\n",$oniom->{ORIGINALLINEProdGjf}[$i] if $debug >= 0;
    print " Original line in $prodonb is\n",$oniom->{ORIGINALLINEProdONB}[$i] if $debug >= 0;
    print "\n\n" if $debug >= 0;
  
   }
   

 } # for ($i=0; $i < $oniom->{ATOMNUMModelGjf}; $i++)


print "Done!\n" if $debug >= 0;


# Check each atom in model system if it moved from its original position in ONB file or not

print "Check each atom in model system if it moved from its original position in ONB file...\n" if $debug >= 0; 

for ($i=0; $i < $oniom->{ATOMNUMModelGjf}; $i++) {
  
  $distance = $oniom->atomdistance(
    $oniom->{CARTCOORDModelGjf}[3*$i],$oniom->{CARTCOORDModelGjf}[3*$i+1],$oniom->{CARTCOORDModelGjf}[3*$i+2],
    $oniom->{CARTCOORDModelONB}[3*$i],$oniom->{CARTCOORDModelONB}[3*$i+1],$oniom->{CARTCOORDModelONB}[3*$i+2]);
  
  if ( $distance > $tol) {
     $oniom->{MovedModelGjf}[$i] = 1;
     print "Atom ",$i+1," in model system moved by $distance.\n" if $debug >= 1; 
     }

 } # for ($i=0; $i < $oniom->{ATOMNUMModelGjf}; $i++)

print "Done!\n" if $debug >= 0;

# Looking for corresponding atom in production system for each atom in model system 
#  if it moved from its original position in ONB file 

print "Looking for corresponding atom in production system for each atom in model system\n" if $debug >= 0; 
print " if it moved from its original position in ONB file...\n" if $debug >= 0; 

for ($i=0; $i < $oniom->{ATOMNUMModelGjf}; $i++) {
  
  if ($oniom->{MovedModelGjf}[$i] == 0) {next;}
  
  $oniom->{ProdIndexModelGjf}[$i] = $oniom->FindModelAtomInProd($i);

  # This corresponding atom in production system is also in model systemsss
  $oniom->{InModelProdGjf}[ $oniom->{ProdIndexModelGjf}[$i] ] = 1;

  # Save atom index in model system for it corresponding atom in production system
  $oniom->{ModelAtmIndexProdGjf}[ $oniom->{ProdIndexModelGjf}[$i] ] = $i;
  
  # If the atom in model system moved during optimization, then the corresponding atom in production
  #  system needs new coordinate. Label it by $oniom->{NeedReplaceProdGjf}
  if ( $oniom->{MovedModelGjf}[$i] ==1 ){
    $oniom->{NeedReplaceProdGjf}[ $oniom->{ProdIndexModelGjf}[$i] ] = 1;
   }
  
  print "Match atom ",$i+1," in model system and atom ",$oniom->{ProdIndexModelGjf}[$i]+1," in production system.\n" if $debug >= 2; 
  

  # Sanity check for matched atom pair

  print "Sanity check for matched atoms\n"  if $debug >= 2; 

  $j = $oniom->{ProdIndexModelGjf}[$i];
  
  $mismatch = 0;
  
  if ( $oniom->{ATOMSModelGjf}[$i] ne  $oniom->{ATOMSProdGjf}[$j]) {
    print "Atom mismatch: model atom ",$i+1," (",$oniom->{ATOMSModelGjf}[$i],
        ") : Production atom ",$j+1," (",$oniom->{ATOMSProdGjf}[$j],")\n" if $debug >= 1;
    $mismatch = 1;
   }
    
  if ( $oniom->{ATOMTYPEModelGjf}[$i] ne  $oniom->{ATOMTYPEProdGjf}[$j]) {
    print "Atom type mismatch: model atom ",$i+1," (",$oniom->{ATOMTYPEModelGjf}[$i],
        ") : Production atom ",$j+1," (",$oniom->{ATOMTYPEProdGjf}[$j],")\n" if $debug >= 1;
    $mismatch = 1;
   }
 
  if ( $oniom->{PARTIALCHARGModelGjf}[$i] !=  $oniom->{PARTIALCHARGProdGjf}[$j]) {
    print "Partial charge mismatch: model atom ",$i+1," (",$oniom->{PARTIALCHARGModelGjf}[$i],
        ") : Production atom ",$j+1," (",$oniom->{PARTIALCHARGProdGjf}[$j],")\n" if $debug >= 1;
    $mismatch = 1;
   }
 
  if ( $oniom->{ATOMFREEZEFLAGModelGjf}[$i] !=  $oniom->{ATOMFREEZEFLAGProdGjf}[$j]) {
    print "Moving flag mismatch: model atom ",$i+1," (",$oniom->{ATOMFREEZEFLAGModelGjf}[$i],
        ") : Production atom ",$j+1," (",$oniom->{ATOMFREEZEFLAGProdGjf}[$j],")\n" if $debug >= 1;
    $mismatch = 1;
   }
 
  if ( $oniom->{AtomLayerModelGjf}[$i] ne  $oniom->{AtomLayerProdGjf}[$j]) {
    print "Layer setup mismatch: model atom ",$i+1," (",$oniom->{AtomLayerModelGjf}[$i],
        ") : Production atom ",$j+1," (",$oniom->{AtomLayerProdGjf}[$j],")\n" if $debug >= 1;
    $mismatch = 1;
   }

  if ( $mismatch ==1 ) {
  
    print " Original line in $oniom->{ModelGjfFile} line(",$oniom->{LineNumModelGjf}[$i],") is \n",
           $oniom->{ORIGINALLINEModelGjf}[$i] if $debug >= 1;
    print " Original line in $oniom->{ProdGjfFile} line(",$oniom->{LineNumProdGjf}[$j],") is \n",
           $oniom->{ORIGINALLINEProdGjf}[$j] if $debug >= 1;
 
   
   }

 } # for ($i=0; $i < $oniom->{ATOMNUMModelGjf}; $i++)

print "Done!\n" if $debug >= 0;


} # sub digest 



#Subroutine FindModelAtomInProd: find corresponding atom for given model atom in production system

sub FindModelAtomInProd{
 my $oniom = shift;
 my $atomindex = shift;
 my $distance;
 my $debug = $oniom->{DEBUG};
 my $i;
 my $modelgjf = $oniom->{ModelGjfFile};
 my $modelonb = $oniom->{ModelONBFile};
 my $prodgjf = $oniom->{ProdGjfFile};
 my $prodonb = $oniom->{ProdONBFile};

 
 for ($i=0; $i < $oniom->{ATOMNUMProdGjf}; $i++){
   
   $distance = $oniom->atomdistance(
    $oniom->{CARTCOORDModelONB}[3*$atomindex],$oniom->{CARTCOORDModelONB}[3*$atomindex+1],
    $oniom->{CARTCOORDModelONB}[3*$atomindex+2],
    $oniom->{CARTCOORDProdONB}[3*$i],$oniom->{CARTCOORDProdONB}[3*$i+1],
    $oniom->{CARTCOORDProdONB}[3*$i+2]);
   
   if ($distance > $oniom->{tolerance}) { next;}
   

   # Sanity check for matched atom pair.
   
   if ($oniom->{ATOMSModelONB}[$atomindex] ne $oniom->{ATOMSProdONB}[$i]) {
    print "Atom ",$atomindex+1," line ",$oniom->{LineNumProdGjf}[$atomindex]," in ",$oniom->{ModelGjfFile},
        " in file $modelgjf is different \n  between $modelonb and $prodonb.\n" if $debug >= 0;
    print " Line number in $prodgjf is ",$oniom->{LineNumProdGjf}[$i],".\n" if $debug >= 0;
    print " Original line in $modelgjf is\n",$oniom->{ORIGINALLINEModelGjf}[$atomindex] if $debug >= 0;
    print " Original line in $modelonb is\n",$oniom->{ORIGINALLINEModelONB}[$atomindex] if $debug >= 0;
    print " Original line in $prodgjf is\n",$oniom->{ORIGINALLINEProdGjf}[$i] if $debug >= 0;
    print " Original line in $prodonb is\n",$oniom->{ORIGINALLINEProdONB}[$i] if $debug >= 0;
    print "\n" if $debug >= 0;
  
   }
   return $i;
 
 
  }
 
 print "Could match any atom in production system to atom ",$atomindex+1," in model system.\n" if $debug >= 0;
 print "This atom is located at line ",$oniom->{LineNumModelGjf}[$atomindex]," in $oniom->{ModelGjfFile}.\n" if $debug >= 0;
 print "Original line in $oniom->{ModelGjfFile} is \n",$oniom->{ORIGINALLINEModelGjf}[$atomindex] if $debug >= 0;
 print "Program ends here.\n" if $debug >= 0;
 die;

} # sub FindModelAtomInProd{


#Subroutine AtomDistance: calculate distance between two atoms

sub atomdistance{

 my $oniom = shift;
 my $X1 = shift;
 my $Y1 = shift;
 my $Z1 = shift;
 my $X2 = shift;
 my $Y2 = shift;
 my $Z2 = shift;
 my $distance;
 
 $distance = sqrt( ( $X1 - $X2 )**2 + ( $Y1 - $Y2 )**2 + ( $Z1 - $Z2 )**2  );
 return $distance;

} # sub atomdistance


#Subroutine WriteProdONIOM: write production ONIOM input file with geometry from model system

sub writeprodoniom{

 my $oniom = shift;
 my $debug = $oniom->{DEBUG};

 my $tol             = $oniom->{tolerance}  ;
 my $atomprodflag    = $oniom->{atomprod}   ;
 my $atmtypprodflag  = $oniom->{atmtypprod} ;
 my $chgprodflag     = $oniom->{chgprod}    ;
 my $moveprodflag    = $oniom->{moveprod}   ;
 my $layerprodflag   = $oniom->{layerprod}  ;
 my $coordprod       = $oniom->{coordprod}  ;
 my $prodnewgjf      = $oniom->{ProdNewInput}; 
 my $i;
 my $currline;
 my ($atom,$atomtyp,$partialcharge,$moveflag,$layer,$x,$y,$z);
 my $modelindex;

 
 open(INPUTFILE,$oniom->{ProdGjfFile}) || die "Could not read $oniom->{ProdGjfFile}\n$!\n";
 
 print "Open $prodnewgjf to write new ONIOM input file...\n" if $debug >= 0;

 open(OUTPUTGJFFILE,">$prodnewgjf") || die "Could not open $prodnewgjf\n$!\n";
 
 while (<INPUTFILE>) {  
     
     # Keep everything before atom list

     print OUTPUTGJFFILE $_ ;

     if (/^\s*-*\d+\s+\d+\s+-*\d+\s+\d+/) { last;} 
     
  } # while (<INPUTFILE>) { 

 for ($i=0; $i< $oniom->{ATOMNUMProdGjf}; $i++) {
   
   $currline = <INPUTFILE>;
   
   # Use original file if do not need to be replaced;
   
   if ( $oniom->{NeedReplaceProdGjf}[$i] == 0 ){
     print OUTPUTGJFFILE $currline ;
     next;
    } 
   
   # If need to be replaced, take info from Model ONIOM input file as default
   $modelindex =  $oniom->{ModelAtmIndexProdGjf}[$i];
   
   $atom = $oniom->{ATOMSModelGjf}[$modelindex];
   $atomtyp = $oniom->{ATOMTYPEModelGjf}[$modelindex];
   $partialcharge = $oniom->{PARTIALCHARGModelGjf}[$modelindex];
   $moveflag = $oniom->{ATOMFREEZEFLAGModelGjf}[$modelindex];
   $layer = $oniom->{AtomLayerModelGjf}[$modelindex];

   $x = $oniom->{CARTCOORDModelGjf}[3*$modelindex];
   $y = $oniom->{CARTCOORDModelGjf}[3*$modelindex+1];
   $z = $oniom->{CARTCOORDModelGjf}[3*$modelindex+2];
   
   # Use setup from Production input file if flags are set
   if ($atomprodflag == 1) { $atom = $oniom->{ATOMSProdGjf}[$i];}
   if ($atmtypprodflag == 1) { $atomtyp = $oniom->{ATOMTYPEProdGjf}[$i];}
   if ($chgprodflag == 1) { $partialcharge = $oniom->{PARTIALCHARGProdGjf}[$i];}
   if ($moveprodflag == 1) { $moveflag = $oniom->{ATOMFREEZEFLAGProdGjf}[$i];}
   if ($layerprodflag == 1) { $layer = $oniom->{AtomLayerProdGjf}[$i];}
   if ($coordprod == 1 ) {
      $x = $oniom->{CARTCOORDProdGjf}[3*$i];
      $y = $oniom->{CARTCOORDProdGjf}[3*$i+1];
      $z = $oniom->{CARTCOORDProdGjf}[3*$i+2];
    }   
   
   print "Before the line is \n",$currline if $debug >=2;
   
   
   $currline =~ s/([a-zA-Z]{1,3})-([a-zA-Z0-9\*]{0,2})-(-*\d*\.\d+)/$atom-$atomtyp-$partialcharge/;
   $currline =~ s/\s(-*[01])\s/ $moveflag /;
   $currline =~ s/\s(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+([HML])/ $x     $y     $z   $layer/;


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

F<ESPT::ONIOMinput>

=head1 AUTHOR

Peng Tao, E<lt>tao.21@osu.eduE<gt>

=head1 COPYRIGHT

Copyright (c) 2009~2010 by Peng Tao


=cut
