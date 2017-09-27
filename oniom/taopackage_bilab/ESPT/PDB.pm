package ESPT::PDB;

=head1 NAME

ESPT::PDB - Read in pdb file, and output ONIOM input file

=head1 SYNOPSIS

    use ESPT::PDB;

=head1 DESCRIPTION

This module reads in PDB (Protein Data Bank) file, and outputs Gaussian ONIOM
input file with atom type and partial charges. This module is designed to deal
with AMBER force field only. Atom type and partial charges are assigned by
AMBERFF package.

=cut

=head1 IMPORTANT

 For ONIOM calculation with protein, you may need connectivity 
 information. However, this program cannot generate those information. 

=cut

=head1 NOTICE

Besides .gjf file, there is also one .onb file as output. This file contains the 
original residue and atom type information. These information might be useful
when coding up ONIOM layer selection part.

=over 15

=item B<$file-E<gt>new()>

Creates a new PDB object

=cut

use base qw(ESPT::ESSfile);
use ESPT::Glib 0.01;
use ESPT::AMBERFF;
use strict;
use warnings;

our $VERSION = '0.04';

### Version History ###
# 0.01	digest pdb file
# 0.02  Deal with atom name starting with number in PDB file. eg "1HH1"
# 0.03  Deal with freeze flag (0 or -1) within certain distance of core residues 
# 0.04  Abstract part of structure within certain distance of core residues 

### To Do List ###

##### NOTE #####
#
################

## the object constructor **

sub new {
	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my $pdb = ESPT::ESSfile->new();

	$pdb->{PROGRAM} = "PDB";
	$pdb->{TYPE} = "pdb";

	# Link 0 & Route commands
	$pdb->{ROUTE} = undef;
	$pdb->{KEYWORDS} = [];


	# molecular info
	$pdb->{ATOMS} = [];
        $pdb->{RESIDUE} = [];
        $pdb->{RESIDUENUM} = [];
        $pdb->{ATOMTYPE} = [];
        $pdb->{ELEMENTNUM} = [];
        $pdb->{ELEMENT} = [];
        $pdb->{PARTIALCHARG} = [];
        $pdb->{CARTCOORD} = [];
        $pdb->{ATOMDEFINE} = [];
        $pdb->{ATOMFREEZEFLAG} = [];
        $pdb->{ATOMINRESIDUE} = [];
        $pdb->{ORIGINALLINE} = [];
        
	# These variables are related to core residues set up
        $pdb->{CORERESIDUE} = [];
        $pdb->{CORERESIDUENAME} = [];
        $pdb->{CORERESIDUENO} = undef;
        $pdb->{RESIDFILE} = undef;
        $pdb->{NEAR} = undef;
        $pdb->{CORESET} = undef;
        $pdb->{RESIDUENEARCORE} = undef;
        
        # These variables are list of core atoms
        $pdb->{COREATOMSID} = [];
        $pdb->{COREATOMSTOT} = undef;
        
        # Variables for labeling the residues to assist freeze flag setup
        $pdb->{RESIDUEIDLIST} = [];
        $pdb->{RESIDUENAMELIST} = [];
        $pdb->{RESIDUEWITHINCORE} = [];
        $pdb->{TOTALRESIDUENUM} = undef;

        # M.Yamada
        # For .top file
        $pdb->{TOPFILENAME} = undef;
        $pdb->{TOPEXIST} = undef;

        # Calculation levels of each atom
        $pdb->{CALCLEVEL} = []; # NATOMS
        
        # Bond Conection
        $pdb->{BOND1} = []; # NBONDS
        $pdb->{BOND2} = []; # NBONDS
	$pdb->{NBONDS} = undef;
        # M.Yamada end

	$pdb->{CHARGE} = undef;
	$pdb->{EELEC} = undef;		# SCF electronic energy
	$pdb->{MULTIPLICITY} = undef;
	$pdb->{NATOMS} = undef;
	bless($pdb, $class);
	return $pdb;
}


## methods ##

# set filename & spin then digest the file
sub analyze : method {
	my $pdb = shift;
	$pdb->prepare(@_);
	$pdb->digest();
	return;
}


## subroutines ##

sub digest {

my $pdb = shift;

# flags & counters
my $atomflag = 0;
my $atomtot = 0;
my $cartflag = 0;
my $carttot = 0;
my $col = 0;
my $counter = 0;
my $geomcount = 0;
my $row = 0;
my $rparsed = 0;
my $Titleflag = 0;
my $Cflag = 0;
my $weightflag = 0;
my $debug = 0;
my $curatmtyp;
my $amberffdum;
my $curatomname;
my $curresidname;
my $curelenum;
my $curparchg;
my $curX;
my $curY;
my $curZ;
my $atomdefined;
my $curresidnum;
my $coreresindex;
my $curcoreresid;
my $curcoreresname;
my $resindex;
my $residref = -999999;
my $resnameref = "";
my $i;
my $coreatomindex;
my $atomcounter = 0; # M.Yamada
my $bondcounter = 0; # M.Yamada
my $moltypecounter = 0; # M.Yamada
my $curtopatomnum; # M.Yamada
my $curtopatomname; # M.Yamada
my $curtopresname; # M.Yamada
my $topread = ""; # M.Yamada

$debug = $pdb->{DEBUG} ;


$amberffdum  = ESPT::AMBERFF->new();

# Check if core residues set up
if ( $pdb->{RESIDFILE} eq "THERE_IS_NO_CORE_RESIDUE_SET_HERE") {
  print "\nNo core residues list provided. All atom are free to move (0) during geometry optimization.\n"  if $debug >= 0;
  $pdb->{CORESET} = 0;
 } else { 
  print "\nCore residues list file $pdb->{RESIDFILE} provided.\n" if $debug >= 0;
  print "All residues within $pdb->{NEAR} angstroms from core region are free to move (0) during geometry optimization.\n"  if $debug >= 0;
  $pdb->{CORESET} = 1; 
 }

# if core residues set up, the read in list of residues numbers in the core

if ( $pdb->{CORESET} == 1 ){
 print "\nOpen $pdb->{RESIDFILE} file for process.\n" if $debug >= 1;
 open(CORERESFILE,$pdb->{RESIDFILE}) || die "Could not read $pdb->{RESIDFILE}\n$!\n";

 $coreresindex = 0;

 while (<CORERESFILE>) {
 
  next if /^$/;
  /\[(.+)\]\s+"(.+)"/;
  $pdb->{CORERESIDUENAME}[$coreresindex] = $1;
  $pdb->{CORERESIDUE}[$coreresindex] = $2;
  
  $curcoreresid = $pdb->{CORERESIDUE}[$coreresindex];
  $curcoreresname = $pdb->{CORERESIDUENAME}[$coreresindex];
  $coreresindex++;

  print "\nRead in residue number $coreresindex with residue name as ($curcoreresname) residue number as ($curcoreresid) \n" if $debug >= 1;
    
  } # while (<CORERESFILE>) {


  $pdb->{CORERESIDUENO} = $coreresindex;
  
  if ( $coreresindex == 0) {
   print "\nThere is no residue number found in file $pdb->{RESIDFILE}, please double check.\n";
   die;
   }
 
  print "\nThere are $coreresindex residues in the core region. \n" if $debug >= 1;

 } # if ( $pdb->{CORESET} == 1 )


# open filename for reading or display error
print "\nOpen $pdb->{FILENAME} file for process.\n" if $debug >= 1;

open(PDBFILE,$pdb->{FILENAME}) || die "Could not read $pdb->{FILENAME}\n$!\n";

# set atom and residues counters as 0
$counter = 0; 
$resindex = -1;

# Process PDB file for ONIOM input file
while (<PDBFILE>) {
  print $_ if $debug >= 3;
  if (/^ATOM\s+\d+/ or /^HETATM\s+\d+/) {
    $pdb->{ORIGINALLINE}[$counter] = $_;
    if (/^ATOM\s+\d+/) {
#       /^ATOM\s+\d+\s+(\S+)\s+(\S+)\s+(.{6})\s+(-?\d+.\d+)\s+(-?\d+.\d+)\s+(-?\d+.\d+)\s+/;
       /^ATOM\s+\d+\s+(\S+)\s+(\S+)\s+(\d+)\s+(-?\d+.\d+)\s+(-?\d+.\d+)\s+(-?\d+.\d+)\s+/;
#       print "Atom is $1, residue is $2, Coordinates are $4, $5, $6,\n" if $debug >= 1;
       $curatomname = $1;
       $curresidname = $2;
       $curresidnum = $3;
       $curX = $4;
       $curY = $5;
       $curZ = $6;
       if ( $curatomname =~ /^[0-9]/) {
         print "\nFound nonstandard atom name $curatomname in residue $curresidname.\n" if $debug >= 1;
         $curatomname =~ s/^([0-9])(\S+)/$2$1/;
         print "\nNew atom name $curatomname to meet AMBER standard name.\n" if $debug >= 1;
         }
         
       

     } else { 
#       /^HETATM\s+\d+\s+(\S+)\s+(\S{3})\s+(.{6})\s+(-?\d+.\d+)\s+(-?\d+.\d+)\s+(-?\d+.\d+)\s+/;
       /^HETATM\s+\d+\s+(\S+)\s+(\S+)\s+(\d+)\s+(-?\d+.\d+)\s+(-?\d+.\d+)\s+(-?\d+.\d+)\s+/;
#       print "Atom is $1, residue is $2, Coordinates are $3, $4, $5,\n" if $debug >= 1;
       $curatomname = $1;
       $curresidname = $2;
       $curresidnum = $3;
       $curX = $4;
       $curY = $5;
       $curZ = $6;
     } #  } else { 
     
     if ( ( $amberffdum->existsresidue($curresidname)) ) {
       # Look up atom type and partial charges for current atom
       $curatmtyp = $amberffdum->atomtype($curresidname,$curatomname);
       $curelenum = $amberffdum->elementnumber($curresidname,$curatomname);
       $curparchg = $amberffdum->partialcharge($curresidname,$curatomname);
       $atomdefined = 1;
       } else {
       # In this case, the residue information is not available from force file library file
        printf "\nResidue %6s does not exist in database. Atom with name %5s may not be defined.\n", $curresidname,$curatomname;
        $curatmtyp = "UDF";
        $curelenum = "UDF";
        $curparchg = 0.0;
        $atomdefined = 0;
        $curelenum = 1;
       } #else
     
#     print "Current atom has atom type $curatmtyp, element number $curelenum, partial charge $curparchg.\n";
#     print "Current residue number is $curresidnum\n";

     $pdb->{ATOMS}[$counter] = $curatomname;
     $pdb->{RESIDUE}[$counter] = $curresidname;
     $pdb->{RESIDUENUM}[$counter] = $curresidnum;
     if ( defined $curatmtyp) {
          $pdb->{ATOMTYPE}[$counter] = $curatmtyp;
         } else {
        print "\nAtom type cannot be assigned to atom $curatomname in residue $curresidname $curresidnum.\n" if $debug >= 0;
          $pdb->{ATOMTYPE}[$counter] = "UDF";
         }

     if ( defined $curparchg) {
        $pdb->{PARTIALCHARG}[$counter] = $curparchg;
       } else {
         print "\nPartial charge cannot be assigned to atom $curatomname in residue $curresidname $curresidnum.\n" if $debug >= 0;
         $pdb->{PARTIALCHARG}[$counter] = 0.0000;
       }

     $pdb->{ELEMENTNUM}[$counter] = $curelenum;
     
     if ($pdb->{CORESET} == 1) {
         $pdb->{ATOMFREEZEFLAG}[$counter] = -1;
        } else {
         $pdb->{ATOMFREEZEFLAG}[$counter] = 0;
        }
     if ( defined $curelenum) {
        $pdb->{ELEMENT}[$counter] = $pdb->atomconvert($curelenum);
       } else {
        print "\nElement type cannot be decided for atom $curatomname in residue $curresidname $curresidnum.\n" if $debug >= 0;
        $pdb->{ELEMENT}[$counter] = "UDF";
       }
     $pdb->{ATOMDEFINE}[$counter] = $atomdefined;
     $pdb->{CARTCOORD}[3*$counter] = $curX;
     $pdb->{CARTCOORD}[3*$counter+1] = $curY;
     $pdb->{CARTCOORD}[3*$counter+2] = $curZ;
     $pdb->{ATOMINRESIDUE}[$counter] = $resindex;

     $counter++;
     
     #Check if encounter a new residue, if yes, store residue number and name. 
     
     if ( ($residref eq $curresidnum) and ( $resnameref eq $curresidname) )
      { next; }
     
     $resindex++;
     $pdb->{RESIDUEIDLIST}[$resindex] = $curresidnum;
     $pdb->{RESIDUENAMELIST}[$resindex] = $curresidname;
     $pdb->{RESIDUEWITHINCORE}[$resindex] = 0;
     $residref = $curresidnum;
     $resnameref = $curresidname;
     $pdb->{ATOMINRESIDUE}[$counter-1] = $resindex;
     print "\nEncounter new residue $resindex, name $curresidname, index $curresidnum.\n" if $debug >= 2;
  
   } # if (/^ATOM\s+\d+/ or /^HETATM\s+\d+/) 
  
 } #  while (<PDBFILE>) {

 $resindex++;
 print "\nThere are $resindex residues in the PDB file.\n" if $debug >= 0;
 $pdb->{TOTALRESIDUENUM} =$resindex; 
 $pdb->{NATOMS} = $counter ;

# M.Yamada
if ( $pdb->{TOPEXIST} == 1 ) {

 print "\nOpen $pdb->{TOPFILENAME} file for process.\n" if $debug >= 1;
 open(TOPFILE,$pdb->{TOPFILENAME}) || die "Could not read $pdb->{TOPFILENAME}\n$!\n";

 # set 0
 $atomcounter = 0; 
 $bondcounter = 0; 
# $moltypecounter = -1; 
# $topread = "mol"; 

 # Process TOP file for ONIOM input file
 while (<TOPFILE>) {
   print $_ if $debug >= 3;
   last if (/^$/ and $topread eq "bonds"); # Read only the first [ bonds ]
   next if /^$/;
   next if /^\s*;/;
   if ( /\[\s*moleculetype\s*\]/ ) {
     $topread = "atom";
     next;
   }
   if ( /\[\s*atoms\s*\]/ and $topread eq "atom" ) {
     $topread = "atoms";
#     $moltypecounter++; 
     next;
   }
#   if ( /^\s*(\d+)\s+\S+\s+\d+\s+(\S+)\s+(\S+)\s+/ and $topread == "atoms" ) {
#     $curtopatomnum = $1;
#     $curtopresname = $2;
#     $curtopatomname = $3;
#     next;
#   }
   if ( /\[\s*bonds\s*\]/ and $topread eq "atoms" ) {
     $topread = "bonds";
     next;
   }
   if ( /^\s*(\d+)\s+(\d+)\s+\.*/ and $topread eq "bonds" ) {
     $pdb->{BOND1}[$bondcounter] = $1;
     $pdb->{BOND2}[$bondcounter] = $2;
     $bondcounter++; 
     next;
   }
  
 } #  while (<TOPFILE>) {
 $pdb->{NBONDS} = $bondcounter;

} # if ( $pdb->{TOPEXIST} == 1 ) {
# M.Yamada END

 # If there is no core residues list provided, no need to set freeze flag
 if ($pdb->{CORESET} == 0) { return; }
 
 # Find list of atoms which belong to core region
 $coreatomindex = 0;
 for ($i=0; $i < $pdb->{NATOMS};$i++) {
   
   if ( $pdb->atombelong2core($i) == 0) {
     $pdb->{CALCLEVEL}[$i] = "L"; # M.Yamada
     next;
   }
    
   $pdb->{COREATOMSID}[$coreatomindex] = $i;
   $pdb->{CALCLEVEL}[$i] = "H"; # M.Yamada
   $coreatomindex++;
   print "\nAtom $i (in residue $pdb->{RESIDUE}[$i] $pdb->{RESIDUENUM}[$i]) belongs to core region.\n" if $debug >= 2;
 
  } # for (my $i=0; $i < $pdb->{NATOMS};$i++)
 
 $pdb->{COREATOMSTOT} = $coreatomindex;
 print "\nThere are total $coreatomindex atoms in the core region.\n" if $debug >= 1;
 
 $pdb->decideflag();

} # sub digest 

#Subroutine that decides the freeze flag of each atom according to its distance
# to core residues.

sub decideflag {
 my $pdb = shift;
 my $curdistance;
 my $curresidue;
 my $debug = $pdb->{DEBUG};
 my $totresiduewithin = 0;
 my $frzbackboneatoms = 0;
 
 
 #  Check if core residues are set
 if ( $pdb->{CORESET} == 0 )
  { return;}

 # go through all the atoms to check if each residue are within certain distance of core or not.
 for (my $i=0; $i < $pdb->{NATOMS};$i++) {
   
   if ($pdb->{RESIDUEWITHINCORE}[$pdb->{ATOMINRESIDUE}[$i]] == 1) {next;}
   
   if ( $pdb->atomwithincore($i) == 1) {
     $pdb->{RESIDUEWITHINCORE}[$pdb->{ATOMINRESIDUE}[$i]] = 1;
     $curresidue = $pdb->{ATOMINRESIDUE}[$i];
     $totresiduewithin++;
     print "\nResidue ($curresidue) is close to core region.\n" if $debug >= 2;
     
    }
   
   
  } # for (my $i=0; $i < $pdb->{NATOMS};$i++)

 print "\nThere are $totresiduewithin residues close to core region (including those residues in the core region). \n" if $debug >= 1;

 $pdb->{RESIDUENEARCORE} = $totresiduewithin;

 # go through all the atoms to setup freeze flag.
 for (my $i=0; $i < $pdb->{NATOMS};$i++) {
  
   if ($pdb->{RESIDUEWITHINCORE}[$pdb->{ATOMINRESIDUE}[$i]] == 0) {next;}
   
   # If this atom belongs to core region, it is free to move.
   if ( $pdb->atombelong2core($i) == 1 ) {
    $pdb->{ATOMFREEZEFLAG}[$i] = 0;
    next;
    }
   
   # If this atom belongs to a residue that is within core, it is free to move.
   
   if ($pdb->{RESIDUEWITHINCORE}[$pdb->{ATOMINRESIDUE}[$i]] == 1){
      $pdb->{ATOMFREEZEFLAG}[$i] = 0;
      next;
    }
   
  }# for (my $i=0; $i < $pdb->{NATOMS};$i++)
  
# print "\nThere are $frzbackboneatoms backbone atoms that are frozen, but in residue close to core region.\n" if $debug >= 1;
 
} # sub decideflag 

# Subroutine  isbackbone: Check if this atom is backbone atom or not.
sub isbackbone{
 my $pdb = shift;
 my $atomid = shift;
 my $debug = $pdb->{DEBUG};
 
 if ($pdb->{ATOMTYPE}[$atomid] eq "C") {return 1;}
 if ($pdb->{ATOMTYPE}[$atomid] eq "O") {return 1;}
 if ($pdb->{ATOMTYPE}[$atomid] eq "N") {return 1;}
 if ($pdb->{ATOMTYPE}[$atomid] eq "H") {
  if ($pdb->{ATOMTYPE}[$atomid-1] eq "N")  {return 1;}
  }
 return 0;


} # sub isbackbone

# Subroutine atomwithincore: Check if one atom is within certain distance of core or not.

sub atomwithincore{

 my $pdb = shift;
 my $atomid = shift;
 my $curcoreatom;
 my $distance;
 my $debug = $pdb->{DEBUG};

 for (my $i=0; $i < $pdb->{COREATOMSTOT};$i++) {
   $curcoreatom = $pdb->{COREATOMSID}[$i];
   $distance = $pdb->atomdistance($atomid,$curcoreatom);
   if ( $distance < $pdb->{NEAR}) {
   
      print "\nAtom $atomid is at least $distance A to core region. \n" if  $debug >= 2;
      return 1;
     } 
   
 
  } # for (my $i=0; $i < $pdb->{NATOMS};$i++)

 print "\nAtom $atomid is NOT close to core region. \n" if  $debug >= 2;

 return 0;
 
} # sub atomwithincore


#Subroutine AtomDistance: calculate distance between two atoms

sub atomdistance{

 my $pdb = shift;
 my $atomid1 = shift;
 my $atomid2 = shift;
 my $distance;
 
 $distance = sqrt( ( $pdb->{CARTCOORD}[3*$atomid1] - $pdb->{CARTCOORD}[3*$atomid2] )**2 + 
                   ( $pdb->{CARTCOORD}[3*$atomid1+1] - $pdb->{CARTCOORD}[3*$atomid2+1] )**2 + 
                   ( $pdb->{CARTCOORD}[3*$atomid1+2] - $pdb->{CARTCOORD}[3*$atomid2+2] )**2  );
 return $distance;


} # sub atomdistance

# Subroutine atombelong2core: check if given atom belong to core residues or not
sub atombelong2core{

 my $pdb = shift;
 my $atomid = shift;
 my $debug = $pdb->{DEBUG};
 
 my $oriname = $pdb->{RESIDUE}[$atomid];
 my $oriid = $pdb->{RESIDUENUM}[$atomid];
 my $curname;
 my $curid;
 
 
 for (my $j=0; $j < $pdb->{CORERESIDUENO}; $j++) {
   $curname  = $pdb->{CORERESIDUENAME}[$j];
   $curid    = $pdb->{CORERESIDUE}[$j];
   
   print "\nCompare (".$oriname.") with (".$curname.") and (".
                    $oriid.") with (".$curid.")\n" if $debug >= 2;

   if ( (($oriname =~ /$curname/) or ($curname =~ /$oriname/)) and 
        (($oriid =~ /$curid/) and ($curid =~ /$oriid/)) ){ 
          print "\nMatch (".$oriname.") with (".$curname.") and (".
                    $oriid.") with (".$curid.")\n" if $debug >= 2;
          return 1;
        }
  }
 
 return 0;

} # sub atombelong2core

# Subroutine writeoniom:  write to output files
sub writeoniom {
 
 my $pdb = shift;
 my $filename = shift;
 my $debug = $pdb->{DEBUG};
 my $oniompdb;
 my $currfrozenflag;
 my $modellevel; # M.Yamada
 my $bondcounter = 0; # M.Yamada
  
 print "\nWrite ONIOM input file $filename from PDB file.\n" if $debug >= 0; 
 
 $oniompdb = join "", $filename, ".onb";
 
 open(ONIOMPDBFILE, ">$oniompdb") || die "Could not create $oniompdb\n$!\n";

 # Open output file
 print "\nGaussian output file name is $filename\n" if $debug >= 1;
 open(ONIOMGJFFILE, ">$filename") || die "Could not create $filename\n$!\n";
 print "\nOpening file ", $filename, " for output ...\n" if $debug >= 0;

 # % commands
 print ONIOMGJFFILE "%chk=$filename.chk\n%mem=3700MB\n%nprocshared=4\n\n";

 # Route
 print ONIOMGJFFILE "#P ONIOM(B3LYP/6-31G(d):AMBER=HardFirst) geom=connectivity nosymm iop(2/15=3) test\n\n";
 if ( $pdb->{TOPEXIST} == 0 ) {
   print ONIOMGJFFILE "ONIOM inputfile generated by pdb2oniom from PDB file $pdb->{FILENAME}. No connectivity generated.\n";
   print ONIOMGJFFILE "Please use GaussView read $filename, and generate connectivity information.\n\n";
 } else {
   print ONIOMGJFFILE "ONIOM inputfile generated by pdb2oniom from PDB file $pdb->{FILENAME} and TOP file $pdb->{TOPFILENAME}.\n\n";
 }
 print ONIOMGJFFILE "0 1 0 1 0 1\n";
 
# print "Notice: For ONIOM calculation with protein, you may need connectivity information. \n";
# print "Notice: However, this program cannot generate those information. \n";
# print "Notice: You may use GaussView to read in original PDB file, and generate those information separately. \n";
# print "Notice: After you have connectivity information, you can copy them to appropriate place in Gaussian input file. \n";
 

 
 for (my $i=0; $i < $pdb->{NATOMS};$i++) {
   $currfrozenflag = $pdb->{ATOMFREEZEFLAG}[$i];
   $modellevel = $pdb->{CALCLEVEL}[$i]; # M.Yamada
#  print "Atom number $i: $pdb->{ATOMS}[$i], $pdb->{RESIDUE}[$i], $pdb->{ATOMTYPE}[$i], $pdb->{PARTIALCHARG}[$i], ";
#  print " $pdb->{ELEMENTNUM}[$i], $pdb->{CARTCOORD}[3*$i], $pdb->{CARTCOORD}[3*$i+1], $pdb->{CARTCOORD}[3*$i+2]\n";
   if ( $pdb->{ATOMDEFINE}[$i] == 1 ) {
     printf ONIOMGJFFILE " %s-%s-%.6f\t %2d  %15.7f %15.7f %15.7f   %s\n", $pdb->{ELEMENT}[$i],$pdb->{ATOMTYPE}[$i],
             $pdb->{PARTIALCHARG}[$i],$currfrozenflag,$pdb->{CARTCOORD}[3*$i], $pdb->{CARTCOORD}[3*$i+1], $pdb->{CARTCOORD}[3*$i+2],$modellevel;
             
     printf ONIOMPDBFILE " %s-%s-%.6f\t %2d  %15.7f %15.7f %15.7f   %s %6s %5s %7s\n", $pdb->{ELEMENT}[$i],$pdb->{ATOMTYPE}[$i],
             $pdb->{PARTIALCHARG}[$i],$currfrozenflag,$pdb->{CARTCOORD}[3*$i], $pdb->{CARTCOORD}[3*$i+1], $pdb->{CARTCOORD}[3*$i+2],$modellevel,
             $pdb->{ATOMS}[$i], $pdb->{RESIDUE}[$i],$pdb->{RESIDUENUM}[$i];
             
   } else {
     printf ONIOMGJFFILE " %s-%s-%.6f\t %2d  %15.7f %15.7f %15.7f   %s\n", $pdb->{ATOMS}[$i],$pdb->{ATOMTYPE}[$i],
             $pdb->{PARTIALCHARG}[$i],$currfrozenflag,$pdb->{CARTCOORD}[3*$i], $pdb->{CARTCOORD}[3*$i+1], $pdb->{CARTCOORD}[3*$i+2],$modellevel;

     printf ONIOMPDBFILE " %s-%s-%.6f\t %2d  %15.7f %15.7f %15.7f   %s %6s %5s %7s\n", $pdb->{ATOMS}[$i],$pdb->{ATOMTYPE}[$i],
             $pdb->{PARTIALCHARG}[$i],$currfrozenflag,$pdb->{CARTCOORD}[3*$i], $pdb->{CARTCOORD}[3*$i+1], $pdb->{CARTCOORD}[3*$i+2],$modellevel,
             $pdb->{ATOMS}[$i], $pdb->{RESIDUE}[$i],$pdb->{RESIDUENUM}[$i];
   }
 }
 # M.Yamada : write conectivity
 if ( $pdb->{TOPEXIST} == 1 ) {
   print ONIOMGJFFILE "\n";
   for (my $i=0; $i < $pdb->{NATOMS};$i++) {
     printf ONIOMGJFFILE " %d", $i+1;
     while () {
       last if $bondcounter >= $pdb->{NBONDS};
       if ( $pdb->{BOND1}[$bondcounter] == $i+1 ) {
         printf ONIOMGJFFILE " %d 1.0", $pdb->{BOND2}[$bondcounter];
         $bondcounter++; 
         next;
       } else {
         last;
       }
     } # while () {
     if ( $bondcounter >= $pdb->{NBONDS} and ($pdb->{RESIDUE}[$i] eq "WAT" or $pdb->{RESIDUE}[$i] eq "TP3" or $pdb->{RESIDUE}[$i] eq "HOH") and $pdb->{ATOMTYPE}[$i] eq "OW" ) {
       printf ONIOMGJFFILE " %d 1.0 %d 1.0", $i+2, $i+3;
     }
     printf ONIOMGJFFILE "\n";
   } # for NATOMS
 } # if ( $pdb->{TOPEXIST} == 1 ) {
 # M.Yamada END
 print ONIOMGJFFILE "\n\n\n";
 close (ONIOMGJFFILE);
 close (ONIOMPDBFILE);
 return;

} # sub writeniom {


# Subroutine writepdb: write PDB file which contains all the residues close to core region
sub writepdb {

 my $pdb = shift;
 my $filename = shift;
 my $debug = $pdb->{DEBUG};
 my $totalatom = 0;

 print "\nWrite PDB file $filename.\n" if $debug >= 0; 
 
 open(NEWPDBFILE, ">$filename") || die "Could not create $filename\n$!\n";
 
 for (my $i=0; $i < $pdb->{NATOMS};$i++) {
   
   # Skip if this residue is not close to core enough.
   if ($pdb->{RESIDUEWITHINCORE}[$pdb->{ATOMINRESIDUE}[$i]] == 0) {next;}
   $totalatom++;
   print NEWPDBFILE $pdb->{ORIGINALLINE}[$i]; 
  
  } 
 
 print "\nThere are $pdb->{RESIDUENEARCORE} residues and $totalatom atoms in output file $filename.\n" if $debug >= 1;

} # sub writepdb


1;

__END__

=back

=head1 VERSION

0.04

=head1 SEE ALSO

F<ESPT::AMBERFF>

=head1 AUTHOR

Peng Tao, E<lt>tao.21@osu.eduE<gt>

=head1 COPYRIGHT

Copyright (c) 2009~2010 by Peng Tao

=cut
