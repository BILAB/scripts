package ESPT::ONIOMtoPDB;

use base qw(ESPT::ESSfile);
use ESPT::Glib 0.01;
use ESPT::AMBERFF;
use strict;
use warnings;

our $VERSION = '0.2';

=head1 NAME

ESPT::ONIOMtoPDB - Use Oniom file and PDB template file object to create new PDB file.

=head1 SYNOPSIS

   use ESPT::ONIOMtoPDB;

   my $oniomtopdb = ONIOMtoPDB->new();

=head1 DESCRIPTION

This module provides methods to put optimized geometry from ONIOM file to PDB file.

=cut

### Version History ###
# 0.1	Treat ONIOM job and log file
# 0.2   Write partial charge into output PDB file


### To Do ###


## the object constructor **

sub new {
	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my $oniomtopdb = ESPT::ESSfile->new();

	$oniomtopdb->{TYPE} = "ONIOMtoPDB";

	# program info
	$oniomtopdb->{PROGRAM}     = "Gaussian";
	$oniomtopdb->{VERSION}     = undef;
	$oniomtopdb->{REVISION}    = undef;
	$oniomtopdb->{COMPILE}     = undef;
	$oniomtopdb->{COMPILEDATE} = undef;
	

        # Operation info
        $oniomtopdb->{G_FILE_NAME}      = undef;
	$oniomtopdb->{PDB_FILE_NAME}    = undef;
	$oniomtopdb->{STEP_NUMBER}      = 0;
        $oniomtopdb->{OUTPUT_FILE_NAME} = undef;
        $oniomtopdb->{IS_LOG_FILE}      = 0;
        $oniomtopdb->{IS_JOB_FILE}     = 0;
        $oniomtopdb->{WHOLE_PDB_FILE}     = 0;
        $oniomtopdb->{OutputCharge}     = 0;
        $oniomtopdb->{ChargeDec}        = 2;
        $oniomtopdb->{OutputOccupancy}  = 0;
        $oniomtopdb->{Occupancy}        = 1.00;


	# PDB file info
	$oniomtopdb->{PDB_ATOMS} = [];
        $oniomtopdb->{PDB_RESIDUE} = [];
        $oniomtopdb->{PDB_RESIDUENUM} = [];
        $oniomtopdb->{PDB_ATOMTYPE} = [];
        $oniomtopdb->{PDB_ELEMENTNUM} = [];
        $oniomtopdb->{PDB_ELEMENT} = [];
        $oniomtopdb->{PDB_PARTIALCHARG} = [];
        $oniomtopdb->{PDB_CARTCOORD} = [];
        $oniomtopdb->{PDB_ATOMDEFINE} = [];
        $oniomtopdb->{PDB_ATOMFREEZEFLAG} = [];
        $oniomtopdb->{PDB_ATOMINRESIDUE} = [];
        $oniomtopdb->{PDB_ORIGINALLINE} = [];
        $oniomtopdb->{PDB_LINENUM} = [];
        $oniomtopdb->{PDB_NATOMS} = undef;
        $oniomtopdb->{PDB_TOTALRESIDUENUM} = undef;
        $oniomtopdb->{PDB_RESIDUEIDLIST} = [];
        $oniomtopdb->{PDB_RESIDUENAMELIST} = [];
        $oniomtopdb->{PDB_TOTALRESIDUENUM} = undef;
        
	# Gaussian ONIOM file info
	$oniomtopdb->{ONIOM_ATOMS}    = [];    # Element symbol from input file
        $oniomtopdb->{ONIOM_RESIDUE}    = [];
        $oniomtopdb->{ONIOM_RESIDUENUM}    = [];
        $oniomtopdb->{ONIOM_ATOMTYPE}   = [];
        $oniomtopdb->{ONIOM_ELEMENTNUM} = [];
        $oniomtopdb->{ONIOM_ELEMENT}  = [];
        $oniomtopdb->{ONIOM_PARTIALCHARG}  = [];
        $oniomtopdb->{ONIOM_CARTCOORD}   = [];
        $oniomtopdb->{ONIOM_ATOMDEFINE}  = [];
        $oniomtopdb->{ONIOM_ATOMFREEZEFLAG} = [];
        $oniomtopdb->{ONIOM_AtomLayer}    = [];
        $oniomtopdb->{ONIOM_ATOMINRESIDUE}  = [];
        $oniomtopdb->{ONIOM_ORIGINALLINE}   = [];
        $oniomtopdb->{ONIOM_LineNum}  = [];
        $oniomtopdb->{ONIOM_ATOMNUM}    = 0;

        
	bless($oniomtopdb, $class);
	return $oniomtopdb;
}

# set filename & spin then digest the file
sub analyze : method {
	my $oniomtopdb = shift;
	$oniomtopdb->digest();
	return;
}

## subroutines ##

sub digest {

my $oniomtopdb = shift;
my $debug = $oniomtopdb->{DEBUG};
my $gfile = $oniomtopdb->{G_FILE_NAME};
my $pdbfile = $oniomtopdb->{PDB_FILE_NAME};
my $stepnumber = $oniomtopdb->{STEP_NUMBER};
my $outputfile = $oniomtopdb->{OUTPUT_FILE_NAME};
my $isLogFile = $oniomtopdb->{IS_LOG_FILE};
my $isJobFile = $oniomtopdb->{IS_JOB_FILE};

# flags & counters from ONIOMGEOM package
my $atomindex;
my $structindicator;
my $lineindex;
my $i;
my $j;
my $distance;
my $mismatch;

# flags & counters from PDB package
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
my $curatmtyp;
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
my $coreatomindex;
my $amberffdum  = ESPT::AMBERFF->new($debug);

# flags & counters from PDB package defined for ONIOMtoPDB package

my $pdb_atom_num;
my $pdb_line_counter;

# Read in PDB file

# open pdb filename for reading or display error
print "Open $pdbfile file for process.\n" if $debug >= 1;

open(PDBFILE,$pdbfile) || die "Could not read $pdbfile\n$!\n";

# set atom and residues counters as 0
$counter = 0; 
$resindex = -1;
$pdb_line_counter = 0;
# Process PDB file for ONIOM input file
while (<PDBFILE>) {
  print $_ if $debug >= 6;
  
  $pdb_line_counter++;
  
  if (/^ATOM\s+\d+/ or /^HETATM\s+\d+/) {
    $oniomtopdb->{PDB_ORIGINALLINE}[$counter] = $_;
    $oniomtopdb->{PDB_LINENUM}[$counter] = $pdb_line_counter;
    
    if (/^ATOM\s+\d+/) {
#       /^ATOM\s+\d+\s+(\S+)\s+(\S+)\s+(.{6})\s+(-?\d+.\d+)\s+(-?\d+.\d+)\s+(-?\d+.\d+)\s+/;
       /^ATOM\s+\d+\s+(\S+)\s+(\S+)\s+(\d+)\s+(-?\d+.\d+)\s+(-?\d+.\d+)\s+(-?\d+.\d+)\s+/;
       print "Atom is $1, residue is $2, Coordinates are $4, $5, $6,\n" if $debug >= 4;
       $curatomname  = $1;
       $curresidname = $2;
       $curresidnum  = $3;
       $curX = $4;
       $curY = $5;
       $curZ = $6;

       if ( $curatomname =~ /^[0-9]/) {
         print "Found nonstandard atom name $curatomname in residue $curresidname.\n" if $debug >= 1;
         $curatomname =~ s/^([0-9])(\S+)/$2$1/;
         print "New atom name $curatomname to meet AMBER standard name.\n" if $debug >= 1;
         }

     } else { 
#       /^HETATM\s+\d+\s+(\S+)\s+(\S{3})\s+(.{6})\s+(-?\d+.\d+)\s+(-?\d+.\d+)\s+(-?\d+.\d+)\s+/;
       /^HETATM\s+\d+\s+(\S+)\s+(\S+)\s+(\d+)\s+(-?\d+.\d+)\s+(-?\d+.\d+)\s+(-?\d+.\d+)\s+/;
       print "Atom is $1, residue is $2, Coordinates are $3, $4, $5,\n" if $debug >= 4;
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
       $atomdefined = 1;
       } else {
       # In this case, the residue information is not available from force file library file
         if ($debug >= 2) {
            printf "Residue %6s does not exist in available force field. Atom with name %5s may not be treated.\n", $curresidname,$curatomname;
           }
        $curatmtyp = "UDF";
        $curelenum = "UDF";
        $atomdefined = 0;
        $curelenum = 1;
       } #else


     $oniomtopdb->{PDB_ATOMS}[$counter] = $curatomname;
     $oniomtopdb->{PDB_RESIDUE}[$counter] = $curresidname;
     $oniomtopdb->{PDB_RESIDUENUM}[$counter] = $curresidnum;
     $oniomtopdb->{PDB_ELEMENTNUM}[$counter] = $curelenum;
     $oniomtopdb->{PDB_ELEMENT}[$counter] = $oniomtopdb->atomconvert($curelenum);
     $oniomtopdb->{PDB_ATOMDEFINE}[$counter] = $atomdefined;
     $oniomtopdb->{PDB_CARTCOORD}[3*$counter] = $curX;
     $oniomtopdb->{PDB_CARTCOORD}[3*$counter+1] = $curY;
     $oniomtopdb->{PDB_CARTCOORD}[3*$counter+2] = $curZ;
     $oniomtopdb->{PDB_ATOMINRESIDUE}[$counter] = $resindex;

     $counter++;
     
     #Check if encounter a new residue, if yes, store residue number and name. 
     
     if ( ($residref eq $curresidnum) and ( $resnameref eq $curresidname) )
      { next; }
     
     $resindex++;
     $oniomtopdb->{PDB_RESIDUEIDLIST}[$resindex] = $curresidnum;
     $oniomtopdb->{PDB_RESIDUENAMELIST}[$resindex] = $curresidname;
     $residref = $curresidnum;
     $resnameref = $curresidname;
     $oniomtopdb->{PDB_ATOMINRESIDUE}[$counter-1] = $resindex;
     print "Encounter new residue $resindex, name $curresidname, index $curresidnum.\n" if $debug >= 5;
  
   } # if (/^ATOM\s+\d+/ or /^HETATM\s+\d+/) 
   
   
  
 } #  while (<PDBFILE>) {
 
 close (PDBFILE);
 $resindex++;
 print "\nThere are $resindex residues in the given PDB files.\n" if $debug >= 0 ;
 print "\nThere are $counter atoms in the given PDB files.\n" if $debug >= 0 ;
 $oniomtopdb->{PDB_TOTALRESIDUENUM} =$resindex; 
 $oniomtopdb->{PDB_NATOMS} = $counter ;
 
 print "\nRead in PDB file done!\n" if $debug >= 0 ;


if ($isJobFile == 1) {
 # Read in Model GJF input file

 print "\nReading ONIOM input file $gfile...\n" if $debug >= 0;

 open(GFILE,$gfile) || die "Could not read $gfile\n$!\n";

 $structindicator = 0;
 $atomindex = 0;
 $lineindex = 0;

 while (<GFILE>) {  #NO1
    $lineindex++;
    if   ($structindicator == 0){
      if (/^\s*-*\d+\s+\d+\s+-*\d+\s+\d+/) {

        print $_ if $debug >= 2;
        $structindicator = 1;
        next;
      } # if (/^\s*-*\d+\s+\d+\s+-*\d+\s+\d+/) 
      next;
    }  # if   ($structindicator ==0){
    else {
      if (/^\s+$/ or /^$/){ last;}
 #     print $_;
     
      $oniomtopdb->{ONIOM_ORIGINALLINE}[$atomindex] = $_;
 #  Read in each atom information
      /^\s*([a-zA-Z]+)-([a-zA-Z0-9\*]+)-(-*\d*\.\d+)\s+(-*[01])\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+([HML])/;
      
       print "$1\t$2\t$3\t$4\t$5\t$6\t$7\t$8\t$atomindex\n" if $debug >= 3; 
       $oniomtopdb->{ONIOM_ATOMS}[$atomindex] = $1;
       $oniomtopdb->{ONIOM_ATOMTYPE}[$atomindex] = $2;
       $oniomtopdb->{ONIOM_PARTIALCHARG}[$atomindex] = $3;
       $oniomtopdb->{ONIOM_ATOMFREEZEFLAG}[$atomindex] = $4;
       $oniomtopdb->{ONIOM_CARTCOORD}[3*$atomindex] = $5;
       $oniomtopdb->{ONIOM_CARTCOORD}[3*$atomindex+1] = $6;
       $oniomtopdb->{ONIOM_CARTCOORD}[3*$atomindex+2] = $7;
       $oniomtopdb->{ONIOM_AtomLayer}[$atomindex] = $8;
       $oniomtopdb->{ONIOM_LineNum}[$atomindex] = $lineindex;
       $atomindex++;
    
    } # else

  } # while (<MODELGJFFILE>) { #NO1

 $oniomtopdb->{ONIOM_ATOMNUM} = $atomindex;

 print "\nThere are $atomindex atoms in file $gfile. \n" if $debug >= 0; 
 print "\nRead in ONIOM job file done!\n" if $debug >= 0;

 close (GFILE);
 } # if ($isJobFile == 1) {
   else {
   # Use subroutine to read in ONIOM log file

   $oniomtopdb->ReadInONIOMLog;
 
 }  # else {
 
# Compare number of atoms from two files
if ($oniomtopdb->{PDB_NATOMS} != $oniomtopdb->{ONIOM_ATOMNUM}){
  print "The number of atoms in $gfile (",$oniomtopdb->{ONIOM_ATOMNUM},
        ") is different from number of atoms in $pdbfile (",
        $oniomtopdb->{PDB_NATOMS},").\n" if $debug >= 0;
  die "Exit \n$!\n";

 }

# Sanity check for ONIOM file

print "\nSanity check for Model system...\n" if $debug >= 3;

for ($i=0; $i < $oniomtopdb->{ONIOM_ATOMNUM}; $i++) {
 
  if ($oniomtopdb->{ONIOM_ATOMS}[$i] ne $oniomtopdb->{PDB_ELEMENT}[$i]) {
    print "Atom ",$i+1," (line ",$oniomtopdb->{ONIOM_LineNum}[$i],
        " in file $gfile) is different \n from $pdbfile. (Line ",
        $oniomtopdb->{PDB_LINENUM}[$i],")\n" if $debug >= 3;
    print " Original line in $gfile is\n",$oniomtopdb->{ONIOM_ORIGINALLINE}[$i] if $debug >= 3;
    print " Original line in $pdbfile is\n",$oniomtopdb->{PDB_ORIGINALLINE}[$i] if $debug >= 3;
    print "\n" if $debug >= 3;
 
   }
  

 } # for ($i=0; $i < $oniom->{ATOMNUMModelGjf}; $i++)

print "Done with sanity check!\n" if $debug >= 3;




} # sub digest 

#Subroutine ReadInONIOMLog: read coordinate from given ONIOM log file

sub ReadInONIOMLog{

 my $oniomtopdb = shift;
 my $debug = $oniomtopdb->{DEBUG};
 my $gfile = $oniomtopdb->{G_FILE_NAME};
 my $totalstructure = 0;
 my $NAtoms = 0;
 my $structureindex;
 my $i;
 my ($logatom,$currentX,$currentY,$currentZ,@currentXYZ);
 my ($logatomsymbol,$logatomsymbolUP, $linecount);
 my ($line1, $line2, $templine);
 
 # open log name for reading or display error
 open(LOGFILE,$gfile) || die "Could not read $gfile\n$!\n";
 
 # grab Total number of atoms and find out how many structures in LOG file
 while (<LOGFILE>){
 	# skip blank lines
	next if /^$/;

        if (/^\s+NAtoms=\s*(\d+)/) {
           $NAtoms=$1;
 #           print "There are $NAtoms atoms.\n";
           }

       if ((/^\s+Z-Matrix\s+orientation:/) || (/^Standard\s+orientation:/)) {
            $totalstructure++;
           }          
        
 } # while (<LOGFILE>){
 
 $oniomtopdb->{ONIOM_ATOMNUM} = $NAtoms;
 
 print "\nThere are $NAtoms atoms in $gfile file.\n" if $debug >= 1;

 print "\nThere are $totalstructure structures in $gfile.\n" if $debug >= 1;

 close (LOGFILE);
 
 
 if ($oniomtopdb->{STEP_NUMBER} == 0) {
      $oniomtopdb->{STEP_NUMBER} = $totalstructure;
    } else {
      print "\nUpon request, geometry No. $oniomtopdb->{STEP_NUMBER} from $gfile. will be extracted.\n"  if $debug >= 0;
    
    }

 # If the user ask for structure that does not exist, complain.
 if ($oniomtopdb->{STEP_NUMBER} > $totalstructure) {
   printf "There are %4d of structures in the log file.\n",$totalstructure  if $debug >= 0;
   printf "You asked for No. %4d of structures!\n",$oniomtopdb->{STEP_NUMBER}  if $debug >= 0;
   die "Please double check!";
  }  # if ($log->{STRUCTURENUMBER} > $totalstructure) {


# Read in partial charges of each atoms from log file
# from section starting as 
#  Symbolic Z-matrix:
# Charge =-11 Multiplicity = 1 for low   level calculation on real  system.
# Charge =  1 Multiplicity = 1 for high  level calculation on model system.
# Charge =  1 Multiplicity = 1 for low   level calculation on model system.
# N-N3-0.096600        -1    41.507    23.202    52.275   L 
# H-H-0.216500         -1    40.56     22.841    52.249   L 
# H-H-0.216500         -1    41.687    23.653    51.384   L 

open(LOGFILE,$gfile) || die "Could not read $gfile\n$!\n";


# Process ONIOM Gaussian output file 

 $line1 = "";
 $line2 = "";
 $templine = "";


while ($templine = <LOGFILE>) {
  
   $line1 = $line2; 
   $line2 = $templine;
   
   # looking for the starting of atom list in given log file.
   if (( $line1 =~ /^\s*Charge\s+=\s+-*\d+\s+Multiplicity/) and (not($line2 =~ /^\s*Charge\s+=\s+-*\d+\s+Multiplicity/))) {
    
    $line1 = $line2;
    #  $line1 =~ /^(\S\S\S)\s+/;
    print "\nStart reading atom partial charges from Gaussian output file $gfile \nat line:\n$line1" if $debug >= 1;
    
    for ($i=0; $i < $NAtoms; $i++)  {
      $line1 =~ /\s*([a-zA-Z]+)-([a-zA-Z0-9\*]+)-(-*\d*.\d*)\s+-*\d\s+-*\d*\.\d*\s+-*\d*\.\d*\s+-*\d*\.\d*\s+([HMLhml])/;
      $oniomtopdb->{ONIOM_PARTIALCHARG}[$i] = $3;
      if ( not(defined $oniomtopdb->{ONIOM_PARTIALCHARG}[$i]) ) {
         $oniomtopdb->{ONIOM_PARTIALCHARG}[$i] = 0.0;
       }
      print "Reading atom number ",$i+1,": Name ($1), type ($2), charge ($3).\n" if $debug >= 2;
      $line1 = <LOGFILE>;
     } #  for ($i=1; $i < $NAtoms; $i++)  {

   last;
   } # if (( $line1 =~ 
 } # while (<GOUTFILE>)


 close (LOGFILE);


 open(LOGFILE,$gfile) || die "Could not read $gfile\n$!\n";

 $structureindex = 0;
 $linecount = 0;

 while (<LOGFILE>) {  #No2
      # Count line number
      $linecount++;
     
      if ((/^\s*Z-Matrix\s+orientation:/) || (/^Standard\s+orientation:/)) { #NO3
            $structureindex++;
            if ( $structureindex > $oniomtopdb->{STEP_NUMBER} ) {last;}
            if ( $structureindex < $oniomtopdb->{STEP_NUMBER} ) { 
              # Skip this structure
              for($i=1;$i<=$NAtoms+4;$i++) {
                 <LOGFILE>;
                # Count line number
                $linecount++;
                 }
              next;
            } else { #NO5
             @currentXYZ = ();

             # Count line number
             $linecount++;
             $_=<LOGFILE>;

             if ( not(/^\s+-+/)) {
              print "Line($_)\n"  if $debug >= 2;
              print "Gaussian log format may change. Pay attention to output structure(s).\n"  if $debug >= 2;
              }

             # Count line number
             $linecount++;
             $_ = <LOGFILE>;

             if ( not(/^\s+Center/) ) {
              print "Line($_)\n"  if $debug >= 2;
              print "Gaussian log format may change. Pay attention to output structure(s).\n"  if $debug >= 2;
              }

             # Count line number
             $linecount++;
             $_ = <LOGFILE>;

             if ( not(/^\s+Number\s+Number\s+Type\s+X\s+Y\s+Z/) ) {
              print "Line($_)\n"  if $debug >= 2;
              print "Gaussian log format may change. Pay attention to output structure(s).\n"  if $debug >= 2;
              }

             # Count line number
             $linecount++;
             $_=<LOGFILE>;
             if ( not(/^\s+-+/)) {
              print "Line($_)\n"  if $debug >= 2;
              print "Gaussian log format may change. Pay attention to output structure(s).\n" if $debug >= 2;
              }
              
             for ($i = 0; $i < $NAtoms;$i++){ #NO6
                $linecount++;
                $_ = <LOGFILE>;
                $oniomtopdb->{ONIOM_ORIGINALLINE}[$i] = $_;
#                print $_;

                /^\s+.+\s+(\d+)\s+\d+\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)/;
                $logatom  = $1;
                $currentX = $2;
                $currentY = $3;
                $currentZ = $4;

                  if ($debug >= 5) {printf "Taking coordinates (%16.10f %16.10f %16.10f) from log file.\n",$2,$3,$4;}

                  $logatomsymbol=$oniomtopdb->atomconvert($logatom);

#		  change all the symbol to the upper case                    

                  $logatomsymbolUP = $logatomsymbol;
                  $logatomsymbolUP =~ tr/a-z/A-Z/;
                  
                  $oniomtopdb->{ONIOM_ATOMS}[$i] = $logatomsymbolUP;
                  $oniomtopdb->{ONIOM_CARTCOORD}[3*$i] = $currentX;
                  $oniomtopdb->{ONIOM_CARTCOORD}[3*$i+1] = $currentY;
                  $oniomtopdb->{ONIOM_CARTCOORD}[3*$i+2] = $currentZ;
                  $oniomtopdb->{ONIOM_LineNum}[$i] = $linecount;
                  

                } #NO6
              
            } #NO5
             

           } #NO3          
      

 } # while (<LOGFILE>) { #No2


} #sub ReadInONIOMLog{


#Subroutine WriteOutputPDB: write PDB file with coordinates from ONIOM job

sub writeoutputpdb{

 my $oniomtopdb = shift;
 my $debug = $oniomtopdb->{DEBUG};
 my $wholepdbfile  = $oniomtopdb->{WHOLE_PDB_FILE}  ;
 my $pdbfile  = $oniomtopdb->{PDB_FILE_NAME}   ;
 my $outpdbfile  = $oniomtopdb->{OUTPUT_FILE_NAME} ;
 my $atomnum  = $oniomtopdb->{PDB_NATOMS}   ;
 my $i;
 my $currline;
 my ($atom,$atomtyp,$x,$y,$z,$roundx,$roundy,$roundz);
 my $atomindex;
 my $OutputCharge = $oniomtopdb->{OutputCharge};
 my $ChargeDec = $oniomtopdb->{ChargeDec};
 my ($charge, $roundcharge, $roundoccup);
 my $OutputOccupancy = $oniomtopdb->{OutputOccupancy};
 my $occup = $oniomtopdb->{Occupancy};
 my $addline;

 
 open(PDBFILE,$pdbfile) || die "Could not read $pdbfile\n$!\n";
 
 print "\nOpen $outpdbfile to write new PDB input file...\n" if $debug >= 0;

 open(OUTPUTPDBFILE,">$outpdbfile") || die "Could not open $outpdbfile\n$!\n";
 
 $atomindex = 0;
 
 while (<PDBFILE>) {  
  
  # If this is coordinat line (begin with ATOM or HETATM)   
  if (/^ATOM\s+\d+/ or /^HETATM\s+\d+/) {
    $currline = $_;
    $x = $oniomtopdb->{ONIOM_CARTCOORD}[3*$atomindex];
    $y = $oniomtopdb->{ONIOM_CARTCOORD}[3*$atomindex+1];
    $z = $oniomtopdb->{ONIOM_CARTCOORD}[3*$atomindex+2];
#    $roundx = sprintf ("%.3f",$x);
#    $roundy = sprintf ("%.3f",$y);
#    $roundz = sprintf ("%.3f",$z);
    
# Update to treat different length of digit better
# now replacing original coordinate by new coordinate
# with different length of digits are treated better:
# For example:
#  9.934 --> 10.1223
# 99.999 --> 100.032
    
    if ( $x >= 100.000 ) {
        $roundx = sprintf ("     %.3f",$x); 
       } elsif ( $x >= 10.000 ) {
        $roundx = sprintf ("      %.3f",$x); 
       } elsif ( $x >= 0.000 ) {
        $roundx = sprintf ("       %.3f",$x); 
       } elsif ( $x > -10.000 ) {
        $roundx = sprintf ("      %.3f",$x); 
       } elsif ( $x > -100.000 ) {
        $roundx = sprintf ("     %.3f",$x); 
       } else {
        $roundx = sprintf ("    %.3f",$x); 
       }
       
    if ( $y >= 100.000 ) {
        $roundy = sprintf (" %.3f",$y); 
       } elsif ( $y >= 10.000 ) {
        $roundy = sprintf ("  %.3f",$y); 
       } elsif ( $y >= 0.000 ) {
        $roundy = sprintf ("   %.3f",$y); 
       } elsif ( $y > -10.000 ) {
        $roundy = sprintf ("  %.3f",$y); 
       } elsif ( $y > -100.000 ) {
        $roundy = sprintf (" %.3f",$y); 
       } else {
        $roundy = sprintf ("%.3f",$y); 
       }
       
       
    if ( $z >= 100.000 ) {
        $roundz = sprintf (" %.3f",$z); 
       } elsif ( $z >= 10.000 ) {
        $roundz = sprintf ("  %.3f",$z); 
       } elsif ( $z >= 0.000 ) {
        $roundz = sprintf ("   %.3f",$z); 
       } elsif ( $z > -10.000 ) {
        $roundz = sprintf ("  %.3f",$z); 
       } elsif ( $z > -100.000 ) {
        $roundz = sprintf (" %.3f",$z); 
       } else {
        $roundz = sprintf ("%.3f",$z); 
       }
       
       
   
    print "Before substitution the line is \n",$currline if $debug >=6 ;

#    $currline =~ s/\s(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)/ $roundx  $roundy  $roundz/;

    $currline =~ s/\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)/$roundx$roundy$roundz/;

    print "After substitution the line is \n",$currline if $debug >=6;
    
    # Output partial charge and/or given occupancy
    if ($OutputCharge == 1)  {
        $charge = $oniomtopdb->{ONIOM_PARTIALCHARG}[$atomindex];
        $roundcharge = sprintf ("%.${ChargeDec}f",$charge); 
#        print "Charge $charge; RoundCharge $roundcharge.\n";
        $roundoccup  = sprintf ("%.2f",$occup); 
#        print "Occupancy $occup; Round Occupancy $roundoccup.\n";
        
        # if original PDB file does not have occupancy or B factor
        if ($currline =~ /\s+\d+\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s*$/) {

           if ($charge >= 0.0 ) {
                $addline = $roundoccup."  ".$roundcharge."\n";
               }
              else {
               $addline = $roundoccup." ".$roundcharge."\n";
              }
             
             if ($currline =~ /\s+\d+\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)$/) {
                  chop($currline);
                  $currline = $currline."  ".$addline;
                } elsif ($currline =~ /\s+\d+\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s$/) {
                  chop($currline);
                  $currline = $currline." ".$addline;
                } elsif ($currline =~ /\s+\d+\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s\s$/) {
                  chop($currline);
                  $currline = $currline.$addline;
                } else {
                 $currline =~ s/(\s+\d+\s+-*\d+\.\d+\s+-*\d+\.\d+\s+-*\d+\.\d+\s\s)\s*/$1$addline/;
                }
       
        # if original PDB file has occupancy but no B factor
         } elsif ($currline =~ /\s+\d+\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d*\.\d+)\s*$/) {
             if ($OutputOccupancy == 1) {
                $currline =~ s/(\s+\d+\s+-*\d+\.\d+\s+-*\d+\.\d+\s+-*\d+\.\d+\s+)\S+\s*/$1$roundoccup/;
              if ($charge >= 0.0 ) {
                 $currline = $currline."  ".$roundcharge."\n";
                 } else {
                 $currline = $currline." ".$roundcharge."\n";
                  }
              } else {
              
              if ($currline =~ /\s+\d+\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d*\.\d+)$/) {
                   chop($currline);
                   if ($charge >= 0.0 ) {
                        $currline = $currline."  ".$roundcharge."\n";
                      } else { 
                      $currline = $currline." ".$roundcharge."\n";
                     }
                 } elsif ($currline =~ /\s+\d+\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d*\.\d+)\s$/) {
                   chop($currline);
                   if ($charge >= 0.0 ) {
                        $currline = $currline." ".$roundcharge."\n";
                      } else { 
                      $currline = $currline.$roundcharge."\n";
                     }
                 } elsif ($currline =~ /\s+\d+\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d*\.\d+)\s\s$/) {
                    if ($charge >= 0.0 ) {
                         chop($currline);
                         $currline = $currline.$roundcharge."\n";
                       } else { 
                         chop($currline);
                         chop($currline);
                         $currline = $currline.$roundcharge."\n";
                      }
                 } else {
                    if ($charge >= 0.0 ) {
                         $currline =~ s/(\s+\d+\s+-*\d+\.\d+\s+-*\d+\.\d+\s+-*\d+\.\d+\s+-*\d*\.\d+\s\s)\s*/$1$roundcharge/;
                         $currline = $currline."\n";
                       } else { 
                         $currline =~ s/(\s+\d+\s+-*\d+\.\d+\s+-*\d+\.\d+\s+-*\d+\.\d+\s+-*\d*\.\d+\s)\s*/$1$roundcharge/;
                         $currline = $currline."\n";
                      }
                 }
              
              }

        # if original PDB file has both occupancy and B factor
          } elsif ($currline =~ /\s+\d+\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d*.\d+)\s+(-*\d*.\d+)\s*$/) {
#             print "I should be here 3\n";
             if ($OutputOccupancy == 1) {
                $currline =~ s/(\s+\d+\s+-*\d+\.\d+\s+-*\d+\.\d+\s+-*\d+\.\d+\s+)(-*\d*.\d+)(\s+-*\d*\.\d+)/$1$roundoccup$3/;
             }
                if ($charge >= 0.0 ) {
                   $currline =~ s/(\s+\d+\s+-*\d+\.\d+\s+-*\d+\.\d+\s+-*\d+\.\d+\s+-*\d*.\d+\s+)-*\d*\.\d+/$1$roundcharge/;
                  } else {
                   $currline =~ s/(\s+\d+\s+-*\d+\.\d+\s+-*\d+\.\d+\s+-*\d+\.\d+\s+-*\d*.\d+\s+)\s-*\d*\.\d+/$1$roundcharge/;
                  }
                
           }
         
      } elsif ($OutputOccupancy == 1) {
      
        if ($currline =~ /\s+\d+\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s*$/) {
           chop($currline);
           $currline = $currline."  ".$roundoccup."\n";
         } elsif ($currline =~ /\s+\d+\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d*.\d+)\s*$/) {
                $currline =~ s/(\s+\d+\s+-*\d+\.\d+\s+-*\d+\.\d+\s+-*\d+\.\d+\s+)(-*\d*.\d+)/$1$roundoccup/;
           } elsif ($currline =~ /\s+\d+\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d*.\d+)\s+(-*\d*.\d+)\s*$/) {
                $currline =~ s/(\s+\d+\s+-*\d+\.\d+\s+-*\d+\.\d+\s+-*\d+\.\d+\s+)(-*\d*.\d+)(\s+-*\d*.\d+)/$1$roundoccup$3/;
            }
      } # elsif ($OutputOccupancy == 1) {
    

    print OUTPUTPDBFILE $currline ;
   
    $atomindex++;
   } # if (/^ATOM\s+\d+/ or /^HETATM\s+\d+/) {
   
    # if this is not coordinat line
    elsif ($wholepdbfile ==1) {
      # Conserve whole PDB file in output PDB file
       print OUTPUTPDBFILE $_;
      }
 
  } 
 
 

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
