package ESPT::PDBCRD2ONIOM;

use base qw(ESPT::ESSfile);
use ESPT::Glib 0.01;
use ESPT::AMBERFF;
use strict;
use warnings;

our $VERSION = '0.1';

=head1 NAME

ESPT::PDBCRD2ONIOM - Use Oniom file and PDB file object to create new ONIOM file
                     with coordinates form PDB file.

=head1 SYNOPSIS

   use ESPT::PDBCRD2ONIOM;

   my $pdbcrd2oniom = PDBCRD2ONIOM->new();

=head1 DESCRIPTION

This module provides methods to put coordinates from PDB file to ONIOM file.

=cut

### Version History ###
# 0.1	Treat ONIOM job and PDB file


### To Do ###


## the object constructor **

sub new {
	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my $pdbcrd2oniom = ESPT::ESSfile->new();

	$pdbcrd2oniom->{TYPE} = "PDBCRD2ONIOM";

	# program info
	$pdbcrd2oniom->{PROGRAM}     = "Gaussian";
	$pdbcrd2oniom->{VERSION}     = undef;
	$pdbcrd2oniom->{REVISION}    = undef;
	$pdbcrd2oniom->{COMPILE}     = undef;
	$pdbcrd2oniom->{COMPILEDATE} = undef;
	

        # Operation info
        $pdbcrd2oniom->{G_FILE_NAME}      = undef;
	$pdbcrd2oniom->{PDB_FILE_NAME}    = undef;
        $pdbcrd2oniom->{OUTPUT_FILE_NAME} = undef;
        $pdbcrd2oniom->{IS_LOG_FILE}      = 0;
        $pdbcrd2oniom->{IS_JOB_FILE}     = 0;


	# PDB file info
	$pdbcrd2oniom->{PDB_ATOMS} = [];
        $pdbcrd2oniom->{PDB_RESIDUE} = [];
        $pdbcrd2oniom->{PDB_RESIDUENUM} = [];
        $pdbcrd2oniom->{PDB_ATOMTYPE} = [];
        $pdbcrd2oniom->{PDB_ELEMENTNUM} = [];
        $pdbcrd2oniom->{PDB_ELEMENT} = [];
        $pdbcrd2oniom->{PDB_PARTIALCHARG} = [];
        $pdbcrd2oniom->{PDB_CARTCOORD} = [];
        $pdbcrd2oniom->{PDB_ATOMDEFINE} = [];
        $pdbcrd2oniom->{PDB_ATOMFREEZEFLAG} = [];
        $pdbcrd2oniom->{PDB_ATOMINRESIDUE} = [];
        $pdbcrd2oniom->{PDB_ORIGINALLINE} = [];
        $pdbcrd2oniom->{PDB_LINENUM} = [];
        $pdbcrd2oniom->{PDB_NATOMS} = undef;
        $pdbcrd2oniom->{PDB_TOTALRESIDUENUM} = undef;
        $pdbcrd2oniom->{PDB_RESIDUEIDLIST} = [];
        $pdbcrd2oniom->{PDB_RESIDUENAMELIST} = [];
        $pdbcrd2oniom->{PDB_TOTALRESIDUENUM} = undef;
        
	# Gaussian ONIOM file info
	$pdbcrd2oniom->{ONIOM_ATOMS}    = [];    # Element symbol from input file
        $pdbcrd2oniom->{ONIOM_RESIDUE}    = [];
        $pdbcrd2oniom->{ONIOM_RESIDUENUM}    = [];
        $pdbcrd2oniom->{ONIOM_ATOMTYPE}   = [];
        $pdbcrd2oniom->{ONIOM_ELEMENTNUM} = [];
        $pdbcrd2oniom->{ONIOM_ELEMENT}  = [];
        $pdbcrd2oniom->{ONIOM_PARTIALCHARG}  = [];
        $pdbcrd2oniom->{ONIOM_CARTCOORD}   = [];
        $pdbcrd2oniom->{ONIOM_ATOMDEFINE}  = [];
        $pdbcrd2oniom->{ONIOM_ATOMFREEZEFLAG} = [];
        $pdbcrd2oniom->{ONIOM_AtomLayer}    = [];
        $pdbcrd2oniom->{ONIOM_ATOMINRESIDUE}  = [];
        $pdbcrd2oniom->{ONIOM_ORIGINALLINE}   = [];
        $pdbcrd2oniom->{ONIOM_LineNum}  = [];
        $pdbcrd2oniom->{ONIOM_ATOMNUM}    = 0;

        
	bless($pdbcrd2oniom, $class);
	return $pdbcrd2oniom;
}

# set filename  then digest the file
sub analyze : method {
	my $pdbcrd2oniom = shift;
	$pdbcrd2oniom->digest();
	return;
}

## subroutines ##

sub digest {

my $pdbcrd2oniom = shift;
my $debug = $pdbcrd2oniom->{DEBUG};
my $gfile = $pdbcrd2oniom->{G_FILE_NAME};
my $pdbfile = $pdbcrd2oniom->{PDB_FILE_NAME};
my $outputfile = $pdbcrd2oniom->{OUTPUT_FILE_NAME};
my $isLogFile = $pdbcrd2oniom->{IS_LOG_FILE};
my $isJobFile = $pdbcrd2oniom->{IS_JOB_FILE};

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

# flags & counters from PDB package defined for pdbcrd2oniom package

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
  print $_ if $debug >= 5;
  
  $pdb_line_counter++;
  
  if (/^ATOM\s+\d+/ or /^HETATM\s+\d+/) {
    $pdbcrd2oniom->{PDB_ORIGINALLINE}[$counter] = $_;
    $pdbcrd2oniom->{PDB_LINENUM}[$counter] = $pdb_line_counter;
    
    if (/^ATOM\s+\d+/) {
#       /^ATOM\s+\d+\s+(\S+)\s+(\S+)\s+(.{6})\s+(-?\d+.\d+)\s+(-?\d+.\d+)\s+(-?\d+.\d+)\s+/;
       /^ATOM\s+\d+\s+(\S+)\s+(\S+)\s+(\d+)\s+(-?\d+.\d+)\s+(-?\d+.\d+)\s+(-?\d+.\d+)\s+/;
       print "Atom is $1, residue is $2, Coordinates are $4, $5, $6,\n" if $debug >= 3;
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
       print "Atom is $1, residue is $2, Coordinates are $3, $4, $5,\n" if $debug >= 3;
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
        if ( $debug >= 1) {
          printf "Residue %6s does not exist. Atom with name %5s cannot be treated.\n", $curresidname,$curatomname;
        }
        $curatmtyp = "UDF";
        $curelenum = "UDF";
        $atomdefined = 0;
        $curelenum = 1;
       } #else


     $pdbcrd2oniom->{PDB_ATOMS}[$counter] = $curatomname;
     $pdbcrd2oniom->{PDB_RESIDUE}[$counter] = $curresidname;
     $pdbcrd2oniom->{PDB_RESIDUENUM}[$counter] = $curresidnum;
     $pdbcrd2oniom->{PDB_ELEMENTNUM}[$counter] = $curelenum;
     $pdbcrd2oniom->{PDB_ELEMENT}[$counter] = $pdbcrd2oniom->atomconvert($curelenum);
     $pdbcrd2oniom->{PDB_ATOMDEFINE}[$counter] = $atomdefined;
     $pdbcrd2oniom->{PDB_CARTCOORD}[3*$counter]   = $curX;
     $pdbcrd2oniom->{PDB_CARTCOORD}[3*$counter+1] = $curY;
     $pdbcrd2oniom->{PDB_CARTCOORD}[3*$counter+2] = $curZ;
     $pdbcrd2oniom->{PDB_ATOMINRESIDUE}[$counter] = $resindex;

     $counter++;
     
     #Check if encounter a new residue, if yes, store residue number and name. 
     
     if ( ($residref eq $curresidnum) and ( $resnameref eq $curresidname) )
      { next; }
     
     $resindex++;
     $pdbcrd2oniom->{PDB_RESIDUEIDLIST}[$resindex] = $curresidnum;
     $pdbcrd2oniom->{PDB_RESIDUENAMELIST}[$resindex] = $curresidname;
     $residref = $curresidnum;
     $resnameref = $curresidname;
     $pdbcrd2oniom->{PDB_ATOMINRESIDUE}[$counter-1] = $resindex;
     print "Encounter new residue $resindex, name $curresidname, index $curresidnum.\n" if $debug >= 4;
  
   } # if (/^ATOM\s+\d+/ or /^HETATM\s+\d+/) 
   
   
  
 } #  while (<PDBFILE>) {
 
 close (PDBFILE);
 $resindex++;
 print "\nThere are $resindex residues in the given PDB files.\n\n" if $debug >= 0 ;
 print "There are $counter atoms in the given PDB files.\n\n" if $debug >= 0 ;
 $pdbcrd2oniom->{PDB_TOTALRESIDUENUM} =$resindex; 
 $pdbcrd2oniom->{PDB_NATOMS} = $counter ;
 
 print "\nRead in PDB file done!\n" if $debug >= 0 ;


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
     
      $pdbcrd2oniom->{ONIOM_ORIGINALLINE}[$atomindex] = $_;
 #  Read in each atom information
      /^\s*([a-zA-Z]+)-([a-zA-Z0-9\*]+)-(-*\d*\.\d+)\s+(-*[01])\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+([HMLhml])/;
      
       print "$1\t$2\t$3\t$4\t$5\t$6\t$7\t$8\t$atomindex\n" if $debug >= 3; 
       $pdbcrd2oniom->{ONIOM_ATOMS}[$atomindex] = $1;
       $pdbcrd2oniom->{ONIOM_ATOMTYPE}[$atomindex] = $2;
       $pdbcrd2oniom->{ONIOM_PARTIALCHARG}[$atomindex] = $3;
       $pdbcrd2oniom->{ONIOM_ATOMFREEZEFLAG}[$atomindex] = $4;
       $pdbcrd2oniom->{ONIOM_CARTCOORD}[3*$atomindex] = $5;
       $pdbcrd2oniom->{ONIOM_CARTCOORD}[3*$atomindex+1] = $6;
       $pdbcrd2oniom->{ONIOM_CARTCOORD}[3*$atomindex+2] = $7;
       $pdbcrd2oniom->{ONIOM_AtomLayer}[$atomindex] = $8;
       $pdbcrd2oniom->{ONIOM_LineNum}[$atomindex] = $lineindex;
       $atomindex++;
    
    } # else

  } # while (<MODELGJFFILE>) { #NO1

 $pdbcrd2oniom->{ONIOM_ATOMNUM} = $atomindex;

 print "\nThere are $atomindex atoms in file $gfile. \n" if $debug >= 0; 
 print "\nRead in ONIOM job file done!\n" if $debug >= 0;

 close (GFILE);
 
# Compare number of atoms from two files
if ($pdbcrd2oniom->{PDB_NATOMS} != $pdbcrd2oniom->{ONIOM_ATOMNUM}){
  print "\nThe number of atoms in $gfile (",$pdbcrd2oniom->{ONIOM_ATOMNUM},
        ") is different from number of atoms in $pdbfile (",
        $pdbcrd2oniom->{PDB_NATOMS},").\n" if $debug >= 0;
  die "Exit \n$!\n";

 }

# Sanity check for ONIOM file

print "\nSanity check for Model system...\n" if $debug >= 3;

for ($i=0; $i < $pdbcrd2oniom->{ONIOM_ATOMNUM}; $i++) {
 
  if ($pdbcrd2oniom->{ONIOM_ATOMS}[$i] ne $pdbcrd2oniom->{PDB_ELEMENT}[$i]) {
    print "Atom ",$i+1," (line ",$pdbcrd2oniom->{ONIOM_LineNum}[$i],
        " in file $gfile) is different \n from $pdbfile. (Line ",
        $pdbcrd2oniom->{PDB_LINENUM}[$i],")\n" if $debug >= 3;
    print " Original line in $gfile is\n",$pdbcrd2oniom->{ONIOM_ORIGINALLINE}[$i] if $debug >= 3;
    print " Original line in $pdbfile is\n",$pdbcrd2oniom->{PDB_ORIGINALLINE}[$i] if $debug >= 3;
    print "\n" if $debug >= 3;
 
   }
  

 } # for ($i=0; $i < $oniom->{ATOMNUMModelGjf}; $i++)

print "\nDone with sanity check!\n" if $debug >= 3;




} # sub digest 



#Subroutine WriteOutputPDB: write PDB file with coordinates from ONIOM job

sub writeoutput{

 my $pdbcrd2oniom = shift;
 my $debug = $pdbcrd2oniom->{DEBUG};
 my $gfile = $pdbcrd2oniom->{G_FILE_NAME};
 my $pdbfile  = $pdbcrd2oniom->{PDB_FILE_NAME}   ;
 my $outgjffile  = $pdbcrd2oniom->{OUTPUT_FILE_NAME} ;
 my $atomnum  = $pdbcrd2oniom->{PDB_NATOMS}   ;
 my $i;
 my $currline;
 my ($atom,$atomtyp,$x,$y,$z,$roundx,$roundy,$roundz);
 my $atomindex;

 print "\nReading ONIOM input file $gfile to create file $outgjffile...\n" if $debug >= 0;

 open(GFILE,$gfile) || die "Could not read $gfile\n$!\n";
 
 open(OUTPUTGJFFILE,">$outgjffile") || die "Could not open $outgjffile\n$!\n";

 $atomindex = $pdbcrd2oniom->{ONIOM_ATOMNUM};

 while (<GFILE>) {  #NO10

      print OUTPUTGJFFILE $_;

      if (/^\s*-*\d+\s+\d+\s+-*\d+\s+\d+/) { last;} 

   }  #NO10

 for ( $i = 0 ; $i < $atomindex; $i++)  { #NO11
      
      $currline = <GFILE>;

      print "Before substitution the line is \n",$currline if $debug >=4 ;
      
      $x = $pdbcrd2oniom->{PDB_CARTCOORD}[3*$i];
      $y = $pdbcrd2oniom->{PDB_CARTCOORD}[3*$i+1];
      $z = $pdbcrd2oniom->{PDB_CARTCOORD}[3*$i+2];
      

      $currline =~ 
           s/\s(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+/ $x    $y    $z   /;
             

      print "After substitution the line is \n",$currline if $debug >=4;

      print OUTPUTGJFFILE $currline ;

       } #NO11  

    while (<GFILE>) { print OUTPUTGJFFILE $_; }
 
 print "\nSuccessfully wrote file $outgjffile!\n" if $debug >= 0;


} #sub writeoutput{




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
