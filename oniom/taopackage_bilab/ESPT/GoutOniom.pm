package ESPT::GoutOniom;

our $VERSION = '0.01';

### Version History ###
# 0.01	Read in ONIOM Gaussian log file for atom list with atom types.
#

=head1 NAME

ESPT::GoutOniom - Read in ONIOM gaussian log file for parameters look up.

=head1 SYNOPSIS

    use ESPT::GoutOniom;


=head1 DESCRIPTION

This module read in ONIOM gaussian log file for parameters look up.
Information obtained includes input stream (atom list with atom type), link
atom type, missing paramters list.

=head1 ATTRIBUTES

=over 15

=item B<$file-E<gt>new()>

Creates a new AMBERFF object

=cut


use base qw(ESPT::ESSfile);
use ESPT::Glib 0.01;
use ESPT::AMBERPARM 0.01;
use strict;
use warnings;


### Version History ###
# 0.01	digest Gaussian ONIOM output file
#

### To Do List ###

##### NOTE #####
#
################

## the object constructor **


sub new {
	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my $goutoniom = ESPT::ESSfile->new();

	$goutoniom->{PROGRAM} = "GAUSSIAN";
	$goutoniom->{TYPE} = "Log";

	# Link 0 & Route commands
	$goutoniom->{ROUTE} = undef;
	$goutoniom->{KEYWORDS} = [];
        $goutoniom->{OUTPUTFILE} = "";


	# molecular info
	$goutoniom->{ATOMLIST} = [];
        $goutoniom->{RESIDUE} = [];
        $goutoniom->{ATOMTYPE} = [];
        $goutoniom->{LAYERLIST} = [];
        $goutoniom->{ELEMENTNUM} = [];
        $goutoniom->{ISLINKATOM} = [];
        $goutoniom->{LINKATOMTYPE} = [];
        $goutoniom->{ELEMENT} = [];
        $goutoniom->{PARTIALCHARG} = [];
        $goutoniom->{CARTCOORD} = [];
        $goutoniom->{ATOMDEFINE} = [];
        $goutoniom->{MISSPARMBONDLIST} = [];
        $goutoniom->{MISSPARMANGLELIST} = [];

        $goutoniom->{FOUNDSTRETCHL1} = [];
        $goutoniom->{FOUNDSTRETCHL2} = [];
        $goutoniom->{FOUNDSTRETCHPARML1} = [];
        $goutoniom->{FOUNDSTRETCHPARML2} = [];
        $goutoniom->{FOUNDSTRETCHNUM} = undef;

        $goutoniom->{FOUNDSTRETCHL1REDUCE} = [];
        $goutoniom->{FOUNDSTRETCHL2REDUCE} = [];
        $goutoniom->{FOUNDSTRETCHPARML1REDUCE} = [];
        $goutoniom->{FOUNDSTRETCHPARML2REDUCE} = [];
        $goutoniom->{FOUNDSTRETCHNUMREDUCE} = undef;
        
        
        $goutoniom->{FOUNDBENDL1} = [];
        $goutoniom->{FOUNDBENDL2} = [];
        $goutoniom->{FOUNDBENDL3} = [];
        $goutoniom->{FOUNDBENDPARML1} = [];
        $goutoniom->{FOUNDBENDPARML2} = [];
        $goutoniom->{FOUNDBENDNUM} = undef;

        $goutoniom->{FOUNDBENDL1REDUCE} = [];
        $goutoniom->{FOUNDBENDL2REDUCE} = [];
        $goutoniom->{FOUNDBENDL3REDUCE} = [];
        $goutoniom->{FOUNDBENDPARML1REDUCE} = [];
        $goutoniom->{FOUNDBENDPARML2REDUCE} = [];
        $goutoniom->{FOUNDBENDNUMREDUCE} = undef;

        $goutoniom->{MISSPARMBONDNUM} = 0;
        $goutoniom->{MISSPARMANGLENUM} = 0;
      
        
	$goutoniom->{CHARGE} = undef;
	$goutoniom->{EELEC} = undef;		# SCF electronic energy
	$goutoniom->{MULTIPLICITY} = undef;
	$goutoniom->{NATOMS} = undef;
        
        $goutoniom->{HCHARGE} = 0.0;
        $goutoniom->{MCHARGE} = 0.0;
        $goutoniom->{LCHARGE} = 0.0;
        $goutoniom->{RCHARGE} = 0.0;

        $goutoniom->{XDIPOLECAL} = 0.0;
        $goutoniom->{YDIPOLECAL} = 0.0;
        $goutoniom->{ZDIPOLECAL} = 0.0;
        $goutoniom->{TOTDIPOLECAL} = 0.0;


        $goutoniom->{RAW} = 0;
        $goutoniom->{PARMLOOKUP} = 0;
        
       
	bless($goutoniom, $class);
	return $goutoniom;
}


## methods ##

# set filename & spin then digest the file
sub analyze : method {
	my $goutoniom = shift;
	$goutoniom->prepare(@_);
	$goutoniom->digest();
	return;
}


sub digest {

my $goutoniom = shift;

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
my ($curatmtyp, $amberffdum, $curatomname, $curresidname, $curelenum, $curparchg, $curX);
my ($curY,  $curZ, $atomdefined, $line1, $line2, $templine, $isendofatomlist, $missbondcounter,$missanglecounter);

$debug = $goutoniom->{DEBUG} ;

# open filename for reading or display error
print "Open $goutoniom->{FILENAME} file for process.\n" if $debug >= 1;

open(GOUTFILE,$goutoniom->{FILENAME}) || die "Could not read $goutoniom->{FILENAME}\n$!\n";

# set atom counter as 0
$counter = 0; 

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
    print "Start reading atom list from Gaussian output file $goutoniom->{FILENAME} \nat line:\n$line2" if $debug >= 1;
    
    $counter = 1;

    
    while ($isendofatomlist == 0) {
    
      $line1 = $line2; 
      $line2 = <GOUTFILE>;
      
    
      $line1 =~ /\s*([a-zA-Z]+)-([a-zA-Z0-9\*]+)-(-*\d*.\d*)\s+-*\d\s+(-*\d*\.\d*)\s+(-*\d*\.\d*)\s+(-*\d*\.\d*)\s+([HMLhml])/;
      $goutoniom->{ATOMLIST}[$counter] = $1;
      $goutoniom->{ATOMTYPE}[$counter] = $2;
      $goutoniom->{LAYERLIST}[$counter] = $7;
      $goutoniom->{ISLINKATOM}[$counter] = 0;
      
      print "Reading atom number $counter: Name ($1), type ($2), layer ($7).\n" if $debug >= 1;
      
      if ( uc($7) eq "L") {$goutoniom->{LCHARGE} += $3;}
      if ( uc($7) eq "H") {$goutoniom->{HCHARGE} += $3;}
      if ( uc($7) eq "M") {$goutoniom->{MCHARGE} += $3;}
      $goutoniom->{RCHARGE} += $3;
      $goutoniom->{XDIPOLECAL} += $4 * $3 * $goutoniom->{DEBYECONVER};
      $goutoniom->{YDIPOLECAL} += $5 * $3 * $goutoniom->{DEBYECONVER};
      $goutoniom->{ZDIPOLECAL} += $6 * $3 * $goutoniom->{DEBYECONVER};
      
      
      # If this is a link atom

      if ($line1 =~ /\s*[a-zA-Z]+-[a-zA-Z0-9\*]+--*\d*.\d*\s+-*\d\s+-*\d*\.\d*\s+-*\d*\.\d*\s+-*\d*\.\d*\s+[HMLhml]\s+([a-zA-Z])-(\S+)/) {
         $goutoniom->{ISLINKATOM}[$counter] = 1;
         $goutoniom->{LINKATOMTYPE}[$counter] = $2;
         print " Atom number $counter is a link atom ($1), and has type ($2).\n" if $debug >= 1;
         
        } # if $line 1
      
      $counter++;
      
      if ( $line2 =~ /^\s*$/) {$isendofatomlist = 1;}

     } # while ($isendofatomlist != 0) 
    
    
   $goutoniom->{NATOMS} = $counter-1;
   print "There are $goutoniom->{NATOMS} atoms in gaussian output file $goutoniom->{FILENAME}.\n" if $debug >= 1;
   last;
   
   } # if (( $line1 =~ 
 
 } # while (<GOUTFILE>)

$goutoniom->{TOTDIPOLECAL} = sqrt($goutoniom->{XDIPOLECAL}**2+$goutoniom->{YDIPOLECAL}**2+$goutoniom->{ZDIPOLECAL}**2);

$missbondcounter = 0;
$missanglecounter = 0;

while ($templine = <GOUTFILE>) {

 if ($templine =~ /\s+Bondstretch\s+undefined\s+between\s+atoms\s+(\d+)\s+(\d+)/) {
 
   $goutoniom->{MISSPARMBONDLIST}[2*$missbondcounter] = $1;
   $goutoniom->{MISSPARMBONDLIST}[2*$missbondcounter+1] = $2;
   $missbondcounter++;
   print "Bond stretch parameter is missing between atom $1 and $2.\n" if $debug >= 1;
    
  } # if for missing bond parameter

 if ($templine =~ /\s+Angle\s+bend\s+undefined\s+between\s+atoms\s+(\d+)\s+(\d+)\s+(\d+)/) {
 
   $goutoniom->{MISSPARMANGLELIST}[3*$missanglecounter] = $1;
   $goutoniom->{MISSPARMANGLELIST}[3*$missanglecounter+1] = $2;
   $goutoniom->{MISSPARMANGLELIST}[3*$missanglecounter+2] = $3;
   $missanglecounter++;
   print "Bend parameter is missing for bond angle defined by atom $1, $2 and $3.\n" if $debug >= 1;
    
  } # if for missing bond parameter

 } # while ($templine = <GOUTFILE>) 
 
 $goutoniom->{MISSPARMBONDNUM} = $missbondcounter;
 $goutoniom->{MISSPARMANGLENUM} = $missanglecounter;

 # Print out missing parameters info when parmlookup is running
 if ($goutoniom->{PARMLOOKUP} == 1) {
    print "\nThere are $missbondcounter bond stretch parameters missing, and $missanglecounter bond angle bend missing.\n" if $debug >= 0;
  }
 
} # sub digest 

# Look up bond stretch and bend paramters

sub lookupparms {
 
 my $goutoniom = shift;
 my $debug = $goutoniom->{DEBUG} ;
 my $i;
 my $NAtoms; 
 my $amberparmdum;
 my $misstrechnum;
 my $misbendnum;
 my ($atomnum1,$atomnum2,$atomnum3);
 my ($atomtyp1,$atomtyp2,$atomtyp3);
 my @curparm;
 my $foundstrechnum = 0;
 my $foundbendnum = 0;
 
 $amberparmdum = ESPT::AMBERPARM->new();
 
 $misstrechnum = $goutoniom->{MISSPARMBONDNUM};

 print "Now look up missing bond stretch parameters...\n" if $debug >= 1;

 for ($i = 0; $i < $misstrechnum; $i++) {
   $atomnum1 = $goutoniom->{MISSPARMBONDLIST}[2*$i];
   $atomnum2 = $goutoniom->{MISSPARMBONDLIST}[2*$i+1];
   $atomtyp1 = $goutoniom->{ATOMTYPE}[$atomnum1];
   $atomtyp2 = $goutoniom->{ATOMTYPE}[$atomnum2];
   
   print "Looking up bond stetech for atom $atomnum1 ($atomtyp1) and $atomnum2 ($atomtyp2).\n" if $debug >=1;
   @curparm = $amberparmdum->bondstretchlookup($atomtyp1,$atomtyp2,$debug);
   
   $goutoniom->{FOUNDSTRETCHL1}[$foundstrechnum] = $atomtyp1;
   $goutoniom->{FOUNDSTRETCHL2}[$foundstrechnum] = $atomtyp2;
   $goutoniom->{FOUNDSTRETCHPARML1}[$foundstrechnum] = $curparm[0];
   $goutoniom->{FOUNDSTRETCHPARML2}[$foundstrechnum] = $curparm[1];
   $foundstrechnum++;
   
   if ( $goutoniom->{ISLINKATOM}[$atomnum1] == 1 ){
      $atomtyp1 = $goutoniom->{LINKATOMTYPE}[$atomnum1];
      
      print "Looking up bond stretch for link atom $atomnum1 ($atomtyp1) and $atomnum2 ($atomtyp2).\n" if $debug >=1;
      @curparm = $amberparmdum->bondstretchlookup($atomtyp1,$atomtyp2,$debug);
      $goutoniom->{FOUNDSTRETCHL1}[$foundstrechnum] = $atomtyp1;
      $goutoniom->{FOUNDSTRETCHL2}[$foundstrechnum] = $atomtyp2;
      $goutoniom->{FOUNDSTRETCHPARML1}[$foundstrechnum] = $curparm[0];
      $goutoniom->{FOUNDSTRETCHPARML2}[$foundstrechnum] = $curparm[1];
      $foundstrechnum++;
    }
   
   if ( $goutoniom->{ISLINKATOM}[$atomnum2] == 1 ){
      $atomtyp2 = $goutoniom->{LINKATOMTYPE}[$atomnum2];
      
      print "Looking up bond stretch for atom $atomnum1 ($atomtyp1) and link atom $atomnum2 ($atomtyp2).\n" if $debug >=1;
      @curparm = $amberparmdum->bondstretchlookup($atomtyp1,$atomtyp2,$debug);
      $goutoniom->{FOUNDSTRETCHL1}[$foundstrechnum] = $atomtyp1;
      $goutoniom->{FOUNDSTRETCHL2}[$foundstrechnum] = $atomtyp2;
      $goutoniom->{FOUNDSTRETCHPARML1}[$foundstrechnum] = $curparm[0];
      $goutoniom->{FOUNDSTRETCHPARML2}[$foundstrechnum] = $curparm[1];
      $foundstrechnum++;
    }
   
 } #  for ($i = 1; $i <= $NAtoms; $i++) 
 
 $goutoniom->{FOUNDSTRETCHNUM} = $foundstrechnum;
 
 print "\nThere are $foundstrechnum bond stretch parameter processed before reduction.\n" if  $debug >=0;


# Look up bond bending parameters

 $misbendnum = $goutoniom->{MISSPARMANGLENUM};

 print "Now look up missing bond angle bend parameters...\n" if $debug >= 1;

 for ($i = 0; $i < $misbendnum; $i++) {
   $atomnum1 = $goutoniom->{MISSPARMANGLELIST}[3*$i];
   $atomnum2 = $goutoniom->{MISSPARMANGLELIST}[3*$i+1];
   $atomnum3 = $goutoniom->{MISSPARMANGLELIST}[3*$i+2];
   $atomtyp1 = $goutoniom->{ATOMTYPE}[$atomnum1];
   $atomtyp2 = $goutoniom->{ATOMTYPE}[$atomnum2];
   $atomtyp3 = $goutoniom->{ATOMTYPE}[$atomnum3];
   
   print "Looking up bond angle bend for atom $atomnum1 ($atomtyp1), $atomnum2 ($atomtyp2) and $atomnum3 ($atomtyp3).\n" if $debug >=1;
   @curparm = $amberparmdum->bondbendlookup($atomtyp1,$atomtyp2,$atomtyp3,$debug);
   
   $goutoniom->{FOUNDBENDL1}[$foundbendnum] = $atomtyp1;
   $goutoniom->{FOUNDBENDL2}[$foundbendnum] = $atomtyp2;
   $goutoniom->{FOUNDBENDL3}[$foundbendnum] = $atomtyp3;
   $goutoniom->{FOUNDBENDPARML1}[$foundbendnum] = $curparm[0];
   $goutoniom->{FOUNDBENDPARML2}[$foundbendnum] = $curparm[1];
   $foundbendnum++;
   
   if ( $goutoniom->{ISLINKATOM}[$atomnum1] == 1 ){
      $atomtyp1 = $goutoniom->{LINKATOMTYPE}[$atomnum1];
      
      print "Looking up bond angle bend for atom $atomnum1 ($atomtyp1), $atomnum2 ($atomtyp2) and $atomnum3 ($atomtyp3).\n" if $debug >=1;
      @curparm = $amberparmdum->bondbendlookup($atomtyp1,$atomtyp2,$atomtyp3,$debug);
      $goutoniom->{FOUNDBENDL1}[$foundbendnum] = $atomtyp1;
      $goutoniom->{FOUNDBENDL2}[$foundbendnum] = $atomtyp2;
      $goutoniom->{FOUNDBENDL3}[$foundbendnum] = $atomtyp3;
      $goutoniom->{FOUNDBENDPARML1}[$foundbendnum] = $curparm[0];
      $goutoniom->{FOUNDBENDPARML2}[$foundbendnum] = $curparm[1];
      $foundbendnum++;
    }
   
   if ( $goutoniom->{ISLINKATOM}[$atomnum3] == 1 ){
      $atomtyp3 = $goutoniom->{LINKATOMTYPE}[$atomnum3];
      
      print "Looking up bond angle bend for atom $atomnum1 ($atomtyp1), $atomnum2 ($atomtyp2) and $atomnum3 ($atomtyp3).\n" if $debug >=1;
      @curparm = $amberparmdum->bondbendlookup($atomtyp1,$atomtyp2,$atomtyp3,$debug);
      $goutoniom->{FOUNDBENDL1}[$foundbendnum] = $atomtyp1;
      $goutoniom->{FOUNDBENDL2}[$foundbendnum] = $atomtyp2;
      $goutoniom->{FOUNDBENDL3}[$foundbendnum] = $atomtyp3;
      $goutoniom->{FOUNDBENDPARML1}[$foundbendnum] = $curparm[0];
      $goutoniom->{FOUNDBENDPARML2}[$foundbendnum] = $curparm[1];
      $foundbendnum++;
    }
   
 } #  for ($i = 1; $i <= $NAtoms; $i++) 
 
 $goutoniom->{FOUNDBENDNUM} = $foundbendnum;
 
 print "\nThere are $foundbendnum bond angle bend parameter processed before reduction.\n" if  $debug >=0;

 
} # sub lookupparms 


# Remove the redundant bond stretch and angle bend parameters.

sub parmreduction {

 my $goutoniom = shift;
 my $debug = $goutoniom->{DEBUG};
 my $amberparmdum;
 my $foundstretchnum;
 my $foundbendnum;
 my $foundstretchnumreduce = 0;
 my $foundbendnumreduce = 0;
 my $i;
 my $j;
 my $needtake;
 my ($atmtyp1,$atmtyp2,$atmtyp3);
 my ($cmpatmtyp1,$cmpatmtyp2,$cmpatmtyp3);
 

 $amberparmdum = ESPT::AMBERPARM->new();
  
 
 $foundstretchnum = $goutoniom->{FOUNDSTRETCHNUM};
 
 # Remove redundant bond stretch parameters.
 
 for ($i=0; $i<$foundstretchnum ;$i++) {
   
   $atmtyp1 = $goutoniom->{FOUNDSTRETCHL1}[$i];
   $atmtyp2 = $goutoniom->{FOUNDSTRETCHL2}[$i];
   
   $needtake = 1;
   for ($j=0; $j<$foundstretchnumreduce;$j++) {
     $cmpatmtyp1 = $goutoniom->{FOUNDSTRETCHL1REDUCE}[$j];
     $cmpatmtyp2 = $goutoniom->{FOUNDSTRETCHL2REDUCE}[$j];
#     print "Compare bond stretch parameter $atmtyp1-$atmtyp2 and $cmpatmtyp1-$cmpatmtyp2 \n" if $debug >= 1;

     if ( (bondmatch($atmtyp1,$atmtyp2,$cmpatmtyp1,$cmpatmtyp2)) == 1)
       {
        print "Bond $atmtyp1-$atmtyp2 and $cmpatmtyp1-$cmpatmtyp2 are the same.\n" if $debug >= 1;
        $needtake = 0;
        last;
       }
     } # for foundstretchnumreduce
   
   if ($needtake == 1) {
     $goutoniom->{FOUNDSTRETCHL1REDUCE}[$foundstretchnumreduce] = $atmtyp1;
     $goutoniom->{FOUNDSTRETCHL2REDUCE}[$foundstretchnumreduce] = $atmtyp2;
     $goutoniom->{FOUNDSTRETCHPARML1REDUCE}[$foundstretchnumreduce] = $goutoniom->{FOUNDSTRETCHPARML1}[$i];
     $goutoniom->{FOUNDSTRETCHPARML2REDUCE}[$foundstretchnumreduce] = $goutoniom->{FOUNDSTRETCHPARML2}[$i];
     $foundstretchnumreduce++;
    } # if
   
  } # for foundstretchnum
 
  $goutoniom->{FOUNDSTRETCHNUMREDUCE} = $foundstretchnumreduce;
  
  print "\nThere are $foundstretchnumreduce set of unique bone stretch parameters.\n" if $debug >=0;



# Remove redundant bond bend parameters.

 $foundbendnum    = $goutoniom->{FOUNDBENDNUM};

  
 for ($i=0; $i<$foundbendnum ;$i++) {
   
   $atmtyp1 = $goutoniom->{FOUNDBENDL1}[$i];
   $atmtyp2 = $goutoniom->{FOUNDBENDL2}[$i];
   $atmtyp3 = $goutoniom->{FOUNDBENDL3}[$i];
   
   $needtake = 1;
   for ($j=0; $j<$foundbendnumreduce;$j++) {
     $cmpatmtyp1 = $goutoniom->{FOUNDBENDL1REDUCE}[$j];
     $cmpatmtyp2 = $goutoniom->{FOUNDBENDL2REDUCE}[$j];
     $cmpatmtyp3 = $goutoniom->{FOUNDBENDL3REDUCE}[$j];
     
#     print "Compare bond bend parameter $atmtyp1-$atmtyp2-$atmtyp3 and $cmpatmtyp1-$cmpatmtyp2-$cmpatmtyp3 \n" if $debug >= 1;

     if ( (anglematch($atmtyp1,$atmtyp2,$atmtyp3,$cmpatmtyp1,$cmpatmtyp2,$cmpatmtyp3)) == 1)
       {
        print "Bond angle $atmtyp1-$atmtyp2-$atmtyp3 and $cmpatmtyp1-$cmpatmtyp2-$cmpatmtyp3 are the same.\n" if $debug >= 1;
        $needtake = 0;
        last;
       }
     } # for foundbendnumreduce
   
   if ($needtake == 1) {
     $goutoniom->{FOUNDBENDL1REDUCE}[$foundbendnumreduce] = $atmtyp1;
     $goutoniom->{FOUNDBENDL2REDUCE}[$foundbendnumreduce] = $atmtyp2;
     $goutoniom->{FOUNDBENDL3REDUCE}[$foundbendnumreduce] = $atmtyp3;
     $goutoniom->{FOUNDBENDPARML1REDUCE}[$foundbendnumreduce] = $goutoniom->{FOUNDBENDPARML1}[$i];
     $goutoniom->{FOUNDBENDPARML2REDUCE}[$foundbendnumreduce] = $goutoniom->{FOUNDBENDPARML2}[$i];
     $foundbendnumreduce++;
    } # if
   
  } # for foundbendnum
 
  $goutoniom->{FOUNDBENDNUMREDUCE} = $foundbendnumreduce;
  
  print "\nThere are $foundbendnumreduce set of unique bone bend parameters.\n" if $debug >=0;
  
  

} # lookupparms


# Write out missing parameter 

sub outputparm{

 my $goutoniom = shift;
 my $filename  = shift;
 my $debug = $goutoniom->{DEBUG};
 my $bondstrechnumfinal;
 my $bondbendnumfinal;
 my $i;
 my $j;
 
 
 open(OUTPUTFILE,">$filename") || die "Could not open file $filename\n$!\n";
 
 $bondstrechnumfinal = $goutoniom->{FOUNDSTRETCHNUMREDUCE};
 
 for ($i=0;$i<$bondstrechnumfinal;$i++) {
  
  printf OUTPUTFILE "Hrmstr1 %4s %4s %10.4f %10.4f\n", $goutoniom->{FOUNDSTRETCHL1REDUCE}[$i],$goutoniom->{FOUNDSTRETCHL2REDUCE}[$i],
       $goutoniom->{FOUNDSTRETCHPARML1REDUCE}[$i],$goutoniom->{FOUNDSTRETCHPARML2REDUCE}[$i];
  
  }
  
 $bondbendnumfinal = $goutoniom->{FOUNDBENDNUMREDUCE};
 
 for ($i=0;$i<$bondbendnumfinal;$i++) {
  
  printf OUTPUTFILE "HrmBnd1 %4s %4s %4s %10.4f %10.4f\n", $goutoniom->{FOUNDBENDL1REDUCE}[$i],$goutoniom->{FOUNDBENDL2REDUCE}[$i],
       $goutoniom->{FOUNDBENDL3REDUCE}[$i],$goutoniom->{FOUNDBENDPARML1REDUCE}[$i],$goutoniom->{FOUNDBENDPARML2REDUCE}[$i];
  
  }
 
 close OUTPUTFILE;

} # sub outputparm{


# Write out total charge

sub writecharges{

 my $goutoniom = shift;
 my $debug = $goutoniom->{DEBUG};
 my $raw = $goutoniom->{RAW};
 
 if ($raw == 0) {
  printf "Total charge of real system is %10.6f.\n",$goutoniom->{RCHARGE};
  printf "Total charge of high layer is %10.6f.\n",$goutoniom->{HCHARGE};
  printf "Total charge of medium layer is %10.6f.\n",$goutoniom->{MCHARGE};
  printf "Total charge of low layer is %10.6f.\n",$goutoniom->{LCHARGE};
  printf "Total charge of high plus medium layer is %10.6f.\n",$goutoniom->{HCHARGE}+$goutoniom->{MCHARGE};
 } else {
  print "Total charge of real system is $goutoniom->{RCHARGE}\n",
  print "Total charge of high layer is $goutoniom->{HCHARGE}\n";
  print "Total charge of medium layer is $goutoniom->{MCHARGE}\n";
  print "Total charge of low layer is $goutoniom->{LCHARGE}\n";
  print "Total charge of high plus medium layer is ", $goutoniom->{HCHARGE}+$goutoniom->{MCHARGE},"\n";
 }
 

} # sub writecharges{


# Calculate dipolemoment

sub dipolemoment{

 my $goutoniom = shift;
 my $debug = $goutoniom->{DEBUG};
 my $raw = $goutoniom->{RAW};
 
 if ($raw == 0) {
  printf "\nDipole moment (Debye) (X, Y, Z) is \n(%15.4f, %15.4f, %15.4f).\n",$goutoniom->{XDIPOLECAL},
                                         $goutoniom->{YDIPOLECAL},$goutoniom->{ZDIPOLECAL};
  printf "Total Dipole moment (Debye) is %15.4f.\n",$goutoniom->{TOTDIPOLECAL};
 } else {
  print "\nDipole moment (Debye) (X, Y, Z) is \n(",$goutoniom->{XDIPOLECAL},", ", $goutoniom->{YDIPOLECAL},", ", $goutoniom->{ZDIPOLECAL},").\n";
  print "Total Dipole moment (Debye) is $goutoniom->{TOTDIPOLECAL}\n";
 }

} # sub dipolemoment{


# Subroutine to compare two bonds are the same or not

sub bondmatch {

  my $bond1atm1;
  my $bond1atm2;
  my $bond2atm1;
  my $bond2atm2;
  my $ucbond1atm1;
  my $ucbond1atm2;
  my $ucbond2atm1;
  my $ucbond2atm2;
  
  my $ismatch = 0;
  
#  shift;
  $bond1atm1 = shift;
  $bond1atm2 = shift;
  $bond2atm1 = shift;
  $bond2atm2 = shift;

  $ucbond1atm1 = uc($bond1atm1);
  $ucbond1atm2 = uc($bond1atm2);
  $ucbond2atm1 = uc($bond2atm1);
  $ucbond2atm2 = uc($bond2atm2);
  
  
#  print "Compare bond $ucbond1atm1-$ucbond1atm2 and $ucbond2atm1-$ucbond2atm2\n";
  if ( ($ucbond1atm1 eq  $ucbond2atm1) and ($ucbond1atm2 eq  $ucbond2atm2) ) {$ismatch = 1;}
  if ( ($ucbond1atm1 eq  $ucbond2atm2) and ($ucbond1atm2 eq  $ucbond2atm1) ) {$ismatch = 1;}
  
  
  return $ismatch;

} #sub bondmatch 

# Subroutine to compare two bond angles are the same or not

sub anglematch {

  my $angle1atm1;
  my $angle1atm2;
  my $angle1atm3;
  
  my $angle2atm1;
  my $angle2atm2;
  my $angle2atm3;
  
  my $ucangle1atm1;
  my $ucangle1atm2;
  my $ucangle1atm3;
  
  my $ucangle2atm1;
  my $ucangle2atm2;
  my $ucangle2atm3;
  
  my $ismatch = 0;
  
#  shift;
  $angle1atm1 = shift;
  $angle1atm2 = shift;
  $angle1atm3 = shift;
  
  $angle2atm1 = shift;
  $angle2atm2 = shift;
  $angle2atm3 = shift;
  

  $ucangle1atm1 = uc($angle1atm1);
  $ucangle1atm2 = uc($angle1atm2);
  $ucangle1atm3 = uc($angle1atm3);
  
  $ucangle2atm1 = uc($angle2atm1);
  $ucangle2atm2 = uc($angle2atm2);
  $ucangle2atm3 = uc($angle2atm3);
  
  if ( $ucangle1atm2 eq $ucangle2atm2) {
  
    if ( ($ucangle1atm1 eq  $ucangle2atm1) and ($ucangle1atm3 eq  $ucangle2atm3) ) {$ismatch = 1;}
    if ( ($ucangle1atm1 eq  $ucangle2atm3) and ($ucangle1atm3 eq  $ucangle2atm1) ) {$ismatch = 1;}
   } # if $ucangle1atm2 eq $ucangle2atm2
   

#  print "Compare bond angle $ucangle1atm1-$ucangle1atm2-$ucangle1atm3 and $ucangle2atm1-$ucangle2atm2-$ucangle2atm3\n";
#  print "Compare bond angle $angle1atm1-$angle1atm2-$angle1atm3 and $angle2atm1-$angle2atm2-$angle2atm3\n";
  
  return $ismatch;

} #sub anglematch 


1;
__END__

=head1 VERSION

0.01

=head1 AUTHOR

Peng Tao, E<lt>tao.21@osu.edu<gt>

=head1 COPYRIGHT

Copyright (c) 2009~2010 by Peng Tao


=cut

