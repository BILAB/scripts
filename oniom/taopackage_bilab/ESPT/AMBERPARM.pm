package ESPT::AMBERPARM;

our $VERSION = '0.01';

### Version History ###
# 0.01	Read in AMBER Parameter file PARM99.dat and gaff.dat
#

=head1 NAME

ESPT::AMBERPARM - AMBER Force Field for bond strectch and bend lookup

=head1 SYNOPSIS

    use ESPT::AMBERPARM;


=head1 DESCRIPTION

This module read in AMBER parameter file parm99.dat and gaff.dat 
for bond stretch and bond bend parameter lookup.

=head2 Bond stretch and bend bend parameter lookup

When passing atom types, this object will be able to return the
paramters for bond stretch or bond bend.

=cut

=head1 ATTRIBUTES

=over 15

=item B<$file-E<gt>new()>

Creates a new AMBERPARM object

=cut

use strict;
use warnings;

our @parm99stretchatmtyp1;
our @parm99stretchatmtyp2;
our @parm99stretchparm1;
our @parm99stretchparm2;

our @parm99bendatmtyp1;
our @parm99bendatmtyp2;
our @parm99bendatmtyp3;
our @parm99bendparm1;
our @parm99bendparm2;

our @gaffstretchatmtyp1;
our @gaffstretchatmtyp2;
our @gaffstretchparm1;
our @gaffstretchparm2;

our @gaffbendatmtyp1;
our @gaffbendatmtyp2;
our @gaffbendatmtyp3;
our @gaffbendparm1;
our @gaffbendparm2;

our $amberparm99file = "/Users/YoshitakaM/apps/taopackage/ESPT/amberff/parm99.dat";
our $ambergafffile   = "/Users/YoshitakaM/apps/taopackage/ESPT/amberff/gaff.dat";

our $parm99strtchnum;
our $parm99bendnum;
our $gaffstrtchnum;
our $gaffbendnum;


my $curresid;
my $endresidue;
my $curatmnam;
my $curatmtyp;
my $curelenum;
my $curparchg;
my $prepfile;
my $prepfilewpath;
my $line1;
my $line2;
my $templine;

#print "Open file $amberparm99file to read in AMBER PARM99 force field.\n";

open(AMBERPARMFILE,$amberparm99file) || die "Could not read AMBEF force file $amberparm99file\n$!\n";

# Reading AMBER PARM 99 file $amberparm99file for bond stretch and bend parameters.

$parm99strtchnum = 0;

$parm99bendnum = 0;

while ($templine = <AMBERPARMFILE>) {

  if ($templine =~ /^\s*([A-Z]\S?)\s*-([A-Z]\S?)\s+(\d+.\d+)\s+(\d+.\d+)/) {
     $parm99stretchatmtyp1[$parm99strtchnum] = $1;
     $parm99stretchatmtyp2[$parm99strtchnum] = $2;
     $parm99stretchparm1[$parm99strtchnum] = $3;
     $parm99stretchparm2[$parm99strtchnum] = $4;

#     print "Read in Parm99 bond stretch parameter $parm99strtchnum, ($1)-($2): ($3) ($4)\n";

     $parm99strtchnum++;
    } # if $templine 
    
  if ($templine =~ /^\s*([A-Z]\S?)\s*-([A-Z]\S?)\s*-\s*([A-Z]\S?)\s+(\d*.\d*)\s+(\d*.\d*)/) {
     $parm99bendatmtyp1[$parm99bendnum] = $1;
     $parm99bendatmtyp2[$parm99bendnum] = $2;
     $parm99bendatmtyp3[$parm99bendnum] = $3;
     $parm99bendparm1[$parm99bendnum] = $4;
     $parm99bendparm2[$parm99bendnum] = $5;

#     print "Read in Parm99 bond angle bend parameter $parm99bendnum, ($1)-($2)-($3): ($4) ($5)\n";

     $parm99bendnum++;
    } # if $templine 

} # while (<AMBERPARMFILE>) 

close (AMBERPARMFILE);

#print "Read $parm99strtchnum set of bond stretch parameters and $parm99bendnum set of bond bend parameters from \n$amberparm99file.\n";

#print "Open file $ambergafffile to read in AMBER GAFF force field.\n";

open(AMBERGAFFFILE,$ambergafffile) || die "Could not read AMBEF force file $ambergafffile\n$!\n";

# Reading AMBER PARM 99 file $ambergafffile for bond stretch and bend parameters.

$gaffstrtchnum = 0;

$gaffbendnum = 0;

while ($templine = <AMBERGAFFFILE>) {

  if ($templine =~ /^\s*([a-z]\S?)\s*-([a-z]\S?)\s+(\d+.\d+)\s+(\d+.\d+)/) {
     $gaffstretchatmtyp1[$gaffstrtchnum] = $1;
     $gaffstretchatmtyp2[$gaffstrtchnum] = $2;
     $gaffstretchparm1[$gaffstrtchnum] = $3;
     $gaffstretchparm2[$gaffstrtchnum] = $4;

#     print "Read in gaff bond stretch parameter $gaffstrtchnum, ($1)-($2): ($3) ($4)\n";

     $gaffstrtchnum++;
    } # if $templine 
    
  if ($templine =~ /^\s*([a-z]\S?)\s*-([a-z]\S?)\s*-\s*([a-z]\S?)\s+(\d*.\d*)\s+(\d*.\d*)/) {
     $gaffbendatmtyp1[$gaffbendnum] = $1;
     $gaffbendatmtyp2[$gaffbendnum] = $2;
     $gaffbendatmtyp3[$gaffbendnum] = $3;
     $gaffbendparm1[$gaffbendnum] = $4;
     $gaffbendparm2[$gaffbendnum] = $5;

#     print "Read in gaff bond angle bend parameter $gaffbendnum, ($1)-($2)-($3): ($4) ($5)\n";

     $gaffbendnum++;
    } # if $templine 

} # while (<AMBERGAFFFILE>) 

close (AMBERGAFFFILE);

# print "Read $gaffstrtchnum set of bond stretch parameters and $gaffbendnum set of bond bend parameters from \n$ambergafffile.\n";

sub new {
 my $AMBERPARM = {};
 
 bless ($AMBERPARM);
 return $AMBERPARM; 
} # new



# Subroutine to look up bond stretch parameter

sub bondstretchlookup {
  my $atomtype1;
  my $atomtype2;
  my $force = 0.0;
  my $bondlen = 0.0;
  my $foundparm = 0;
  my $debug;
  my $i;
  
  shift;
  $atomtype1 = shift;
  $atomtype2 = shift;
  $debug     = shift; 
#  print "  Looking up bond stretch parameter for bond between $atomtype1 and $atomtype2.\n" if $debug >= 1;
  
  for ($i = 0; $i < $parm99strtchnum; $i++) {
    
    if ( bondmatch($atomtype1,$atomtype2,$parm99stretchatmtyp1[$i],$parm99stretchatmtyp2[$i]) ) {
       $force = $parm99stretchparm1[$i];
       $bondlen = $parm99stretchparm2[$i];
       $foundparm = 1;
       if ($debug>=1) {print "Found missing parameter for $atomtype1 and $atomtype2 in Parm99.dat as $force  $bondlen\n";}
       last;
      } # if bondmatch
   } # for parm99strtchnum
  
  if ( $foundparm == 0) {
   for ($i = 0; $i < $gaffstrtchnum; $i++) {
    
     if ( bondmatch($atomtype1,$atomtype2,$gaffstretchatmtyp1[$i],$gaffstretchatmtyp2[$i]) ) {
        $force = $gaffstretchparm1[$i];
        $bondlen = $gaffstretchparm2[$i];
        $foundparm = 1;
        if ($debug>=1) {print "Found missing parameter for $atomtype1 and $atomtype2 in gaff.dat as $force  $bondlen\n";}
        last;
       } # if bondmatch
    } # for parm99strtchnum
   
   } #  if ( $foundparm == 0)
  
  
  if ($foundparm == 0)  {
    print "Could not find bond stretch paramters for bond $atomtype1-$atomtype2.\n";
    print "Dummy value (0.0) for bond stretch froce and bond length are used.\n";
   } # if $foundparm
  
   
  return $force, $bondlen;

} # bondstretchlookup

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


# Subroutine to look up bond angle bend parameter

sub bondbendlookup {
  my $atomtype1;
  my $atomtype2;
  my $atomtype3;
  
  my $force = 0.0;
  my $bondlen = 0.0;
  my $foundparm = 0;
  my $debug;
  my $i;
  
  shift;
  $atomtype1 = shift;
  $atomtype2 = shift;
  $atomtype3 = shift;
  $debug     = shift; 
#  print "  Looking up bond stretch parameter for bond between $atomtype1 and $atomtype2.\n" if $debug >= 1;
  
  for ($i = 0; $i < $parm99bendnum; $i++) {
    
    if ( anglematch($atomtype1,$atomtype2,$atomtype3,$parm99bendatmtyp1[$i],$parm99bendatmtyp2[$i],$parm99bendatmtyp3[$i]) ) {
       $force   = $parm99bendparm1[$i];
       $bondlen = $parm99bendparm2[$i];
       $foundparm = 1;
       if ($debug>=1) {print "Found missing parameter for $atomtype1-$atomtype2-$atomtype3 in Parm99.dat as $force  $bondlen\n";}
       last;
      } # if anglematch
   } # for parm99bendnum
  
  if ( $foundparm == 0) {
   for ($i = 0; $i < $gaffbendnum; $i++) {
    
     if ( anglematch($atomtype1,$atomtype2,$atomtype3,$gaffbendatmtyp1[$i],$gaffbendatmtyp2[$i],$gaffbendatmtyp3[$i]) ) {
        $force   = $gaffbendparm1[$i];
        $bondlen = $gaffbendparm2[$i];
        $foundparm = 1;
        if ($debug>=1) {print "Found missing parameter for $atomtype1-$atomtype2-$atomtype3 in gaff.dat as $force  $bondlen\n";}
#        if ($debug>=1) {print "Found missing parameter for $atomtype1 and $atomtype2 in gaff.dat as $force  $bondlen\n";}
        last;
       } # if anglematch
    } # for gaffbendhnum
   
   } #  if ( $foundparm == 0)
  
  
  if ( ($foundparm == 0) && ($debug >= 1) )  {
    print "Could not find bond stretch paramters for bond $atomtype1-$atomtype2.\n";
    print "Dummy value (0.0) for bond stretch froce and bond length are used.\n";
   } # if $foundparm
  
   
  return $force, $bondlen;

} # bondbendlookup


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

=back

=head1 VERSION

0.01

=head1 SEE ALSO

F<ESPT::GoutOniom>

=head1 AUTHOR

Peng Tao, E<lt>tao.21@osu.eduE<gt>

=head1 COPYRIGHT

Copyright (c) 2009~2010 by Peng Tao


=cut
