package ESPT::GTable;

use ESPT::ONIOMLog 0.01;
use strict;
use warnings;

our $VERSION = '0.01';

=head1 NAME

ESPT::GTable - Create table based on gaussian energies.

=head1 SYNOPSIS

   use ESPT::GTable;

   my $gtable = GTable->new();

=head1 DESCRIPTION

This module creates table based on gaussian energies.

=cut

### Version History ###
# 0.01	digest ONIOM gtable files

### To Do ###
# Treat general gaussian log files

=head1 ATTRIBUTES

=over 15

=cut

## the object constructor **

sub new {
	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my $gtable = ();

	$gtable->{TYPE} = "gtable";

	# program info
	$gtable->{PROGRAM} = "Gaussian";
        
        $gtable->{DEBUG} = 0;
        $gtable->{RAW}   = 0;
	
	# Link 0 & Route commands

	
	# Arrays for ONIOM Energies
	$gtable->{ONIOMtotal} = [];     #Final ONIOM Energy
	$gtable->{ONIOMRealLow} = [];   # Real System Low Level ONIOM Energy
	$gtable->{ONIOMModelLow} = [];  # Model System Low Level ONIOM Energy
	$gtable->{ONIOMModelHigh} = []; # Real System Low Level ONIOM Energy
        $gtable->{ONIOMModelMed} = [];  # Model System Medium Level ONIOM Energy
        $gtable->{ONIOMMidLow} = [];    # Middle System Low Level ONIOM Energy
        $gtable->{ONIOMMidMed} = [];    # Middle System Medium Level ONIOM Energy

        $gtable->{ONIOMReal_HihgLow} = [];   # Real System Low Level ONIOM Energy - Model System Low Level ONIOM Energy 
        
        $gtable->{LOGFILELIST} = [];    # All the log file names (within complete path) in current block
        $gtable->{LOGFILENAMELIST} = [];    # All the log file names (no path) in current block

        $gtable->{InputFilename} = undef;    # ONIOM input filename
        $gtable->{OutputFilename} = undef;    # ONIOM output filename
        
        $gtable->{ISTHREELAYER} = undef; # Flag for three-layer ONIOM calculation
        
        
	bless($gtable, $class);
	return $gtable;
}


## methods ##

# set filename & spin then digest the file
sub analyze : method {
	my $gtable = shift;
	$gtable->digest();
	return;
}


## subroutines ##

sub digest {

my $gtable = shift;

my $curoniom = ESPT::ONIOMLog->new();

# flags & counters
my $debug = $gtable->{DEBUG};
my $raw   = $gtable->{RAW};

my ($blocknum,$blocktag);

my ($currline, @currgfile, $pathdepth, $currfilename,$filewpath, $isthreelayer, $logindex,$i);

my ($conversionfactor);

$blocknum = 0;
$conversionfactor = 627.5095;

# open filename for reading 

print "Open input file ", $gtable->{InputFilename}, " for process...\n\n"  if $debug >= 0;

open(INPUTFILE,$gtable->{InputFilename}) || die "Could not read $gtable->{InputFilename}\n$!\n";

open(OUTPUTFILE,">$gtable->{OutputFilename}") || die "Could not open $gtable->{OutputFilename} for output.\n$!\n";


# grab everything which may be useful
while (<INPUTFILE>){
        
	# skip blank lines
	next if /^\s*$/;
        
        # skip reference        
        next if /^\s*#/;

#        print $_;

	# Find tag of one block
	if ( /^\s*\<\s*LABEL\s*\>\s*(.*)\s*\<\/LABEL\>/i ) {
                $blocktag = $1;
		next;
	}
        
	# Find beginning of one ONIOM block
	if ( /^\s*\<\s*ONIOM\s*\>/i ) {
		$blocknum++;
                if ( not (defined $blocktag) ) {$blocktag = "ONIOM jobs output"};
                print "The label of block $blocknum is $blocktag.\n" if $debug >= 1;
                print "Since ONIOMSTART is used for this block, all log files will be treated as ONIOM jobs.\n" if $debug >= 2;
                
                printf OUTPUTFILE "=================================================\n";
                printf OUTPUTFILE "Output for block $blocknum: $blocktag\n";
                
                $logindex = 0;
                
                $gtable->{ONIOMtotal} = ();    
                $gtable->{ONIOMRealLow} = ();  
                $gtable->{ONIOMModelLow} = (); 
                $gtable->{ONIOMModelHigh} = ();
                $gtable->{ONIOMModelMed} = (); 
                $gtable->{ONIOMMidLow} = ();   
                $gtable->{ONIOMMidMed} = ();   
                $gtable->{LOGFILELIST} = ();     
                $gtable->{LOGFILENAMELIST} = ();
                $gtable->{ONIOMReal_HihgLow} = ();
                
                while (<INPUTFILE>){        

                  $currline = $_;
     	          next if /^$/;
                  next if /^\s*#/;
                  
                  # if encounter end of one block. exit this loop
                  
                  if (/^\s*\<\s*\/ONIOM\s*\>/i) {

                     printf OUTPUTFILE "-------------------------------------------------\n";
                     printf OUTPUTFILE  "Relative values in kcal/mol:\n\n";

                     if ( $gtable->{ISTHREELAYER} ==0 ) {
                        printf OUTPUTFILE  "%23s\t %23s\t %23s\t %23s\t %23s\n\n", "ONIOM Total","Model High",
                                          "Model Low","Real Low","Real-Model Low";
                       } else {
                        printf OUTPUTFILE  "%23s\t %23s\t %23s\t %23s\n\n", "ONIOM Total","Model High","Middle Medium","Real Low";
                      
                       }                
                    
                    for ($i=0;$i < $logindex; $i++) {
                    
                       # Output relative energy with respect to first log file in kcal/mol.
                       # Treat 2-layer and 3-layer ONIOM job differently.
                       
                       # Output file name for every line
                       printf  OUTPUTFILE $gtable->{LOGFILENAMELIST}[$i]."\n";
                       
#                       print  "Processing file ", $gtable->{LOGFILENAMELIST}[$i]." ...\n\n" if $debug >= 0;

                       if ( $gtable->{ISTHREELAYER} ==0 ) {
                       # 2-Layer ONIOM job
                       
                       if ($raw == 1) {
                                         printf  OUTPUTFILE "%23.12f\t %23.12f\t %23.12f\t %23.12f\t %23.12f\n\n", 
                                          ($gtable->{ONIOMtotal}[$i] - $gtable->{ONIOMtotal}[0]) * $conversionfactor,
                                          ($gtable->{ONIOMModelHigh}[$i] - $gtable->{ONIOMModelHigh}[0]) * $conversionfactor,
                                          ($gtable->{ONIOMModelLow}[$i] - $gtable->{ONIOMModelLow}[0]) * $conversionfactor,
                                          ($gtable->{ONIOMRealLow}[$i] - $gtable->{ONIOMRealLow}[0]) * $conversionfactor,
                                          ($gtable->{ONIOMReal_HihgLow}[$i] - $gtable->{ONIOMReal_HihgLow}[0]) * $conversionfactor;
                           } else {
                                         printf  OUTPUTFILE "%23.6f\t %23.6f\t %23.6f\t %23.6f\t %23.6f\n\n", 
                                          ($gtable->{ONIOMtotal}[$i] - $gtable->{ONIOMtotal}[0]) * $conversionfactor,
                                          ($gtable->{ONIOMModelHigh}[$i] - $gtable->{ONIOMModelHigh}[0]) * $conversionfactor,
                                          ($gtable->{ONIOMModelLow}[$i] - $gtable->{ONIOMModelLow}[0]) * $conversionfactor,
                                          ($gtable->{ONIOMRealLow}[$i] - $gtable->{ONIOMRealLow}[0]) * $conversionfactor,
                                          ($gtable->{ONIOMReal_HihgLow}[$i] - $gtable->{ONIOMReal_HihgLow}[0]) * $conversionfactor;
                           }
                              
                       } else {
                       # 3-Layer ONIOM job

                       if ($raw == 1) {
                                        printf  OUTPUTFILE "%23.12f\t %23.12f\t %23.12f\t %23.12f\t %23.12f\n\n", 
                                          ($gtable->{ONIOMtotal}[$i] - $gtable->{ONIOMtotal}[0]) * $conversionfactor,
                                          ($gtable->{ONIOMModelHigh}[$i] - $gtable->{ONIOMModelHigh}[0]) * $conversionfactor,
                                          ($gtable->{ONIOMMidMed}[$i] - $gtable->{ONIOMMidMed}[0]) * $conversionfactor,
                                          ($gtable->{ONIOMRealLow}[$i] - $gtable->{ONIOMRealLow}[0]) * $conversionfactor;
                            } else {
                                        printf  OUTPUTFILE "%23.6f\t %23.6f\t %23.6f\t %23.6f\t %23.6f\n\n", 
                                          ($gtable->{ONIOMtotal}[$i] - $gtable->{ONIOMtotal}[0]) * $conversionfactor,
                                          ($gtable->{ONIOMModelHigh}[$i] - $gtable->{ONIOMModelHigh}[0]) * $conversionfactor,
                                          ($gtable->{ONIOMMidMed}[$i] - $gtable->{ONIOMMidMed}[0]) * $conversionfactor,
                                          ($gtable->{ONIOMRealLow}[$i] - $gtable->{ONIOMRealLow}[0]) * $conversionfactor;
                            }
                       
                       } 
                        
                     
                     } # for ($i=0;$i < $logindex; $i++) {

                    printf OUTPUTFILE "=================================================\n\n\n";

                    last;
                   } # if (/^\.{3}END\.{3}/) {

                  # if encounter end of file. without end of block, print warning sign
                  
                  if ( eof(INPUTFILE) == 1) {  print "End of block $blocknum is missing.\n" if $debug >= 0;  }
                  
                  if ( /^\s*\<\s*ONIOM\s*\>/i ) {
                      print "Encouter a start of another ONIOM block before reaching end of the current block $blocknum.\n" if $debug >= 0;
                      die "Exit \n";
                     }
   
                  if ( /^\s*\<\s*Gaussian\s*\>/i ) {
                      print "Encouter a start of another normal Gaussian block before reaching end of the current block $blocknum.\n" if $debug >= 0;
                      die "Exit \n";
                     }
                  
                  # Read in each log file with full path
                  @currgfile = split(/\//,$currline);
                  $pathdepth = @currgfile;
                  $currfilename = $currgfile[$pathdepth-1];
                  $filewpath = $currline;
                  chop ($filewpath);
                  chop ($currfilename);
#                  print "@currgfile";
                  
                  # Store these file names w/o path
                  $gtable->{LOGFILELIST}[$logindex] = $filewpath;
                  $gtable->{LOGFILENAMELIST}[$logindex] = $currfilename;
                  
                  print "Process ONIOM output file $currfilename ...\n\n" if $debug >= 0;

                  $curoniom->analyze($filewpath, $debug);
                  
                  # If this log file job exit abnormally. print warning sign.
                  if ($curoniom->{ISNORMAL} == 0){
                    print "File $currfilename exited abnormally. Its results may not be used for energy comparison.\n" if $debug >= 0;
                   }
                  
                  
                  # Check if it is three layer ONIOM job
                  if ( $logindex == 0 ) {
                     $gtable->{ISTHREELAYER} = $curoniom->{ISTHREELAYER};

                    if ($gtable->{ISTHREELAYER} == 0) {
                       print "The first log file is not three layer ONIOM calculation.\n" if $debug >= 3;
                      } else {
                       print "The first log file is three layer ONIOM calculation.\n" if $debug >= 3;
                       
                       print OUTPUTFILE "First file of block $blocktag is three layer ONIOM calculation.\n";
                      }
                     printf OUTPUTFILE  "Absolute values in Hartee:\n";
 
                     if ( $gtable->{ISTHREELAYER} ==0 ) {
                        printf OUTPUTFILE  "%23s\t %23s\t %23s\t %23s\t %23s\n\n", "ONIOM Total","Model High",
                                          "Model Low","Real Low","Real-Model Low";
                       } else {
                        printf OUTPUTFILE  "%23s\t %23s\t %23s\t %23s\n\n", "ONIOM Total","Model High","Middle Medium","Real Low";
                      
                       }                

                    }  # if ( $logindex == 0 ) {

                  
                  # If any log file is different from the first log file (3-layer vs. 2-layer), exit with warning.
                  if ( $logindex > 0 ) {
                    if ( $gtable->{ISTHREELAYER} != $curoniom->{ISTHREELAYER}) {
                         if ( $gtable->{ISTHREELAYER} == 0 ){
                          print "The first file \n",$gtable->{LOGFILENAMELIST}[0],
                                " \n is not three layer ONIOM job. But file ", $logindex+1,
                                "\n",$gtable->{LOGFILENAMELIST}[$logindex],"\n is a three layer ONIOM job.\n" if $debug >= 0;
                          die "Exit here.\n";
                         } else {
                          print "The first file \n",$gtable->{LOGFILENAMELIST}[0],
                                " \n is a three layer ONIOM job. But file ", $logindex+1,
                                "\n",$gtable->{LOGFILENAMELIST}[$logindex],"\n is not a three layer ONIOM job.\n" if $debug >= 0;
                          die "Exit here.\n";
                         }
                      }
                    }
                   
                   # Treat 2-layer and 3-layer ONIOM job differently.
                   
                   if ( $gtable->{ISTHREELAYER} ==0 ) {
                    # 2-Layer ONIOM job
                    $gtable->{ONIOMtotal}[$logindex]        = $curoniom->{ONIOMENERGY};
                    $gtable->{ONIOMRealLow}[$logindex]      = $curoniom->{LOWREAL};
                    $gtable->{ONIOMModelLow}[$logindex]     = $curoniom->{LOWMODEL};
                    $gtable->{ONIOMModelHigh}[$logindex]    = $curoniom->{HIGHMODEL};
                    $gtable->{ONIOMReal_HihgLow}[$logindex] = $curoniom->{LOWREAL_MODEL};
                    
                    printf OUTPUTFILE $currfilename."\n";
                    if ($raw == 1) {
                              printf  OUTPUTFILE "%23.12f\t %23.12f\t %23.12f\t %23.12f\t %23.12f\n\n", $curoniom->{ONIOMENERGY},$curoniom->{HIGHMODEL},
                              $curoniom->{LOWMODEL}, $curoniom->{LOWREAL},$curoniom->{LOWREAL_MODEL};
                          } else {
                              printf  OUTPUTFILE "%23.6f\t %23.6f\t %23.6f\t %23.6f\t %23.6f\n\n", $curoniom->{ONIOMENERGY},$curoniom->{HIGHMODEL},
                              $curoniom->{LOWMODEL}, $curoniom->{LOWREAL},$curoniom->{LOWREAL_MODEL};
                          }
                    } else {
                    # 3-Layer ONIOM job
                    $gtable->{ONIOMtotal}[$logindex]        = $curoniom->{ONIOMENERGY};
                    $gtable->{ONIOMRealLow}[$logindex]      = $curoniom->{LOWREAL};
                    $gtable->{ONIOMModelLow}[$logindex]     = $curoniom->{LOWMODEL};
                    $gtable->{ONIOMModelHigh}[$logindex]    = $curoniom->{HIGHMODEL};
                    $gtable->{ONIOMModelMed}[$logindex]     = $curoniom->{MEDMODEL};
                    $gtable->{ONIOMMidLow}[$logindex]       = $curoniom->{LOWMID};
                    $gtable->{ONIOMMidMed}[$logindex]       = $curoniom->{MEDMID};

                    printf OUTPUTFILE $currfilename."\n";
                    if ($raw == 1) 
                       {
                          printf  OUTPUTFILE "%23.12f\t %23.12f\t %23.12f\t %23.12f\n\n", $curoniom->{ONIOMENERGY},$curoniom->{HIGHMODEL},
                             $curoniom->{MEDMID}, $curoniom->{LOWREAL};
                        } else {
                          printf  OUTPUTFILE "%23.6f\t %23.6f\t %23.6f\t %23.6f\n\n", $curoniom->{ONIOMENERGY},$curoniom->{HIGHMODEL},
                             $curoniom->{MEDMID}, $curoniom->{LOWREAL};
                        }
                    } 

                   $logindex++;
                 } #   while (<INPUTFILE>)  #  for ONIOM block
                 
		next;
	}
        
        
	# Find beginning of one normal Gaussian job block
	if ( /^\s*\<\s*Gaussian\s*\>/i ) {
		$blocknum++;
                print "The label of block $blocknum is $blocktag.\n" if $debug >= 0;
                print "Since Gaussian tag is used for this block, all log files are normal Gaussian job log files.\n" if $debug >= 0;
                print "Current version of gaussiantable cannot treat normal Gaussian log files. This block will be skipped.\n\n" if $debug >= 0;
                
		next;
	} # if ( /^\s*\<\s*Gaussian\s*\>/i ) {
        
        
} # overall while (<INPUTFILE>){

print "Finished output file ", $gtable->{OutputFilename}, ".\n\n"  if $debug >= 0;

} # digest

1;

__END__

=head1 VERSION

0.07

=head1 AUTHOR

Peng Tao, E<lt>tao.21@osu.edu<gt>

=head1 COPYRIGHT

Copyright (c) 2009~2010 by Peng Tao


=cut

