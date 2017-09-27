#!/usr/bin/perl -w

use warnings;
use strict;

my @chain = qw(A B C D E F G H I J K);
my $count = 0;
my $prev_resid = "";
my $resid = "";
my $prev_wato_resid="";
my $wato_resid = "";

open(IN,$ARGV[0]);
while(<IN>){
  next if ($_ =~ /^TER/ or /^END/);
  next if ($_ !~ /^ATOM/ or /^HETATM/);
  if ($_ =~ /O   WAT/){
    if ($wato_resid eq ""){
      $wato_resid = substr($_,22,4);
    }else{
      $wato_resid++;
      if ($wato_resid == 10000){
	$wato_resid = 0;
      }
    }
  }
  if ($_ =~ /WAT/){
    substr($_,22,4) = sprintf("%4s",$wato_resid);
  }
  $resid = substr($_,22,4);
  $resid = trim($resid);
  $prev_resid = $resid if ($prev_resid eq "");
  if (($resid eq "0" ) and ($resid ne $prev_resid)){
    $count++;
  }
  substr ($_,21,1) = $chain[$count];
  print;
  $prev_resid = $resid;
}
close(IN);

sub trim {
  my $val = shift;
  $val =~ s/^ *(.*?) *$/$1/;
  return $val;
}
