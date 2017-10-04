#!/usr/bin/perl
# vim:set foldmethod=marker:
use strict;
use 5.10.0;

my $list_n = shift;
open my $list, $list_n or die "$list_n can't open.\n";

my $MAX_CYCLE = shift;
die "Specify MAX_CYCLE...\n" if $MAX_CYCLE <= 0;

my %init_strcture_of;
my $cyc = 1;
while (<$list>){
    my @cols = split;
    $init_strcture_of{$cyc} = \@cols;
    $cyc++;
}
close $list;

my @trace;

# last init structures.
for (my $j; $j<8; $j++) {
    push @{$trace[$j]}, $init_strcture_of{$MAX_CYCLE+1}->[$j];
}
# backtrace
for (my $i=$MAX_CYCLE; $i>=0; $i--) {
    my $inits = $init_strcture_of{$i};
    for (my $j; $j<8; $j++) {
        my $next = $trace[$j][0];
        my ($cy, $th, $st) = split "-",$next;
        unshift @{$trace[$j]}, $inits->[$th];
    }
}

foreach  (@trace){
    my @tmp = @$_;
    $" = "\n";
    print "@tmp";
    print "\n//\n";
}

