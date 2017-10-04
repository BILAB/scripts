#!/usr/bin/perl
# vim:set foldmethod=marker:
use strict;
use 5.10.0;

my $file_n = shift;
open my $list, $file_n or die "File $file_n can't open.\n";

my $total_time = 0;
my $traj_cnt = 0;
while (<$list>){
    chomp;
    my ($c, $th, $last) = split "-";
    system "gmx trjconv  -f ../cyc$c/${c}-${th}.trr ".
                        "-s ../cyc$c/${c}-${th}.tpr ".
                        "-n ../../../grotop/index.ndx ".
                        "-o $c.cut.trr ".
                        "-t0 $total_time -b 0 -e $last << zzz\n".
                            "0\n".
                        "zzz";
    $total_time += $last;
    $traj_cnt++;
}

close $list;

my @trrs = glob "*.cut.trr";
if (@trrs != $traj_cnt){
    die "Some trajectory would'nt be faild to conv.\n";
}
else {
    system "gmx trjcat -f *.cut.trr -o full_traj.trr";
} 
