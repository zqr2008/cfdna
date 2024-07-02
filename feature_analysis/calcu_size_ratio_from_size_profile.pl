#!/usr/bin/perl -w
use strict;

#This script is to calculate size ratio between short(60~155bp) and long(170~255bp) DNA fragments.

while (@ARGV<2){
   warn "perl $0 <sample> <*size.freq>\n";
   exit 1;
   }

my $sample=shift;

my $short_cnt=0;
my $long_cnt=0;

my $short_ratio=0;
my $long_ratio=0;
while (<>)
{
        chomp;
   #104     84248   0.0438926874432612      55      1.62914691943128
   my @F=split/\s+/;
   if ($F[0]>=60 && $F[0]<=155)
   {
   	   $short_cnt+=$F[1];
   	   $short_ratio+=$F[2];
   }
   elsif($F[0]>=170 && $F[0]<=250)
   {
   	   $long_cnt+=$F[1];
   	   $long_ratio+=$F[2];
   }
}

my $size_ratio=0;
if ($long_cnt>0)
{
	$size_ratio=$short_cnt/$long_cnt;
}

print join("\t", $sample, $short_ratio, $long_ratio, $size_ratio)."\n";


