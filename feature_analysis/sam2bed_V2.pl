use strict;
use warnings;
use List::Util qw /min max/;

#This version use the alignment record of the downstream read to calculate the DNA fragment length, which is more accurate.

if(@ARGV!=1)
{
        warn "usage: <*.filtered.sam>\n";
        exit 1;
}

while(<>)
{
	chomp;
	#A00869:328:HL5MHDSXY:4:1553:27968:31955	99	chr22	16050237	30	70M	=	16050304	142
	my @F=split;
        next if ($F[6] ne "=");
	if($F[8]<0 && $F[3]>=$F[7])
	{
		my $s=min($F[3], $F[7])-1;
		my $len=abs($F[3]-$F[7])+length($F[9]);
		my $e=$s+$len;
		print join("\t", $F[2], $s, $e, $len), "\n";
	}
}
