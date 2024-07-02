#!/usr/bin/perl -w
use strict;

while (@ARGV<3){
   warn "perl $0 <sample_type> <sample> <*motif.freq>\n";
   exit 1;
   }

my $sample_type=shift;
my $sample=shift;
my $divScore=0;
my $short_divScore=0;
my $peak_divScore=0;
my $long_divScore=0;
my $SL_Ratio=0;
while(<>){
   chomp;
   next if (/Type/);
   my $line=$_;
   my @F=split/\s+/;
   if ($F[1]==0){
	$F[1]=10**-10;
   }
   my $motif_score=-$F[1]*log($F[1])/log(256);
   $divScore+=$motif_score;

   if ($F[2]==0){
        $F[2]=10**-10;
   }
   my $short_motif_score=-$F[2]*log($F[2])/log(256);
   $short_divScore+=$short_motif_score;

   if ($F[3]==0){
        $F[3]=10**-10;
   }
   my $peak_motif_score=-$F[3]*log($F[3])/log(256);
   $peak_divScore+=$peak_motif_score;

   if ($F[4]==0){
        $F[4]=10**-10;
   }
   my $long_motif_score=-$F[4]*log($F[4])/log(256);
   $long_divScore+=$long_motif_score;

}
$SL_Ratio=$short_divScore/$long_divScore;
print "$sample_type\t$sample\t$divScore\t$short_divScore\t$peak_divScore\t$long_divScore\t$SL_Ratio\n";

