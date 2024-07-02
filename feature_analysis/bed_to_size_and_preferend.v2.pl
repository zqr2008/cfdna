
use strict;
use warnings;

# In this version, the fragment length was calculated by end-start, not directly used the length included in bed file
# Both 4mer motif and size will be calculated by this script

if(@ARGV!=2)
{
        warn "usage: <*.rmdup.bed> <sample_prefix>\n";
        exit 1;
}
my $bed_file=shift;
my $sample_prefix=shift;

my @END=("AA","AG","AC","AT","CC","CG","CT","GG","GT","TT");
my $sum_site=0;
my $short_site=0;
my $peak_site=0;
my $long_site=0;
my %hash;
my %short;
my %peak;
my %long;

if ($bed_file=~/\.gz/){
	open FILE, "gzip -dc $bed_file|";
}else{
	open FILE, "$bed_file";
}
while(<FILE>)
{
        chomp;
        #chr1    13043   13164   121     AGTCCATAGG      GGAGAAGGGG
        my @F=split;
	next if ($F[0]!~/chr\d+$/);
        $sum_site+=1;
        my $motif_up_end=substr($F[4],0,1);       #4mer motif of upstream end
        my $motif_down_end=revcom($F[5]);         #reverse and complement
        $motif_down_end=substr($motif_down_end,0,1);  #4mer motif of downstream end
	my $combine=$motif_up_end."".$motif_down_end;
	$combine="AG" if ($combine eq "GA");
	$combine="AC" if ($combine eq "CA");
	$combine="AT" if ($combine eq "TA");
	$combine="CG" if ($combine eq "GC");
	$combine="CT" if ($combine eq "TC");
	$combine="GT" if ($combine eq "TG");
	$hash{$combine}++;
	my $len=$F[2]-$F[1];
	if ($len<=150){
        	$short{$combine}++;
		$short_site++;
     	}elsif($len>=160 && $len<=170){
        	$peak{$combine}++;
		$peak_site++;
     	}elsif($len>=250){
        	$long{$combine}++;
		$long_site++;
     	}

}
close FILE;


open OUT2, ">$sample_prefix\.end.freq";
print OUT2 "Preferend\tAll\n";
foreach my $i(@END){
    if (exists $hash{$i}){
            my $ratio=$hash{$i}/$sum_site;
            print OUT2 "$i\t$ratio\n";
                                        }
        else{
            print OUT2 "$i\t0\n";
        }
}

open OUT3, ">$sample_prefix\.short.end.freq";
print OUT3 "Preferend\tShort\n";
foreach my $i(@END){
    if (exists $short{$i}){
            my $short_ratio=$short{$i}/$short_site;
            print OUT3 "$i\t$short_ratio\n";
                                        }
        else{
            print OUT3 "$i\t0\n";
        }
}
open OUT4, ">$sample_prefix\.peak.end.freq";
print OUT4 "Preferend\tPeak\n";
foreach my $i(@END){
    if (exists $peak{$i}){
            my $peak_ratio=$peak{$i}/$peak_site;
            print OUT4 "$i\t$peak_ratio\n";
                                        }
        else{
            print OUT4 "$i\t0\n";
        }
}
open OUT5, ">$sample_prefix\.long.end.freq";
print OUT5 "Preferend\tLong\n";
foreach my $i(@END){
    if (exists $long{$i}){
            my $long_ratio=$long{$i}/$long_site;
            print OUT5 "$i\t$long_ratio\n";
                                        }
        else{
            print OUT5 "$i\t0\n";
        }
}
`paste $sample_prefix\.end.freq $sample_prefix\.short.end.freq $sample_prefix\.peak.end.freq $sample_prefix\.long.end.freq > $sample_prefix\.merge.end.freq`;
`less $sample_prefix\.merge.end.freq |tail -n +2 | awk 'BEGIN{print "Type\tAll\tShort\tPeak\tLong"}{print \$1"\t"\$2"\t"\$4"\t"\$6"\t"\$8}'> $sample_prefix\.merge.end.freq2`;
`mv  $sample_prefix\.merge.end.freq2 $sample_prefix\.merge.end.freq`;
`less $sample_prefix\.merge.end.freq 	|tail -n +2 | awk -F '\t' 'BEGIN{print "Type\tAll\tShort\tPeak\tLong\tSLRatio"}{if (\$5==0) {print \$1"\t"\$2"\t"\$3"\t"\$4"\t"\$5"\t""0"}  else{print \$1"\t"\$2"\t"\$3"\t"\$4"\t"\$5"\t"(\$3/\$5)}}'  > $sample_prefix\.merge.end.freq2`;
`mv $sample_prefix\.merge.end.freq2  $sample_prefix\.merge.end.freq`;















sub revcom
{
        # A subroutine to compute the reverse complement of DNA sequence 
        my ($tmp)=@_;
        my @dna=split//,$tmp;
        my @revcom=reverse(@dna);
        my $rev=join("",@revcom);
        $rev=~tr/ACGTacgt/TGCAtgca/;
        return $rev;
}

