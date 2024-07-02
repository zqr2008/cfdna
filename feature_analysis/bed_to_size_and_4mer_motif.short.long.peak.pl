
use strict;
use warnings;

# In this version, the fragment length was calculated by end-start, not directly used the length included in bed file
# Both 4mer motif and size will be calculated by this script

if(@ARGV!=2)
{
        warn "usage: <*.rmdup.bed> <sample_prefix>\n";
        exit 1;
}
my %h;
my %tot;
my %auto;
my $auto_n;

my $bed_file=shift;
my $sample_prefix=shift;

my $motif_list="4mer.motif.list";
my @MER;
open IN, $motif_list;
while(<IN>)
{
    chomp;
        next if (/Type/);
        push @MER, $_;
}
close IN;

my $sum_site=0;
my %hash;
my $sum_long_site=0;
my $sum_short_site=0;
my $sum_peak_site=0;
my %hash_long;
my %hash_short;
my %hash_peak;
if ($bed_file=~/\.gz$/){
	open FILE, "gzip -dc $bed_file|";
}else{
	open FILE, "$bed_file";
}
while(<FILE>)
{
        chomp;
        #chr1    13043   13164   121     AGTCCATAGG      GGAGAAGGGG
        my @F=split;
        my $length=$F[2]-$F[1];
        next if ($length>600);  # only focus on DNA fragments below 600bp
        $h{$F[0]}{$length}++;
        $tot{$F[0]}++;
        if($F[0]=~/chr\d+$/)
        {
                $h{chrA}{$length}++;
                $tot{chrA}++;
		$sum_site+=2;
        	my $motif_up=substr($F[4],0,4);       #4mer motif of upstream end
        	my $motif_down=revcom($F[5]);         #reverse and complement
        	$motif_down=substr($motif_down,0,4);  #4mer motif of downstream end

		$hash{$motif_up}++;
	        $hash{$motif_down}++;
		if ($length<=150){
                	$sum_short_site+=2;
                	my $short_motif_up=substr($F[4],0,4);       #4mer motif of upstream end
                	my $short_motif_down=revcom($F[5]);         #reverse and complement
                	$short_motif_down=substr($motif_down,0,4);  #4mer motif of downstream end
                	$hash_short{$short_motif_up}++;
                	$hash_short{$short_motif_down}++;
        	}
        	if ($length>=250){
                	$sum_long_site+=2;
                	my $long_motif_up=substr($F[4],0,4);       #4mer motif of upstream end
                	my $long_motif_down=revcom($F[5]);         #reverse and complement
                	$long_motif_down=substr($motif_down,0,4);  #4mer motif of downstream end
                	$hash_long{$long_motif_up}++;
                	$hash_long{$long_motif_down}++;
        	}
        	if ($length<=170 && $length>=160){
                	$sum_peak_site+=2;
                	my $peak_motif_up=substr($F[4],0,4);       #4mer motif of upstream end
                	my $peak_motif_down=revcom($F[5]);         #reverse and complement
                	$peak_motif_down=substr($motif_down,0,4);  #4mer motif of downstream end
                	$hash_peak{$peak_motif_up}++;
                	$hash_peak{$peak_motif_down}++;
        	}
        }
}
close FILE;

open OUT1, ">$sample_prefix\.size";
for my $len (0..600)
{
        my @res;
        for my $chr ('chrA', 'chrM')
        {
                $h{$chr}{$len}=defined $h{$chr}{$len} ? $h{$chr}{$len} : 0;
                $tot{$chr}=$tot{$chr} || 1;
                my $pct=$h{$chr}{$len}/$tot{$chr}*100;
                push @res, $h{$chr}{$len}, $pct;
        }
        print OUT1 join("\t", $len, @res), "\n";
}

open OUT2, ">$sample_prefix\.4mer.motif.freq";
print OUT2 "Type\t$sample_prefix\n";
foreach my $i(@MER){
    if (exists $hash{$i}){
            my $ratio=$hash{$i}/$sum_site;
            print OUT2 "$i\t$ratio\n";
                                        }
        else{
            print OUT2 "$i\t0\n";
        }
}

open OUT3, ">$sample_prefix\.long.4mer.motif.freq";
print OUT3 "Type\t$sample_prefix.long\n";
foreach my $i(@MER){
    if (exists $hash_long{$i}){
            my $ratio=$hash_long{$i}/$sum_long_site;
            print OUT3 "$i\t$ratio\n";
                                        }
        else{
            print OUT3 "$i\t0\n";
        }
}

open OUT4, ">$sample_prefix\.short.4mer.motif.freq";
print OUT4 "Type\t$sample_prefix.short\n";
foreach my $i(@MER){
    if (exists $hash_short{$i}){
            my $ratio=$hash_short{$i}/$sum_short_site;
            print OUT4 "$i\t$ratio\n";
                                        }
        else{
            print OUT4 "$i\t0\n";
        }
}

open OUT5, ">$sample_prefix\.peak.4mer.motif.freq";
print OUT5 "Type\t$sample_prefix.peak\n";
foreach my $i(@MER){
    if (exists $hash_peak{$i}){
            my $ratio=$hash_peak{$i}/$sum_peak_site;
            print OUT5 "$i\t$ratio\n";
                                        }
        else{
            print OUT5 "$i\t0\n";
        }
}

`paste $sample_prefix\.4mer.motif.freq $sample_prefix\.short.4mer.motif.freq $sample_prefix\.peak.4mer.motif.freq $sample_prefix\.long.4mer.motif.freq >$sample_prefix\.merge.4mer.motif.freq`;
`less $sample_prefix\.merge.4mer.motif.freq |tail -n +2 | awk 'BEGIN{print "Type\tAll\tShort\tPeak\tLong"}{print \$1"\t"\$2"\t"\$4"\t"\$6"\t"\$8}' >$sample_prefix\.merge.4mer.motif.freq2`;
`mv $sample_prefix\.merge.4mer.motif.freq2 $sample_prefix\.merge.4mer.motif.freq`;
`less $sample_prefix\.merge.4mer.motif.freq |tail -n +2 | awk -F '\t' 'BEGIN{print "Type\tAll\tShort\tPeak\tLong\tSLRatio"}{if (\$5==0) {print \$1"\t"\$2"\t"\$3"\t"\$4"\t"\$5"\t""0"}  else{print \$1"\t"\$2"\t"\$3"\t"\$4"\t"\$5"\t"(\$3/\$5)}}'  > $sample_prefix\.merge.4mer.motif.freq2`;
`mv $sample_prefix\.merge.4mer.motif.freq2 $sample_prefix\.merge.4mer.motif.freq`;

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

