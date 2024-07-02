use strict;
use warnings;

#This verison will use less memory, this version will directly filter the raw bam to the clean bam

if(@ARGV!=3)
{
        warn "usage: <*.bam | sam > <mini_map_qual> <max_mismatch>\n";
        exit 1;
}

my $file=shift;
my $mini_map_qual=shift;
my $mismatch=shift;

my %hash;
my %hash_pair;

#A00402:131:HFVHYDSXY:2:1235:14995:29559 163  chrM  5195  48  150M  =   5211   166  seq   qual  NM:i:2  MD:Z:109T22C17   AS:i:140        XS:i:21 RG:Z:F321

if( $file =~ /bam$/ ) {
        open IN, "samtools view $file |" or die("$!");
} 
else 
{
        open IN, "$file" or die( "$!" );
}

while(<IN>)
{
    chomp;
    my $line=$_;
    my @F=split/\s+/,$line;
    next if ($F[6] ne "=");                  #  the two reads should map to the sample chr
    next if ($F[4] < $mini_map_qual);        #  mapping quality cutoff
    next if ($F[1] & 0x0400 );               #  the read is either a PCR duplicate or an optical duplicate 
#   next unless $F[1] & 0x02;                #  the fragment is propoerly mapped to the reference genome
    next if $F[1] & 0x4;                     #  new: filter unmapped seqment
    next if $F[1] & 0x8;                     #  new: filter "next segment in the template unmapped"
    next if ($F[1] & 0x10 && $F[1] & 0x20);   #  new: filter "both reads are reverse complemented"
    next unless ($F[1] & 0x10 || $F[1] & 0x20); #new:  keep "only one read is reverse complemented"
    next if $F[1] & 0x100;                   #  filter secondary alignment
    next if ($F[8]>600 || $F[8]<-600);       #  filter fragments with length more than 600bp
    next if (exists $hash{$F[0]});           #  filter DNA fragment with many mismatches
    $hash_pair{$F[0]}++;
    if ($hash_pair{$F[0]}==2)          
    {
       delete $hash_pair{$F[0]};       #only store the upaired read to hash_pair for futher filtering
    }
    my $distance=10;   #initial the value
    if ($F[5]!~/^(\d+)M$/)             #  only allow mismatch
    {
        $hash{$F[0]}++;
    }
    if ($line=~/\s+NM:i:(\d+)\s+/)
    {
    $distance=$1;
       if  ($distance > $mismatch)               # filter the fragments with more than the max mismatches num in reads
       {
        $hash{$F[0]}++;
       }
    }
}
close IN;


if( $file =~ /bam$/ ) {
        open IN, "/hwfssz1/ST_HEALTH/P17Z10200N0306/yeqingshi/software/anaconda3/envs/py37/bin/samtools view $file |" or die("$!");
} 
else 
{
        open IN, "$file" or die( "$!" );
}
while(<IN>){
    chomp;
    my $line=$_;
    my @F=split/\s+/,$line;
    next if ($F[6] ne "=");                                        #  the two reads should map to the sample chr
    next if ($F[4] < $mini_map_qual);                              #  mapping quality cutoff
    next if ($F[1] & 0x0400 );                                     #  the read is either a PCR duplicate or an optical duplicate 
#   next unless $F[1] & 0x02;                                      #  the fragment is propoerly mapped to the reference genome
    next if $F[1] & 0x4;                                           #  new: filter unmapped seqment
    next if $F[1] & 0x8;                                           #  new: filter "next segment in the template unmapped"
    next if ($F[1] & 0x10 && $F[1] & 0x20);                        #  new: filter "both reads are reverse complemented"
    next unless ($F[1] & 0x10 || $F[1] & 0x20);                    #  new:  keep "only one read is reverse complemented"
    next if $F[1] & 0x100;                                         #  filter secondary alignment
    next if ($F[8]>600 || $F[8]<-600);                             #  filter fragments with length more than 600bp   
    next if (exists $hash{$F[0]} || exists $hash_pair{$F[0]});     #  filter DNA fragment with many mismatches
    print "$line\n";
}
close IN;


