#!/usr/bin/perl -w
use strict;

#This script is to extract the context seq all from ref.

while (@ARGV<2){
   warn "perl $0 <motif_length> <*bed>\n";
   exit 1;
   }

my $extent_size=shift;
my $fa="Homo_sapiens_assembly38.fa"; #human genome
my %hg;
open (my $HG,'<', $fa) or die "Could not open file '$fa' $!";
my $S="";
while (<$HG>)
{
        chomp;
        if ( /^>/ )
        {
                $_ =~ s/>//g;
                my @fields=split;
                $S= $fields[0];
                $hg{$S}="";
                #print "processing $S\t";
        }
        else
        {
                my $str=uc $_;   #transform all lowcase to uppercase
                $hg{$S} .= $str ;
        }
}


#chr1   564760  564881  121     0       +       TA      CT
#chr1   567478  567543  65      0       +       TA  CT


my $for_count=0;
my $rev_count=0;
while(<>){
     chomp;
         my $line=$_;
         my @F=split/\s+/,$line;
#        next if ($F[0]!~/chr[0-9M]/);     #only use autosome
#        next if ($F[0]=~/\_/);
           my $start_forward=$F[1]-$extent_size;
           next if ($F[0] eq "chrM" && $start_forward<0);
           $for_count++;
           my $len=$extent_size;
           my $E5=substr($hg{$F[0]},$F[1],$extent_size);     #4 bp downstream motif
           
           my $x=$F[2]-$extent_size; 
           my $E3=substr($hg{$F[0]},$x,$extent_size);
           
           my $com_seq=$E5.$E3;
           if ($com_seq!~/N/){
              print join("\t", $line, $E5, $E3)."\n";
           }
}


