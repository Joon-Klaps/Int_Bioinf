#!/usr/bin/perl
use strict;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);

my $error_sentence = "USAGE : perl $0\n --prom promotor positions\n --out output file\n --genome genome.fasta";
my $prom;
my $genome;
my $OUT;

#================================
#getting the options :
GetOptions ("prom=s" => \$prom,
        "out=s" => \$OUT,
        "genome=s" => \$genome
    ) or die "$error_sentence";



#================================
#checking if everything is fine :
if (!$prom || !$OUT || !$genome) {die "$error_sentence\n";}
if(-e $OUT) { die "File $OUT Exists, please remove old file or rename output file (--out)"};
#=================================

open (FASTA, $genome);
my @seq;

my $firstline = 1;
foreach my $line (<FASTA>)
{
    if ($firstline == 1){
        $firstline = 0;
        next;
    }
    else{
        chomp $line;
        my @tmp = split("", $line);
        push @seq, @tmp;
    }
}
close FASTA;


open (OUT, ">$OUT");
open (PROMS, $prom);

foreach my $line (<PROMS>)
{
    my @tmp = split /\t/, $line;
    
    my $startpos = $tmp[3];
    my $endpos = $tmp[4];
    my $abscoverage = $tmp[5];
    my $CPM = $tmp[10];
    
    my $promseq = join("",@seq[$startpos..$endpos]);
    
    print OUT
        ">pos[$startpos _ $endpos]\n$promseq\n";
}

close PROMS;
close OUT;
