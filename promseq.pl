#!/usr/bin/perl
use strict;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);

my $error_sentence = "USAGE : perl $0\n --prom promotor positions\n --out output file\n --genome genome.fasta \n OPTIONAL:\n --id";
my $prom;
my $genome;
my $OUT;
my $ID = "";

#================================
#getting the options :
GetOptions ("prom=s" => \$prom,
        "out=s" => \$OUT,
        "genome=s" => \$genome,
        "id=s" => \$ID,
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

my $r = 0;
foreach my $line (<PROMS>)
{
    my @tmp = split /\t/, $line;
    
    my $startpos = $tmp[3];
    my $endpos = $tmp[4];
    my $abscoverage = $tmp[5];
    my $CPM = $tmp[10];
    my $strand = $tmp[6];
    my $len = @seq[$startpos..$endpos];
    my $promseq;
    
    if ($strand eq "+"){
        my $promseq = join("",@seq[$startpos..$endpos]);
        print OUT
            ">$ID\[$startpos\_$endpos]$r\n$promseq\n";
    }
    elsif ($strand eq "-"){
        my @comp = reverse @seq[$startpos..$endpos];
        my @original = @comp;
        for (my $i = 0; $i<$len; $i++){
            if ($comp[$i] eq "T"){
                $original[$i] = "A";
            }
            elsif ($comp[$i] eq "A"){
                $original[$i] = "T";
            }
            elsif ($comp[$i] eq "G"){
                $original[$i] = "C";
            }
            elsif ($comp[$i] eq "C"){
                $original[$i] = "G";
            }
        }
        my $promseq = join("",@original);
        print OUT
            ">$ID\[$startpos\_$endpos]$r\n$promseq\n";

    }
    $r = $r+1;
}

close PROMS;
close OUT;
