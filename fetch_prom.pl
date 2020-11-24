#!/usr/bin/perl
use strict;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);

my $error_sentence = "USAGE : perl $0\n --tss tss_file.gtf\n --out output file name\nOPTIONAL :\n --bp 40 (in bp DEFAULT 40 bp upstream)";
my $tssfile;
my $bp = 40;
my $OUT;

#================================
#getting the options :
GetOptions ("tss=s" => \$tssfile,    # file from  filter_tss.pl
        "out=s" => \$OUT,
        "bp=s" => \$bp #output file
    ) or die "$error_sentence";



#================================
#checking if everything is fine :
if (!$tssfile || !$OUT) {die "$error_sentence\n";}
if(-e $OUT) { die "File $OUT Exists, please remove old file or rename output file (--out)"};
#=================================

parse_tss($tssfile);

#unlink $generic;

sub parse_tss{
    my ($file)=@_;
    open (OUT, ">$OUT") or die "can't save the output to $OUT\n";
    open (FILE, $file) or die;
    foreach my $line (<FILE>)
    {
        chomp $line;
        my @tmp = split /\t/, $line;
        my $chr = $tmp[0];
        my $tss = $tmp[3];
        my $strand = $tmp[6];
        my $absreads = $tmp[8];
        my $totalreads = $tmp[9];
        my $relreads = $tmp[10];
        my $a = 0;
        my $b = 0;
        
        if ($tss < $bp){
            print OUT
            "";
        }
        elsif ($strand eq "+"){
            $a = $tss - $bp;
            $b = $tss;
            print OUT
                "$chr\tCAPPABLE_SEQ\tTSS\t$a\t$b\t$absreads\t$strand\t.\t$absreads\t$totalreads\t$relreads\n";
        }
        else{
            $a = $tss;
            $b = $tss + $bp;
            print OUT
                "$chr\tCAPPABLE_SEQ\tTSS\t$a\t$b\t$absreads\t$strand\t.\t$absreads\t$totalreads\t$relreads\n";
        }
    }
    close FILE;
    close OUT;

}
