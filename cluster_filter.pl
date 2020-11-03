#!/usr/bin/perl
use strict;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);

my $error_sentence = "USAGE : perl $0\n --tss tss_file.gtf (the output from firstbase.pl)\n --out output file name\nOPTIONAL :\n --combine 20 (in bp DEFAULT 20 bp merging)\n --filter 10 (minimal absolute coverage DEFAULT 10)\n --rpm 5 (minimal reads per million DEFAULT 5)";
my $tssfile;
my $COMBINE = 20;
my $FILTER = 10;
my $minRPM = 5;
my $OUT;

#================================
#getting the options :
GetOptions ("tss=s" => \$tssfile,
	    "combine=s" => \$COMBINE,
	    "out=s" => \$OUT,
        "filter=s" => \$FILTER,
        "rmp=s" => \$minRPM
    ) or die "$error_sentence";

#================================
#checking if everything is fine :
if (!$tssfile || !$OUT) {die "$error_sentence\n";}
if(-e $OUT) { die "File $OUT Exists, please remove old file or rename output file (--out)"};
#=================================                                                                                                                         

my $generic = new File::Temp( UNLINK => 1 );
my $generic = "generic";
my $command = "bedtools merge -s -c 3 -o collapse -d $COMBINE -i $tssfile > $generic";
system($command);

my $totalreads = -9; # this parses the total number of reads
open (TSS, $tssfile) or die;
foreach my $entry (<TSS>)
{
    chomp $entry;
    my @tmp = split /\t/, $entry;
    if ($totalreads == -9){
        $totalreads = -8;
        next;
    }
    $totalreads = $tmp[9];
    
    last if ($totalreads != -9);
}
close TSS;

parse_merged($generic);
unlink($generic);

#unlink $generic;

sub parse_merged{
    my ($file)=@_;
    open (OUT, ">$OUT") or die "can't save the output to $OUT\n";
    open (FILE, $file) or die;
    foreach my $line (<FILE>)
    {
        chomp $line;
        my @tmp = split /\t/, $line;
        my @positions = split /\,/, $tmp[-1];
        my $merge_size = @positions;
        my @result=[];
        my $total_reads_line = 0;
        for (my $i=0; $i<$merge_size; $i++)
        {
            my $pos = $positions[$i];
            my @coords = split/\_/, $pos;
            my $score = $coords[2];
            $result[$i][0]=$score;
            $result[$i][1]=$pos;
            $total_reads_line = $total_reads_line + $score;
        }           
        my @sorted = sort { $b->[0] <=> $a->[0] } @result;
        my $highest_pos = $sorted[0][1];
        
    my $rpm = $total_reads_line*1000000/$totalreads;
    my $chr = $tmp[0];
    
    $highest_pos =~ s/\Q$chr//;
    $highest_pos =~ s/^\_//;
   
    my @col = split /\_/, $highest_pos;
        
        if ($total_reads_line >= $FILTER && $rpm >= $minRPM){
            print OUT
                "$chr\tCAPPABLE_SEQ\tTSS\t$col[0]\t$col[0]\t$total_reads_line\t$col[2]\t.\t$total_reads_line\t$totalreads\t$rpm\n";
        }
    }
    close FILE;
    close OUT;

}
