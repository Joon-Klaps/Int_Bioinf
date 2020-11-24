#!/usr/bin/perl
use strict;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);

my $error_sentence = "USAGE : perl $0\n --in sorted.bam\n --out output.gtf\nOPTIONAL:\n --combine bp to merge TSS\n --filter minimal absolute coverage of TSS\n --rpm minimal reads per million of TSS" ;

#OPTIONS :
my $IN;
my $OUT;
my $COMBINE = 20;
my $FILTER = 5;
my $RPM = 5;
#================================= 
#get options :
GetOptions ("in=s" => \$IN,
            "out=s" => \$OUT,
            "combine=s" => \$COMBINE,
            "filter=s" => \$FILTER,
            "rpm=s" => \$RPM
    ) or die "USAGE : perl $0 $error_sentence";

#=================================
#if something went wrong, notify :
if (!$IN || !$OUT) {die "$error_sentence\n"};

#do not override a file 
if(-e $OUT) { die "File $OUT Exists, please remove old file or rename output file (--out)"};
#================================= 
#start the main program :
my $generic =  clean_name($IN);
my $resulting_bam = $IN;
my $total_reads = `samtools view -c -F4 $resulting_bam`;
chomp $total_reads;

my $file_tmp = new File::Temp( UNLINK => 1 );
my $command = "bedtools bamtobed -cigar  -i $resulting_bam > $file_tmp";
system($command);
my $result = parse_bed($file_tmp);
unlink($file_tmp);
my $total_count = 0;

my $firstbp = new File::Temp( UNLINK => 1 );
open(FIRSTBP, ">$firstbp") or die "can't save result into $OUT\n";

foreach my $chr (sort keys %$result)
{
    my $tmp = $$result{$chr};
    foreach my $keys (sort {$a<=>$b} keys %$tmp)
    {
    if ($$result{$chr}{$keys}{"-"}{"count"})
        {
        my $neg;
        my $relative_neg = sprintf('%.3f',($$result{$chr}{$keys}{"-"}{"count"}/$total_reads)*1000000);
        my $abs_neg = $$result{$chr}{$keys}{"-"}{"count"};
        $neg = $relative_neg;

        if ($neg >= 0)
            {
                my $start = $keys;
                $total_count ++;
                my $id = "TSS"."_".$start."_".$abs_neg."_-_";

                print FIRSTBP "$chr\t$generic\t$id\t$start\t$start\t$relative_neg\t-\t.\tnumber_of_read=$abs_neg;total_number_of_read=$total_reads\n";
            }
        }
    
    if ($$result{$chr}{$keys}{"+"}{"count"})
    {
        my $pos;
        my $relative_pos = sprintf('%.3f',($$result{$chr}{$keys}{"+"}{"count"}/$total_reads)*1000000);
        my $abs_pos = $$result{$chr}{$keys}{"+"}{"count"};
        $pos= $relative_pos;
        if ($pos >= 0)
        {
            my $start = $keys +1;
            $total_count ++;
            my $id = "TSS"."_".$start."_".$abs_pos."_+_";
            print FIRSTBP "$chr\t$generic\t$id\t$start\t$start\t$relative_pos\t+\t.\tnumber_of_read=$abs_pos;total_number_of_read=$total_reads\n";
        }
    }
    }
}
close FIRSTBP;
#unlink($generic);

#==================================

my $generic = new File::Temp( UNLINK => 1 );
my $generic = "generic";
my $command = "bedtools merge -s -c 3 -o collapse -d $COMBINE -i $firstbp > $generic";
system($command);

parse_merged($generic);
unlink($firstbp);
unlink($generic);

#=================================
sub clean_name {
    my ($file)=@_;
    my $generic = $file;
    $generic =~ s/\.tss//;
    $generic =~ s/.*\///g;
    return $generic;
}

sub parse_bed {
    my ($file)=@_;
    open (FILE, $file) or die;
    my %result;
    foreach my $line (<FILE>)
    {
    chomp $line;
    my ($start, $orientation)=undef;
    my @tmp = split /\t/, $line;
    my $chr = $tmp[0];
    $orientation = $tmp[5];
    if ($orientation eq "+")
    {
        $start = $tmp[1];
    }
    elsif ($orientation eq "-"){
        $start = $tmp[2];
    }
    $result{$chr}{$start}{$orientation}{"count"}++;
    }
    return \%result;
}

#====================================

sub parse_merged{
    my ($file)=@_;
    open (OUT, ">$OUT") or die "can't save the output to $OUT\n";
    open (FILE, $file) or die;
    my $chr;
    foreach my $line (<FILE>)
    {
        chomp $line;
        my @tmp = split /\t/, $line;
        $chr = $tmp[0];
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
        
        my $rpms = $total_reads_line*1000000/$total_reads;
        my $chr = $tmp[0];
    
        $highest_pos =~ s/\Q$chr//;
        $highest_pos =~ s/^\_//;
   
        my @col = split /\_/, $highest_pos;
        
        if ($total_reads_line >= $FILTER && $rpms >= $RPM){
            print OUT
                "$chr\tCAPPABLE_SEQ\tTSS\t$col[0]\t$col[0]\t$total_reads_line\t$col[2]\t.\t$total_reads_line\t$total_reads\t$rpms\n";
        }
    }
    close FILE;
    close OUT;

}









