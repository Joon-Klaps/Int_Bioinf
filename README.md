# Int_Bioinf

Welcome to our github respository of the integrated bioinformatics project. In here, you will be able to keep track of our scripts and modules used. 

## Data cleaning and assembly 

#### Trimming of reads
##### Software and installation 
Porechop software was used; downloaded and installed from https://github.com/rrwick/Porechop via the command line:
```bash 
git clone https://github.com/rrwick/Porechop.git
cd Porechop
python3 setup.py install
porechop -h
```
Next, we can edit the .bash_profile and add the alias line, so we can call upon it with porechop and not the whole path "alias porechop="~/Porechop/porechop-runner.py". 
```bash 
nano ./bash_profile
alias porechop="~/Porechop/porechop-runner.py"
```
##### Protocol
The basic version of the script was used, as follows:
```bash 
porechop -i input_reads.fastq.gz -o output_reads.fastq.gz
```
This script looks for adaptors at the 3’ and 5’ end, as well as known adaptors not completely at the end of the read (can be skipped via: --no_split but I didn’t do this).

#### Genome assembly
##### Software : NGLMR
CoNvex Gap-cost alignMents for Long Reads (ngmlr) is a long-read mapper designed to sensitively align PacBilo or Oxford Nanopore to (large) reference genomes. Ngmlr uses an SV aware k-mer search to find approximate mapping locations for a read and then a banded Smith-Waterman alignment algorithm to compute the final alignment. Ngmlr uses a convex gap cost model that penalizes gap extensions for longer gaps less than for shorter ones to compute precise alignments. The gap model allows ngmlr to account for both the sequencing error and real genomic variations at the same time and makes it especially effective at more precisely identifying the position of breakpoints stemming from structural variations. The k-mer search helps to detect and split reads that cannot be aligned linearly, enabling ngmlr to reliably align reads to a wide range of different structural variations including nested SVs (e.g. inversions flanked by deletions).
(https://www.nature.com/articles/s41592-018-0001-7)

to download a precompiled NGMLR version on linux operator

wget https://github.com/philres/ngmlr/releases/download/v0.2.7/ngmlr-0.2.7-linux-x86_64.tar.gz
tar xvzf ngmlr-0.2.7-linux-x86_64.tar.gz
cd ngmlr-0.2.7/

##### Protocol
For Oxford Nanopore run:

ngmlr -t 4 -r reference.fasta -q reads.fastq -o test.sam -x ont

where
reference.fasta is the reference genome we would like to map to 
reads.fasta is the fastq sequence file we would like to map to the reference genome 
test.sam is the output sam file storing mapping information 


## TSS determination
Determination of TSS was done using the perl scripts modified from Ettwiller et al. (https://github.com/Ettwiller).
The input for this analysis are sorted mapped reads (.bam).
1) The first script indicates the 1st bp of every read from the .bam file.
2) The second script clusters starting sites together within a specified distance, then filters for coverage (both absolute an relative).
3) The third script gives the promotor region for every TSS, defined as a region op X bp directly upstream of the TSS.

Note: all the output files from all these scripts can be viewed as a track in IGV, which I recommend :)

#### Protocol
##### Setting up
```bash
cd Desktop/
```
##### 1) Extracting the first bp positions for every read
```bash
perl firstbase.pl --bam sorted.bam --out firstbasepos.gtf
```
REQUIRED: <br>
--bam: this is the input, a sorted bam file <br>
--out: this is the prefered output file containing the positions of the first bp of every read <br>
<br>
The output of this command is the input for the following:<br>

##### 2) Clustering start sites and filtering for coverage
```bash
perl cluster_filter.pl --tss firstbasepos.gtf --out tss.gtf --combine 20 --filter 10 --rmp 5
```
REQUIRED: <br>
--tss: this is the input, which is the output from the firstbase.pl, a .gtf file with all first bp positions <br>
--out: this is the prefered output file containig the positions of all TSS, including coverage (absolute and in reads per million) <br>
OPTIONAL: <br>
--combine: the distance in bp where start positions are merged; default = 20 <br>
--filter: the minimal absolute coverage for a start position to be included; default = 10 <br>
--rpm: the minimal relative coverage for a start position to be included (in counts per million); default = 5 <br>
NOTE:<br>
Combining restrictions for both absolute and relative coverage assures consistence and accuracy.<br>
<br>
The output of this command is the input for the following:

##### 3) Extract promotor regions
```bash
perl fetch_prom.pl --tss tss.gtf --out promotors.gtf --bp 40
```
REQUIRED:<br>
--tss: this is the input, which is the output from cluster_filter.pl, a .gtf file with all the TSS sites<br>
--out: this is the prefered output file containing the positions of all promotor regions of found TSS<br>
OPTIONAL:<br>
--bp: the region directly upstream of the TSS in bp; default = 40<br>

##### 4) Lookup promotor regions in genome
```bash
perl promseq.pl --prom promotors.gtf --out promotor_seqs.fasta --genome genome.fasta
```
REQUIRED:<br>
--prom: this is the input, which is the output from fetch_prom.pl, a gtf file with all the promotor regions<br>
--out: this is the prefered output file containing the sequences of all promotor regions for each identified TSS in fasta format<br>
--genome: this is the reference genome in fasta format<br>

## Secquence analysis 
### SMRT-cappable-seq
[analysis scripts of smrt](https://github.com/elitaone/SMRT-cappable-seq)
### Cappable-seq
[analysis scripts of Ettwiller](https://github.com/Ettwiller/TSS/)



