# Int_Bioinf

Welcome to our github respository of the integrated bioinformatics project. In here, you will be able to find all our scripts and modules used.

Here are some noteworthy github repositories of similar studies:
- SMRT-cappable-seq:
[analysis scripts of Yan et al., 2018](https://github.com/elitaone/SMRT-cappable-seq)
- Cappable-seq:
[analysis scripts of Ettwiller et al., 2016](https://github.com/Ettwiller/TSS/)


## Data cleaning and assembly

### Trimming of reads
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

### Genome assembly
##### Software : NGLMR
CoNvex Gap-cost alignMents for Long Reads (ngmlr) is a long-read mapper designed to sensitively align PacBilo or Oxford Nanopore to (large) reference genomes. Ngmlr uses an SV aware k-mer search to find approximate mapping locations for a read and then a banded Smith-Waterman alignment algorithm to compute the final alignment. Ngmlr uses a convex gap cost model that penalizes gap extensions for longer gaps less than for shorter ones to compute precise alignments. The gap model allows ngmlr to account for both the sequencing error and real genomic variations at the same time and makes it especially effective at more precisely identifying the position of breakpoints stemming from structural variations. The k-mer search helps to detect and split reads that cannot be aligned linearly, enabling ngmlr to reliably align reads to a wide range of different structural variations including nested SVs (e.g. inversions flanked by deletions).
(https://www.nature.com/articles/s41592-018-0001-7)

to download a precompiled NGMLR version on linux operator

```
wget https://github.com/philres/ngmlr/releases/download/v0.2.7/ngmlr-0.2.7-linux-x86_64.tar.gz

tar xvzf ngmlr-0.2.7-linux-x86_64.tar.gz

cd ngmlr-0.2.7/
```

##### Protocol
For Oxford Nanopore we run:
```
ngmlr -t 4 -r reference.fasta -q reads.fastq -o test.sam -x ont
```
where

- reference.fasta is the reference genome we would like to map to

- reads.fasta is the fastq sequence file we would like to map to the reference genome

- test.sam is the output sam file storing mapping information

We transform the SAM files to BAM files, sort the sequences and create an index
```
samtools view -S -b ./SAM_files/trimmed_15_pseudo.sam > trimmed_15_pseudo.bam

samtools sort trimmed_15_pseudo.bam -o trimmed_15_pseudo.sorted.bam

samtools index trimmed_15_pseudo.sorted.bam
```
## TSS and TTS determination
### TSS
Determination of TSS was done using the perl scripts modified from Ettwiller et al. (https://github.com/Ettwiller).
The input for this analysis are sorted mapped reads (.bam).
1) The first script clusters starting sites of reads together within a specified distance, then filters for coverage (both absolute an relative).
2) The second script gives the promotor region for every TSS, defined as a region op X bp directly upstream of the TSS.
3) The third script looks the promotor regions up in the genome and returns a fasta file with the sequences.

Note: all the output files from all these scripts can be viewed as a track in IGV, which I recommend :)

#### Protocol
##### Setting up
```bash
cd Desktop/
```
##### 1) Predicting TSS
```bash
perl TSS.pl --in sorted.bam --out TSS.gtf --combine 20 --filter 5 --rpm 5
```
REQUIRED: <br>
--in: this is the input, a sorted bam file <br>
--out: his is the prefered output file containig the positions of all TSS, including coverage (absolute and in reads per million) <br>
OPTIONAL: <br>
--combine: the distance in bp where start positions are merged; default = 20 <br>
--filter: the minimal absolute coverage for a start position to be included; default = 5 <br>
--rpm: the minimal relative coverage for a start position to be included (in counts per million); default = 5 <br>
NOTE:<br>
Combining restrictions for both absolute and relative coverage assures consistence and accuracy.<br>
The new version of the script has an additional filter: the total coverage must be at least 1.5x higher downstream of the TSS compared to upstream. This difference is averaged over a window of 15 bp up- and downstream of the predicted TSS. <br>
A version that does not use this extra filter is available as 'TSS_quick.pl", which is a bit faster and almost equally accurate on bacterial RNAseq libraries. <br>
<br>

##### 2) Extract promotor regions
```bash
perl prompos.pl --tss tss.gtf --out promotors.gtf --bp 40
```
REQUIRED:<br>
--tss: this is the input, which is the output from TSS.pl, a .gtf file with all the TSS sites<br>
--out: this is the prefered output file containing the positions of all promotor regions of found TSS<br>
OPTIONAL:<br>
--bp: the region directly upstream of the TSS in bp; default = 40<br>

##### 3) Lookup promotor regions in genome
```bash
perl promseq.pl --prom promotors.gtf --out promotor_seqs.fasta --genome genome.fasta --id "t0"
```
REQUIRED:<br>
--prom: this is the input, which is the output from prompos.pl, a gtf file with all the promotor regions<br>
--out: this is the prefered output file containing the sequences of all promotor regions for each identified TSS in fasta format<br>
--genome: this is the reference genome in fasta format<br>
OPTIONAL:<br>
--id: adds and ID (string) to the front of the fasta identifier

### TTS
#### Protocol (retired)

First the BAM files are converted to BED files using the tool bamtobed from the Bedtools software.

```
Installing Bedtools
mv bedtools.static.binary bedtools
chmod a+x bedtools
```

Converting BAM to BED files

```
bedtools bamtobed [OPTIONS] -i <reads.BAM> > reads.bed
```
#### Protocol

#### 1) Predicting TTS for Pseudomonas
```bash
perl TTS_pseudo.pl --in sorted.bam --out TTS.gtf --ratio 1.5 (#the following are TSS parameters:) --combine 20 --filter 5 --rpm 5
```
REQUIRED: <br>
--in: this is the input, a sorted bam file <br>
--out: his is the prefered output file containig the positions of all TTS, including coverage (absolute and in reads per million) <br>
OPTIONAL: <br>
--ratio: minimal ratio of coverage before and after TTS (over a window of 15bp)
TSS PARAMETERS: # Because the scripts use TSS in the identification of TTS <br>
--combine: the distance in bp where start positions are merged; default = 20 <br>
--filter: the minimal absolute coverage for a start position to be included; default = 5 <br>
--rpm: the minimal relative coverage for a start position to be included (in counts per million); default = 5 <br>
NOTE:<br>
1) The script works well only on bacterial libraries.
2) Due to a bug, the script currently writes the same entry (TTS) multiple times to the output file. This can cause downstream problems. Hence, a simple solution is to sort the output file using the following bash command line script: <br>
```bash
sort -u -n -k4,4 TTS.gtf > TTS_unique.gtf  
```
<br>
#### 2) Predicting TTS for LUZ7
```bash
perl TTS_LUZ.pl --in sorted.bam --out TTS.gtf --ratio 1.5 (#the following are TSS parameters:) --combine 20 --filter 5 --rpm 5
```
REQUIRED: <br>
--in: this is the input, a sorted bam file <br>
--out: his is the prefered output file containig the positions of all TTS, including coverage (absolute and in reads per million) <br>
OPTIONAL: <br>
--ratio: minimal ratio of coverage before and after TTS (over a window of 15bp)
TSS PARAMETERS: # Because the scripts use TSS in the identification of TTS <br>
--combine: the distance in bp where start positions are merged; default = 20 <br>
--filter: the minimal absolute coverage for a start position to be included; default = 5 <br>
--rpm: the minimal relative coverage for a start position to be included (in counts per million); default = 5 <br>
NOTE:<br>
The script works well only on dense, viral genomes. <br>
#### 3) Getting terminator site positions
To get the terminator regions, the same scripts to get the promotor regions can be used.
```bash
perl prompos.pl --tss tts.gtf --out terminators.gtf --bp 40
```
REQUIRED:<br>
--tss: this is the input, which is the output from TTS.pl, a .gtf file with all the TTS sites<br>
--out: this is the prefered output file containing the positions of all terminator regions of found TTS<br>
OPTIONAL:<br>
--bp: the region directly upstream of the TTS in bp; default = 40, but recommended to be changed<br>
#### 4) Looking up terminator regions in genome
```bash
perl promseq.pl --prom terminators.gtf --out terminator_seqs.fasta --genome genome.fasta --id "t0"
```
REQUIRED:<br>
--prom: this is the input, which is the output from prompos.pl, a gtf file with all the terminator regions<br>
--out: this is the prefered output file containing the sequences of all terminator regions for each identified TTS in fasta format<br>
--genome: this is the reference genome in fasta format<br>
OPTIONAL:<br>
--id: adds and ID (string) to the front of the fasta identifier









