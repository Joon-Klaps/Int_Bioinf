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
##### Software
If Bowtie2 is not yet installed (`module available` if so load it `module load bowtie2`), Bowtie2 can be installed using conda from the command line via:
```bash 
conda install bowtie2
```

##### Protocol
In order to align to a reference genome, an indexed genome must be built first as follows:
```bash 
Bowtie2-build reference_genome.fasta reference_genome_name
```
To align use the following command:
```bash 
Bowtie2 -x indexed_reference_genome -U reads.fastq.gz
```

Info about [Job-arrays](https://rc.dartmouth.edu/index.php/using-discovery/scheduling-jobs/using-job-arrays/)


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
perl cluster_filter.pl --tss firstbasepos.gtf --out tss.gtf --combine 20 --filter 10 --rpm 5
```
REQUIRED: <br>
--tss: this is the input, which is the output from the firstbase.pl, a .gtf file with all first bp positions <br>
--out: this is the prefered output file containig the positions of all TSS, including coverage (absolute and in reads per million) <br>
OPTIONAL: <br>
--combine: the distance in bp where start positions are merged; default = 20
--filter: the minimal absolute coverage for a start position to be included; default = 10
--rpm: the minimal relative coverage for a start position to be included (in counts per million); default = 5
NOTE:
Combining restrictions for both absolute and relative coverage assures consistence and accuracy.


The output of this command is the input for the following:

##### 3) Extract promotor regions
```bash
perl fetch_prom.pl --tss tss.gtf --out promotors.gtf --bp 40
```
REQUIRED:
--tss: this is the input, which is the output from cluster_filter.pl, a .gtf file with all the TSS sites
--out: this is the prefered output file containing the positions of all promotor regions of found TSS
OPTIONAL:
--bp: the region directly upstream of the TSS in bp; default = 40



## Secquence analysis 
### SMRT-cappable-seq
[analysis scripts of smrt](https://github.com/elitaone/SMRT-cappable-seq)
### Cappable-seq
[analysis scripts of Ettwiller](https://github.com/Ettwiller/TSS/)



