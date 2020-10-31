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
Determination of TSS was done using the perl scripts by Ettwiller et al. (https://github.com/Ettwiller)
A first script looks for enriched start sites, where a second script combines close (< several bp) into 1 site.
The input for this analysis is a bam file of mapped reads.

#### Protocol
##### Setting up
```bash
cd Desktop/
git clone https://github.com/Ettwiller/TSS.git
```
##### Analysis
```bash
perl ./TSS/bam2firstbasegtf.pl --bam trimmed_5_pseudo.sorted.bam --cutoff 10 --out enriched_cutoff10.gtf
```
A cutoff of 10 - 20 seems to capture most TSS.

The output of this command is the input for the following:

```bash
perl ./TSS/cluster_tss.pl  --tss enriched_cutoff10.gtf --cutoff  5 --out enriched_cutoff_10_cluster_5.gtf
```
This script combines TSS within 5bp into 1 TSS. A cutoff of 5 was used in the paper and seems accurate.



## Secquence analysis 
### SMRT-cappable-seq
[analysis scripts of smrt](https://github.com/elitaone/SMRT-cappable-seq)
### Cappable-seq
[analysis scripts of Ettwiller](https://github.com/Ettwiller/TSS/)



