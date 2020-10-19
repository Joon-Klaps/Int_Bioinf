# Int_Bioinf

Welcome to our github respository of the integrated bioinformatics project. In here, you will be able to keep track of our scripts and modules used. 

## Data cleaning and assembly 

#### trimming of reads
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

#### genome assembly
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
## Secquence analysis 
#### ....

