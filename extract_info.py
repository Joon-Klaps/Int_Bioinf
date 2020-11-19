# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 20:00:21 2020

@author: int_bioinf_team


How to use:
     python extract_info.py TTS -i input_file.bed -o output_file.gtf -b pasepairs_to_add -r RPM_cutoff
     python extract_info.py TSS -i input_file.bed -o output_file.gtf -b pasepairs_to_add -r RPM_cutoff
     python extract_info.py SEQ -i input_file.gtf -o output_file.fasta -g genome.fasta
"""
try:
    import sys
    import csv
    import argparse
    import numpy as np
    import pandas as pd
    from Bio import SeqIO
except ImportError:
    print("Something went wrong, maybe this helps (check if the package is installed <pip list>):")
    raise
    
# Define my functions

""" 
Layout of bed_file: 
   [0] genome version, [1]start, [2]end, [3]end of the longest read with this TSS, [4]. ,[5] strand,[6] number of TSS, [7] number of the same reads (same TSS and 3'end), [8]distance, [9]pvalue, [10]oneside95, [11]oneside99   
   Note that if its a pos strand, it started 40bp earlielr (-40bp=start)
             if its a neg strand, it started 40bp later (+40bp=start)                                                                                                                         

"""
   
def main_TSS(input_file, output_file, bpsPromotor, RPM_cutoff):
    bed_file = open(input_file)
    df_bed = csv.reader(bed_file, delimiter= "\t")
    prom_gtf = convertToTSS_GTF(df_bed, int(bpsPromotor))
    bed_file.close()
    
    prom_gtf=filterTNS(prom_gtf, RPM_cutoff)
    writeTSV(prom_gtf,output_file)
    return None

def main_TTS(input_file, output_file,bpsTerminator, RPM_cutoff):
    bed_file = open(input_file)
    df_bed = csv.reader(bed_file, delimiter= "\t")
    term_gtf = convertToTTS_GTF(df_bed, int(bpsTerminator))
    bed_file.close()
    
    term_gtf=filterTNS(term_gtf, RPM_cutoff)
    writeTSV(term_gtf,output_file)
    return None


def main_SEQ(input_file,output_file, genome):
    
    fasta_sequence= SeqIO.read(open(genome), 'fasta')
    with open(output_file, "w") as out, open(input_file, "r") as inp:
        df_gtf=csv.reader(inp, delimiter="\t")
        for i in df_gtf:
            start=int(i[3])
            end =int(i[4])
            out.write(">location[{}_{}] \n]".format( start, end ))
            out.write("{} \n".format(str(fasta_sequence.seq[start:end])))
            
    out.close()
    inp.close()
    return None
    
    
## Extract the TSS from the bed file to a GTF file. 
def convertToTSS_GTF(df_bed, bps) :
    # Define array with locations that have already been searched
    searched_locations_plus=[]
    searched_locations_min=[]
    
    #Define our preliminary GTF file
    result = []
    
    for i in df_bed :
        if i[5] =="+":
            if i[1] not in searched_locations_plus :
                genome=i[0]
                cap="CAPPABLE_SEQ"
                TNS="TSS"
                start=int(i[1])
                end=int(i[1])-bps
                Nreads=int(i[6])
                strand="+"
                values=[genome,cap, TNS, int(end), int(start), int(Nreads), strand, ".", int(Nreads), 0,0]
                
                result.append(values)               
                
                searched_locations_plus.append(i[1])
                
            else : 
                continue
                
                
        else : # if strand is neg
            if i[2] not in searched_locations_min :
                genome=i[0]
                cap="CAPPABLE_SEQ"
                TNS="TSS"
                start=int(i[2])
                end=int(i[2])+bps
                Nreads=int(i[6])
                strand="-"
                values=[genome,cap, TNS, int(start), int(end), int(Nreads), strand, ".", int(Nreads), 0,0]
                
                result.append(values)               
                
                searched_locations_min.append(i[2])
                
            else : 
                continue
            
    df=createPandaDF(result)
    
    Tot=df['Nreads'].sum()
    df=df.assign(TOTREADS=Tot)
    
    df['rpm']=df.apply(lambda x:calculateRPM(x['Nreads'], x['TOTREADS']), axis=1)
    return df

## Extract the TTS from the bed file to a GTF file. 
def convertToTTS_GTF(df_bed, bps):
    
    #Define dictionary to keep track of the amount of reads that end at that location
    dict_TTS_plus={}
    dict_TTS_neg={}
    #Define our preliminary GTF file
    result = []
    
    for i in df_bed: 
        if i[5] =="+":
            
            if i[2] not in dict_TTS_plus.keys():
              
                genome=i[0]
                cap="CAPPABLE_SEQ"
                TNS="TTS"
                start=int(i[2])
                end=int(i[2])-bps
                Nreads=int(i[2]) # For now this is zero, we will update it later with the dictionary value. 
                strand="+"
                values=[genome,cap, TNS, int(end), int(start), int(Nreads), strand, ".", int(Nreads), 0,0]
                
                dict_TTS_plus[i[2]]=int(i[7])
                result.append(values)
                
            else:
               dict_TTS_plus[i[2]] += int(i[7])
               continue
             
        #If strand is negative 
        else :
             if i[1] not in dict_TTS_neg.keys():
                 genome=i[0]
                 cap="CAPPABLE_SEQ"
                 TNS="TTS"
                 start=int(i[1])
                 end=int(i[1])+bps
                 Nreads=int(i[1]) # We set this to the dict value so we can easily replace it later on. 
                 strand="-"
                 values=[genome,cap, TNS, int(start), int(end), int(Nreads), strand, ".", int(Nreads), 0,0]
                    
                 result.append(values)               
                 dict_TTS_neg[i[1]]=int(i[7])
                 
             else :
                 dict_TTS_neg[i[1]]+=int(i[7])
                 continue
             
    df=createPandaDF(result)
    

    df['Nreads']=df.apply(lambda x:equal(x['end'],x['start'], x['strand'],dict_TTS_plus, dict_TTS_neg), axis=1)
    df['Nreads_2']= df['Nreads']
    
    Tot=df['Nreads'].sum()
    df=df.assign(TOTREADS=Tot)
    
    df['rpm']=df.apply(lambda x:calculateRPM(x['Nreads'], x['TOTREADS']), axis=1)
    return df
    
def equal(end,start, strand, dic_plus, dic_neg):
    if strand == "+":
        return dic_plus[str(end)]
    if strand == "-":
        return dic_neg[str(start)]


# Creates the desired GTF format for the perl script to extract the sequences. 
def createPandaDF(matrix):
    df_result=pd.DataFrame(data=np.array(matrix))
    df_result.columns=['genome', 'cap', 'TNS', 'start', 'end', 'Nreads', 'strand', '.', 'Nreads_2', 'TOTREADS', 'rpm']
    df_result[['start', 'end', 'Nreads', 'Nreads_2', 'TOTREADS', 'rpm']]=df_result[['start', 'end', 'Nreads', 'Nreads_2', 'TOTREADS', 'rpm']].astype(int)
    return df_result
    
def calculateRPM(Nreads,TOTREADS):
    return float(Nreads)*1000000/TOTREADS

def filterTNS(df, RPMcutoff):
    return df.loc[df['rpm'] >=float(RPMcutoff)].copy()

def writeTSV(df, output):
    df.to_csv(output, sep='\t', header=False, index= False)
    return None

##--------------Parser
parser = argparse.ArgumentParser()

subparsers = parser.add_subparsers(help='sub-command help', dest = 'mode')

##--------------1.1. main_TSS(input_bed, output_bed, rpm)
parser_1 = subparsers.add_parser('TSS', help='extract the TSS locations from bed file')
parser_1.add_argument('--input','-i', help='input bed file from binomial', dest='input_file')
parser_1.add_argument('--output','-o', help='output GTF file, filtered based on RPM', dest='output_file')
parser_1.add_argument('--BPSP','-b', help='define the amount of bps you wish to extract default is 40',type=int, default=40, dest='bpsPromotor')
parser_1.add_argument('--RPM','-r', help='define the RPM cutoff value default is 5',type=int, default=5, dest='RPM_cutoff')

##--------------1.2. main_TTS(input_file, output_file, output_control, pvalue)
parser_2 = subparsers.add_parser('TTS', help='extract the TSS locations from bed file binomial')
parser_2.add_argument('--input','-i', help='input bed file from binomial', dest='input_file')
parser_2.add_argument('--output','-o', help='output bed file', dest='output_file')
parser_2.add_argument('--BPST','-b', help='define the amount of bps you wish to extract default is 40',type=int, default=40, dest='bpsTerminator')
parser_2.add_argument('--RPM','-r', help='cutoff for selecting the RPM', type=int, default=5, dest='RPM_cutoff')


parser_3 = subparsers.add_parser('SEQ', help='extract the fast sequence  from given gtf file')
parser_3.add_argument('--input','-i', help='input gtf file', dest='input_file')
parser_3.add_argument('--output','-o', help='output fasta file', dest='output_file')
parser_3.add_argument('--genome','-g', help='give the genome fasta file from where to extract sequence', dest='genome')


if __name__ == '__main__':
    args = parser.parse_args()

    if args.mode =='TSS':
        main_TSS(args.input_file, args.output_file,args.bpsPromotor, args.RPM_cutoff)
    elif args.mode == 'TTS':
        main_TTS(args.input_file, args.output_file,args.bpsTerminator, args.RPM_cutoff)
    elif args.mode == 'SEQ':
        main_SEQ(args.input_file, args.output_file,args.genome)