"""
=========================================================================================================================
Scan_gene_coverage: Algorithm for calculating the coverage depth for each coding gene and the percentage of each gene covered at 10X depth or more using a per-base coverage depth file.
 
Coppée et al.
5WBF: A low-cost and straightforward whole blood filtration method suitable for whole-genome sequencing of Plasmodium falciparum clinical isolates
 
Required: 	- Python 3.5.3
			- BioPython 1.79
            - Bcbio-gff 0.6.6
			
REQUIRED OPTIONS:
Per-base coverage depth file obtained with bedtools genomecov. Path to the input per-base coverage depth file.
Genome sequence reference file (fasta). Path to the input fasta reference file.
Gene coordinate file (gff). Path to the input gff reference file.

Version: 1.0.0 (april. 2021)
Author: Romain Coppée
Contact: romain.coppee@gmail.com
=========================================================================================================================
"""

#!/usr/bin/env python
import argparse
import ast
import logging
import multiprocessing as mp
import subprocess
import os
import sys
import builtins
#import numpy
from Bio import SeqIO
import pprint
from BCBio import GFF
from BCBio.GFF import GFFExaminer


#############################################
#############     CLASS     #################
#############################################

#Definition of Chromosomes
class Chromosome:
    def __init__(self, name, length, seq):
        self.name = name
        self.seq = seq
        self.coverage = [0] * (length+1)
        self.length = length
        self.genes = {}
    def set_genes(self, genes):
        self.genes = genes
    def set_coverage(self, coverage):
        self.coverage = coverage


#############################################
#############     FUNCTIONS     #############
#############################################

#Help and argument parser
def ArgumentsParser():
  parser=argparse.ArgumentParser(description = '''Scan_gene_coverage: Algorithm for calculating the coverage depth for each coding gene and the percentage of each gene covered at 10X depth or more using a per-base coverage depth file.''',
                                 epilog = '''The reads must be aligned on the same version of the genome as the gff provided.''')
  parser.add_argument('-p', '--pileup', type = str, required = True,
                      help = 'Path to a valid per-base coverage depth file. It must have been done on the same reference as the reference file.')
  parser.add_argument('-o', '--output', type = str, required = True,
                      help = 'Path to the output file. Output a coverage depth by coding genes.')
  parser.add_argument('-g', '--gff', type = str, required = True,
                      help = 'Path to a valid .gff reference file.')
  parser.add_argument('-f', '--fasta', type = str, required = True,
                      help = 'Path to a valid .fasta file. Should be the same reference as the .gff file.')
  logger_args = parser.add_argument_group('Logger configuration')
  logger_args.add_argument( '--log_level',
                         type = str,
                         nargs = '?',
                         default = 'INFO',
                         help = 'log level',
                         choices = ['ERROR', 'error', 'WARNING', 'warning', 'INFO', 'info', 'DEBUG', 'debug'],
                         required = False )
  logger_args.add_argument( '--log_file',
                         type = str,
                         nargs = '?',
                         help = 'log file (use the stderr by default)',
                         required = False )
  return(parser.parse_args())

#Read a fasta and gff file for building a dict of Chromosome class
def extract_all_from_gff(gff_file_name, fasta_file_name):
    chromosomes = dict()
    #Record the sequence of chromosome from fasta
    for seq_record in SeqIO.parse(fasta_file_name, "fasta"):
        chromosomes[seq_record.id] = Chromosome(seq_record.id, len(seq_record), seq_record.seq)
    #Record the gene positions from gff
    gff_handle = open(gff_file_name)
    for rec in GFF.parse(gff_handle):
        print(rec.id)
        if rec.id not in chromosomes:
            logger.error(rec.id+" not in fasta. Please check the names of chromosomes: they must be the same in fasta and gff.")
            continue
        genes = dict()
        for gene in rec.features:
            for mRNA in gene.sub_features:
                if mRNA.type != "mRNA":
                    continue
                for exon in mRNA.sub_features:
                    if exon.type == "exon":
                        if gene.id in genes:
                            genes[gene.id] = genes[gene.id]+'@'+str(exon.location.nofuzzy_start)+'_'+str(exon.location.nofuzzy_end)
                        else:
                            genes[gene.id] = str(exon.location.nofuzzy_start)+'_'+str(exon.location.nofuzzy_end)
        chromosomes[rec.id].set_genes(genes)
    gff_handle.close()
    return(chromosomes)

#Record coverage in all the genome
def extract_coverage_from_pileup(pileup_file_name, chromosomes):
    coverage = {}
    with open(pileup_file_name, 'r') as pileup:
        for line in pileup:
            #Record coverage
            record = line.split("\t")
            chromosomes[record[0]].coverage[int(record[1])]=int(record[2])
    return(chromosomes)

def coverage_scan(chromosomes, output_name):
    output = open(output_name, 'w')
    for chromosome in chromosomes:
        sequence = list(chromosomes[chromosome].seq)
        coverage = chromosomes[chromosome].coverage
        genes = chromosomes[chromosome].genes
        for gene in genes:
            ratio = 0
            gene_size = 0
            prop = 0
            exons = genes[gene].split('@')
            for exon in exons:
                list_exon = exon.split("_")
                if len(list_exon) > 2:
                    logger.warning("Error: end, start and another number in exon from "+gene+" Keeping just the 2 firsts.")
                    list_exon = list_exon[0:1]
                for i in range(int(list_exon[0])+1,int(list_exon[1])+1):
                    ratio = ratio+coverage[i]
                    gene_size=gene_size+1
                    if coverage[i] >= 10:
                        prop = prop+1
            prop_gene_cov=prop/gene_size*100
            final_ratio=ratio/gene_size
            #print(final_ratio)
            output.write(gene+"\t"+str(round(final_ratio, 2))+"\t"+str(round(prop_gene_cov, 2))+"\n")


#############################################
#############     MAIN     ##################
#############################################
if __name__ == "__main__":
    args=ArgumentsParser()
    # Logger config
    logging_std_format = '[%(levelname)s] %(message)s'
    logging_debug_format = '%(asctime)s [%(levelname)s] [%(threadName)s - %(name)s] %(message)s'
    log_level = args.log_level.upper()
    if ( log_level == 'DEBUG' ):
        logging_std_format = logging_debug_format
    logging_datefmt = '%d/%m/%Y - %H:%M:%S'
    if ( args.log_file != None ):
        logging.basicConfig( format = logging_std_format,
                             datefmt = logging_datefmt,
                             filename = args.log_file,
                             filemode = 'w',
                             level = log_level )
    else:
        logging.basicConfig( format = logging_std_format,
                             datefmt = logging_datefmt,
                             level = log_level )
    logger = logging.getLogger( 'main' )
    
    chromosomes = dict()
    if args.gff:
        if args.fasta:
            logger.info("Input GFF file: "+args.gff+". Input FASTA file: "+args.fasta)
            chromosomes = extract_all_from_gff(args.gff, args.fasta)
        else:
            logger.error("Error: missing fasta reference input. With -g option you should have a valid -f fasta file.")
            sys.exit(1)
    else:
        logger.error("Error: missing reference input.")
        sys.exit(1)
    #Parse per-base coverage depth file
    chromosomes = extract_coverage_from_pileup(args.pileup, chromosomes)
    #Compute coverage for each gene
    coverage_scan(chromosomes, args.output)
    logger.info("End Program")
