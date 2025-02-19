import sys
import argparse
import host_profiler as hp
import os
import pathogenprofiler as pp
from uuid import uuid4
import json
import glob
import atexit
import logging
from rich_argparse import ArgumentDefaultsRichHelpFormatter
from rich.logging import RichHandler
from packaging.version import Version
from pydantic import BaseModel
from typing import List
from .profile import profile
from .collate import collate

__softwarename__ = "human-profiler"
__default_data_dir__ = f'{sys.base_prefix}/share/{__softwarename__}/'

def main():
    #### Argument Parsing ####

    parser = argparse.ArgumentParser(description='Host-Profiler pipeline',formatter_class=ArgumentDefaultsRichHelpFormatter)
    parser.add_argument('--version', action='version', version="Host-Profiler version %s" % hp.__version__)
    parser.add_argument('--no_cleanup', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--logging',type=str.upper,default="INFO",choices=["DEBUG","INFO","WARNING","ERROR","CRITICAL"],help='Logging level')
    subparsers = parser.add_subparsers(help="Task to perform")


    ###### Profile #######
    parser_sub = subparsers.add_parser('profile', help='Run whole profiling pipeline', formatter_class=ArgumentDefaultsRichHelpFormatter)
    input=parser_sub.add_argument_group("Input options")
    group = input.add_mutually_exclusive_group(required=True)
    group.add_argument('--read1','-1',help='First read file')
    input.add_argument('--read2','-2',help='Second read file')
    group.add_argument('--bam','-a',help='BAM file. Make sure it has been generated using the H37Rv genome (GCA_000195955.2)')
    group.add_argument('--fasta','-f',help='Fasta file')
    group.add_argument('--vcf','-v',help='VCF file')
    input.add_argument('--platform','-m',choices=["illumina","nanopore"],default="illumina",help='NGS Platform used to generate data')
    input.add_argument('--db',help='Mutation panel name')

    output=parser_sub.add_argument_group("Output options")
    output.add_argument('--prefix','-p',default="hprofiler",help='Sample prefix for all results generated')
    output.add_argument('--dir','-d',default=".",help='Storage directory')
    output.add_argument('--csv',action="store_true",help="Add CSV output")
    output.add_argument('--txt',action="store_true",help="Add text output")
    output.add_argument('--add_columns',default=None,type=str,help="Add additional columns found in the mutation database to the text and csv results")
    output.add_argument('--add_mutation_metadata',action="store_true",help=argparse.SUPPRESS)
    output.add_argument('--call_whole_genome',action="store_true",help="Call whole genome")

    filtering=parser_sub.add_argument_group("Variant filtering options")
    filtering.add_argument('--depth',default="0,10",type=str,help='Minimum depth hard and soft cutoff specified as comma separated values')
    filtering.add_argument('--af',default="0,0.1",type=str,help='Minimum allele frequency hard and soft cutoff specified as comma separated values')
    filtering.add_argument('--strand',default="0,3",type=str,help='Minimum read number per strand hard and soft cutoff specified as comma separated values')
    filtering.add_argument('--sv_depth',default="0,10",type=str,help='Structural variant minimum depth hard and soft cutoff specified as comma separated values')
    filtering.add_argument('--sv_af',default="0.5,0.9",type=str,help='Structural variant minimum allele frequency hard cutoff specified as comma separated values')
    filtering.add_argument('--sv_len',default="100000,50000",type=str,help='Structural variant maximum size hard and soft cutoff specified as comma separated values')


    algorithm=parser_sub.add_argument_group("Algorithm options")
    algorithm.add_argument('--mapper',default="bwa", choices=["bwa","minimap2","bowtie2","bwa-mem2"],help="Mapping tool to use. If you are using nanopore data it will default to minimap2",type=str)
    algorithm.add_argument('--caller',default="freebayes", choices=["bcftools","gatk","freebayes"],help="Variant calling tool to use",type=str)
    algorithm.add_argument('--barcode_caller',default="mpileup", choices=["mpileup","bcftools","gatk","freebayes"],help="Variant calling tool to use",type=str)
    algorithm.add_argument('--kmer_counter',default="kmc", choices=["kmc","dsk"],help="Kmer counting tool to use",type=str)
    algorithm.add_argument('--coverage_tool',default="samtools", choices=["bedtools","samtools"],help="Coverage  tool to use",type=str)
    algorithm.add_argument('--calling_params',type=str,help='Override default parameters for variant calling')
    algorithm.add_argument('--species_only',action="store_true",help="Predict species and quit")
    algorithm.add_argument('--no_species',action="store_false",dest="run_species",help="Skip species prediction")
    algorithm.add_argument('--no_trim',action="store_true",help="Don't trim files using trimmomatic")
    algorithm.add_argument('--no_coverage_qc',action="store_true",help="Don't collect coverage statistics")
    algorithm.add_argument('--no_clip',action="store_false",help="Don't clip reads")
    algorithm.add_argument('--no_delly',action="store_true",help="Don't run delly")
    algorithm.add_argument('--no_mash',action="store_true",help="Don't run mash if kmers speciation fails")
    algorithm.add_argument('--no_samclip',action="store_true",help="Don't run mash if kmers speciation fails")
    algorithm.add_argument('--output_kmer_counts',action="store_true",help=argparse.SUPPRESS)
    algorithm.add_argument('--add_variant_annotations',action="store_true",help=argparse.SUPPRESS)
    algorithm.add_argument('--threads','-t',default=1,help='Threads to use',type=int)
    algorithm.add_argument('--ram',default=8,help='Max memory to use',type=int)

    other=parser_sub.add_argument_group("Other options")
    other.add_argument('--verbose',default=0, choices=[0,1,2],help="Verbosity increases from 0 to 2",type=int)
    other.add_argument('--temp',help="Temp firectory to process all files",type=str,default=".")
    other.add_argument('--version', action='version', version="Host-Profiler version %s" % hp.__version__)
    other.add_argument('--no_cleanup',action="store_true",help="Don't remove temporary files on error")
    other.add_argument('--supplementary_bam',help=argparse.SUPPRESS)
    other.add_argument('--delly_vcf',help=argparse.SUPPRESS)
    other.add_argument('--barcode_stdev',type=float,default=0.15,help=argparse.SUPPRESS)
    other.add_argument('--barcode_snps','--barcode-snps',help='Dump barcoding mutations to a file')
    other.add_argument('--db_dir',type=os.path.abspath,default=__default_data_dir__,help='Storage directory')
    other.add_argument('--logging',type=str.upper,default="INFO",choices=["DEBUG","INFO","WARNING","ERROR","CRITICAL"],help='Logging level')
    other.add_argument('--debug',action="store_true",help="Debug mode")

    parser_sub.set_defaults(func=profile)


    ### new subparser for collate ###

    parser_sub = subparsers.add_parser('collate', help='Collate results from multiple samples', formatter_class=ArgumentDefaultsRichHelpFormatter)
    parser_sub.add_argument('--out',type=str,help='File with samples',required = True)
    parser_sub.add_argument('--samples',type=str,help='File with samples')
    parser_sub.add_argument('--dir',default=".",type=str,help='Directory containing results')
    parser_sub.add_argument('--suffix',default=".results.json",type=str,help='File suffix')
    parser_sub.add_argument('--db_dir',type=os.path.abspath,default=__default_data_dir__,help='Storage directory')
    parser_sub.add_argument('--db',default='malaria',help='Mutation panel name')
    parser_sub.add_argument('--temp',help="Temp firectory to process all files",type=str,default=".")
    parser_sub.add_argument('--min-depth',default=10,type=int,help='Minimum depth to pass')

    parser_sub.set_defaults(func=collate)

    args = parser.parse_args()

    logging.basicConfig(
        level=args.logging, format="%(message)s", datefmt="[%X]", handlers=[RichHandler()]
    )

    if hasattr(args, 'func'):
        args.software_version = hp.__version__
        args.tmp_prefix = 'xxxxx'
        args.files_prefix = f"{args.temp}/{args.tmp_prefix}"
        args.func(args)
    else:
        parser.print_help()
