import math
import sys
import os
import argparse
import numpy as np
from os.path import join as opj
from Bio import Seq
bases = "ATCG" # base pair: A-T C-G
def parse_args():
    """
    Parses command line arguments.
    """
    parser = argparse.ArgumentParser('DNA Concantatente')
    model_settings = parser.add_argument_group('model seetting')
    model_settings.add_argument('--model', choices=['DBG'], default='DBG',
                                help='choose the algorithm to use')
    model_settings.add_argument('--mode',default="default")
    DBG_settings = parser.add_argument_group('DBG model settings')
    DBG_settings.add_argument('--k', type=int, default=50,
                                help='the kmers-k')
    DBG_settings.add_argument('--step', type=int, default=1,
                                help='step between adjacent kmer.')
    DBG_settings.add_argument('--limit', type=int, default=1,
                                help='min num of the words in count dict.')
    DBG_settings.add_argument('--top', type=int, default=0,
                                help='min num of the words in count dict.')

    path_settings = parser.add_argument_group('dir settings')
    path_settings.add_argument('--data_dir', default='../data/data1',
                               help='the data source dir.')
    path_settings.add_argument('--file_type', choices=['short','long'], default='short',
                               help='short or long file.')
    path_settings.add_argument('--result_dir', default='./',
                               help='the data source dir.')
    return parser.parse_args()
    

def twin(kmer):
    return Seq.reverse_complement(kmer)

def front_seq(kmer):
    for x in bases:
        yield kmer[1:]+x

def back_seq(kmer):
    for x in bases:
        yield x + kmer[:-1]

def step_front(kmer,base):
    return kmer[1:] + base

def step_back(base,kmer):
    return base + kmer[:-1]

def contig2str(contig):
    return contig[0] + "".join(x[-1] for x in contig[1:])