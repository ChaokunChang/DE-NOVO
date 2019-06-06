import math
import sys
import os
import numpy as np
from os.path import join as opj
from matplotlib import pyplot as plt
from collections import defaultdict
from Bio import Seq, SeqIO, SeqRecord

class Node():
    def __init__(self, gene):
        self.gene = gene
        self.in_degree = 0
        self.out_degree = 0

class DBG():
    """ The DBG Algorithm to implement the de novo problem. """
    def __init__(self):
        print("Initializing...")

    def load_data(self,data_dir='../data/data1',file_type='short'):
        filenames = os.listdir(data_dir)
        for fid,fname in enumerate(filenames):
            file_prefix = fname.split('.')[0].split('_')[0]
            if file_prefix != file_type:
                continue
            file_path = opj(data_dir,fname)
            reads = SeqIO.parse(file_path,'fasta')
            self.build_graph(reads)

    def build_graph(self, reads_gen):
        assert(isinstance(reads_gen,generator)
        for read in reads_gen: