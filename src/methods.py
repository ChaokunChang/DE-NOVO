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
        filenames = 
        for fid,fname in enumerate(filenames):
            if filenames.split('_','.')[0] != file_type:
                continue
            file_path = opj(data_dir,fname)
            reads = SeqIO.parse(file_path,'fasta')





if __name__=="__main__":
    dbg = DBG()
    print("Test Done.")
    
    lst = [1,2,3,4]
    print(lst)