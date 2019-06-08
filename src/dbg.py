import math
import sys
import os
import argparse
import numpy as np
from os.path import join as opj
from matplotlib import pyplot as plt
from collections import defaultdict
from Bio import Seq, SeqIO, SeqRecord
from itertools import chain
from utils import *
# from DNA.src.utils import *
class DBG():
    """ The DBG Algorithm to implement the de novo problem. """
    def __init__(self,k=31,step=1,limit=1):
        print("Initializing...")
        self.k = k
        self.step = step
        self.limit = limit
        self.count_dict = None
        self.nodes = None
        self.graph = None
        self.in_degree = None
        self.out_degree = None

    def load_data(self,data_dir='../data/data1',file_type='short'):
        filenames = os.listdir(data_dir)
        rlist = []
        for fid,fname in enumerate(filenames):
            file_prefix = fname.split('.')[0].split('_')[0]
            if file_prefix != file_type:
                continue
            file_path = opj(data_dir,fname)
            reads = SeqIO.parse(file_path,'fasta')
            rlist.append(reads)
        print("Load the sequences in {} done".format(data_dir))
        print("Building counting dict for them...")
        self.build_graph(rlist)

    def build_graph(self,reads_list):
        assert(isinstance(reads_list,list))
        graph={}
        k = self.k
        count_dict = defaultdict(int)
        in_degree = defaultdict(int)
        out_degree = defaultdict(int)
        for r_id,reads in enumerate(reads_list):
            for read in reads:
                seq = str(read.seq)
                if r_id == 1:
                    seq = twin(seq)
                # print(seq)
                for i in range(len(seq)-k +1 -1,self.step):
                    kmer = seq[i:i+k]
                    nxtkmer = seq[i+1:i+1+k]
                    if count_dict[kmer] == 0:
                        graph[kmer] = [[],[]] # 0 is in and 1 is out
                    if count_dict[nxtkmer] == 0:
                        graph[nxtkmer] = [[],[]]
                    graph[kmer][1].append(nxtkmer)
                    out_degree[kmer] += 1
                    graph[nxtkmer][0].append(kmer)
                    in_degree[nxtkmer] += 1
                    count_dict[kmer] += 1
                    count_dict[nxtkmer] += 1
                
                seq = twin(seq)
                for i in range(len(seq)-k +1 -1):
                    kmer = seq[i:i+k]
                    nxtkmer = seq[i+1:i+1+k]
                    if count_dict[kmer] == 0:
                        graph[kmer] = [[],[]] # 0 is in and 1 is out
                    if count_dict[nxtkmer] == 0:
                        graph[nxtkmer] = [[],[]]
                    graph[kmer][1].append(nxtkmer)
                    out_degree[kmer] += 1
                    graph[nxtkmer][0].append(kmer)
                    in_degree[nxtkmer] += 1
                    count_dict[kmer] += 1
                    count_dict[nxtkmer] += 1
                count_dict[seq[-k:]] += 1
                count_dict[seq[:k]] += 1
                # all count is 2 times
        
        # print("The low frequency kmers (limited up to {}) will be remove from dict.".format(self.limit))
        # low_freq_kmers = [x for x in count_dict if count_dict[x] <= self.limit]
        # for x in low_freq_kmers:
        #     del count_dict[x]
        self.graph = graph
        self.in_degree = in_degree
        self.out_degree = out_degree
        self.count_dict = count_dict
    
    def fit(self,data_dir='../data/data1',file_type='short',mode="default"):
        self.load_data(data_dir,file_type)
        return self.get_answers()

    def get_answers(self):
        G = self.graph
        dic = self.count_dict
        # for kmer in G:
        #     G[kmer][1] = sorted(G[kmer][1], key=lambda x: self.in_degree[x])
        candidate = sorted(self.in_degree.items(),key= lambda x : x[1])
        done = set()
        output = []
        ans = ""
        for kmer,degree in candidate:
            if kmer in done:
                continue
            done.add(kmer)
            ans = kmer
            terminate = False
            while not terminate:
                terminate = True
                for nxtkmer in G[kmer][1]:
                    if nxtkmer not in done:
                        ans += nxtkmer[-1]
                        terminate = False
                        kmer = nxtkmer
                        done.add(nxtkmer)
                        break
            output.append(ans)
        print("Number of results:{}".format(len(output)))
        # print(output)
        return output


if __name__ == "__main__":
    args = parse_args()
    dbg = DBG(k=args.k,step=args.step,limit=args.limit)
    res = dbg.fit(args.data_dir,args.file_type,args.mode)

    if args.top == 0:
        n = len(res)
    else:
        n = args.top
    seq_len = []
    for i,seq in enumerate(res):
        seq_len.append((i,len(seq)))
    seq_len = sorted(seq_len,key=lambda x: -x[1])
    ans = [res[x[0]] for x in seq_len[:n]]
    # print(ans)
    with open(opj(args.result_dir ,"results100.fasta") ,'w') as fout:
        for i,seq in enumerate(ans):
            fout.write(">short_read_{}/1".format(i))
            fout.write("\n")
            fout.write(seq)
            fout.write("\n")
