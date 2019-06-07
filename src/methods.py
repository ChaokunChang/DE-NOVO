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

class DBG():
    """ The DBG Algorithm to implement the de novo problem. """
    def __init__(self,k=31,step=1,limit=1,mode="default"):
        print("Initializing...")
        self.k = k
        self.step = step
        self.limit = limit
        self.count_dict = None
        self.results = None
        self.graph = None

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
        self.count_dict = self.build_dict(rlist)

    def build_dict(self,reads_list):
        assert(isinstance(reads_list,list))
        count_dict = defaultdict(int)
        for reads in reads_list:
            for read in reads:
                seq = str(read.seq)
                # print(seq)
                k = self.k
                for i in range(len(seq)-k +1):
                    kmer = seq[i:i+k]
                    count_dict[kmer] += 1
                seq = twin(seq)
                for i in range(len(seq)-k +1):
                    count_dict[seq[i:i+k]] += 1
        print("The kmers which count is less than limit({}) will be remove from dict.".format(self.limit))
        # count_dict = [x for x in count_dict if count_dict[x] > self.limit ]
        del_ele = [x for x in count_dict if count_dict[x] <= self.limit]
        for x in del_ele:
            del count_dict[x]
        
        return count_dict

    def get_forward_contig(self,kmer):
        contig = [kmer]
        while True:
            # if sum((step_front(contig[-1],base) in self.count_dict) for base in bases ) != 1:
            if sum( x in self.count_dict for x in front_seq(contig[-1]) ) != 1:
                break # continue when exits only
            # next_kmers = [step_front(contig[-1],base) for base in bases if step_front(contig[-1],base) in self.count_dict]
            next_kmers = [x for x in front_seq(contig[-1]) if x in self.count_dict]
            next_kmer = next_kmers[0]

            if next_kmer == kmer or next_kmer == twin(kmer):
                break # break out when meat circle and mobius contig
            if next_kmer == twin(contig[-1]):
                break # break out when meat hairpins
            # if sum(step_back(base,next_kmer) in self.count_dict for base in bases ) != 1:
            if sum(x in self.count_dict for x in back_seq(next_kmer)) != 1:
                break # continue when exits only
            contig.append(next_kmer)
        return contig

    def get_contig(self,kmer):
        forward_contig = self.get_forward_contig(kmer)
        backword_contig = self.get_forward_contig(twin(kmer))
        # if kmer[:-1] == forward_contig[-1][1:]:
        if kmer in front_seq(forward_contig[-1]):
            # The kmer is the forward_contig's next contig, means that the while loop
            # break out at the second one(meat circle).
            # The then contig only contain the forward one, add backword one will caused circle.
            contig = forward_contig
        else:
            # There won't form circle, so we must add the backword contig, however,
            # we must reverse the backword contig before.
            revers_back_contig = [twin(x) for x in backword_contig[-1:0:-1]] # the first one is ignored to avoid circle.
            contig = revers_back_contig + forward_contig
        return contig

    def get_contigs(self):
        visited = set()
        self.results = []
        for kmer in self.count_dict:
            if kmer not in visited:
                contig = self.get_contig(kmer)
                for km in contig:
                    visited.add(km)
                    visited.add(twin(km))
                self.results.append(contig2str(contig))
        print("All results are generated, total number:{}".format(len(self.results)))
        return self.results

    def build_graph(self):
        """ Build Graph with the long sequence generated from get_contigs. """
        graph = {}
        heads = {}
        tails = {}
        k = self.k
        for i,sequence in enumerate(self.results):
            graph[i] = ([],[])
            heads[sequence[:k]] = (i,'+')
            tails[twin(sequence[-k:])] = (i,'-')

        for node_id in graph:
            sequence = self.results[node_id]
            seq_tail = sequence[-k:]
            seq_head = sequence[:k]
            for next_seq in front_seq(seq_tail):
                if next_seq in heads:
                    graph[i][0].append(heads[next_seq])
                if next_seq in tails:
                    graph[i][0].append(tails[next_seq])
            for before_seq in front_seq(twin(seq_head)):
                if before_seq in heads:
                    graph[i][1].append(heads[before_seq])
                if before_seq in tails:
                    graph[i][1].append(tails[before_seq])
        self.graph = graph
        return self.graph
    
    def fit(self,data_dir='../data/data1',file_type='short'):
        self.load_data(data_dir,file_type)
        self.get_contigs()
        self.build_graph()
        return self.show_graph()

    def get_longest_path(self):
        pass

    def generate_GFA(self,store_path=None):
        k = self.k
        output = ""
        print("H  VN:Z:1.0")
        output += "H  VN:Z:1.0" + '\n'
        for i,seq in enumerate(self.results):
            print("S\t%d\t%s\t*"%(i,seq))
            output += "S\t%d\t%s\t*"%(i,seq) + '\n'
        for i in self.graph:
            for j,o in self.graph[i][0]:
                print("L\t%d\t+\t%d\t%s\t%dM"%(i,j,o,k-1))
                output += "L\t%d\t+\t%d\t%s\t%dM"%(i,j,o,k-1) +'\n'
            for j,o in self.graph[i][1]:
                print("L\t%d\t-\t%d\t%s\t%dM"%(i,j,o,k-1))
                output += "L\t%d\t-\t%d\t%s\t%dM"%(i,j,o,k-1) + '\n'
        if store_path is not None:
            with open(store_path,'w') as fout:
                fout.write(output)

    def show_graph(self,console_print=False,store_path=None):
        lines = ["digraph G {", "graph [rankdir=LR, fontname=\"Courier\"];", "node [shape=record];"]
        for i in self.graph:
            lines.append("%s[label=\"<F> %s+ | <R> %s-\"];"%(i,i,i))
        t = dict(('+F','-R'))
        for i in self.graph:
            for j in self.graph[i]:
                o = self.graph[i][j]
                lines.append("%s:%s -> %s:%s;"%(i,t[o[0]],j,t[o[1]]))
        lines.append("}")
        character_graph = '\n'.join(lines)
        if store_path is not None:
            with open(store_path,'w') as f:
                f.write(character_graph)
        if console_print:
            print(character_graph)
        return character_graph



if __name__ == "__main__":
    # for x in fw('CCTTAG'):
    #     print(x)
    args = parse_args()
    dbg = DBG(k=args.k,step=args.step,limit=args.limit,mode=args.mode)
    dbg.load_data()
    res = dbg.get_contigs()
    if args.top == 0:
        n = len(res)
    else:
        n = args.top
    seq_len = []
    for i,seq in enumerate(res):
        seq_len.append((i,len(seq)))
    seq_len = sorted(seq_len,key=lambda x: -x[1])
    ans = [res[x[0]] for x in seq_len[:n]]
    print(ans)
    with open(opj(args.result_dir ,"results100.fasta") ,'w') as fout:
        for i,seq in enumerate(ans):
            fout.write(">short_read_{}/1".format(i))
            fout.write("\n")
            fout.write(seq)
            fout.write("\n")
