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
        print("The low frequency kmers (limited up to {}) will be remove from dict.".format(self.limit))
        low_freq_kmers = [x for x in count_dict if count_dict[x] <= self.limit]
        for x in low_freq_kmers:
            del count_dict[x]
        
        return count_dict

    def get_forward_contig(self,kmer):
        """ This function expand every kmer to a forward contig, which every kmer of each contig
            must be 1-in-1-out kmer, and there will be no possible circle or mobius or hairpins
            Causion: We must ensure that k is big enough to avoid the duplicate occurance of a kmer
        """
        contig = [kmer]
        while True:
            if sum( x in self.count_dict for x in front_seq(contig[-1]) ) != 1:
                break # continue when exits only
            next_kmers = [x for x in front_seq(contig[-1]) if x in self.count_dict]
            next_kmer = next_kmers[0]

            if next_kmer == kmer or next_kmer == twin(kmer):
                break # break out when meat circle and mobius contig
            if next_kmer == twin(contig[-1]):
                break # break out when meat hairpins
            if sum(x in self.count_dict for x in back_seq(next_kmer)) != 1:
                break # continue when exits only
            contig.append(next_kmer)
        return contig

    def get_contig(self,kmer):
        """ 
        The contig(expanded kmer) can not form circle itself, so if the fc's(forward_contig) next is kmer, 
        The the bc'next must be twin(kmer), and can form same circle,
        So, the final contig = fc in this case
        While, if kmer is not the next of fc, then we must add the reversed bc before it.
         """
        forward_contig = self.get_forward_contig(kmer)
        backword_contig = self.get_forward_contig(twin(kmer))
        contig = forward_contig
        # if kmer[:-1] == forward_contig[-1][1:]:
        if kmer not in front_seq(forward_contig[-1]):
            revers_back_contig = [twin(x) for x in backword_contig[-1:0:-1]] # the first one is ignored to avoid circle.
            contig = revers_back_contig + contig
        return contig

    def compress_path(self):
        """ Compress the path in graph(Though the graph is not generated till now).
        All contigs in graph will be compressed preliminary(maybe we will compress it later.), 
        each new contig(compressed contigs) if formed by 1-in-1-out contigs,
        and each new_contig's head and tail is not 1-in-1-out.
        The new contigs can be seen as new reads
        It seems that we must consider the new contigs relationship with original reads.
         """
        visited = set()
        self.nodes = []
        for kmer in self.count_dict:
            if kmer not in visited:
                contig = self.get_contig(kmer)
                for km in contig:
                    # avoid duplicated kmers in contigs. 
                    visited.add(km)
                    visited.add(twin(km))
                self.nodes.append(contig2str(contig))
        print(" All contigs are processed(group by the contiguous 1-in-1-out contigs ), \
                New contigs number:{}".format(len(self.nodes)) )
        return self.nodes

    def simple_nodes(self):
        self.nodes = []
        visited = set()
        for id,kmer in enumerate(self.count_dict):
            if kmer not in visited:
                self.nodes.append(kmer)
                self.nodes.append(twin(kmer))
                visited.add(kmer)
                visited.add(twin(kmer))


    def build_graph(self):
        """ Build Graph with the long sequence generated from get_contigs. """
        graph = {}
        heads = {}
        tails = {}
        k = self.k
        for i,sequence in enumerate(self.nodes):
            graph[i] = ([],[])
            heads[sequence[:k]] = (i,'+')
            tails[twin(sequence[-k:])] = (i,'-')

        for i in graph:
            sequence = self.nodes[i]
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
    
    def fit(self,data_dir='../data/data1',file_type='short',mode="default"):
        self.load_data(data_dir,file_type)
        if mode=="default":
            self.compress_path()
        elif mode == "simple":
            self.simple_nodes()
        else:
            pass
        return self.build_graph()

    def get_answers(self):
        G = self.graph
        k = self.k
        nodes = self.nodes
        degree = defaultdict(int)
        for i in G:
            degree[i] += len(G[i][0]) + len(G[i][1])
        candidate = sorted(degree.items(),key= lambda x : x[1])
        done = set()
        output = []
        for i,seq_len in candidate:
            if i in done:
                continue
            done.add(i)
            visited = set()
            visited.add(i)
            ans = nodes[i]
            terminate = False
            while not terminate:
                terminate = True
                for node_id,tag in G[i][0]:
                    if node_id not in done:
                        #done.add(node_id)
                        if node_id in visited:
                            continue
                        visited.add(node_id)
                        terminate = False
                        if tag == '+':
                            ans += nodes[node_id][k-1:]
                        else:
                            ans += twin(nodes[node_id])[k-1:]
                        break
            output.append(ans)
        print("Number of results:{}".format(len(output)))
        # print(output)
        return output


    # def generate_GFA(self,store_path=None):
    #     k = self.k
    #     output = ""
    #     print("H  VN:Z:1.0")
    #     output += "H  VN:Z:1.0" + '\n'
    #     for i,seq in enumerate(self.nodes):
    #         print("S\t%d\t%s\t*"%(i,seq))
    #         output += "S\t%d\t%s\t*"%(i,seq) + '\n'
    #     for i in self.graph:
    #         for j,o in self.graph[i][0]:
    #             print("L\t%d\t+\t%d\t%s\t%dM"%(i,j,o,k-1))
    #             output += "L\t%d\t+\t%d\t%s\t%dM"%(i,j,o,k-1) +'\n'
    #         for j,o in self.graph[i][1]:
    #             print("L\t%d\t-\t%d\t%s\t%dM"%(i,j,o,k-1))
    #             output += "L\t%d\t-\t%d\t%s\t%dM"%(i,j,o,k-1) + '\n'
    #     if store_path is not None:
    #         with open(store_path,'w') as fout:
    #             fout.write(output)

    # def show_graph(self,console_print=False,store_path=None):
    #     lines = ["digraph G {", "graph [rankdir=LR, fontname=\"Courier\"];", "node [shape=record];"]
    #     lines = []
    #     for i in self.graph:
    #         lines.append("%s[label=\"<F> %s+ | <R> %s-\"];"%(i,i,i))
    #     t = dict(('+F','-R'))
    #     for i in self.graph:
    #         for j in self.graph[i]:
    #             o = self.graph[i][j]
    #             lines.append("%s:%s -> %s:%s;"%(i,t[o[0]],j,t[o[1]]))
    #     lines.append("}")
    #     character_graph = '\n'.join(lines)
    #     if store_path is not None:
    #         with open(store_path,'w') as f:
    #             f.write(character_graph)
    #     if console_print:
    #         print(character_graph)
    #     return character_graph



if __name__ == "__main__":
    args = parse_args()
    dbg = DBG(k=args.k,step=args.step,limit=args.limit)
    dbg.fit(args.data_dir,args.file_type,args.mode)
    res = dbg.get_answers()

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
