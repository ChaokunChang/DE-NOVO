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
sys.setrecursionlimit(100000) # set the maximum depth as 1500

def longest_path(graph,cur_p,visited,count):
    """ 
    graph:
    cur_p: a kmer in graph
     """
    best_ans = ""
    best_sp = ""
    for sp in graph[cur_p][1]:
        if sp not in visited:
            # print(len(visited))
            if count >= 3600 :
                break
            visited.add(sp)
            next_ans = sp[-1]+longest_path(graph,sp,visited,count+1)            
            if len(next_ans) > len(best_ans):
                best_ans = next_ans
                best_sp = sp
            visited.remove(sp)
    if best_sp != "":
        visited.add(best_sp)
    return best_ans

def graph_mining(graph, discription:dict, mode="simple"):
    candidate = sorted(discription['in_degree'].items(),key= lambda x : x[1])
    done = set()
    output = []
    ans = ""
    k = discription['k']
    print("Mining graph in {} mode".format(mode))
    print("Candidate size:{}".format(len(candidate)))
    count = 0
    if mode == "simple" or mode == "default":
        for kmer in graph:
            graph[kmer][1] = sorted(graph[kmer][1], key=lambda x: -discription['count_dict'][x])
        for kmer,degree in candidate:
            if degree == 0:
                count+=1
            if kmer in done:
                continue
            done.add(kmer)
            ans = kmer
            terminate = False
            onedone = set()
            onedone.add(kmer)
            while not terminate:
                terminate = True
                for nxtkmer in graph[kmer][1]:
                    if nxtkmer not in done and nxtkmer not in onedone:
                        ans += nxtkmer[-1]
                        terminate = False
                        kmer = nxtkmer
                        if degree > 0:
                            done.add(nxtkmer)
                        onedone.add(nxtkmer)
                        break
            output.append(ans)
    elif mode == "longest":
        for kmer,degree in candidate:
            if kmer in done:
                continue
            done.add(kmer)
            longest = kmer
            k = len(kmer)
            while(len(longest)<10000):
                tmp = longest_path(graph,longest[-k:],done,0)
                longest = longest + tmp
                if len(tmp) < 1:
                    # print("break:{}".format(len(tmp)))
                    break
                else:
                    print("OK:{}".format(len(tmp)))
            print(len(longest))
            output.append(longest)
    elif mode=="compress":
        # graph is a collection of class Node
        nodes = discription['nodes']
        for kmer,degree in candidate:
            if degree == 0:
                count += 1
            if kmer in done:
                continue
            done.add(kmer)
            ans = kmer
            terminate = False
            onedone = set()
            onedone.add(kmer)
            while not terminate:
                terminate = True
                for nxtkmer in graph[kmer].childs:
                    nxtkmer = nodes[nxtkmer]
                    if nxtkmer not in done and nxtkmer not in onedone:
                        ans += nxtkmer[k-1:]
                        terminate = False
                        kmer = nxtkmer
                        if degree > 0:
                            done.add(nxtkmer)
                        onedone.add(nxtkmer)
                        break
            output.append(ans)
    else:
        pass
    print("Count degree 0 {}".format(count))
    return output
    
def kmerize(seq,k,count_dict, graph, in_degree, out_degree):
    for i in range(len(seq)-k +1 -1):
        kmer = seq[i:i+k]
        nxtkmer = seq[i+1:i+1+k]
        if kmer not in count_dict:
            graph[kmer] = [[],[]] # 0 is in and 1 is out
            out_degree[kmer] = 0
            in_degree[kmer] = 0
            count_dict[kmer] = 0
        if nxtkmer not in count_dict:
            graph[nxtkmer] = [[],[]]
            out_degree[nxtkmer] = 0
            in_degree[nxtkmer] = 0
            count_dict[nxtkmer] = 0
        
        if nxtkmer not in graph[kmer][1]:
            graph[kmer][1].append(nxtkmer)
        if kmer not in graph[nxtkmer][0]:
            graph[nxtkmer][0].append(kmer)

        out_degree[kmer] += 1
        in_degree[nxtkmer] += 1
        
        count_dict[kmer] += 1
    count_dict[seq[-k:]] += 1

    

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
                kmerize(seq,k,count_dict,graph,in_degree,out_degree)
                kmerize(twin(seq),k,count_dict,graph,in_degree,out_degree)
        # print("The low frequency kmers (limited up to {}) will be remove from dict.".format(self.limit))
        # low_freq_kmers = [x for x in count_dict if count_dict[x] <= self.limit]
        # for x in low_freq_kmers:
        #     del count_dict[x]
        print("Size of count dict:{}".format(len(count_dict)))
        print("Number of the graph nodes:{}".format(len(graph)))
        self.graph = graph
        self.in_degree = in_degree
        self.out_degree = out_degree
        self.count_dict = count_dict
        # count_list = self.count_dict.values()
        # plt.hist(count_list,bins=100)
        # plt.show()
    
    def graph_simplify(self):
        old_graph = self.graph
        new_graph = {}
        new_in = defaultdict(int)
        new_out = defaultdict(int)
        new_dict = defaultdict(int)
        done = set()
        k = self.k
        node_num = 0
        nodes = []
        count = 0
        for kmer in self.count_dict:
            if kmer in done:
                continue
            count += 1
            new_node = kmer
            size = 0
            node_num += 1
            done.add(kmer)
            cur_kmer = kmer
            # back_ans = []
            # for_ans = []
            while len(old_graph[kmer][1]) == 1:
                nxtkmer = old_graph[kmer][1][0]
                if len(old_graph[nxtkmer][0]) == 1:
                    done.add(nxtkmer)
                    # for_ans += [nxtkmer]
                    new_node += nxtkmer[-1]
                    kmer = nxtkmer
                    size += 1
                else:
                    break
            print(size)
            kmer = cur_kmer
            while len(old_graph[kmer][0]) == 1:
                befkmer = old_graph[kmer][0][0]
                if len(old_graph[befkmer][1]) == 1:
                    done.add(befkmer)
                    # back_ans = [befkmer] + back_ans
                    new_node = befkmer[0] + new_node
                    kmer = befkmer
                    size += 1
                else:
                    break
            # ans = back_ans + [cur_kmer] + for_ans
            # ans = contig2str(ans)
            # if (ans != new_node):
            #     for i in range(len(ans)):
            #         if ans[i] != new_node[i]:
            #             print("{},{}:{}".format(count,i,ans[i]))
            # assert(ans == new_node)
            # new_node = ans
            if new_node not in nodes:
                nodes.append(new_node)
                new_graph[new_node] = [[],[]]
                new_dict[new_node] = 0
                new_out[new_node] = 0
                new_in[new_node] = 0
            new_dict[new_node] += size
            # print(len(new_node))
        graph_sz = len(new_graph)

        for i in range(graph_sz):
            for j in range(graph_sz):
                if i == j:
                    continue
                node1 = nodes[i]
                node2 = nodes[j]
                # print("node1:",node1[-k:])
                # print("node2:",node2[:k])
                if node1[k-1:] == node2[:k-1]:
                    new_graph[node1][1].append(node2)
                    new_out[node1] += 1
                    new_graph[node2][0].append(node1)
                    new_in[node2] += 1
        self.graph = new_graph
        self.in_degree = new_in
        self.out_degree = new_out
        self.count_dict = new_dict
        print("Compressed Graph size:{}".format(len(self.graph)))
        return nodes

    def fit(self,data_dir='../data/data1',file_type='short',mode="simple"):
        self.load_data(data_dir,file_type)
        return self.get_answers(mode)

    def get_answers(self,mode):
        nodes = None
        if mode == "compress":
            nodes,self.graph,self.in_degree,self.out_degree = graph_compress(self.k,self.count_dict,self.graph)
        discription = { 'in_degree':self.in_degree,'count_dict':self.count_dict,
                        'out_degree':self.out_degree,'nodes':nodes,'k':self.k }
        output = graph_mining(self.graph,discription,mode)
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
