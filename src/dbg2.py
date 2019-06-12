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

class PairEnd:
    def __init__(self,length,reads1,reads2,pair_dis=500):
        self.length = length
        self.pairs = {}
        self.pair_dis = pair_dis
        for i in range(len(reads1)):
            self.pairs[reads1[i]] = twin(reads2[i])
            self.pairs[twin(reads2[i])] = reads1[i]

            self.pairs[twin(reads1[i]) ] = reads2[i]
            self.pairs[twin(reads2[i])] = twin(reads1[i])
    
    def contain_pair(self,contig):
        ans = []
        seqset = set()
        locmap = {}
        c_size = len(contig)
        for i in range(c_size - self.length+1):
            subc = contig[i : i+self.length]
            seqset.add(subc)
            locmap[subc] = i
        for read in self.pairs:
            if read in seqset:
                ans.append((read,locmap[read]))
        return ans
    
    def get_pair(self,p1):
        return self.pairs[p1]

def graph_mining(graph, discription:dict, mode="simple"):
    candidate = sorted(discription['in_degree'].items(),key= lambda x : x[1])
    done = set()
    output = []
    ans = ""
    k = discription['k']
    PES = discription['pair_ends']
    print("Mining graph in {} mode".format(mode))
    print("Candidate size:{}".format(len(candidate)))
    count = 0
    if mode == "simple" or mode == "default" or mode=="compress":
        # for kmer in graph:
        #     graph[kmer][1] = sorted(graph[kmer][1], key=lambda x: -discription['count_dict'][x])
        for kmer,degree in candidate:
            if len(graph[kmer][0]) == 1 and len(graph[kmer][1])==1:
                nxt = graph[kmer][1][0]
                bef = graph[kmer][0][0]
                assert(len(graph[nxt][0])>1 and len(graph[bef][1])>1 )
                count += 1
            if degree == 2:
                n1 = graph[kmer][1][0]
                n2 = graph[kmer][1][1]
                # assert(graph[n1][0] == 1 and graph[n1][1] == 1)
                # assert(graph[n2][0] == 1 and graph[n2][1] == 1)
                # assert(graph[n1][1][0] == graph[n2][1][0])
            if degree > 2:
                print("!!!!!!!!!!!!!!!!!!!")
            if kmer in done:
                continue
            done.add(kmer)
            ans = kmer
            terminate = False
            while not terminate:
                terminate = True
                single_pairs = PES.contain_pair(kmer)

                check_pos = []
                for p1,loc1 in single_pairs:
                    endpos = (loc1 + PES.pair_dis) - len(kmer)
                    loc2 = endpos - PES.length
                    if loc2 >= 0:
                        start = loc2 + k-1
                        remain_len = PES.length
                        p2 = PES.get_pair(p1)
                        check_pos.append((p2,start,remain_len))
                    elif endpos >= 0:
                        start = k-1
                        remain_len = -loc2
                        p2 = PES.get_pair(p1)
                        check_pos.append((p2,k-1,remain_len))

                # print(single_pairs)
                selection = ""
                directions = graph[kmer][1]
                # assert(0==1)
                if len(directions) == 1:
                    if directions[0] not in done:
                        selection = directions[0]
                elif len(directions) > 1:
                    for nxtkmer in directions:
                        if nxtkmer in done:
                            continue
                        for p2,start,rl in check_pos:
                            if rl < PES.length:
                                if start+rl <= len(nxtkmer):
                                    tp2 = nxtkmer[start:start+rl]
                                    rl = 0
                                    addlist = [nxtkmer]
                                else:
                                    rl = start + rl - len(nxtkmer)
                                    start = k-1
                            elif rl > 0:
                                if start+rl <= len(nxtkmer):
                                    tp2 = nxtkmer[start:start+rl]
                                    addlist = [nxtkmer]
                                else:
                                    tmp = nxtkmer[start:]
                                    rl = rl - (len(nxtkmer) - start)
                                    start = k-1
                                    nxtnxtkmers = graph[nxtkmer][1]
                                    assert(len(nxtnxtkmers) == 1)
                                    nxtnxtkmer = nxtnxtkmers[0]
                                    tp2 = tmp + nxtnxtkmer[start:start+rl]

                                    addlist = [nxtkmer,nxtnxtkmer]
                                    
                            if tp2 == p2[-rl:]:
                                selection = "".join(addlist)
                                for ele in addlist:
                                    done.add(ele)
                                break
                        if selection != "":
                            print(selection)
                            break
                    # assert(selection != "")

                if selection != "":
                    terminate = False
                    kmer = selection
                    done.add(selection)
                    ans += selection[k-1:]
                
            output.append(ans)
        print("Count 1-in-1-out {}".format(count))
    else:
        pass
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
        self.pair_ends = None

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
        cache1 = []
        cache2 = []
        for r_id,reads in enumerate(reads_list):
            for read in reads:
                seq = str(read.seq)
                kmerize(seq,k,count_dict,graph,in_degree,out_degree)
                kmerize(twin(seq),k,count_dict,graph,in_degree,out_degree)
                if r_id == 0:
                    cache1.append(seq)
                else:
                    cache2.append(seq)
        self.pair_ends = PairEnd(len(seq),cache1,cache2)

        print("Size of count dict:{}".format(len(count_dict)))
        print("Number of the graph nodes:{}".format(len(graph)))
        self.graph = graph
        self.in_degree = in_degree
        self.out_degree = out_degree
        self.count_dict = count_dict
        count = 0
        for kmer in self.in_degree:
            if self.in_degree[kmer] == 0:
                count += 1
        print("0 Degree:{}".format(count))
        self.graph_simplify()
        # node_len_lst = []
        # for node in self.graph:
        #     node_len_lst.append(len(node))
        # plt.hist(node_len_lst,bins=100)
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
        for kmer in self.count_dict:
            if kmer in done:
                continue
            size = 0
            node_num += 1
            new_node = kmer
            done.add(kmer)
            cur_kmer = kmer
            while len(old_graph[kmer][1]) == 1:
                nxtkmer = old_graph[kmer][1][0]
                if len(old_graph[nxtkmer][0]) == 1:
                    done.add(nxtkmer)
                    new_node += nxtkmer[-1]
                    kmer = nxtkmer
                    size += 1
                else:
                    break
            kmer = cur_kmer
            while len(old_graph[kmer][0]) == 1:
                befkmer = old_graph[kmer][0][0]
                if len(old_graph[befkmer][1]) == 1:
                    done.add(befkmer)
                    new_node = befkmer[0] + new_node
                    kmer = befkmer
                    size += 1
                else:
                    break
            if new_node not in nodes:
                nodes.append(new_node)
                new_graph[new_node] = [[],[]]
                new_dict[new_node] = 0
                new_out[new_node] = 0
                new_in[new_node] = 0
            new_dict[new_node] += size

        for i in range(node_num):
            for j in range(node_num):
                if i == j:
                    continue
                node1 = nodes[i]
                node2 = nodes[j]
                if node1[-k+1:] == node2[:k-1]:
                    new_graph[node1][1].append(node2)
                    new_out[node1] += 1
                    new_graph[node2][0].append(node1)
                    new_in[node2] += 1
        self.graph = new_graph
        self.in_degree = new_in
        self.out_degree = new_out
        self.count_dict = new_dict
        print("Compressed Graph size:{}".format(len(self.graph)))

    def fit(self,data_dir='../data/data1',file_type='short',mode="simple"):
        self.load_data(data_dir,file_type)
        return self.get_answers(mode)

    def get_answers(self,mode):
        discription = { 'in_degree':self.in_degree,'count_dict':self.count_dict,
                        'out_degree':self.out_degree,'k':self.k,'pair_ends':self.pair_ends }
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
