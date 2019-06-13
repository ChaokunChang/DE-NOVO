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
    
    def contain_pairs(self,contig):
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
    for node, degree in candidate:
        assert(graph[node][1] is not None)
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
            if len(graph[kmer][1]) == 2:
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
                single_pairs = PES.contain_pairs(kmer)

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
                addlist = None
                if len(directions) == 1:
                    if directions[0] not in done:
                        selection = directions[0]
                elif len(directions) > 1:
                    for nxtkmer in directions:
                        if nxtkmer in done:
                            continue
                        for p2,start,rl in check_pos:
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
                            else:
                                addlist = None
                        if selection != "":
                            print(selection)
                            break
                    # assert(selection != "")

                if selection != "":
                    terminate = False
                    if addlist is not None:
                        kmer = addlist[-1]
                    else:
                        kmer = selection
                    done.add(selection)
                    ans += selection[k-1:]
                
            output.append(ans)
        print("Count 1-in-1-out {}".format(count))
    else:
        pass
    return output

def longest_path(graph, root,visited:set,depth):
    if len(graph[root][1])==0:
        return [root]
    best = [root]
    cur = [root]
    # print("cur depth:{}".format(depth))
    for child in graph[root][1]:
        if child not in visited:
            visited.add(child)
            cur = [root] + longest_path(graph,child,visited,depth+1)
            visited.remove(child)
            if len("".join(cur)) > len("".join(best)) :
                best = cur
    if len(graph[root][1])>0 and best==[root]:
        # circle occur caused terminate
        print("!!!!")
    return best

def nodes_combine(nodes,k):
    assert(len(nodes) > 0)
    return nodes[0] + "".join(n[k-1:] for n in nodes[1:])

def choose_path(cur_path,cur_node,graph,k,target_pairs,pair_length):
    count = 0
    for p,loc in target_pairs:
        if len(cur_path)>=loc:
            if p == cur_path[loc-pair_length:loc]:
                return cur_path
            else:
                target_pairs.remove( (p,loc) )
            count += 1
    if len(target_pairs) == 0:
        return ""
    else:
        for node in graph[cur_node][1]:
            next_path = cur_path + node[k-1:]
            ans = choose_path( next_path, node, graph, k, target_pairs, pair_length)
            return ans

def compressed_graph_mining(graph, discription):
    k = discription['k']
    in_degree = discription['in_degree']
    out_degree = discription['out_degree']
    PES = discription['pair_ends']
    count_dict = discription['count_dict']
    results = []
    start_nodes = []
    end_nodes = []
    single_nodes = []
    split_nodes = []
    merge_nodes = []
    hub_nodes = []
    total_length = 0
    for node in graph:
        if len(graph[node][0]) == 0 : # in degree is 0
            start_nodes.append(node)
        if len(graph[node][1]) == 0 : # out degree is 0
            end_nodes.append(node) # both in-out degree is 0
        if len(graph[node][0]) == 0 and len(graph[node][1]) == 0:
            single_nodes.append(node)
        if len(graph[node][1]) > 1:
            assert(len(graph[node][1]) == 2)
            split_nodes.append(node)
            if len(graph[node][0]) > 1:
                assert(len(graph[node][0]) == 2)
                hub_nodes.append(node)
        if len(graph[node][0]) > 1:
            assert(len(graph[node][0]) == 2)
            merge_nodes.append(node)
        total_length += len(node)
    print("All asserts passed.")
    print("Split nodes number:{}".format(len(split_nodes)))
    print("Merge nodes number:{}".format(len(merge_nodes)))
    print("Hub nodes number:{}".format(len(hub_nodes)))
    print("Start nodes number:{}".format(len(start_nodes)))
    print("End nodes number:{}".format(len(end_nodes)))
    print("Sigle nodes number:{}".format(len(single_nodes)))
    print("The Start node and End nodes maybe same(single node case),Let's remove them.")
    for node in single_nodes:
        print(len(node))
        results.append(node)
        total_length -= len(node)
        start_nodes.remove(node)
        end_nodes.remove(node)
    print("Length:{}".format(total_length))
    assert(len(start_nodes) == 4)
    visited = set()
    for node in split_nodes:
        single_pairs = PES.contain_pairs(node)
        target_pairs = []
        for pair,loc in single_pairs:
            target_pairs.append(( PES.get_pair(pair),loc+PES.pair_dis ))
        right_path = choose_path(node,node,graph,k,target_pairs,PES.length)
        if right_path != "":
            print(right_path)
        else:
            print("DONE")
        checked_len = len(node) - (k-1)
        cur_check = node
        nxt_check = node
        while(checked_len < len(right_path)- (k-1)):
            choice = 0
            if len(graph[cur_check][1]) == 0:
                break
            # nxt_check = graph[cur_check][1][0]
            for i in range( len(graph[cur_check][1]) ):
                cur = graph[cur_check][1][i]
                if cur == right_path[checked_len:checked_len+len(cur)]:
                    choice = i
                    nxt_check = cur
                    break
            checked_len += len(graph[cur_check][1][choice]) - (k-1)
            graph[cur_check][1] = [ graph[cur_check][1][choice] ]
            cur_check = nxt_check
        
        # af_pairs = []
        # new_af = None
        # for af in graph[node][1]:
        #     pair = PES.contain_pairs(af)
        #     af_pairs.append(pair)
        #     for p1,loc1 in single_pairs:
        #         p2 = PES.get_pair(p1)
        #         for tp2,loc2 in pair:
        #             if tp2 == p2:
        #                 gap_dis = (len(node)-loc1 + loc2+PES.length+1 - (k-1) )
        #                 if gap_dis == PES.pair_dis:
        #                     new_af = [af]
        #                     break
        #         if new_af is not None:
        #             break
        #     if new_af is not None:
        #         break
        # if new_af is not None:
        #     print(new_af)
        #     graph[node][1] = new_af

        # bef_pairs = []
        # new_bef = None
        # for bef in graph[node][0]:
        #     pair = PES.contain_pairs(bef)
        #     bef_pairs.append(pair)
        #     for p1,loc1 in single_pairs:
        #         p2 = PES.get_pair(p1)
        #         for tp2,loc2 in pair:
        #             if tp2 == p2:
        #                 gap_dis = (loc1+PES.length+1 + (len(bef)-loc2) - (k-1) )
        #                 if gap_dis == PES.pair_dis:
        #                     new_bef = [bef]
        #                     break
        #         if new_bef is not None:
        #             break
        # if new_bef is not None:
        #     break
        # graph[node][0] = new_bef

    for node in start_nodes[2:]:
        visited.add(node)
        path = longest_path(graph,node,visited,0)
        for p in path:
            visited.add(p)
        ans = nodes_combine(path,k)
        print(ans)
        results.append(ans)

    # for node in start_nodes:
    #     visited = set()
    #     ans = node
    #     selection = node
    #     candidates = graph[selection][1]
    #     while(len(candidates) > 0):
    #         selection = ""
    #         for contig in candidates:
    #             if contig not in visited:
    #                 if len(contig) > len(contig):
    #                     selection = contig
    #         if selection != "":
    #             ans += selection[k-1:]
    #             visited.add(selection)
    #             candidates = graph[selection][1]
    #         else:
    #             results.append(ans)
    #             break
    # for node in end_nodes:
    #     visited = set()
    #     ans = node
    #     selection = node
    #     candidates = graph[selection][0]
    #     while(len(candidates) > 0):
    #         selection = ""
    #         for contig in candidates:
    #             if contig not in visited:
    #                 if len(contig) > len(contig):
    #                     selection = contig
    #         if selection != "":
    #             ans = selection[:-k+1] + ans
    #             visited.add(selection)
    #             candidates = graph[selection][0]
    #         else:
    #             results.append(ans)
    #             break
    
    return results


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
    def __init__(self,k=29,step=1,limit=1):
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
            print("Get data from {}".format(fname))
            file_path = opj(data_dir,fname)
            reads = SeqIO.parse(file_path,'fasta')
            rlist.append(reads)
        print("Load the sequences in {} done".format(data_dir))
        print("Building counting dict for them...")
        self.build_graph(rlist)
        self.graph_simplify()
        # node_len_lst = []
        # for node in self.graph:
        #     node_len_lst.append(len(node))
        # plt.hist(node_len_lst,bins=100)
        # plt.show()

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
        print("Before Compress, 0 Degree count:{}".format(count))
    
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
        # output = graph_mining(self.graph,discription,mode)
        output = compressed_graph_mining(self.graph,discription)
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
