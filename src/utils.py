import math
import sys
import os
import argparse
import numpy as np
from os.path import join as opj
from Bio import Seq
from collections import defaultdict
bases = "ATCG" # base pair: A-T C-G
def parse_args():
    """
    Parses command line arguments.
    """
    parser = argparse.ArgumentParser('DNA Concantatente')
    model_settings = parser.add_argument_group('model seetting')
    model_settings.add_argument('--model', choices=['DBG'], default='DBG',
                                help='choose the algorithm to use')
    model_settings.add_argument('--mode',default="compress")
    DBG_settings = parser.add_argument_group('DBG model settings')
    DBG_settings.add_argument('--k', type=int, default=29,
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
                               help='the result dir.')
    path_settings.add_argument('--result_name', default="results100.fasta",
                               help='the result file name.')
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

class Node:
	def __init__(self,contig_id,childs=[],fathers=[]):
		self.id = contig_id
		self.childs = childs
		self.fathers = fathers
	def add_child(self,child_id):
		if isinstance(child_id,list):
			self.childs += child_id
		else:
			self.childs.append(child_id)

	def add_father(self,father_id):
		if isinstance(father_id,list):
			self.fathers += father_id
		else:
			self.fathers.append(father_id)
	
	def get_in_num(self):
		return len(self.fathers)
	def get_out_num(self):
		return len(self.childs)

def get_forward_contig(count_dict,kmer,graph):
	ans = []
	nxtkmers = graph[kmer][1]
	while len(nxtkmers) == 1:
		kmer = nxtkmers[0]
		if len(graph[kmer][0]) > 1:
			break
		ans = ans + [kmer]
		nxtkmers = graph[kmer][1]
	return ans

def get_backward_contig(count_dict,kmer,graph):
	ans = []
	befkmers = graph[kmer][0]
	while len(befkmers) == 1:
		kmer = befkmers[0]
		if len(graph[kmer][1]) > 1:
			break
		ans = [kmer] + ans
		befkmers = graph[kmer][0]
	return ans

def get_contig(count_dict,kmer,graph):
	forward_contig = get_forward_contig(count_dict,kmer,graph)
	backward_contig = get_backward_contig(count_dict,kmer,graph)
	contig = backward_contig + [kmer] + forward_contig
	return contig

def graph_compress(k,count_dict,graph):
	nodes = []
	count = 0
	visited = set()
	new_graph = {}
	new_in = defaultdict(int)
	new_out = defaultdict(int)
	for kmer in count_dict:
		if kmer not in visited:
			contig = get_contig(count_dict,kmer,graph)
			for km in contig:
				visited.add(km)
				# visited.add(twin(km))
			contigstr = contig2str(contig)
			nodes.append(contigstr)
			new_graph[contigstr] = Node(count,[],[])
			new_in[contigstr] = 0
			new_out[contigstr] = 0
			count += 1
	for i in range(len(new_graph)):
		for j in range(len(new_graph)):
			if i==j :
				continue
			if nodes[j][:k] in front_seq(nodes[i][-k:]):
				new_graph[nodes[i]].add_child(j)
				new_out[nodes[i]] += 1
				new_graph[nodes[j]].add_father(i)
				new_in[nodes[j]] += 1
	count1 = 0
	for c in new_graph:
		if new_graph[c].get_out_num() == 1:
			count1 += 1
	print("After compress:{},{}".format(count,count1) )
	return nodes,new_graph,new_in,new_out
