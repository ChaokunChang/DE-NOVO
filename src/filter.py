from utils import *
from Bio import Seq, SeqIO, SeqRecord

file_path = './data4_final.fasta'
reads = SeqIO.parse(file_path,'fasta')

answers = []

for read in reads:
    seq = str(read.seq)
    answers.append((seq,len(seq)) )

results = sorted(answers,key=lambda x: -x[1])

with open('./data4_ans.fasta','w') as f:
    for i,re in enumerate(results):
        if i<5 or i % 2 == 0:
            f.write('>short_read_{}/0 \n'.format(i))
            f.write(re[0][:])
            print(re[1])
            f.write('\n')

