#!/usr/bin/env python3

from math import ceil
from collections import defaultdict

with open("sample_sex_info.txt") as f:
    store_male_files, store_female_files = [], []
    for line in f:
        file_name, sex = line.strip().split('\t')
        if sex == "male":
            store_male_files.append(file_name)
        else:
            store_female_files.append(file_name)
            
def encode_seq(seq, store_map_rel=('A', 'T', 'G', 'C')):
    v = 0
    for c in seq:
        v <<= 2 
        v |= store_map_rel.index(c)
    return v

def decode_seq(v, seq_len, store_map_rel=('A', 'T', 'G', 'C')):
    seq = ''
    for i in range(seq_len):
        seq += store_map_rel[(v & 3)]
        v >>= 2
    return seq[::-1]
            
store_male_kmer = defaultdict(int)
for n, file in enumerate(store_male_files):
    with open(file) as f:
        for i, line in enumerate(f):
            if i % 2:
                v = encode_seq(line.strip())
                store_male_kmer[v] += 1
                
male_kmer_threshold = ceil(len(store_male_files) * 0.8)
male_kmer = dict((k, v) for k, v in store_male_kmer.items() if v >= male_kmer_threshold)


store_female_kmer = defaultdict(int)
for n, file in enumerate(store_female_files):
    with open(file) as f:
        for i, line in enumerate(f):
            if i % 2:
                v = encode_seq(line.strip())
                store_female_kmer[v] += 1

female_kmer_threshold = ceil(len(store_female_files) * 0.8)
female_kmer = dict((k, v) for k, v in store_female_kmer.items() if v >= female_kmer_threshold)

with open("male_specific.fa", 'w') as f1:
    for i, k in enumerate((male_kmer.keys() - female_kmer.keys()), 1):
        f1.write(f'>{i}_{male_kmer[k]}\n{decode_seq(k, 19)}\n')
        
        
with open("female_specific.fa", 'w') as f1:
    for i, k in enumerate((female_kmer.keys() - male_kmer.keys()), 1):
        f1.write(f'>{i}_{female_kmer[k]}\n{decode_seq(k, 19)}\n')