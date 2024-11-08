import re
import pysam 
from collections import defaultdict
from operator import itemgetter
import pandas as pd 

reg = re.compile("([\d]+)([MIDNSHP=X])")
with pysam.FastaFile("ref.fa") as fa, pysam.AlignmentFile("V91.sort.bam", "rb") as f, open("read_map_info.tsv", 'w') as f1:
    store_read_info = defaultdict(set)
    for read in f:
        if read.mapq >= 60 and not read.flag & 0xb00: #过滤掉比对质量值大于等于60且比对为主要比对的reads比对信息
            left_clip, right_clip, ref_len, qry_len = 0, 0, 0, 0
            store_ins_info, store_del_info, store_mut_info = [], [], []
            ref_consensus_seq, qry_consensus_seq = '', ''
            for i, info in enumerate(reg.finditer(read.cigarstring)):
                n, op = info.groups()
                n = int(n)
                if i == 0 and op in ('S', 'P'):
                    qry_consensus_seq += read.seq[left_clip:left_clip + n].lower()
                    ref_consensus_seq += fa.fetch(read.reference_name, read.pos - n, read.pos)
                    left_clip = n
                elif op in ('M', '=', 'X'):
                    #寻找SNP变异
                    ref_seq = fa.fetch(read.reference_name, read.pos + ref_len, read.pos + ref_len + n)
                    qry_seq = read.seq[left_clip + qry_len:left_clip + qry_len + n]
                    for i, (ref_base, qry_base) in enumerate(zip(ref_seq, qry_seq), 1):
                        if ref_base.upper() != qry_base.upper():
                            store_mut_info.append((read.pos + ref_len + i, ref_base, qry_base))
                    qry_consensus_seq += qry_seq
                    ref_consensus_seq += ref_seq
                    ref_len += n
                    qry_len += n
                elif op == 'I':
                    qry_consensus_seq += read.seq[left_clip + qry_len:left_clip + qry_len + n]
                    ref_consensus_seq += '-' * n 
                    qry_len += n
                    store_ins_info.append((read.pos + ref_len, n))
                elif op in ('D', 'N'):
                    ref_consensus_seq += fa.fetch(read.reference_name, read.pos + ref_len, read.pos + ref_len + n)
                    qry_consensus_seq += '-' * n
                    if op == 'D':
                        store_del_info.append((read.pos + ref_len, n))
                    ref_len += n
            if op in ('S', 'P'):
                qry_consensus_seq += read.seq[-n:].lower()
                ref_consensus_seq += fa.fetch(read.reference_name, read.pos + ref_len, read.pos + ref_len + n)
                right_clip = n
            #
            ref_start = read.pos - left_clip + 1
            ref_end = ref_start + ref_len + right_clip - 1
            store_read_info[(read.reference_name, ref_start, ref_end, tuple(store_ins_info), tuple(store_del_info), tuple(store_mut_info), ref_consensus_seq, qry_consensus_seq)].add(read.query_name)
    #
    f1.write('\t'.join(('Aligned_Sequence', 'Reference_Sequence', 'Reference_Name', 'Read_Status', 'n_deleted', 'n_inserted', 'n_mutated', '#Reads', '%Reads')) + '\n')
    #大于1条reads支持该结果
    total_n = sum(len(v) for v in store_read_info.values())
    store_read_info = sorted(((k, len(v)) for k, v in store_read_info.items() if len(v) > 1), key=itemgetter(1), reverse=True)
    for k, n in store_read_info:
        ref_name, ref_start, ref_end, store_ins_info, store_del_info, store_mut_info, ref_seq, seq = k
        status = 'MODIFIED' if store_ins_info or store_del_info or store_mut_info else 'UNMODIFIED'
        f1.write('\t'.join((seq, ref_seq, f'{ref_name}:{ref_start}-{ref_end}', status, 
                                *map(lambda x: str(len(x)), (store_del_info, store_ins_info, store_mut_info)), str(n), f"{n * 100 / total_n:.5f}")) + '\n')
    data = pd.read_table("read_map_info.tsv", sep='\t', header=0)
    data.to_excel("V91.read_map_info.xlsx", header=True, index=False)