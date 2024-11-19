import re 
import pysam 
import pathlib
import pandas as pd  

def parse_cigar(cigar, reg=re.compile("([\d]+)([MIDNSHP=X])")):
    cigar_info = reg.findall(cigar)
    clip_left = int(cigar_info[0][0]) if cigar_info[0][1] in ("S", "H") else 0
    clip_right = int(cigar_info[-1][0]) if cigar_info[-1][1] in ("S", "H") else 0
    query_len, ref_len = 0, 0
    for (n, op) in cigar_info:
        n = int(n)
        if op in ("M", "=", "X"):
            query_len += n
            ref_len += n
        elif op == "I":
            query_len += n
        elif op in ("D", "N"):
            ref_len += n
    return clip_left, clip_right, ref_len, query_len

base_map_rel = str.maketrans({'A':'T', 'T':'A', 'G': 'C', 'C':'G'})
#490-4027
for sam_file in pathlib.Path("../02.mapping").glob("*.sort.bam"):
    name = sam_file.name.split('.', 1)[0]
    with pysam.AlignmentFile(sam_file, "rb") as f, open(f"{name}.breakpoint_reads_info.tsv", 'w') as f1:
        f1.write('\t'.join(("read_id", "seq", "reference_name", "map_info", "residual_base", "total_reads")) + '\n')
        store_stat_info = []
        break_point_reads = [0, 0]
        for i, break_point in enumerate((490, 4027)):
            for read in f.fetch("CmAPPR2-hAT-3538", break_point - 1, break_point + 1):
                if not read.flag & 0x904:
                    break_point_reads[i] += 1
                    if read.is_reverse:
                        strand = "-"
                        start, end = read.query_length - read.qend + 1, read.query_length - read.qstart
                    else:
                        strand = "+"
                        start, end = read.qstart + 1, read.qend
                    read_map_info = [(start, end, strand, read.reference_start + 1, read.reference_end)]
                    if read.has_tag("SA"):
                        sa = read.get_tag("SA").strip(';')
                        for info in sa.split(';'):
                            sa_name, sa_pos, sa_strand, sa_cigar, sa_mapq, sa_nm = info.split(',')
                            clip_left, clip_right, ref_len, qry_len = parse_cigar(sa_cigar)
                            sa_pos = int(sa_pos)
                            sa_end = sa_pos + ref_len - 1
                            if (break_point == 490 and 4027 >= sa_pos and sa_end >= 4027) or (break_point == 4027 and 490 >= sa_pos and sa_end >= 490):
                                if sa_strand == '-':
                                    clip_left, clip_right = clip_right, clip_left
                                read_map_info.append((clip_left + 1, clip_left + qry_len, sa_strand, sa_pos, sa_pos + ref_len - 1))
                        if len(read_map_info) != 2:
                            continue
                        read_map_info.sort(key=lambda x: x[:2])
                        a_info, b_info = read_map_info
                        print(read.query_name, read.mapq, read.reference_name, read.seq[a_info[1]:b_info[0] - 1] if strand == "+" else read.seq[::-1][a_info[1]:b_info[0] - 1].translate(base_map_rel), read_map_info)
                        store_stat_info.append('\t'.join((read.query_name, read.seq if strand == "+" else read.seq[::-1].translate(base_map_rel),
                                            read.reference_name, str(read_map_info), read.seq[a_info[1]:b_info[0] - 1] if strand == "+" else read.seq[::-1][a_info[1]:b_info[0] - 1].translate(base_map_rel))))
        total_reads = sum(break_point_reads)
        f1.writelines(stat_info + f'\t{total_reads}\n' for stat_info in store_stat_info )
    data = pd.read_table(f"{name}.breakpoint_reads_info.tsv", header=0, sep='\t')
    data.to_excel(f"{name}.breakpoint_reads_info.xlsx", index=False, header=True)
    pathlib.Path(f"{name}.breakpoint_reads_info.tsv").unlink()
