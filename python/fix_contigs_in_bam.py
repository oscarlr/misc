#!/bin/env python
import sys
import pysam

samfile = pysam.AlignmentFile(sys.argv[1])
instructionsfn = sys.argv[2]

def get_query_pos(ref_name,ref_pos_to_look_for,bamfile,contig_name):
    query_pos_to_get = None
    samfile = pysam.AlignmentFile(bamfile)
    for read in samfile:
        if read.query_name != contig_name:
            continue
        for query_pos, ref_pos in read.get_aligned_pairs():
            #print contig_name,query_pos,ref_pos,ref_pos_to_look_for
            if ref_pos == ref_pos_to_look_for:
                query_pos_to_get = query_pos
    assert query_pos_to_get != None
    return query_pos_to_get

def get_sequence(ref_name,ref_start,ref_end,bamfile,contig_name):
    samfile = pysam.AlignmentFile(bamfile)
    sequence = None
    for read in samfile:
        if read.query_name != contig_name:
            continue
        sequence = read.query_sequence
    assert sequence != None
    if ref_start != -1:
        ref_start_query_pos = get_query_pos(ref_name,ref_start,bamfile,contig_name)
    else:
        ref_start_query_pos = 0
    if ref_end != -1:
        ref_end_query_pos = get_query_pos(ref_name,ref_end,bamfile,contig_name)
    else:
        ref_end_query_pos = -1
    return sequence[ref_start_query_pos:ref_end_query_pos]

contigs_to_overwrite = set()
contigs_to_keep = set()

with open(instructionsfn,'r') as instructionsfh:
    header = instructionsfh.readline()
    for line in instructionsfh:
        line = line.rstrip().split("\t")
        contigs_to_overwrite.add(line[1])

for contig in samfile:
    if contig.query_name in contigs_to_overwrite:
        continue
    contigs_to_keep.add((contig.query_name,contig.query_sequence))

out_sequences = {}
with open(instructionsfn,'r') as instructionsfh:
    header = instructionsfh.readline()
    for line in instructionsfh:
        line = line.rstrip().split("\t")
        out_sequence_name = line[0]
        contig_name = line[1]        
        ref_name = line[2]
        ref_start = int(line[3])
        ref_end = int(line[4])
        bamfile = line[5]
        if out_sequence_name not in out_sequences:
            out_sequences[out_sequence_name] = []
        sequence = get_sequence(ref_name,ref_start,ref_end,bamfile,contig_name)
        out_sequences[out_sequence_name].append(sequence)

for contig_name, contig_sequence in contigs_to_keep:
    print ">%s\n%s" % (contig_name,contig_sequence)

for contig_name in out_sequences:
    print ">%s" % contig_name
    print "".join(out_sequences[contig_name])

