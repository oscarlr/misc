#!/bin/env python
import sys
import pysam

svs_bed = sys.argv[1]
hap1_reads_to_svs_bam = sys.argv[2]
hap2_reads_to_svs_bam = sys.argv[3]
hap1_reads_to_ref_bam = sys.argv[4]
hap2_reads_to_ref_bam = sys.argv[5]
padding = int(sys.argv[6])

def load_svs(bed):
    svs = []
    with open(bed,'r') as bedfh:
        for line in bedfh:
            line = line.rstrip().split('\t')
            svs.append([(line[0],line[1],line[2]),line])
    return svs
        
def score(read,start,end):
    score = None
    padding = 50
    start -= padding
    end += padding
    if read.reference_start > start:
        return score
    if read.reference_end < end:
        return score
    seq_start = None
    seq_end = None
    length = 0.0
    errors = 0.0
    start_scoring = False
    for query,ref in read.get_aligned_pairs():
        if start_scoring == False:
            if ref == None:
                continue
            if query == None:
                continue
            if ref >= start:
                start_scoring = True
            continue
        if ref == None or query == None:
            errors += 1.0
        length += 1.0
        if ref >= end:
            break
    if length > 0:
        score = (length - errors)/(length)
    return score

def score_svs(bam,padding,ref=False):
    svs = {}
    samfile = pysam.AlignmentFile(bam)
    for read in samfile.fetch():
        if read.is_secondary:
            continue
        if read.is_supplementary:
            continue
        if read.is_unmapped:
            continue
        ref_name = samfile.get_reference_name(read.reference_id)
        sv_chrom,coord = ref_name.split(":")
        sv_start,sv_end = coord.split("-")
        sv_coord = (sv_chrom,sv_start,sv_end)
        if sv_coord not in svs:
            svs[sv_coord] = {}
        end_of_padding = padding
        if ref:
            end_of_padding = padding + (int(sv_end) - int(sv_start))
        alignment_score = score(read,padding,end_of_padding)
        if alignment_score != None:
            svs[sv_coord][read.query_name] = float(alignment_score)
    return svs

svs = load_svs(svs_bed)
#print "hap1_svs"
hap1_svs = score_svs(hap1_reads_to_svs_bam,padding)
#print "hap2_svs"
hap2_svs = score_svs(hap2_reads_to_svs_bam,padding)
#print "hap1_ref"
hap1_ref = score_svs(hap1_reads_to_ref_bam,padding,True)
#print "hap2_ref"
hap2_ref = score_svs(hap2_reads_to_ref_bam,padding,True)

def print_winner(reads,hap_svs,ref,sv):    
    number_ref_reads = 0
    number_hap_reads = 0
    total_reads = set()
    ref_winner = 0
    hap_winner = 0
    for read in reads:
        score = 0
        ref_score = 0
        if sv in hap_svs:
            if read in hap_svs[sv]:
                score = hap_svs[sv][read]
                number_hap_reads += 1
                total_reads.add(read)
        if sv in ref:
            if read in ref[sv]:
                ref_score = ref[sv][read]
                number_ref_reads += 1
                total_reads.add(read)
        if score > ref_score:
            hap_winner += 1
        if ref_score > score:
            ref_winner += 1
    return (ref_winner,hap_winner,number_ref_reads,number_hap_reads,len(total_reads))

out = ["chrom","start","end","h1_ref","h1","h1_ref_reads","h1_reads","h2_ref","h2","h2_ref_reads","h2_reads","h","r","h1_total_reads","h2_total_reads"]
print "\t".join(map(str,out))
for sv,line in svs:
    # hap 1
    hap1_reads = set()
    if sv in hap1_svs:
        for read in hap1_svs[sv]:
            hap1_reads.add(read)
    if sv in hap1_ref:
        for read in hap1_ref[sv]:
            hap1_reads.add(read)
    # hap 2
    hap2_reads = set()
    if sv in hap2_svs:
        for read in hap2_svs[sv]:
            hap2_reads.add(read)
    if sv in hap2_ref:
        for read in hap2_ref[sv]:
            hap2_reads.add(read)
    # (ref_winner,hap_winner,number_ref_reads,number_hap_reads)
    h1_ref,h1,h1_ref_reads,h1_reads,h1_total_reads = print_winner(hap1_reads,hap1_svs,hap1_ref,sv)
    h2_ref,h2,h2_ref_reads,h2_reads,h2_total_reads = print_winner(hap2_reads,hap2_svs,hap2_ref,sv)
    h = h1 + h2
    r = h1_ref + h2_ref
    out = list(sv)
    out += [h1_ref,h1,h1_ref_reads,h1_reads,h2_ref,h2,h2_ref_reads,h2_reads,h,r,h1_total_reads,h2_total_reads]
    out += line
    # if min(h,r) == 0:
    #     if h > r:
    #         out.append("hap")
    #     elif r > h:
    #         out.append("ref")
    #     else:
    #         out.append("None")
    # else:
    #     out.append("None")
    # if min([h1_ref_reads,h1_reads,h2_ref_reads,h2_reads]) > 1:
    #     if h > 1:
    #         out.append("hap")
    #     elif r > 0:
    #         out.append("ref")
    #     else:
    #         out.append("None")
    # else:
    #     out.append("None")
    # out += line
    print "\t".join(map(str,out))
