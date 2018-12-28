#!/bin/env python
import sys
import pysam
import datetime
from collections import Counter

samfile = pysam.AlignmentFile(sys.argv[1])
fastafile = pysam.FastaFile(sys.argv[2])

def vcf_header():
    i = datetime.datetime.now()
    line = [ "##fileformat=VCFv4.2",
             "##fileDate=%s%s%s" % (i.year,i.month,i.day),
             "##source=Quiver",
             "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
             "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
             "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">",
             "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">",
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"]
    return line

def get_snps(samfile,fastafile):
    het_snps = []
    for pileupcolumn in samfile.pileup():
        bases = []
        for pileupread in pileupcolumn.pileups:
            if pileupread.is_del:
                continue 
            if pileupread.is_refskip:
                continue
            bases.append(pileupread.alignment.query_sequence[pileupread.query_position])
        if len(bases) == 0:
            continue
        count = Counter(bases)
        if len(count.items()) != 2:
            continue    
        ref_name = samfile.get_reference_name(pileupcolumn.reference_id)
        ref_base= fastafile.fetch(ref_name,pileupcolumn.pos,pileupcolumn.pos+1)
        alt_base = None
        for element,count in count.items():
            if element == ref_base:
                continue
            alt_base = element
        out = [ref_name,pileupcolumn.pos + 1,".",ref_base,alt_base]
        het_snps.append(out)
    return het_snps

header = vcf_header()
print "\n".join(header)
het_snps = get_snps(samfile,fastafile)
for het_snp in het_snps:
    out = het_snp
    extra = [60,"PASS","SVTYPE=SNP","GT:GQ:DP","0/1:60:100"]
    out += extra
    print "\t".join(map(str,out))
    


