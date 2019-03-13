#!/bin/env
import sys
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def fix_reads_for_blasr(file_to_fix,fixed_file,type_):
    records = []
    for record in SeqIO.parse(file_to_fix,type_):
        record.id = "%s/0/0_8" % record.id
        record.description = ""
        records.append(record)
    SeqIO.write(records,fixed_file,type_)


file_to_fix = sys.argv[1]
fixed_file = sys.argv[2]
type_ = sys.argv[1].split(".")[-1]

fix_reads_for_blasr(file_to_fix,fixed_file,type_)
