#!/usr/bin/python

#clip a number (x) of bases from ends of bam reads based on fasta sequence of primers.

import pysam
from sys import argv
from cigar import Cigar
import re
import regex #regex with fuzzy matching
from read2_dict import R2_dict

def RC(seq):
    # return reverse complement of DNA string
    # source: http://crazyhottommy.blogspot.ca/2013/10/python-code-for-getting-reverse.html
    seq1 = 'ATCGTAGCatcgtagc'
    seq_dict = { seq1[i]:seq1[i+4] for i in range(16) if i < 4 or 8<=i<12 }
    return "".join([seq_dict[base] for base in reversed(seq)])

ibam_file = argv[1]
#obam_file = argv[2]
n_errors = argv[3]
obam_file = ibam_file[:-4] + '_clipped_d'+str(n_errors)+'.bam'
bed_file = argv[2]
reads = {}
ibam = pysam.AlignmentFile(ibam_file, 'rb')
obam = pysam.AlignmentFile(obam_file, 'wb', template=ibam)
bf=open(bed_file)
r_primers={}
f_primers={}
bf_line_number=0
#read counters
read1s=0
total_read1=0
read2s=0
total_read2=0
total_reads=0

#reads = {read.query_name + '_' + str(read.flag): read for read in ibam if read.is_proper_pair}
read_qn_pairs = []
read_qn = [r.query_name for r in ibam]
read_qn_pairs = [qn for qn in read_qn if read_qn.count(qn)==2]
read_qn_pairs=list(set(read_qn_pairs)) #get unique query names
#print 'done getting query names'
#get query names for pairs only
# for qn in read_qn:
#     if read_qn.count(qn)==2:
#         read_qn_pairs.append(qn)
#print len(read_qn)
#print len(read_qn_pairs)
#exit(0)
#exit(0)
#print len(reads)

ibam_indexed = pysam.IndexedReads(ibam)
ibam_indexed.build()
# try:

#create dictionary of primers
for line in bf.readlines():
    bf_line_number += 1 #   in case primer sequences are duplicated
    line_parts = line.strip().split('\t')
    chromosome = line_parts[0]
    start = int(line_parts[1])
    end = int(line_parts[2])
    f_primer = str(line_parts[3])
    f_primer_l = len(f_primer)
    r_primer = str(line_parts[4])
    r_primer_l = len(r_primer)
    f_primer_rgx = regex.compile('(' + f_primer + '){d<='+n_errors+'}', regex.IGNORECASE)  # (f_primer_rgx = regex.compile(f_primer, regex.IGNORECASE)  allow up to 5 errors in match '(' + f_primer + '){e<=5}', regex.IGNORECASE
    r_primer_rgx = regex.compile('(' + r_primer + '){d<='+n_errors+'}', regex.IGNORECASE)  #r_primer_rgx = regex.compile(r_primer, regex.IGNORECASE)
    f_primers[f_primer+'_'+str(chromosome)+'_'+str(start)+'_'+str(end)] = f_primer_rgx
    r_primers[r_primer+'_'+str(chromosome)+'_'+str(start)+'_'+str(end)] = r_primer_rgx
print len(f_primers)
print len(r_primers)

for qn in read_qn_pairs:
    t=ibam_indexed.find(qn)
    for r in t:
        if r.is_proper_pair:
            if r.is_read1:
                total_read1 += 1
                r1 = r
                #print 'read 1: ' + str(r1.query_alignment_sequence)
                for fp in f_primers:
                    len_fp = len(fp.split('_')[0])
                    csome = fp.split('_')[1]
                    start_pos = int(fp.split('_')[2])
                    end_pos = int(fp.split('_')[3])
                    if f_primers[fp].match(r1.query_alignment_sequence,0,len_fp-1) and 'chr' + str(r1.rname) == csome:  # and r1.pos>=s and r1.pos<=e and r1.rname and r1.pos >= start_pos
                        read1s += 1
                        len_match_r1 = f_primers[fp].match(r1.query_alignment_sequence,0,len_fp-1)  # length of matching sequence
                        len_match_r1 = len(len_match_r1.group(1)) #only one match is expected
                        r1_cigar = Cigar(r1.cigarstring)
                        r1_cigar = r1_cigar.mask_left(len_match_r1)
                        r1.cigarstring = str(r1_cigar)
                        r1.pos += len_match_r1  # offset start position by length of forward primer
                        obam.write(r1)
                        print str(r1.pos)
                        print 'no. clipped read 1: ' + str(read1s)
                        print 'no. read 1: ' + str(total_read1)
                        print 'no. reads: ' + str(total_reads)
                        print 'length of match: ' + str(len_match_r1)
                        print 'length of primer: ' + str(len_fp)
                        break
            if r.is_read2:
                total_read2 += 1
                r2 = r
                #print 'read 2: ' + str(r2.query_alignment_sequence)
                for rp in r_primers:
                    len_rp = len(rp.split('_')[0])
                    csome = rp.split('_')[1]
                    s = int(rp.split('_')[2])
                    e = int(rp.split('_')[3])
                    if r_primers[rp].match(RC(r2.query_alignment_sequence),0,len_rp-1) and 'chr' + str(r2.rname) == csome:  # missing reverse complement # and r2.pos>=s and r2.pos<=e and r2.pos>=s
                        read2s += 1
                        len_match_r2 = r_primers[rp].match(RC(r2.query_alignment_sequence),0,len_rp-1)
                        len_match_r2 = len(len_match_r2.group(1))
                        r2_cigar = Cigar(r2.cigarstring)
                        r2_cigar = r2_cigar.mask_right(len_match_r2) #test clip len of primer
                        r2.cigarstring = str(r2_cigar)
                        #r2.pos += len_match_r2
                        obam.write(r2)
                        print str(r2.pos)
                        print 'no. clipped read 2: ' + str(read2s)
                        print 'no. read 2: ' + str(total_read2)
                        print 'no. reads: ' + str(total_reads)
                        print 'length of match: ' + str(len_match_r2)
                        print 'length of primer: ' + str(len_rp)
                        break
        else:
            #read does not align as a pair or is an orphan read (singleton)
            continue

