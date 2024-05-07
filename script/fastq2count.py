#!/usr/bin/python
import sys
from Bio import SeqIO
from collections import Counter

def translation(seq):
  dnamap = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"_", "TAG":"_",
    "TGT":"C", "TGC":"C", "TGA":"_", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
    "---":"-"}
  pep = []
  i = 0
  while i < len(seq):
    codon = seq[i:i+3]
    aa = dnamap[codon]
    pep.append(aa)
    i = i + 3
  pep = ''.join(pep)
  return pep

def sum_mut(aa1,aa2):
    return sum ( aa1[i] != aa2[i] for i in range(len(aa1)) )


def call_mutid(mutpep,refseq,shift):
  mut_id_ls = []
  assert(len(mutpep)==len(refseq))
  for n in range(len(mutpep)):
    pos = n+shift
    if refseq[n]!=mutpep[n]:
       mut_id_ls.append(refseq[n]+str(pos)+mutpep[n])
  return mut_id_ls

def cal_fastq_dic(fastq, ref_dna):
  print ("reading %s" % fastq)
  Rrecords = SeqIO.parse(fastq,"fastq")
  mut_id_ls = []
  len_error = 0
  shift = 1
  for record in Rrecords:
      ID = str(record.id)
      seq = str(record.seq)
      if len(seq) != 363:
          len_error += 1
          continue
      mut_aa = translation(seq)
      ref_aa = translation(ref_dna)
      if sum_mut(seq, ref_dna) == 0:
        mut_id = 'WT'
        mut_id_ls.append(mut_id)
      else:
          if sum_mut(mut_aa, ref_aa) == 0:
            call_mut_id = call_mutid(seq, ref_dna, shift)
            mut_ids = "-".join(sorted(call_mut_id, key=lambda x: int(x[1:-1])))
            mut_id = 'silent' + '_' + mut_ids
            mut_id_ls.append(mut_id)
          else:
            call_mut_id = call_mutid(mut_aa, ref_aa,shift)
            mut_id = "-".join(sorted(call_mut_id, key=lambda x:int(x[1:-1])))
            mut_id_ls.append(mut_id)
  AA_dict = Counter(mut_id_ls)
  return AA_dict

def write_mut_table(mut_dic,outfilename):
  outfile = open(outfilename,'w')
  outfile.write("\t".join(['NEP_Mutation', 'count'])+"\n")
  for mut in mut_dic.keys():
    NEP_Mutation = mut
    count = mut_dic[mut]
    outfile.write("\t".join(map(str,[NEP_Mutation, count]))+"\n")
  outfile.close()
  print('Written %s' %outfilename)


def main():
    if len(sys.argv) != 3:
        sys.exit('[usage] python %s <fasq file> <mutation count output filename>')
    fastqfile = sys.argv[1]
    outfilename = sys.argv[2]
    ref_dna = 'ATGGATCCAAACACTGTGTCAAGCTTTCAGGACATACTGATGAGGATGTCAAAAATGCAGTTGGGGTCCTCATCGGAGGACTTGAATGGAATAATAACACAGTTCGAGTCTCTGAAACTCTACAGAGATTCGCTTGGAGAAGCAGTAATGAGAATGGGAGACCTCCACTCACTCCAAAACAGAAACGGAAAATGGCGGGAACAATTAGGTCAGAAGTTTGAAGAAATAAGGTGGTTGATTGAAGAAGTGAGACACAGACTGAAGATAACAGAGAATAGTTTTGAGCAAATAACATTTATGCAAGCCTTACAACTATTGCTTGAAGTGGAGCAAGAGATAAGAACTTTCTCGTTTCAGCTTATT'
    mutation_dic = cal_fastq_dic(fastqfile,ref_dna)
    write_mut_table(mutation_dic, outfilename)


if __name__ == "__main__":
  main()

