#!/usr/bin/python
import string

def design_cassettles_Fprimer(output_file,input,cassettles):
    outfile = open(output_file,'w')
    for n in cassettles:
        seq = input[n*24:n*24+60]
        id = '>Cassettles'+str(n+1)
        for aa in range(0,8,1):
            outfile.write(id+ '_'+ str(aa+1)+ "\n"+ seq[0:21+aa*3] + 'NNK' + seq[aa*3+21+3::]+ "\n")
    outfile.close()


def design_cassettles_Rprimer(output_file,input, cassettles):
    outfile = open(output_file,'w')
    for n in cassettles:
        seq = input[n*24:n*24+21]
        complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
        rcseq = seq.translate(complements)[::-1]
        id = '>Cassettles'+str(n+1)
        outfile.write(id+ '_Rprimer'+  "\n"+ rcseq + "\n")
    outfile.close()


def main():
    input = 'GAGGAGAATCCCGGGCCCATGGATCCAAACACTGTGTCAAGCTTTCAGGACATACTGATGAGGATGTCAAAAATGCAGTTGGGGTCCTCATCGGAGGACTTGAATGGAATAATAACACAGTTCGAGTCTCTGAAACTCTACAGAGATTCGCTTGGAGAAGCAGTAATGAGAATGGGAGACCTCCACTCACTCCAAAACAGAAACGGAAAATGGCGGGAACAATTAGGTCAGAAGTTTGAAGAAATAAGGTGGTTGATTGAAGAAGTGAGACACAGACTGAAGATAACAGAGAATAGTTTTGAGCAAATAACATTTATGCAAGCCTTACAACTATTGCTTGAAGTGGAGCAAGAGATAAGAACTTTCTCGTTTCAGCTTATTTAATAATAAAAAACA'
    cassettles = range(0,15,1)
    output_file = 'result/NEP_lib_Fprimer.fas'
    design_cassettles_Rprimer('result/NEP_lib_Rprimer.fas',input,cassettles)
    design_cassettles_Fprimer(output_file,input,cassettles)

if __name__ == "__main__":
  main()
