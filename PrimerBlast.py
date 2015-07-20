# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 12:54:14 2015

@author: Dave
"""

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
from Bio import SeqIO
from Bio import Seq

#This initializes primer1 and primer2
primer1=[]
primer2=[]

#This portion receives the input to fill primer1 and primer2
primer1.append(input("Enter first primer sequence: "))
primer1.append(input("Enter first primer name: "))
primer2.append(input("Enter second primer sequence: "))
primer2.append(input("Enter second primer name: "))

#This opens up a blast dialog with NCBI to perform a blast search
result_handle = NCBIWWW.qblast("blastn", "nt", primer1[0])

#This opens an .xml file for saving the blast results for primer 1
save_file = open("Blast_results_%s.xml" % primer1[1], "w")
#This writes the results into the file
save_file.write(result_handle.read())
save_file.close()
result_handle.close()
#This opens up the file that we just wrote
result_handle = open("Blast_results_%s.xml" % primer1[1])
#This assigns the data to the variable "result_handle"
blast_record = NCBIXML.read(result_handle)

# E_VALUE_THRESH = 0.04
#
#for alignment in blast_record.alignments:
#    for hsp in alignment.hsps:
#        if hsp.expect < E_VALUE_THRESH:
#            print('****Alignment****')
#            print('sequence:', alignment.title)
#            print('length:', alignment.length)
#            print('e value:', hsp.expect)
#            print(hsp.query[0:75] + '...')
#            print(hsp.match[0:75] + '...')
#            print(hsp.sbjct[0:75] + '...')

#This saves the GenInfo Identifier number from the first blast hit
gi=blast_record.alignments[0].title.split('|')[1]
#This saves the starting base number from the gi
gi_start=blast_record.alignments[0].hsps[0].sbjct_start


#This block of code is the equivalent of above
result_handle = NCBIWWW.qblast("blastn", "nt", primer2[0])
save_file = open("Blast_results_%s.xml" % primer2[1], "w")
save_file.write(result_handle.read())
save_file.close()
result_handle.close()
result_handle2 = open("Blast_results_%s.xml" % primer2[1])
blast_record2 = NCBIXML.read(result_handle2)


gis = []
#This cycles through the blast record and makes a list of all of the gis
for aligns in range(0,len(blast_record2.alignments)):
    gis.append(blast_record2.alignments[aligns].title.split('|')[1])

#This finds the first blast record with the same gi from the first blast search
#It then finds the starting base.
gi_end=blast_record2.alignments[gis.index(gi)].hsps[0].sbjct_start

    
#This portion takes the gi number, the start nucleotide and the end
#nucleotide and retrieves the full sequence between from Entrez
#Insert your e-mail address below.
Entrez.email = "******@***.com"     # Always tell NCBI who you are
handle = Entrez.efetch(db="nucleotide", 
                       id=gi, 
                       rettype="fasta", 
                       strand=1, 
                       seq_start=gi_start, 
                       seq_stop=gi_end)
record = SeqIO.read(handle, "fasta")
handle.close()
if gi_end < gi_start:
    record=record.seq.reverse_complement()
output_handle = open("Sequence_%s_to_%s.fasta" % (primer1[1],primer2[1]),"w")
SeqIO.write(record, output_handle, "fasta")
output_handle.close()
print(record.seq)
    