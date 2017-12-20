# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 12:30:25 2016

@author: mada_
"""


#gravar a sequencia em formato genbank
from Bio import Entrez, SeqIO
Entrez.email = "madalenaa.castro@gmail.com" 
handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id="52840256", seq_start="998651", seq_stop="1275600")
seq_record = SeqIO.read(handle, "genbank")
SeqIO.write(seq_record, 'sequenceGenbank.gb', "genbank") #Guarda em formato genbank
handle.close()

#gravar a sequencia em formato txt
Dfile= open('sequenceTxt.txt', 'w')
Dfile.write(str(seq_record.seq)) 
Dfile.close()

#gravar a sequencia em formato fasta
Entrez.email = "madalenaa.castro@gmail.com" 
handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id="52840256", seq_start="998651", seq_stop="1275600")
seq_record = SeqIO.read(handle, "fasta")
SeqIO.write(seq_record, 'sequenceFasta.fasta', "fasta") #Guarda em formato genbank
handle.close()
