from Bio import Entrez
from Bio import SeqIO
import os
import re

# aceder ao NCBI e guardar o ficheiro correspondente 
# ao genoma do organismo filtrando uma lista de genes 
# de interesse para a sua via 
#

def fetch_sequence():
    genesStarteStop = [(1579611,1580459),(3264100,3264669),(3317965,3319431),(88210,89373),(1893413,1894105),(87011,88213),(3285417,3286169)]
    #print (genesStarteStop[1])
    Entrez.email="a71150@alunos.uminho.pt"
    i=0
    # Ficheiro do genoma
    handle=Entrez.efetch(db="nucleotide",rettype="gbwithparts",retmode="text",id="NC_002942.5")
    sequence_filename="genoma.gb"
    seq_record=SeqIO.read(handle,"gb")
    SeqIO.write(seq_record,sequence_filename,"genbank")
    features = []
    #Ficheiros dos genes
    for tuplo in genesStarteStop:
        #Criar o ficheiro GenBank de cada Gene
        handle=Entrez.efetch(db="nucleotide",rettype="gbwithparts",retmode="text",id="NC_002942.5",seq_start=tuplo[0],seq_stop=tuplo[1])
        sequence_filename="gene "
        indice = str(i)
        name = sequence_filename + indice + ".gb"
        seq_record=SeqIO.read(handle,"gb")
        SeqIO.write(seq_record,name,"genbank")
        i+=1
        handle.close()

#Identificação do gene (GeneID NCBI,Accession number NCBI, 
# locus tag, nome do gene, strand)
def geneIdentification():
    dicionarios[] =  
    for i in range(7):
        index = str(i)
        name="gene "+index+".gb" 
        #Ler o ficheiro
        record=SeqIO.read(name,"genbank")
        for feat in record.features:
            #print (feat.qualifiers)
            try:
                #GeneID_NCBI
                geneID=feat.qualifiers["db_xref"][0].split(":")[1]
            except Exception:
                geneID = "-"
            try:
                locustag = feat.qualifiers["locus_tag"][0]
            except Exception:
                locustag = "-"
            try:
                genename = feat.qualifiers["gene"][0]
            except Exception:
                genename = "-"    
            try:
                strand = str(feat.location)
                if "+" in strand:
                    strand="(+)"
                else:
                    strand="(-)"
            except Exception:
                    strand = "EXCEPTION" 
        
        
        print("--------------------")            
        print("GeneID_NCBI:"+geneID)
        print("locustag:"+ locustag)
        print("gene name:"+ genename)
        print("Strand: "+ strand)
        
           
        
def run():
    fetch_sequence()
    geneIdentification() 

run()
