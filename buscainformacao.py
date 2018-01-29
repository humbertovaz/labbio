from Bio import Entrez
from Bio import SeqIO
from Record import Record
from io import StringIO
from RecordUniProt import RecordUniProt
import re
import requests
import anal_prot
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML



def fetch_sequence():
    genesStarteStop = [(87011,88213),(88210,89373),(1589969,1590817),(1903771,1904463),(3274457,3275026),(3295774,3296526),(3328322,3329788)]
    Entrez.email = "a71150@alunos.uminho.pt" #Email do André Germano
    i = 0
    # Ficheiros dos genes
    for tuplo in genesStarteStop:
        # Criar o ficheiro GenBank de cada Gene
        handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id="NZ_CP013742",
                               seq_start=tuplo[0], seq_stop=tuplo[1])
        sequence_filename = "gene"
        indice = str(i)
        name = sequence_filename + indice + ".gb"
        seq_record = SeqIO.read(handle, "gb")
        SeqIO.write(seq_record, name, "genbank")
        i += 1
        handle.close()


# UniProt localização celular, a função da proteína, o ID da Proteína
def uniprot(proteinids): # devolve dicionario de RecordUniprot
    proteinsRecords = {}
    for protein_id in proteinids:
        name, id, locale, status, fmol, bio, function2, leng, ec, seq, evalu, score = "-", "-", "-", "-", "-", "-", "-",\
                                                                                      "-", "-", "-", "-", "-"
        try:
            params = {"query": protein_id, "format": "fasta"}
            frecord = requests.get("http://www.uniprot.org/uniprot/", params)
            idfasta = ""
            for record in SeqIO.parse(StringIO(frecord.text), "fasta"):
                if "pneumophila (strain Philadelphia 1 / ATCC 33152 / DSM 7513" in str(record):
                    aux = record.id
                    seq = str(record.seq)
                    idfasta = str(aux).split("|")[1]
                    name, id, locale, status, fmol, bio, function2, leng, ec = anal_prot.getDataFromProt(idfasta)
                    rec = RecordUniProt(name, id, locale, status, ''.join(fmol), ''.join(bio), function2, leng, ec, seq)
                    proteinsRecords[idfasta] = rec


        except Exception:
            print("Try going online, no internet connection found.")
    return proteinsRecords

# Identificação do gene (GeneID NCBI,Accession number NCBI
# locus tag, nome do gene, strand)
def geneIdentification():
    dicionario = {}
    proteinids = []
    flag=0
    for i in range(7):
        index = str(i)
        name = "gene" + index + ".gb"
        # Ler o ficheiro
        record = SeqIO.read(name, "genbank")
        flag=0
        for feat in record.features:
            # print (feat.qualifiers)
            if (flag == 0):
                try:
                    pos = str(feat.location).replace("(+)","").replace("(-)","")
                except Exception:
                    strand = "EXCEPTION"
                flag=1
            else:
                try:
                    strand = str(feat.location)
                except Exception:
                    strand = "EXCEPTION"
                if pos in strand:
                    if "+" in strand:
                       strand = "(+)"
                    else:
                        strand = "(-)"
                    try:
                        locustag = feat.qualifiers["locus_tag"][0]
                    except Exception:
                        locustag = "-"
                    try:
                        genename = feat.qualifiers["gene"][0]
                    except Exception:
                        genename = "-"
                    try:
                        proteinID = feat.qualifiers["protein_id"][0]
                    except Exception:
                        proteinID = "EXCEPTION PROTEIN ID"
                    rec = Record(locustag, genename, strand)
                    dicionario[proteinID] = rec
    return dicionario


def csv_exportGenes(genedict):
    sep = ";"
    head = ""
    dataCSV = ""
    # Cria ficheiro .csv em modo de escrita
    filecsv = open("GenesTable.csv", "w")
    dataCSV += "Locus Tag" + sep + "Gene Name" + sep + "Strand" + sep + "Protein" + "\n"
    for key in genedict:
        rec = genedict[key]
        dataCSV += rec.getlocustag() + sep + rec.getgenename() + sep + rec.getstrand() + sep + key + "\n"
    filecsv.write(dataCSV + "\n")
    filecsv.close()
    return "GenesTable.csv"

def csv_exportProt(protdict):
    sep = ";"
    head = ""
    dataCSV = ""
    # Cria ficheiro .csv em modo de escrita
    filecsv = open("protTable.csv", "w")
    dataCSV += "Protein Name" + sep + "Id" + sep + "Funcao Molecular"+ sep + "Biological process" + sep + "Localizacao celular" + sep + "Status" + sep + "Function2" + sep + "Length" + sep + "EC" + sep + "Sequencia" + "\n"
    for key in protdict:
        rec = protdict[key]
        dataCSV += rec.getName() + sep + rec.getId() + sep + rec.getfmol() + sep + rec.getbio() +sep + rec.getlocale() + sep + rec.getstatus() + sep + rec.getfunction2() + sep + rec.getleng() + sep + rec.getEC() + sep + rec.getSeq() + "\n"
    filecsv.write(dataCSV + "\n")
    filecsv.close()
    return "ProtTable.csv"

def csv_exportProtBlast(protdict):
    sep = ";"
    head = ""
    dataCSV = ""
    # Cria ficheiro .csv em modo de escrita
    filecsv = open("BlastprotTable.csv", "w")
    dataCSV += "Protein Name" + sep + "Id" + sep + "Funcao Molecular"+ sep + "Biological process" + sep + "Localizacao celular" + sep + "Status" + sep + "Function2" + sep + "Length" + sep + "EC" + sep + "Evalue" + sep + "Score" + "\n"
    for key in protdict:
        rec = protdict[key]
        dataCSV += rec.getName() + sep + rec.getId() + sep + rec.getfmol() + sep + rec.getbio() +sep + rec.getlocale() + sep + rec.getstatus() + sep + rec.getfunction2() + sep + rec.getleng() + sep + rec.getEC() + sep + rec.getEvalu() + sep + rec.getScore() +  "\n"
    filecsv.write(dataCSV + "\n")
    filecsv.close()
    return "ProtTable.csv"

def uniprotBlast(proteinids): # devolve dicionario de RecordUniprot
    proteinsRecords = []
    evalu, score = "-", "-"
    for (protein_id, evalu, score) in proteinids:
        name, id, locale, status, fmol, bio, function2, leng, ec, seq = "-", "-", "-", "-", "-", "-", "-",\
                                                                                      "-", "-", "-"
        try:
            params = {"query": protein_id, "format": "fasta"}
            frecord = requests.get("http://www.uniprot.org/uniprot/", params)
            idfasta = ""
            for record in SeqIO.parse(StringIO(frecord.text), "fasta"):
                #if "pneumophila (strain Philadelphia 1 / ATCC 33152 / DSM 7513" in str(record):
                aux = record.id
                seq = str(record.seq)
                idfasta = str(aux).split("|")[1]
                name, id, locale, status, fmol, bio, function2, leng, ec = anal_prot.getDataFromProt(idfasta)
                rec = RecordUniProt(name, id, locale, status, ''.join(fmol), ''.join(bio), function2, leng, ec, seq, str(evalu), str(score))
                proteinsRecords.append(rec)


        except Exception:
            print("Try going online, no internet connection found.")
    return proteinsRecords



def blast(seq, database="swissprot", maxMatch=5): #retorna tuplos (id,evalu,score) de proteinas
    matches = 0
    idsGoodBlast = []
    print("blasting....")
    result = NCBIWWW.qblast("blastp", database, seq)
    print("blast done")
    save_file = open("blastteste.xml", "w")
    save_file.write(result.read())
    save_file.close()
    result.close()

    resultxml = open("blastteste.xml")
    blast = NCBIXML.parse(resultxml)

    E_VALUE_THRESH = 0.5
    SCORE = 95
    for blast_record in blast:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < E_VALUE_THRESH and hsp.score > SCORE:
                    data = (alignment.title.split("|")[3], hsp.expect, hsp.score)
                    idsGoodBlast.append(data)
                    matches += 1
                    if (matches >= maxMatch):
                        return idsGoodBlast
    return idsGoodBlast

# Apos recolher os 5 ids de cada proteina do blast, utilizar a funcao getInfoFromBlast que me vai buscar a Informacao à Uniprot acerca das melhores possiveis proteinas
def blastTodas(proteinsRecords):
    dicProtAposBlast = {}
    for idfasta in proteinsRecords:
        goodIds = blast(proteinsRecords[idfasta].getSeq(), "swissprot", 5)
        # Busca informacao dos resultados do Blast
        recordsUniprot = uniprotBlast(goodIds) # cada posicao em goodsIds é um tuplo c/ (Id, evalu, score) resultantes do blast

        dicProtAposBlast[idfasta] = recordsUniprot  # recordsUniprot é um array de records com a info de um recordUniprot

    return dicProtAposBlast # cada posicao do dicionario contem record com info relativa às proteinas geradas pelo blast apartir da fornecida (idfasta)

def csv_exportblast(protdict):
    sep = ";"
    head = ""
    dataCSV = ""
    # Cria ficheiro .csv em modo de escrita
    filecsv = open("blast.csv", "w")
    dataCSV += "Id" + sep + "Blast Results" + sep + "Score" + sep + "E_value" + sep + "Function" + "\n"
    for key in protdict:
        dataCSV += key + sep
        record = protdict[key]
        for rec in record:
            name,id,locale,status,fmol,function2,leng,ec,seq,evalu,score=rec.getAll()
            dataCSV += id + sep + score + sep + evalu + sep + fmol + "\n" + " "+ sep
    filecsv.write(dataCSV + "\n")
    filecsv.close()
    return "blast.csv"

def main():
    #Vai buscar a informação dos 7 genes
    fetch_sequence()
    #Cria um dicionario com Records dos genes
    genesdic = geneIdentification()
    #Cria CSV com tabela de genes e sua info. apartir do dicionario obtido acima
    csv_exportGenes(genesdic)
    #Cria um dicionario com "RecordsUniprot" de proteinas apartir do dicionario de genes 
    proteinsRecords = uniprot(genesdic)
    #Cria CSV com tabela com informação importante resultante da pesquisa na Uniprot
    csv_exportProt(proteinsRecords)
    #Executa o blast em todas as proteinas
    dicProtAposBlast = blastTodas(proteinsRecords)
    #Cria CSV com tabela resultante do blast
    csv_exportblast(dicProtAposBlast)
    print("Done!")

main()
