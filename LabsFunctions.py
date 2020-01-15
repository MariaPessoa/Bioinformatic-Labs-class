# -*- coding: utf-8 -*-
"""
Created on Sat Jan 11 11:34:58 2020

@author: liamo
"""
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
Entrez.email = "INSERT VALID EMAIL HERE"
#from Bio import ExPASy, SwissProt


#NG_011631.1 RefSeqGene MITF
#NG_027728.1 RefSeqGene IRF4
#NG_012026.1 RefSeqGene MC1R - 2 from RefSeq
#MC1R NM_002386.4; MITF NM_198159.3; IRF4 NM_001195286.2                             
def genefetch(filename,query,start=None,stop=None): 
    '''
    The genefetch function relies on the Biopython Entrez package to download
    and save GenBank files from the nucleotide database. The function takes
    the name of the file to use and the DNA acession number to retrieve.
    If provided with start and stop arguments, the function retrieves only
    that part of the sequence.
    '''
    if start is None and stop is None:
        fetch_handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text",
                                     id=query)
    else: 
        fetch_handle = Entrez.efetch(db="nucleotide", rettype="gbwithparts",
                                     retmode="text", id=query, seq_start=start,
                                     seq_stop=stop)
    dna_fetch_record = fetch_handle.read()
    with open(filename, "w") as out_handle: out_handle.write(dna_fetch_record)
    return None

def geneinfo(filename):
    '''
    The geneinfo function take the name of file with one or more GenBank records.
    The function prints out some information in the file, such as the the 
    description, acession number and sequence length. It also prints out 
    information specific to "CDS" features in the GenBank file.
    '''
    with open(filename) as file: 
        file_genes = list(SeqIO.parse(file, "gb"))
    for i in range(len(file_genes)):
        #'annotations', 'dbxrefs', 'description', 'features', 'format', 'id',
        #'letter_annotations', 'name', 'reverse_complement', 'seq','translate'
        print(file_genes[i].description)
        print("Acession Number:",file_genes[i].id)
        print("Sequence length:",len(file_genes[i].seq))
        for feat in file_genes[i].features: #source,gene,mRNA,exon,CDS
            if feat.type == "CDS": 
                #'gene', 'gene_synonym', 'note', 'codon_start', 'product',
                #'protein_id', 'db_xref', 'translation'
                print("Gene:",feat.qualifiers['gene'])
                print("Protein:",feat.qualifiers['protein_id'])
                print("X-refs:",feat.qualifiers['db_xref'],end="\n\n")
    return None

#genefetch("GeneMITF.txt","NG_011631.1",4954,233856) #still too big!
#find protein > find smaller sequence (codes same protein) > NM_000248.3

def blast(filename,query): 
    '''
    The blast function relies on the Biopython Entrez package to perform a
    blastn (DNA vs DNA) BLAST for the DNA query sequence proveided, in the 
    nucleotide database, with set parameters. After running said blast, the 
    function saves the results in an xml file with the name provided.
    '''
    blast_handle = NCBIWWW.qblast("blastn", "nucleotide", query,
                                  megablast=True,gapcosts="0 3",word_size=28,
                                  nucl_reward=1, nucl_penalty=-2,filter=True,
                                  hitlist_size=10, expect=0.05)
    with open(filename, "w") as out_handle:
        out_handle.write(blast_handle.read())

def blastinfo(filename):
    '''
    The blastinfo function takes the name of a xml file with the results of
    a BLAST and prints out some of the parameters used, the max alignment
    length and the max score, along with their respective sequence accession
    numbers.
    '''
    with open(filename) as file:
        blast_record = NCBIXML.read(file)
    print("BLAST Parameters:")
    print("Query ID:",blast_record.query_id)
    print("Database:",blast_record.database)
    print("E-value threshold:",blast_record.expect)
    print("Match score:",blast_record.sc_match)
    print("Mismatch score:",blast_record.sc_mismatch)
    max_score = -999
    max_length = -999
    result = {}
    for align in blast_record.alignments:
        temp = align.title.split("|")[3]
        if align.length > max_length: 
            max_length = align.length
            result["Max length"] = (max_length,temp)
        for hsp in align.hsps:
            if hsp.expect < 0.05:
                if hsp.score > max_score:
                    max_score = hsp.score
                    result["Max score"] = (max_score,temp)
    print(result)
    return None

def swissprot(filename):
    '''
    The swissprot function takes the name of a swissprot file in text format
    as an argument and prints the "ID","GN","PROSITE" & "SQ" fields. 
    It is meant to be a simple stand-in for the Biopython version, which has 
    incompatiblity problems with the new swissprot format (15/1/2020).
    '''
    #handle = ExPASy.get_sprot_raw(query)
    #record = SwissProt.read(open(filename))
    #record = SeqIO.read(open(filename),"uniprot-xml")
    with open(filename) as file: file_list = file.readlines()
    for i in range(len(file_list)):
        if "ID" in file_list[i][:5]: print(file_list[i].replace("    ",""),end="")
        if "GN" in file_list[i][:5]: print(file_list[i],end="")
        if "PROSITE" in file_list[i][:20]: print(file_list[i],end="")
        if "SQ" in file_list[i][:5]: 
            print(" ".join(file_list[i].split("  ")[:-2]))
    return None

def menu():
    while True:
        print("""
        ======== Functions ============
        0. Quit
        1.genefetch
        2.geneinfo
        3.blast
        4.blastinfo
        5.swissprot
        """)
        op = input("Choose an option: ")
        if op == "1": 
            file = str(input("Filename? "))
            query = str(input("Query accession number? "))
            start = input("Start position? (Enter if none) ")
            stop = input("Stop position? (Enter if none) ")
            if len(file) > 0 and len(query) > 0:
                if len(start) > 0 and len(stop) > 0:
                    genefetch(file,query,int(start),int(stop))
                else: genefetch(file,query)
            else: print("Invalid option")
        elif op == "2": 
            user = int(input("IRF4 & MC1R(1) or MITF(2)? "))
            if user == 1: geneinfo("GenesMC1R_IRF4.txt")
            elif user == 2: geneinfo("GeneMITF.txt")
            else: print("Invalid option")
        elif op == "3": 
            file = str(input("Filename? "))
            query = str(input("Query accession number? "))
            if len(file) > 0 and len(query) > 0:
                blast(file,query)
            else: print("Invalid option")
        elif op == "4": 
            user = int(input("IRF4(1), MC1R(2) or MITF(3)? "))
            if user == 1: blastinfo("Blast_irf4_on.xml")
            elif user == 2: blastinfo("Blast_mc1r_on.xml")
            elif user == 3: blastinfo("Blast_mitf_on.xml")
            else: print("Invalid option")
        elif op == "5": 
            user = int(input("IRF4(1), MC1R(2) or MITF(3)? "))
            if user == 1: swissprot("Q15306.txt")
            elif user == 2: swissprot("Q01726.txt")
            elif user == 3: swissprot("O75030.txt")
            else: print("Invalid option")
        elif op == "0":
            break
        else:
            print("\n Invalid option.")
    print("\n Goodbye")


if __name__ == "__main__":
    menu()
    pass