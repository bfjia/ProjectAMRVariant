from Bio import SeqIO
import os
from GenerateReqFiles import GenerateReqFiles
import ParseCardData
import random
import numpy as np
import pandas as pd
import re

class Codon:
    def __init__(self) -> None:
        self.__codonTable = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'STOP', 'TAG':'STOP',
            'TGC':'C', 'TGT':'C', 'TGA':'STOP', 'TGG':'W',
        }

    def DNAtoAA(self, DNA):
        return self.__codonTable[DNA]
    
    def AAtoDNA(self,AA):
        key_list = list(self.__codonTable.keys())
        val_list = list(self.__codonTable.values())   
        position = val_list.index(AA)
        return(key_list[position])
    
    def GetRandomAA(self):
        val_list = list(self.__codonTable.values())   
        return(random.sample(val_list,1)[0])

def GenerateFasta(proteinVariants, rnaVariants):
        
    codonTable = Codon()
    types = ['wt', 'mut', 'other']
    fasta = []
    for i in range(len(types)):
        random.seed(25041+i)
        chosenSet = random.sample(proteinVariants, int(len(proteinVariants)/5))

        for seq in chosenSet:
            dna = seq.aro.dna
            mutType = seq.mutType
            header = ">" + str(seq.aro.aro) + "|protein|" + types[i] + "|" + str(mutType) + "|"

            for snp in seq.snp:
                if (types[i] == 'wt'):
                    replacement = codonTable.AAtoDNA(snp.wt)
                elif (types[i] == 'mut'):
                    if (snp.mut == "fs"):
                        replacement = codonTable.AAtoDNA(snp.wt)
                    else:
                        replacement = codonTable.AAtoDNA(snp.mut)
                elif (types[i] == 'other'):
                    aa = codonTable.GetRandomAA()
                    while (aa == snp.mut or aa == snp.wt):
                        aa = codonTable.GetRandomAA()
                    replacement = codonTable.AAtoDNA(aa)
                dna = mutType.GenerateSyntheticDNAForProtein(dna, replacement, snp.position)
                header = header + str(snp.position) + ":" + codonTable.DNAtoAA(replacement) + ";"
            fasta.append(header + "\n" + dna) 
    with open ("SyntheticProteinVariants.fasta", "w") as f:
        f.write("\n".join(fasta))

    types = ['wt', 'mut', 'other']
    fasta = []
    for i in range(len(types)):
        random.seed(25041+i)
        chosenSet = random.sample(rnaVariants, int(len(rnaVariants)/5))

        for seq in chosenSet:
            dna = seq.aro.dna
            mutType = seq.mutType
            header = ">" + str(seq.aro.aro) + "|rna|" + types[i] + "|" + str(mutType) + "|"

            for snp in seq.snp:
                if (types[i] == 'wt'):
                    replacement = snp.wt
                elif (types[i] == 'mut'):
                    replacement = snp.mut
                elif (types[i] == 'other'):
                    replacement = "A"
                    while (replacement == snp.mut or replacement == snp.wt):
                        replacement = random.sample(["A","T","G","C"],1)[0]
                header = header + str(snp.position) + ":" +replacement + ";"
                dna = mutType.GenerateSyntheticDNAForRNA(dna, replacement, snp.position)
        
                #codon = codonTable.AAtoDNA(aa)
                #start = (snp.position -1) *3
                #end = (snp.position -1) *3 + 3
                #originalDNA = dna[start:end]
                #originalAA = codonTable.DNAtoAA(originalDNA)
                #dna = dna[:start] + codon + dna[end:]
            fasta.append(header + "\n" + dna)
    with open ("SyntheticRNAVariants.fasta", "w") as f:
        f.write("\n".join(fasta))

def CheckAccuracy(resultPath, fastaPath):
    with open(resultPath, "r") as f:
        variants = f.readlines()
    
    with open(fastaPath, 'r') as f:
        fasta = f.readlines()
    fasta = [x for x in fasta if x.startswith(">")]
    
    referenceVariants = {}
    for line in fasta:
        l = line.strip().replace(">","").split("|")
        aro = l[0]
        mutationType = l[2]
        mutationClass = l[3]
        snp = l[4]
        if aro in referenceVariants.keys():
            referenceVariants[aro].append([mutationType, mutationClass, snp])
        else:
            referenceVariants[aro] = []
            referenceVariants[aro].append([mutationType, mutationClass, snp])

    detectedVariants = {}
    for line in variants:
        l = line.strip().split("\t")
        aro = l[0]
        mutationType = l[2]
        mutationClass = l[1]
        try:
            snp = l[3]
            #additional = l[5]
        except:
            snp = "None"
            #additional = "None"
        if aro in detectedVariants.keys():
            detectedVariants[aro].append([mutationType, mutationClass, snp])
        else:
            detectedVariants[aro] = []
            detectedVariants[aro].append([mutationType, mutationClass, snp])
    
    keyDictionary = {"wt": 0, "mut": 1, "other":2, 
                    "Wildtype" : 0, "Resistant Variant" : 1, "Other Variant":2, "Not Found":3, "Partial":4}
    keyDictionary = {"wt": 0, "mut": 1, "other":2, 
                "Wildtype" : 0, "Resistant Variant" : 1, "Other Variant":2, "Not Found":3, "Partial":4}
    matrix = np.zeros ((3, 5))
    matrixAro = np.empty((3,5), dtype=object)
    for refARO in referenceVariants.keys():
        for refVariant in referenceVariants[refARO]:
            for detectedVariant in detectedVariants[refARO]:
                refNumbers = re.findall( '(\d+)', refVariant[2])
                detectedNumbers = re.findall( '(\d+)', detectedVariant[2])

                if (refNumbers == detectedNumbers):
                    matrix[keyDictionary[refVariant[0]],keyDictionary[detectedVariant[0]]] += 1 
                    label = refARO + "|" + refVariant[2] + "|" + detectedVariant[2] #+ "|" + detectedVariant[3]
                    print(label)
                    matrixAro[keyDictionary[refVariant[0]],keyDictionary[detectedVariant[0]]] = str(matrixAro[keyDictionary[refVariant[0]],keyDictionary[detectedVariant[0]]]) + str(label)
                # print("test")

    df = pd.DataFrame(matrix)
    df.columns = ["Wildtype", "ResistantVariant", "OtherVariant", "NotFound", "Partial"]
    df["Truth\\Pred"] = ["Wildtype", "ResistantVariant", "OtherVariant"]
    df = df.set_index("Truth\\Pred")
    df.to_csv("accuracy.tsv", sep ='\t', header = True , index = True)

    df = pd.DataFrame(matrixAro)
    df.columns = ["Wildtype", "ResistantVariant", "OtherVariant", "NotFound","Partial"]
    df["Truth\\Pred"] = ["Wildtype", "ResistantVariant", "OtherVariant"]
    df = df.set_index("Truth\\Pred")
    df.to_csv("accuracy_aro.tsv", sep ='\t', header = True , index = True)


#cardJsonPath = "/mnt/f/OneDrive/ProjectAMRSAGE/AMR_Metagenome_Simulator-master/full/card/card.json"
#cardSnpPath = "/mnt/f/OneDrive/ProjectAMRSAGE/AMR_Metagenome_Simulator-master/full/card/snps.txt"

#CARD = ParseCardData.ParseJson(cardJsonPath) 
#proteinVariants, rnaVariants = ParseCardData.ParseSNP(cardSnpPath, CARD)

#GenerateFasta()
resultPath = "./out.tsv"
fastaPath = "./SyntheticProteinVariants.fasta"
CheckAccuracy(resultPath, fastaPath)
#50 protein sequences of each of the following:
# 1) wt
# 2) mut
# 3) other

