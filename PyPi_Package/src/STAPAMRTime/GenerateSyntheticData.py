from . import DataTypes
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

def GenerateProteinSyntheticFasta(variants):
    #with open ("data.tsv", 'w') as f:
    #    for v in proteinVariants:
    #        f.write(v.Export())
    #    for v in rnaVariants:
    #        f.write(v.Export())

    #what do we want?
    #lets say we spike in total 99 AMR genes. 33wt, 33res, 33other
    #randomly choose 99 unique AROs.
    #for wt: make sure the base is wt then insert it into fasta.
    #for mut: make sure the base is mut then insert it
    #   additional: ensure all mutation categories (i.e. fs, single, multi) are in there.  
    #for other: make sure the base is not wt or mut then insert it. 


    codonTable = Codon()
    fasta = []
    uniqueARO = list(set([x.aro.GetARO() for x in variants]))
    if ("3003686" in uniqueARO):
        uniqueARO.remove("3003686")
    random.seed(25041)

    variantMutationTypes = list(set([x.GetMutType() for x in variants]))

    variantSets = {}
    n = 40
    random.shuffle(uniqueARO)
    missing = 0

    for type in variantMutationTypes:
        variantSets[type] = []
        for i in range(missing + int(n/len(variantMutationTypes))):
            added = False
            for aro in uniqueARO:
                seqList = [x for x in variants if x.aro.GetARO() == aro]
                #print(seq.GetMutType() + " : " + type)
                for seq in seqList:
                    if (seq.GetMutType() == type):
                        variantSets[type].append(aro)
                        uniqueARO.remove(aro)
                        added = True
                        break
                if (added):
                    break
        if (len(variantSets[type]) != int(n/len(variantMutationTypes))):
            for i in range(int(n/len(variantMutationTypes)) - len(variantSets[type])):
                added = False
                for aro in uniqueARO:
                    seqList = [x for x in variants if x.aro.GetARO() == aro]
                    #print(seq.GetMutType() + " : " + type)
                    for seq in seqList:
                        if (seq.GetMutType() == "Single"):
                            variantSets[type].append(aro)
                            uniqueARO.remove(aro)
                            added = True
                            break
                    if (added):
                        break
    variantSets["wt"] = uniqueARO[-n:]
    del uniqueARO[-n:]
    variantSets["other"] = uniqueARO[-n:]
    del uniqueARO[-n:]

    for type in variantSets.keys():
        for aro in variantSets[type]:
            variantObjects = [x for x in variants if x.aro.GetARO() == aro]
            if (type == "wt"):
                #take the first object with given ARO
                seq = variantObjects[0]
            elif (type == "other"):
                #again, take the first object with given ARO
                seq = variantObjects[0]
            else:
                #mutants, scan through objects until the muttype matches, return that match
                for vo in variantObjects:
                    if (vo.GetMutType() == type):
                        seq = vo
                        break
            
            dna = seq.aro.dna
            mutType = seq.mutType
            if (type == "wt"):
                label = "wt"
            elif (type == "other"):
                label = "other"
            else:
                label = "mut"
            header = ">" + str(seq.aro.aro) + "|protein|" + label + "|" + str(mutType) + "|"

            for snp in seq.snp:
                if (type == 'wt'):
                    replacement = codonTable.AAtoDNA(snp.wt)
                elif (type == 'other'):
                    aa = codonTable.GetRandomAA()
                    while (aa == snp.mut or aa == snp.wt):
                        aa = codonTable.GetRandomAA()
                    replacement = codonTable.AAtoDNA(aa)
                else: # (type == 'mut'):
                    if (snp.mut == "fs"):
                        replacement = codonTable.AAtoDNA(snp.wt)
                    elif (snp.mut == "*"):
                        replacement = codonTable.AAtoDNA("STOP")
                    else:
                        replacement = codonTable.AAtoDNA(snp.mut)

                dna = mutType.GenerateSyntheticDNAForProtein(dna, replacement, snp.position)
                header = header + str(snp.position) + ":" + codonTable.DNAtoAA(replacement) + ";"
            fasta.append(header + "\n" + dna) 
    with open ("SyntheticProteinVariants.fasta", "w") as f:
        f.write("\n".join(fasta))

def GenerateRRNASyntheticFasta(variants):


    #what do we want?
    #lets say we spike in total 99 AMR genes. 33wt, 33res, 33other
    #randomly choose 99 unique AROs.
    #for wt: make sure the base is wt then insert it into fasta.
    #for mut: make sure the base is mut then insert it
    #   additional: ensure all mutation categories (i.e. fs, single, multi) are in there.  
    #for other: make sure the base is not wt or mut then insert it. 

    codonTable = Codon()
    fasta = []
    uniqueARO = list(set([x.aro.GetARO() for x in variants]))
    if ("3003686" in uniqueARO):
        uniqueARO.remove("3003686")
    random.seed(25041)

    variantMutationTypes = list(set([x.GetMutType() for x in variants]))

    variantSets = {}
    n = 40
    random.shuffle(uniqueARO)
    missing = 0

    for type in variantMutationTypes:
        variantSets[type] = []
        for i in range(missing + int(n/len(variantMutationTypes))):
            added = False
            for aro in uniqueARO:
                seqList = [x for x in variants if x.aro.GetARO() == aro]
                #print(seq.GetMutType() + " : " + type)
                for seq in seqList:
                    if (seq.GetMutType() == type):
                        variantSets[type].append(aro)
                        uniqueARO.remove(aro)
                        added = True
                        break
                if (added):
                    break
        if (len(variantSets[type]) != int(n/len(variantMutationTypes))):
            for i in range(int(n/len(variantMutationTypes)) - len(variantSets[type])):
                added = False
                for aro in uniqueARO:
                    seqList = [x for x in variants if x.aro.GetARO() == aro]
                    #print(seq.GetMutType() + " : " + type)
                    for seq in seqList:
                        if (seq.GetMutType() == "Single"):
                            variantSets[type].append(aro)
                            uniqueARO.remove(aro)
                            added = True
                            break
                    if (added):
                        break
    variantSets["wt"] = uniqueARO[-n:]
    del uniqueARO[-n:]
    variantSets["other"] = uniqueARO[-n:]
    del uniqueARO[-n:]

    for type in variantSets.keys():
        for aro in variantSets[type]:
            variantObjects = [x for x in variants if x.aro.GetARO() == aro]
            if (type == "wt"):
                #take the first object with given ARO
                seq = variantObjects[0]
            elif (type == "other"):
                #again, take the first object with given ARO
                seq = variantObjects[0]
            else:
                #mutants, scan through objects until the muttype matches, return that match
                for vo in variantObjects:
                    if (vo.GetMutType() == type):
                        seq = vo
                        break
            
            dna = seq.aro.dna
            mutType = seq.mutType
            if (type == "wt"):
                label = "wt"
            elif (type == "other"):
                label = "other"
            else:
                label = "mut"
            header = ">" + str(seq.aro.aro) + "|rna|" + label + "|" + str(mutType) + "|"

            for snp in seq.snp:
                if (type == 'wt'):
                    replacement = snp.wt
                elif (type == 'other'):
                    replacement = "A"
                    while (replacement == snp.mut or replacement == snp.wt):
                        replacement = random.sample(["A","T","G","C"],1)[0]
                else:
                    replacement = snp.mut
                header = header + str(snp.position) + ":" +replacement + ";"
                dna = mutType.GenerateSyntheticDNAForRNA(dna, replacement, snp.position)
        
                #codon = codonTable.AAtoDNA(aa)
                #start = (snp.position -1) *3
                #end = (snp.position -1) *3 + 3
                #originalDNA = dna[start:end]
                #originalAA = codonTable.DNAtoAA(originalDNA)
                #dna = dna[:start] + codon + dna[end:]
            fasta.append(header + "\n" + dna)

    with open ("SyntheticRRnaVariants.fasta", "w") as f:
        f.write("\n".join(fasta))

def CheckAccuracy(resultTable, fastaPath):
    variants = resultTable
    
    with open(fastaPath, 'r') as f:
        fasta = f.readlines()
    fasta = [x for x in fasta if (x.startswith(">") and x.find("other") == -1)]
    #variants = [x for x in variants if (x.find("Resistant") > -1)]

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
        mutationClassification = l[3]
        snp = l[4]
        if aro in detectedVariants.keys():
            detectedVariants[aro].append([mutationType, mutationClassification, snp])
        else:
            detectedVariants[aro] = []
            detectedVariants[aro].append([mutationType, mutationClassification, snp])
    
    keyDictionary = {"wt": 0, "mut": 1, "other":0, 
                    "Wildtype" : 0, "Resistant Variant" : 1, "Other Variant":0, "Not Found":0, "Partial":0, "True":1, "False":0}

    matrix = np.zeros ((3, 5))
    matrixAro = np.empty((3,5), dtype=object)

    matrixBinary = np.zeros((2,2))
    matrixBinaryAro =  np.empty((2,2), dtype=object)
    matrixBinaryAro[0,0] = ""
    matrixBinaryAro[0,1] = ""
    matrixBinaryAro[1,0] = ""
    matrixBinaryAro[1,1] = ""
    sum = 0
    for refARO in referenceVariants.keys():
        #refARO is a list
        referenceVariants[refARO] = [list(x) for x in set(tuple(x) for x in referenceVariants[refARO] )]
        for refVariant in referenceVariants[refARO]:
            refPosition = list(set(re.findall( '(\d+)', refVariant[2])))
            if ("531" in refPosition):
                print("test")
            if (refARO in detectedVariants.keys()):
                detected = [x for x in detectedVariants[refARO] if list(set(re.findall( '(\d+)', x[2]))) == refPosition]
            else:
                detected = []

            if (len(detected) > 1):
                expected = [x.split(":")[1] for x in refVariant[2].strip().split(";") if x]
                expected = [x.replace("U", "T") for x in expected]
                expected.sort()
                found = False
                for i in detected:
                    bases = [x[-1].upper() for x in i[2].split(";")]
                    bases = [x.replace("*", "STOP") for x in bases]
                    bases.sort()
                    if expected == bases:
                        matrixBinary[keyDictionary[refVariant[0]],keyDictionary[i[1]]] += 1 
                        label = refARO + "|" + refVariant[2] + "|" + detectedVariant[2] #+ "|" + detectedVariant[3]
                        matrixBinaryAro[keyDictionary[refVariant[0]],keyDictionary[i[1]]] = matrixBinaryAro[keyDictionary[refVariant[0]],keyDictionary[i[1]]] + str(label) + ";"
                        found = True
                if(found == False):
                    #likely a wildtype because the mutant base is not found.
                    if (refVariant[0]=="wt"):
                        matrixBinary[keyDictionary[refVariant[0]],0] += 1 
                        label = refARO + "|" + refVariant[2] + "|" + detectedVariant[2] #+ "|" + detectedVariant[3]
                        matrixBinaryAro[keyDictionary[refVariant[0]],0] = matrixBinaryAro[keyDictionary[refVariant[0]],0] + str(label) + ";"
                    else:
                        #a resist type not ID'd
                        matrixBinary[keyDictionary[refVariant[0]],0] += 1 
                        label = refARO + "|" + refVariant[2] + "|" + detectedVariant[2] #+ "|" + detectedVariant[3]
                        matrixBinaryAro[keyDictionary[refVariant[0]],0] = matrixBinaryAro[keyDictionary[refVariant[0]],0] + str(label) + ";"
        
                        #print(refVariant)
                        #print(detected)
            elif(len(detected) == 1):
                detectedVariant = detected[0]
                matrixBinary[keyDictionary[refVariant[0]],keyDictionary[detectedVariant[1]]] += 1 
                label = refARO + "|" + refVariant[2] + "|" + detectedVariant[2] #+ "|" + detectedVariant[3]
                matrixBinaryAro[keyDictionary[refVariant[0]],keyDictionary[detectedVariant[1]]] = matrixBinaryAro[keyDictionary[refVariant[0]],keyDictionary[detectedVariant[1]]] + str(label) + ";"
            else:
                print("notfound")
            sum += 1
    totalRes = len([x for x in referenceVariants.values() if x[0][0] == "mut"])
    totalGroundTruth = sum
    
    #df = pd.DataFrame(matrixBinary)
    #df.columns = ["Non-Resistant", "Resistant"]
    #df["Truth\\Pred"] = ["Non-Resistant", "Resistant"]
    #df = df.set_index("Truth\\Pred")
    #df.to_csv("accuracy.tsv", sep ='\t', header = True , index = True)

    df = pd.DataFrame(matrixBinaryAro)
    df.columns = ["Non-Resistant", "Resistant"]
    df["Truth\\Pred"] = ["Non-Resistant", "Resistant"]
    df = df.set_index("Truth\\Pred")
    #df.to_csv("accuracy_aro.tsv", sep ='\t', header = True , index = True)

    matrix = AccuracyMetrics(matrixBinary[1,1], matrixBinary[0,0], matrixBinary[0,1], matrixBinary[1,0])

    return matrix, totalRes, totalGroundTruth

#cardJsonPath = "/mnt/d/OneDrive/ProjectAMRSAGE/AMR_Metagenome_Simulator-master/full/card/card.json"
#cardSnpPath = "/mnt/d/OneDrive/ProjectAMRSAGE/AMR_Metagenome_Simulator-master/full/card/snps.txt"

#CARD = ParseCardData.ParseJson(cardJsonPath) 
#proteinVariants, rnaVariants = ParseCardData.ParseSNP(cardSnpPath, CARD)


#GenerateProteinSyntheticFasta(proteinVariants)
#GenerateRRNASyntheticFasta(rnaVariants)

#resultPath = "./syntheticBalanced_DetectedVariants.tsv"
#with open (resultPath , 'r') as f:
#    variants = f.readlines()[1:]
#fastaPath = "./sim/SyntheticProteinVariants.fasta"
#test = CheckAccuracy(variants, fastaPath)
#print(test)
#50 protein sequences of each of the following:
# 1) wt
# 2) mut
# 3) other
#matrix, totalRes, total = CheckAccuracy("synthetic.tsv", "SyntheticVariants.fasta")
