"""
Provides underlying datastructure and functions for variant storage and detection

Structs for storing CARD database, SNP information, Variant information, Mutation Types etc.
"""

import itertools
import re
from Bio.Seq import Seq
import math
import numpy


class AccuracyMetrics():
    def __init__(self, TP, TN, FP, FN) -> None:
        self.tp = numpy.float64(TP)
        self.tn = numpy.float64(TN)
        self.fp = numpy.float64(FP)
        self.fn = numpy.float64(FN)

    def GetTPR(self):
        #sensitivity, recall, hit rate, true positive rate
        return self.tp/(self.tp + self.fn)
    
    def GetTNR(self):
        #specificity, true negative rate
        return self.tn/(self.tn+self.fp)
    
    def GetPPV(self):
        #precision, positive predictive value
        return self.tp/(self.tp+self.fp)

    def GetNPV(self):
        #negative predictive value
        return self.tn/(self.tn+self.fn)
    
    def GetFNR(self):
        #miss rate, false negative rate
        return self.fn/(self.fn+self.tp)
    
    def GetFPR(self):
        #fall out rate or false positive rate
        return self.fp/(self.fp+self.tn)
    
    def GetFDR(self):
        #false discovery rate
        return self.fp/(self.fp/self.tp)
    
    def GetFOR(self):
        #false omission rate
        return self.fn/(self.fn+self.tn)
    
    def GetAccuracy(self):
        return (self.tp+self.tn) / (self.tp+self.tn+self.fp+self.fn)
    
    def GetF1(self):
        return  (2 * self.tp) / (2*self.tp + self.fp+self.fn)
    
    def GetMCC(self):
        return (self.tp*self.tn - self.fp*self.fn)/math.sqrt((self.tp+self.fp) + (self.tp+self.fn) + (self.tn+self.fp) + (self.tn+self.fn))

    def GetFM(self):
        #fm index
        return math.sqrt(self.GetPPV() / self.GetTPR())
    
    def GetAllMetrics(self):
        header = ["TP", "TN", "FP", "FN", "TPR", "TNR", "PPV", "NPV", "FNR", "FPR","FDR","FOR","Accuracy", "F1", "FM", "MCC" ]
        values = [self.tp, self.tn, self.fp, self.fn, self.GetTPR(), self.GetTNR(), self.GetPPV(), self.GetNPV(), self.GetFNR(), self.GetFPR(), self.GetFDR(), self.GetFOR(), self.GetAccuracy(), self.GetF1(), self.GetFM(), self.GetMCC()]
        return values, header

    def WriteConfusionMatrix(self, filename, col=["Positive", "Negative"], row=["Positive", "Negative"]):
        output = [
            "Truth\\Pred\t" + "\t".join(col),
            "\t".join([row[0], self.tp, self.fn]),
            "\t".join([row[1], self.fp, self.tn])
        ]
        with open (filename, 'w') as f:
            f.writelines(output)

#structure for storing useful information from CARD.json
class ARO:
    def __init__(self, aro, name, cvterm, species,dna,protein) :
        self.aro = aro
        self.name = name
        self.cvterm = cvterm
        self.species = species
        self.dna = dna
        self.protein = protein
    
    def GetARO(self):
        return str(self.aro)

#structure for sotring snps
class SNP:
    def __init__(self, wt, mut, position, depth = None, totalDepth = None, codon = None):
        self.wt = wt
        self.mut = mut
        self.position = position #1 based index of base location
        self.depth = depth
        self.totalDepth = totalDepth
        self.codon = codon
    
    def IsIndel(self):
        return False
    
    def __str__(self) -> str:
        return ("{}{}{}".format(self.wt, str(self.position), self.mut))

#structure for sotring Indels
class Indel(SNP):
    def __init__(self, wt, position, sequence, type, depth=None, totalDepth = None, codon=None):
        SNP.__init__(self, wt = wt, mut = "Indel", position = position, depth=depth, totalDepth=totalDepth, codon=codon)
        self.sequences = sequence
        self.type = type

    def IsIndel(self):
        return True

#Base class for storing Variant information.
class Variant:
    def __init__(self, aro, snp, mutType, geneType, detectedSnp = None, detectedIndel = None):
        self.aro = aro
        self.mutType = mutType
        self.geneType = geneType
        self.snp = snp
        self.detectedSnp = detectedSnp
        self.detectedIndel = detectedIndel
    
    def Export(self) -> str:
        return("{}\t{}\t{}\t{}\t{}\n".format(str(self.aro.aro), self.geneType, self.mutType, ",".join([str(x) for x in self.snp]), self.aro.dna))

    def SetDetectedSnp(self, snp):
        self.detectedSnp = snp

    def SetDetectedIndel(self, snp):
        self.detectedIndel = snp

    def FindVariantInVCF(self, vcf):
        raise ("Default Variant Class does not have a function for FindVariantInVCF. Use one of the child classes")
        
    def ParsePileupReadBases(self, refBase, baseString):
        #shamelessly stolen from https://github.com/Niknafs/NGSTools/blob/master/baseParser.py
        # remove end of read character
        ref = refBase.upper()
        string = baseString.upper()
        string = string.replace('$','')
        types = {'A':0,'G':0,'C':0,'T':0,'-':0,'+':0,'X':[], "ins":[], "del":[]}
        

        while string != '':
            if string[0] == '^':
                # skip two characters when encountering '^' as it indicates
                # a read start mark and the read mapping quality
                string = string[2:]
            elif string[0] == '*':
                types['-'] += 1
                types['del'].append(ref)
                # skip to next character
                string = string[1:]
            
            elif string[0] in ['.',',']:
                if (len(string)== 1) or (string[1] not in ['+','-']):
                    # a reference base
                    types[ref] += 1
                    string = string[1:]
                elif string[1] == '+': 
                    types[ref] += 1
                    insertionLength = int(re.match("([0-9]*)",string[2:]).groups()[0])#int(list(filter(str.isdigit, (string[2:5])))[0]) #hack for double+ digit insertion length
                    insertionSeq = string[3:3+ insertionLength]
                    types['ins'].append(insertionSeq)
                    types['+'] += 1
                    string = string[3+insertionLength:]
                elif string[1] == '-':
                    types[ref] += 1
                    deletionLength = int(re.match("([0-9]*)",string[2:]).groups()[0])
                    deletionSeq = string[3:3+deletionLength]
                    types['del'].append(deletionSeq)
                    types['-'] += 1
                    string = string[3+deletionLength:]
                    
            elif (string[0] in types.keys()) and ((len(string)==1) or (string[1] not in ['-','+'])):
                # one of the four bases
                types[string[0]] += 1
                string = string[1:]
            else:
                # unrecognized character
                # or a read that reports a substitition followed by an insertion/deletion
                if (string[1] == "-"):
                    deletionLength = int(re.match("([0-9]*)",string[2:]).groups()[0])
                    deletionSeq = string[3:3+deletionLength]
                    types['del'].append(deletionSeq)
                    types['-'] += 1
                    string = string[3+deletionLength:]
                elif (string[1] == "+"):
                    insertionLength = int(re.match("([0-9]*)",string[2:]).groups()[0])#int(list(filter(str.isdigit, (string[2:5])))[0]) #hack for double+ digit insertion length
                    insertionSeq = string[3:3+ insertionLength]
                    types['ins'].append(insertionSeq)
                    types['+'] += 1
                    string = string[3+insertionLength:]
                else:
                    types['X'].append(string[0])
                    string = string[1:]
        return types
    
    #returns a list of class:SNPs at position X
    def ParsePileupForSNPs(self, relevantSNP, snp):
        detectedSNPs = []
        if(len(relevantSNP) < 1):
            exist = False
        else:
            for s in relevantSNP:
                data = s.strip("\n").split('\t')
                wt = data[2]
                totalDepth = data[3]
                bases = data[4]
                baseTypes = self.ParsePileupReadBases(wt, bases)
                #for key in baseTypes.keys():
                    #types = {'A':0,'G':0,'C':0,'T':0,'-':[],'*':0,'+':[],'X':[]}
                if baseTypes[wt.upper()] > 0:
                    detectedSNPs.append(SNP(wt, wt, snp.position, depth = baseTypes[wt.upper()], totalDepth=totalDepth))
                for b in "ATGC".replace(wt.upper(),""):
                    if baseTypes[b] > 0:
                       detectedSNPs.append(SNP(wt, b, snp.position, depth = baseTypes[b], totalDepth=totalDepth))
                if baseTypes["+"] > 0:
                    detectedSNPs.append(Indel(wt, snp.position, baseTypes["ins"], "Insertion", depth = baseTypes['+'], totalDepth=totalDepth ))
                if baseTypes["-"] > 0:
                    detectedSNPs.append(Indel(wt, snp.position, baseTypes["del"], "Deletion", depth = baseTypes['-'],totalDepth=totalDepth ))
                if len(baseTypes["X"]) > 0:
                    return "???"
        return detectedSNPs
    """
    def FindIndelInVCF(self, vcf):
        detectedIndel = []
        relevantVCF = [x for x in vcf if x.startswith("{}\t".format(str(self.aro.aro)))]
        
        indels = [x for x in relevantVCF if 'INDEL' in x]
        for indel in indels:
            data = self.SplitVCFLine(indel)
            wt = "INDEL"
            position = data[1]
            if (int(position) <=12 or int(position) >len(relevantVCF) - 12):
                continue #ignores the first and last 12 bases, theres a lot of indels in there.
            else:
                mut = data[4]
                if (mut == "."):
                    continue
                quality = float(data[5])
                info = data[7].replace(".","0").split(";")
                #INDEL;IDV=1;IMF=0.0128205;DP=78;VDB=5.91414e-32;SGB=-0.379885;MQSB=6.22462e-06;MQ0F=0.153846;AN=1;DP4=45,32,1,0;MQ=21
                dp = int([x for x in info if x.startswith("DP=")][0].replace("DP=",""))
                dp4 = list(map(int, list(([x for x in info if x.startswith("DP4=")][0].replace("DP4=","").split(",")))))
                mq = int([x for x in info if x.startswith("MQ=")][0].replace("MQ=","")) 
                detectedIndel.append(SNP(wt,mut,position,quality,dp, dp4,mq))
                indels = [x for x in relevantVCF if 'INDEL' in x]

        self.detectedIndel = detectedIndel
        return self
    """
    def GetType(self):
        return self.geneType
    
    def GetMutType(self):
        return str(self.mutType)
        
    def __str__(self) -> str:
        return(str(self.geneType).capitalize() + ":" + str(self.mutType).capitalize())

#child of Variant class specific for RNA variants
class RNAVariant(Variant):
    def __init__(self, aro, snp, mutType):
        for s in snp:
            if (s.mut.lower() == "u"):
                s.mut = "T"
            if (s.wt.lower() == "u"):
                s.wt = "T"
        Variant.__init__(self, aro, snp, mutType, "rna variant")


    def FindVariantInPileup(self, pileup):
        detectedSNPs = []
        detectedIndels = []
        pileup = [x for x in pileup if x.startswith("{}\t".format(str(self.aro.aro)))]

        for snp in self.snp:
            relevantSNP = [x for x in pileup if x.startswith("{}\t{}\t".format(str(self.aro.aro), str(int(snp.position))))]
            indelSnp = [x for x in pileup if x.startswith("{}\t{}\t".format(str(self.aro.aro), str(int(snp.position)-1)))]
            
            base0indel = []
            base1indel = []
            
            if (relevantSNP):
                for s in self.ParsePileupForSNPs(relevantSNP, snp):
                    if (s.IsIndel()):
                        base1indel.append(s)
                    else:
                        detectedSNPs.append(s)
                for s in self.ParsePileupForSNPs(indelSnp, snp):
                    if (s.IsIndel()):
                        base0indel.append(s)
            
            detectedIndels.append([base0indel, base1indel])

        self.SetDetectedSnp(detectedSNPs)
        self.SetDetectedIndel(detectedIndels)
        return self

#Child of Variant class specific for Protein Variants
class ProteinVariant(Variant):
    def __init__(self, aro, snp, mutType, genetype = "protein variant"):
        Variant.__init__(self, aro, snp, mutType, genetype)
    
    def FindVariantInPileup(self, pileup):
        detectedSNPs = []
        detectedIndels = []
        pileup = [x for x in pileup if x.startswith("{}\t".format(str(self.aro.aro)))]

        for snp in self.snp:
            bases = [[],[],[],[]]
            bases[0] = [x for x in pileup if x.startswith("{}\t{}\t".format(str(self.aro.aro), str(int(snp.position)*3-3)))]
            bases[1] = [x for x in pileup if x.startswith("{}\t{}\t".format(str(self.aro.aro), str(int(snp.position)*3-2)))]
            bases[2] = [x for x in pileup if x.startswith("{}\t{}\t".format(str(self.aro.aro), str(int(snp.position)*3-1)))]
            bases[3] = [x for x in pileup if x.startswith("{}\t{}\t".format(str(self.aro.aro), str(int(snp.position)*3)))]

            if (bases[0], bases[1] and bases[2] and bases[3]):
                baseSnp = [{},{},{},{}]
                baseIndel = [[],[],[],[]]
                idx = 0 #counter in case of inserts
                for position in range(4):
                    baseVariant = self.ParsePileupForSNPs(bases[position], snp)
                    for b in baseVariant:
                        if b.IsIndel():
                            baseIndel[position].append(b)
                            #single base inserts followed by a deletion might be actual variants
                            insertBase = [x for x in set(b.sequences) if len(x) == 1]
                            if b.type == "Insertion":
                                idx = idx + 1
                                for s in insertBase:
                                    if idx < 4:
                                        baseSnp[idx][s] = b
                            #if b.type == "Deletion":
                             #   idx = idx -1
                        else:
                            if idx < 4:
                                baseSnp[idx][b.mut] = b
                    idx += 1

                if (baseIndel[0] or baseIndel[1] or baseIndel[2] or baseIndel[3]):
                    detectedIndels.append([baseIndel[0] , baseIndel[1] , baseIndel[2] , baseIndel[3]])

                #sometimes, a snp can be recorded as an indel. so we first gotta check if there are indels at position

                combinations = list(itertools.product(*[baseSnp[1].keys(), baseSnp[2].keys(), baseSnp[3].keys()]))
                for combo in combinations:
                    wt = snp.wt
                    codon = "".join(list(combo))
                    mut = str(Seq(codon).translate())
                    lowestDepth = min((baseSnp[1][codon[0]].depth, baseSnp[2][codon[1]].depth, baseSnp[3][codon[2]].depth))
                    maxTotalDepth = max (baseSnp[1][codon[0]].totalDepth, baseSnp[2][codon[1]].totalDepth, baseSnp[3][codon[2]].totalDepth)
                    detectedSNPs.append(SNP(wt, mut, snp.position, lowestDepth, maxTotalDepth, codon))

        self.SetDetectedSnp(detectedSNPs)
        self.SetDetectedIndel(detectedIndels)
        return self

#child of ProteinVariant class specific for Protein overexpression
class ProteinOverexpressionVariant(ProteinVariant):
    def __init__(self, aro, snp, mutType):
        ProteinVariant.__init__(self, aro, snp, mutType, "protein overexpression variant")



#BaseClass for different mutation types. Used to store rules for parsing and filtering SNPs
class MutationType:
    def __init__(self, type) :
        self.__type = type
        self.__snps = []

    def GetType(self):
        return self.__type
    
    def ParseSNP(self, aro, mutations):
        return self.OnFailure("Error: BaseClass has no default ParseSNP function, use one of the child classes" )
    
    #Synthetic Seq Generator: return a DNA sequence containing the desied amino acid mutation
    def GenerateSyntheticDNAForProtein(self, dna, replacement, position):
        _dna = dna
        start = (position -1) *3
        end = (position -1) *3 + 3
        _dna = dna[:start] + replacement + dna[end:]
        return _dna

    #Synthetic Seq Generator: return DNA sequence with desired nucleotide mutation
    def GenerateSyntheticDNAForRNA(self, dna, replacement, position):
        _dna = dna
        start = (position -1) 
        _dna = dna[:start] + replacement + dna[start+1:]
        return _dna

    def Classify(self, variant):
        raise Exception("Base MutationType has no Classify function, use one of the childrens")

    #generic functions for onsuccess and on failure logics
    def OnSuccess(self, snps):
        return True, snps
    def OnFailure(self, aro, mutations, e):
        return False, ("{} error for {}:{}, {}".format(str(self.GetType()), aro, mutations, str(e)))
        
    def __str__(self) -> str:
        return(self.__type.capitalize())

#single mutations. i.e. A123L
class SingleMutationType(MutationType):
    def __init__(self) :
        MutationType.__init__(self, "single")
        #self.__type="single"
    
    #logic for parsing snp.txt for single mutations
    def ParseSNP(self, aro, mutations):
        try:
            #e.g. L527V
            mutations = mutations.strip()
            wt = mutations[0]
            mut = mutations[-1]
            position = int(mutations[1:-1])
            snp = [SNP(wt, mut, position)]
            return self.OnSuccess(snp)
        except Exception as e:
            return self.OnFailure(aro, mutations, e)
    
    def ClassifySNP(self, variant):
        #we want to compare self.snp versus self.detectedsnp. 
        #single mutation is easy, just check for mut. 
        if len(variant.detectedSnp) ==0:
            return "Not Found", None

        for cardSnp in variant.snp:
            snpCollection = []
            other = False
            wildtype = False
            resistant = False
            for detectedSnp in variant.detectedSnp:
                if detectedSnp.position == cardSnp.position:
                    if detectedSnp.mut == cardSnp.mut:
                        resistant = True
                        snpCollection.append(detectedSnp)
                    elif detectedSnp.mut == cardSnp.wt:
                        wildtype = True
                        snpCollection.append(detectedSnp)
                    else:
                        other = True
                        snpCollection.append(detectedSnp)
            if (resistant):
                return "Resistant Variant", snpCollection
            elif (wildtype):
                return "Wildtype", snpCollection
            elif (other):
                return "Other Variant", snpCollection
            else:
                return "???", None

#multiclass mutations. i.e. A123L, L233Q, K250I
class MultipleMutationType(MutationType):
    def __init__(self) :
        MutationType.__init__(self,"multi")
        #self.__type="multi"
    
    #logic for parsing snp.txt for multi mutations
    def ParseSNP(self, aro, mutations):
        try:
            #e.g. G452C,R659L
            snps = []
            for m in mutations.split(","):
                mutation = m.strip()
                wt = mutation[0]
                mut = mutation[-1]
                position = int(mutation[1:-1])
                snps.append(SNP(wt,mut,position))
            return self.OnSuccess(snps)
        except Exception as e:
            return self.OnFailure(aro, mutations, e)
    
    def ClassifySNP(self, variant):
        #we want to compare self.snp versus self.detectedsnp. 
        #single mutation is easy, just check for mut. 
        if len(variant.detectedSnp) ==0:
            return "Not Found", None

        #sometimes, a snp can be recorded as an indel. so we first gotta check if there are indels at position

        flags = [0] * len(variant.snp)
        detected = [None] * len(variant.snp)
        for i in range(len(variant.snp)):
            cardSnp = variant.snp[i]
            wildtype = False
            other = False
            mutant = False
            snpCollection = []
            for detectedSnp in variant.detectedSnp:
                if detectedSnp.position == cardSnp.position:
                    if detectedSnp.mut == cardSnp.mut:
                        mutant = True
                        snpCollection.append(detectedSnp)
                        break
                    elif detectedSnp.mut == cardSnp.wt:
                        wildtype = True
                        snpCollection.append(detectedSnp)
                    else:
                        other = True
                        snpCollection.append(detectedSnp)
            if (mutant):
                flags[i] = 2
                detected[i] = snpCollection
            elif (other):
                flags[i] = 3
                detected[i] = snpCollection
            elif (wildtype):
                flags[i] = 1
                detected[i] = snpCollection
            else:
                return "Not Found", None

        f = list(set(flags))
        if (len(f) == 1):
            if f[0] == 0:
                return "Not Found", None
            elif f[0] == 1:
                return "Wildtype", detected
            elif f[0] == 2:
                return "Resistant Variant", detected
            elif f[0] == 3:
                return "Other Variant", detected
        else:
            if 2 in f:
                return "Partial", detected
            elif 1 in f:
                return "Wildtype", detected
            elif 3 in f:
                return "Other Variants", detected
            else:
                return "????", detected


#nonsense mutations. i.e. A123STOP
class NonSenseMutationType(MutationType):
    def __init__(self) :
        MutationType.__init__(self,"nonsense")
        #self.__type="nonsense"

    #logic for parsing snp.txt for NS mutations
    def ParseSNP(self, aro, mutations):
        try:
            #e.g. R279STOP
            mutations = mutations.strip()
            wt = mutations[0]
            mut = mutations[-4:]
            if (mut.lower() != "stop"):
                raise Exception("Mut is not STOP for Object:NonSenseMutationType")
            else:
                mut = "*"
            position = int(mutations[1:-4])
            snp = [SNP(wt, mut, position)]
            return self.OnSuccess(snp)
        except Exception as e:
            return self.OnFailure(aro, mutations, e)
    
    def ClassifySNP(self, variant):
        #we want to compare self.snp versus self.detectedsnp. 
        #single mutation is easy, just check for mut. 
        if len(variant.detectedSnp) ==0:
            return "Not Found", None

        for cardSnp in variant.snp:
            other = False
            wildtype = False
            resistant = False
            snpCollection = []
            for detectedSnp in variant.detectedSnp:
                if detectedSnp.position == cardSnp.position:
                    if detectedSnp.mut == cardSnp.mut:
                        resistant = True
                        snpCollection.append(detectedSnp)   
                    elif detectedSnp.mut == cardSnp.wt:
                        wildtype = True
                        snpCollection.append(detectedSnp)
                    else:
                        other = True
                        snpCollection.append(detectedSnp)
            if(resistant):
                return "Resistant Variant", snpCollection
            elif (wildtype):
                return "Wildtype", snpCollection
            elif (other):
                return "Other Variant", snpCollection
            else:
                return "???", None


#frameshift mutations. i.e. A123FS
class FrameshiftMutationType(MutationType):
    def __init__(self): 
        MutationType.__init__(self,"frameshift")
        #Aself.__type="frameshift"
    
    #logic for parsing snp.txt for fs mutations
    def ParseSNP(self, aro, mutations):
        try:
            #e.g. R279STOP
            mutations = mutations.strip()
            wt = mutations[0]
            mut = mutations[-2:]
            if (mut.lower() != "fs"):
                raise Exception("Mut is not FS for Object:FrameshiftMutationType")
            position = int(mutations[1:-2])
            snp = [SNP(wt, mut, position)]
            return self.OnSuccess(snp)
        except Exception as e:
            return self.OnFailure(aro, mutations, e)    
    
    def ClassifySNP(self, variant):
        #we want to compare self.snp versus self.detectedsnp. 
        #single mutation is easy, just check for mut. 
        if len(variant.detectedIndel) ==0 and len(variant.detectedSnp) ==0:
            return "Not Found", None
            

        for cardSnp in variant.snp:
            if (variant.detectedIndel):
                for indel in variant.detectedIndel:
                    for base in indel:
                        if ((len(base) % 3) != 0):
                            return "Resistant Variant", indel
            else:
                return "Wildtype", variant.detectedSnp


    
    #functions to replace baseclass Synthetic Seq Generators
    def GenerateSyntheticDNAForProtein(self, dna, replacement, position):
        _dna = dna
        start = (position -1) *3
        end = (position -1) *3 + 3
        _dna = dna[:start] + "A" + replacement + dna[end:]
        return _dna
    
    def GenerateSyntheticDNAForRNA(self, dna, replacement, position):
        _dna = dna
        start = (position -1) *3
        end = (position -1) *3 + 3
        _dna = dna[:start] + "A" + replacement + dna[end:]
        return _dna

#frameshift mutations. i.e. ARO300333:A123L+ARO300444:L321A
class CodependentMutationType(MutationType):
    def __init__(self) :
        MutationType.__init__(self,"co-dependent")
        #self.__type="co-dependent"
    
    #unsurported for parsing
    def ParseSNP(self, aro, mutations):
        return self.OnFailure(aro, mutations, "Parsing for CodependentMutationType is unsupported")

#indel mutations.
class IndelMutationType(MutationType):
    def __init__(self) :
        MutationType.__init__(self,"indel")
        #self.__type="indel"
    
    #unsurported for parsing
    def ParseSNP(self, aro, mutations):
        return self.OnFailure(aro, mutations, "Parsing for IndelMutationType is unsupported")
    
    def ClassifySNP(self, variant):
        #we want to compare self.snp versus self.detectedsnp. 
        #single mutation is easy, just check for mut. 
        if len(variant.detectedIndel) ==0:
            return "Not Found", None
        other = False
        wildtype = False
        for cardSnp in variant.snp:
            for detectedSnp in variant.detectedIndel:
                if detectedSnp.position == cardSnp.position:
                    if (((len(detectedSnp.mut)-1) % 3) == 0):
                        return "Insertion", [detectedSnp]
                    else:
                        return "Resistant Variant", [detectedSnp]

            for detectedSnp in variant.detectedSnp:
                if detectedSnp.position == cardSnp.position:
                    if detectedSnp.mut == cardSnp.mut:
                        raise Exception ("didnt think this was possible. detectedsnp.mut = fs")
                    elif detectedSnp.mut == cardSnp.wt:
                        return "Wildtype", [detectedSnp]
                    else:
                        return "Other Variant", [detectedSnp]
