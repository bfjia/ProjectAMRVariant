"""
Provides underlying datastructure and functions for variant storage and detection

Structs for storing CARD database, SNP information, Variant information, Mutation Types etc.
"""

import itertools
from Bio.Seq import Seq
from numpy.lib.utils import deprecate_with_doc

#structure for storing useful information from CARD.json
class ARO:
    def __init__(self, aro, name, cvterm, species,dna,protein) :
        self.aro = aro
        self.name = name
        self.cvterm = cvterm
        self.species = species
        self.dna = dna
        self.protein = protein

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
                types['del'].append("*")
                # skip to next character
                string = string[1:]
            
            elif string[0] in ['.',',']:
                if (len(string)== 1) or (string[1] not in ['+','-']):
                    # a reference base
                    types[ref] += 1
                    string = string[1:]
                elif string[1] == '+': 
                    insertionLength = int(list(filter(str.isdigit, (string[2:5])))[0]) #hack for double+ digit insertion length
                    insertionSeq = string[3:3+ insertionLength]
                    types['ins'].append(insertionSeq)
                    types['+'] += 1
                    string = string[3+insertionLength:]
                elif string[1] == '-':
                    deletionLength = int(list(filter(str.isdigit, (string[2:5])))[0])
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
        
    def __str__(self) -> str:
        return(self.geneType.capitalize() + ":" + self.mutType.capitalize())

#child of Variant class specific for RNA variants
class RNAVariant(Variant):
    def __init__(self, aro, snp, mutType):
        Variant.__init__(self, aro, snp, mutType, "rna variant")


    def FindVariantInPileup(self, pileup):
        detectedSNPs = []
        detectedIndels = []
        pileup = [x for x in pileup if x.startswith("{}\t".format(str(self.aro.aro)))]

        for snp in self.snp:
            relevantSNP = [x for x in pileup if x.startswith("{}\t{}\t".format(str(self.aro.aro), str(int(snp.position))))]
            s = self.ParsePileupForSNPs(relevantSNP, snp)
            if (s.IsIndel):
                detectedIndels.extend(s)
            else:
                detectedSNPs.extend(s)

        self.SetDetectedSnp(detectedSNPs)
        self.SetDetectedIndel(detectedIndels)
        return self
    """ 
    def FindVariantInVCF(self, vcf):
        detectSNPs = []
        relevantVCF = [x for x in vcf if x.startswith("{}\t".format(str(self.aro.aro)))]
        SNPs = [x for x in relevantVCF if not 'INDEL' in x]
        for snp in self.snp:
            ##CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	temp/variants.sorted.bam
            relevantSNP = [x for x in SNPs if x.startswith("{}\t{}\t".format(str(self.aro.aro), str(int(snp.position))))]
            
            if(len(relevantSNP) < 1):
                exist = False
            else:
                for s in relevantSNP:
                    data = self.SplitVCFLine(relevantSNP)
                    wt = snp.wt
                    if (data[4] == "."):
                        mut = data[3]
                    else:
                        mut = data[4]
                    quality = float(data[5])
                    info = data[7].replace(".","0").split(";")
                    dp = int([x for x in info if x.startswith("DP=")][0].replace("DP=",""))
                    dp4 = list(map(int, list(([x for x in info if x.startswith("DP4=")][0].replace("DP4=","").split(",")))))
                    mq = int([x for x in info if x.startswith("MQ=")][0].replace("MQ=",""))
                    detectSNPs.append(SNP(wt, mut, snp.position, quality, dp, dp4, mq))
        self.detectedSnp = detectSNPs
        return self
        """
#Child of Variant class specific for Protein Variants
class ProteinVariant(Variant):
    def __init__(self, aro, snp, mutType, genetype = "protein variant"):
        Variant.__init__(self, aro, snp, mutType, genetype)
    
    def FindVariantInPileup(self, pileup):
        detectedSNPs = []
        detectedIndels = []
        pileup = [x for x in pileup if x.startswith("{}\t".format(str(self.aro.aro)))]

        for snp in self.snp:
            base1 = [x for x in pileup if x.startswith("{}\t{}\t".format(str(self.aro.aro), str(int(snp.position)*3-2)))]
            base2 = [x for x in pileup if x.startswith("{}\t{}\t".format(str(self.aro.aro), str(int(snp.position)*3-1)))]
            base3 = [x for x in pileup if x.startswith("{}\t{}\t".format(str(self.aro.aro), str(int(snp.position)*3)))]
            
            base1Snp = {}
            base2Snp = {}
            base3Snp = {}
            base1Indel = []
            base2Indel = []
            base3Indel = []

            for b in self.ParsePileupForSNPs(base1, snp):
                if not b.IsIndel():
                    base1Snp[b.mut] = b
                else:
                    base1Indel.append(b)
            for b in self.ParsePileupForSNPs(base2, snp):
                if not b.IsIndel():
                    base2Snp[b.mut] = b
                else:
                    base2Indel.append(b)
            for b in self.ParsePileupForSNPs(base3, snp):
                if not b.IsIndel():
                    base3Snp[b.mut] = b
                else:
                    base3Indel.append(b)

            detectedIndels.append([base1Indel, base2Indel, base3Indel])

            combinations = list(itertools.product(*[base1Snp.keys(), base2Snp.keys(), base3Snp.keys()]))
            for combo in combinations:
                wt = snp.wt
                codon = "".join(list(combo))
                mut = str(Seq(codon).translate())
                lowestDepth = int((base1Snp[codon[0]].depth+ base2Snp[codon[1]].depth+ base3Snp[codon[2]].depth)/3)
                maxTotalDepth = max (base1Snp[codon[0]].totalDepth, base2Snp[codon[1]].totalDepth, base3Snp[codon[2]].totalDepth)
                detectedSNPs.append(SNP(wt, mut, snp.position, lowestDepth, maxTotalDepth, codon))

        self.SetDetectedSnp(detectedSNPs)
        self.SetDetectedIndel(detectedIndels)
        return self
    """ 
    def FindVariantInVCF(self, vcf):
        detectSNPs = []
        relevantVCF = [x for x in vcf if x.startswith("{}\t".format(str(self.aro.aro)))]
        
        def ExtractDataFromVCFLine(vcfLine):
            data = self.SplitVCFLine(vcfLine)
            wt = None
            if (data[4] == "."):
                mut = data[3]
            else:
                mut = data[4]
            quality = data[5]
            info = data[7].replace(".","0").split(";")

            dp = int([x for x in info if x.startswith("DP=")][0].replace("DP=",""))
            dp4 = list(map(int, list(([x for x in info if x.startswith("DP4=")][0].replace("DP4=","").split(",")))))
            mq = int([x for x in info if x.startswith("MQ=")][0].replace("MQ=",""))
            exist = True
            return mut, SNP(wt, mut, data[1], quality, dp, dp4, mq)

        SNPs = [x for x in relevantVCF if not 'INDEL' in x]
        for snp in self.snp:
            ##CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	temp/variants.sorted.bam
            base1 = [x for x in SNPs if x.startswith("{}\t{}\t".format(str(self.aro.aro), str(int(snp.position)*3-2)))]
            base2 = [x for x in SNPs if x.startswith("{}\t{}\t".format(str(self.aro.aro), str(int(snp.position)*3-1)))]
            base3 = [x for x in SNPs if x.startswith("{}\t{}\t".format(str(self.aro.aro), str(int(snp.position)*3)))]

            base1Snp = {}
            base2Snp = {}
            base3Snp = {}

            if(len(base1) < 1 or len(base2) < 1 or len(base3) < 1 ):
                exist = False
            else:
                
                for b in base1:
                    m,s = ExtractDataFromVCFLine(b)
                    if (m.find(",") > -1):
                        for mBase in m.split(","):
                            base1Snp[mBase] = SNP(s.wt, s.mut, s.position, s.quality, s.DP, s.DP4, s.MQ)
                    else:
                        base1Snp[m] = s
                for b in base2:
                    m,s = ExtractDataFromVCFLine(b)
                    if (m.find(",") > -1):
                        for mBase in m.split(","):
                            base2Snp[mBase] = SNP(s.wt, s.mut, s.position, s.quality, s.DP, s.DP4, s.MQ)
                    else:
                        base2Snp[m] = s                
                for b in base3:
                    m,s = ExtractDataFromVCFLine(b)
                    if (m.find(",") > -1):
                        for mBase in m.split(","):
                            base3Snp[mBase] = SNP(s.wt, s.mut, s.position, s.quality, s.DP, s.DP4, s.MQ)
                    else:
                        base3Snp[m] = s
                combinations = list(itertools.product(*[list(base1Snp.keys()), list(base2Snp.keys()), list(base3Snp.keys())]))
                for combo in combinations:
                    wt = snp.wt
                    codon = "".join(list(combo))
                    mut = str(Seq(codon).translate())

                    quality = [base1Snp[codon[0]].quality, base2Snp[codon[1]].quality, base3Snp[codon[2]].quality]
                    dp = [base1Snp[codon[0]].DP, base2Snp[codon[1]].DP, base3Snp[codon[2]].DP]
                    mq = [(base1Snp[codon[0]].MQ, base2Snp[codon[1]].MQ, base3Snp[codon[2]].MQ)]

                    dp4 = [(base1Snp[codon[0]].DP4), ( base2Snp[codon[1]].DP4), (base3Snp[codon[2]].DP4)]

                    exist = True
                    detectSNPs.append(SNP(wt, mut, snp.position, quality, dp, dp4, mq, codon = codon))
        
            self.__SetDetectedIndel() = detectSNPs
        return self
    """

#child of ProteinVariant class specific for Protein overexpression
class ProteinOverexpressionVariant(ProteinVariant):
    def __init__(self, aro, snp, mutType):
        ProteinVariant.__init__(self, aro, snp, mutType, "protein overexpression")



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
            if(resistant):
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
                return "???", None

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
            other = False
            wildtype = False

            for indel in variant.detectedIndel:
                for i in range(3):
                    if indel[i]:
                        for detectedSnp in indel[i]:
                            if detectedSnp.position == cardSnp.position:
                                for sequence in detectedSnp.sequences:
                                    if ((len(sequence) % 3) != 0):
                                        return "Resistant Variant", indel

            return "Wildtype", []


    
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
