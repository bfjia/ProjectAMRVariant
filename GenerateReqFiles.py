import os
from subprocess import Popen, PIPE, TimeoutExpired
import DataTypes
from Bio import SeqIO
from Bio.Alphabet import generic_dna
#samtools mpileup -AQ 0 -Evur gb|AL123456.3|-|2153888-2156111|ARO:3003392|Mycobacterium:1-10 -f ../../card/card/nucleotide_fasta_protein_variant_model.fasta protein.sorted.bam

#fastaPath = "/mnt/f/OneDrive/ProjectAMRSAGE/AMR_Metagenome_Simulator-master/full/card/nucleotide_fasta_protein_variant_model.fasta"
#snpPath = "/mnt/f/OneDrive/ProjectAMRSAGE/AMR_Metagenome_Simulator-master/full/card/snps.tsv"


#samfile = "/mnt/f/OneDrive/ProjectAMRSAGE/AMR_Metagenome_Simulator-master/filter/ChrErrProtein.sam"
#sortedBamFile =  "/mnt/f/OneDrive/ProjectAMRSAGE/AMR_Metagenome_Simulator-master/filter/ChrErrProtein.sorted.bam"
class GenerateReqFiles:
    def __init__(self,  forwardPath, reversePath, OverWrite = False, outputDir = "temp") -> None:
        self.__forwardPath = forwardPath
        self.__reversePath = reversePath
        self.__OverWrite = OverWrite
        self.__outputDir = outputDir

        self.__referenceFastaPath = self.__outputDir + "/ref.fasta"
        self.__referenceGenbankPath = self.__outputDir + "/ref.gb"
        self.__referenceIndex = self.__outputDir + "/refIndex"
        self.__SAMPath = self.__outputDir + "/variants.sam"
        self.__SortedBAMPath = self.__outputDir + "/variants.sorted.bam"
        self.__mpileupPath = self.__outputDir + "/variants.pileup"
        self.__vcfPath = self.__outputDir + "/variants.vcf"

        #self.__pileupPath = self.__outputDir + "/pileups"
        #self.__depthPath = self.__outputDir + "/depths"
        #self.__hitsPath = self.__outputDir + "/hits"

        if not os.path.exists(self.__outputDir):
            os.makedirs(self.__outputDir)
        #if not os.path.exists(self.__pileupPath):
        #    os.makedirs(self.__pileupPath)
        #if not os.path.exists(self.__depthPath):
        #    os.makedirs(self.__depthPath)
        #if not os.path.exists(self.__hitsPath):
        #    os.makedirs(self.__hitsPath)

    def __str__(self) -> str:
        pass

    def __RunCMD(self, cmd):
        print("Running: {}".format(cmd))
        process = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True, universal_newlines=True)
        (output,err) = process.communicate()
        p_status = process.wait()
        print(output)
        print(err)
        return

    def GenerateReferenceFasta(self, variantList):
        if not os.path.isfile(self.__referenceFastaPath) or self.__OverWrite == True:   
            fasta = []
            for variant in variantList:
                fasta.append(">" + str(variant.aro.aro) + "\n" + variant.aro.dna)
            fasta = list(set(fasta))
            with open(self.__referenceFastaPath, "a") as f:
                f.write("\n".join(fasta))
            
            #genbank file
            records = SeqIO.parse(self.__referenceFastaPath, 'fasta')
            with open(self.__referenceFastaPath, "r") as input_handle:
                with open(self.__referenceGenbankPath, 'w') as output_handle:
                    sequences = list(SeqIO.parse(input_handle, "fasta"))
                    for seq in sequences:
                        seq.seq.alphabet = generic_dna
                    SeqIO.write(sequences, output_handle, "genbank")
        else:
            print("using pre-existing reference file at {}".format(self.__referenceFastaPath))

    def Bowtie2Align(self):
        if not os.path.isfile(self.__SAMPath) or self.__OverWrite == True:   
            __CMD_index = "bowtie2-build {} {} --seed 25041".format(self.__referenceFastaPath, self.__referenceIndex)
            self.__RunCMD(__CMD_index)
            __CMD_align = "bowtie2 -x {} -1 {} -2 {} --no-unal --seed 25041 --very-sensitive > {}".format(self.__referenceIndex, self.__forwardPath, self.__reversePath, self.__SAMPath)
            self.__RunCMD(__CMD_align)
        else:
            print("using pre-existing SAM file at {}".format(self.__SAMPath))

    def SAMtoSortedBAM(self):
        if not os.path.isfile(self.__SortedBAMPath) or self.__OverWrite == True:   
            self.__RunCMD("samtools view -bS {} | samtools sort > {}".format(self.__SAMPath, self.__SortedBAMPath))
            self.__RunCMD("samtools index {} ".format(self.__SortedBAMPath))
        else:
            print("using pre-existing sorted BAM file at {}".format(self.__SortedBAMPath))
    
    def SortedBAMtoVCF(self):
        if not os.path.isfile(self.__mpileupPath) or self.__OverWrite == True:
            __CMD_Pileup = "bcftools mpileup -Ov -d 9999 -EAQ 15 -a AD,ADF,ADR,INFO/AD,INFO/ADF,INFO/ADR -f {} {} > {}".format(self.__referenceFastaPath, self.__SortedBAMPath, self.__mpileupPath)
            self.__RunCMD(__CMD_Pileup)
            __CMD_BCFCall = "bcftools call -Ov --ploidy 1 -A -m {}  >  {}".format(self.__mpileupPath, self.__vcfPath)
            self.__RunCMD(__CMD_BCFCall)
        else:
            print("using pre-existing sorted BAM file at {}".format(self.__SortedBAMPath))

    def GenerateVariantFiles(self):    
        if not os.path.isfile(self.__mpileupPath) or self.__OverWrite == True:
            __CMD_Pileup = "samtools mpileup -AE -Q13 -d 99999 -f {} {} > {}".format(self.__referenceFastaPath, self.__SortedBAMPath, self.__mpileupPath)
            self.__RunCMD(__CMD_Pileup)
        
            #for i in range(start, stop+1):
            #    __CMD_Depth = "samtools depth {} | grep '\<{}\>' | grep '\<{}\>' >> {}/{}_{}.depth".format(self.__SortedBAMPath, aro, i, self.__depthPath, aro, snp)
            #    self.__RunCMD(__CMD_Depth)
        
            #__CMD_alignment = "samtools view {} '{}:{}-{}' > {}/{}_{}.reads".format(self.__SortedBAMPath, aro, start, stop, self.__hitsPath, aro, snp)
            #self.__RunCMD(__CMD_alignment)
        
        else:
            print("using pre-existing variant temp files at {}".format(self.__mpileupPath))



