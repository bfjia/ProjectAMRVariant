import ParseCardData
import GenerateReqFiles
import ParseVCF
import argparse



def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-f",
        "--forward",
        type=str,
        help="forward fastq",
    )
    parser.add_argument(
        "-r",
        "--reverse",
        type=str,
        help="reverse fastq path",
    )
    parser.add_argument(
        "--card-snp",
        default="./snps.txt",
        type=str,
        help="path to card snp",
    )
    parser.add_argument(
        "--card-json",
        default="./card.json",
        type=str,
        help="path to card json",
    )
    parser.add_argument(
        "--temp",
        default="./temp",
        type=str,
        help="output path for temp files",
    )

    args = parser.parse_args()

    
    cardJsonPath = args.card_json #"/mnt/f/OneDrive/ProjectAMRSAGE/AMR_Metagenome_Simulator-master/full/card/card.json"
    cardSnpPath = args.card_snp #/mnt/f/OneDrive/ProjectAMRSAGE/AMR_Metagenome_Simulator-master/full/card/snps.txt"
    forward = args.forward
    reverse = args.reverse
    
    #forward = "/mnt/f/OneDrive/ProjectAMRSAGE/synthetic_1.fq.gz"
    #reverse = "/mnt/f/OneDrive/ProjectAMRSAGE/synthetic_2.fq.gz"

    CARD = ParseCardData.ParseJson(cardJsonPath) 
    proteinVariants, rnaVariants = ParseCardData.ParseSNP(cardSnpPath, CARD)
    variantsCollection = proteinVariants + (rnaVariants)
    #bowtie2 step.

    generator = GenerateReqFiles.GenerateReqFiles(forward, reverse, outputDir=args.temp)

    generator.GenerateReferenceFasta(proteinVariants)
    generator.Bowtie2Align()
    generator.SAMtoSortedBAM()
    generator.GenerateVariantFiles()


    vcfPath = "{}/variants.pileup".format(args.temp)
    with open(vcfPath, "r") as f:
        vcfFile = f.readlines()

    output = ParseVCF.ParseVcfForVariants(variantsCollection, vcfFile)
    with open ("./variants.tsv", 'w') as f:
        f.write("\n".join(output))

main()