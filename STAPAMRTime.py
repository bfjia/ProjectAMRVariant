from . import ParseCardData
from . import GenerateReqFiles
import argparse

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-f",
        "--forward",
        type=str,
        default=None,
        help="forward fastq",
    )
    parser.add_argument(
        "-r",
        "--reverse",
        type=str,
        default=None,
        help="reverse fastq path",
    )
    parser.add_argument(
        "-i",
        "--interleaved",
        type=str,
        default = None,
        help="interleaved fastq path",
    )    
    parser.add_argument(
        "-u",
        "--unpaired",
        type=str,
        default = None,
        help="interleaved fastq path",
    )
    parser.add_argument(
        "--card-snp",
        default=None,
        type=str,
        help="path to card snp",
    )
    parser.add_argument(
        "--card-json",
        default=None,
        type=str,
        help="path to card json",
    )
    parser.add_argument(
        "--card-version",
        default=None,
        type=str,
        help="The URL or version string of CARD database to download.",
    )
    parser.add_argument(
        "--temp",
        default="./temp",
        type=str,
        help="output path for temp files",
    )
    parser.add_argument(
        "--protein-filter",
        default=6,
        type=float,
        help="filter used for detecting protein variants",
    )
    parser.add_argument(
        "--rna-filter",
        default=6,
        type=float,
        help="filter used for detecting protein variants",
    )
    parser.add_argument(
        "--mapping-quality-filter",
        default=13,
        type=int,
        help="samtools filter used for finding variants",
    )
    parser.add_argument(
        "-t",
        "--threads",
        default=8,
        type=int,
        help="output path for temp files",
    )
    parser.add_argument("--overwrite", action="store_true", help="overwrite existing files.")
    parser.add_argument("--max-precision", action="store_true", help="Parameters used to maximize precision")
    parser.add_argument("--max-recall", action="store_true", help="Parameters used to maximize recall.")

    args = parser.parse_args()
    
    if (args.forward != None and args.reverse!=None):
        forward = args.forward
        reverse = args.reverse
    elif (args.interleaved != None):
        forward = args.interleaved
        reverse = None       
    elif (args.unpaired != None):
        reverse = args.unpaired
        forward = None
    else:
        print ("FASTQ filepath error")
        return
    
    downloadCARD = False
    if (args.card_json == None):
        cardJsonPath = "./localDB/card.json"
        downloadCARD = True
    else:
        cardJsonPath = args.card_json #"/mnt/f/OneDrive/ProjectAMRSAGE/AMR_Metagenome_Simulator-master/full/card/card.json"
    
    if (args.card_snp == None):
        cardSnpPath = "./localDB/snps.txt"
        downloadCARD = True
    else:
        cardSnpPath = args.card_snp #"/mnt/f/OneDrive/ProjectAMRSAGE/AMR_Metagenome_Simulator-master/full/card/card.json"

    if (args.max_precision):
        filterDict = {"filterMAPQProtein":31, 
                        "filterAbsoluteProtein":6, 
                        "filterRelativeProtein":0,
                        "filterMAPQRNA": 13,
                         "filterAbsoluteRNA":0, 
                         "filterRelativeRNA": 1}
    elif (args.max_recall):
        filterDict = {"filterMAPQProtein":0, 
                    "filterAbsoluteProtein":0, 
                    "filterRelativeProtein":0,
                    "filterMAPQRNA": 19,
                    "filterAbsoluteRNA":0, 
                    "filterRelativeRNA":0}
    else:
        filterDict = {"filterMAPQProtein":5, 
            "filterAbsoluteProtein":0, 
            "filterRelativeProtein":0,
            "filterMAPQRNA": 19,
            "filterAbsoluteRNA":0, 
            "filterRelativeRNA":0}

    if (downloadCARD):
        if (args.card_version == None):
            ParseCardData.GetCardJsonFromWeb()
        else:
            if args.card_version.find("card.mcmaster.ca"):
                ParseCardData.GetCardJsonFromWeb(downloadURL=args.card_version)
            else:
                ParseCardData.GetCardJsonFromWeb(version=args.card_version)

  
    CARD = ParseCardData.ParseJson(cardJsonPath) 
    proteinVariants, rnaVariants = ParseCardData.ParseSNP(cardSnpPath, CARD)
    variantsCollection = (rnaVariants) + proteinVariants
    #bowtie2 step.

    generator = GenerateReqFiles.GenerateReqFiles(forward, reverse, outputDir=args.temp)

    generator.GenerateReferenceFasta(variantsCollection)
    generator.Bowtie2Align(args.threads)
    generator.SAMtoSortedBAM()
    generator.GenerateVariantFiles(quality=filterDict["filterMAPQProtein"])
    generator.GenerateUnfilteredVariantSummary(variantsCollection)
    
    with open (args.temp + "/results.tsv", 'r') as f:
        unfilteredResult = f.readlines()

    def FilterProteinVariants(unfilteredResult, filterDict):
        for i in range(len(unfilteredResult)):
            line = unfilteredResult[i].strip().split("\t")
            if (line[1].find("rna") > -1):
                if line[3] == "Resistant Variant":
                    #apply logic here
                    abs = min([float(x) for x in line[6].strip().split(";")])
                    rel = min([float(x) for x in line[7].strip().replace("%","").split(";")])

                    if (abs < filterDict["filterAbsoluteRNA"] or rel < filterDict["filterRelativeRNA"] ):
                        line[3] = "False"
                    else:
                        line[3] = "True"
                else:
                    line[3] = "False"
            elif (line[1].find("protein") > -1):
                if line[3] == "Resistant Variant":
                    #apply logic here
                    abs = min([float(x) for x in line[6].strip().split(";")])
                    rel = min([float(x) for x in line[7].strip().replace("%","").split(";")])

                    if (abs < filterDict["filterAbsoluteProtein"] or rel < filterDict["filterRelativeProtein"] ):
                        line[3] = "False"
                    else:
                        line[3] = "True"
                else:
                    line[3] = "False"
            unfilteredResult[i] = "\t".join(line)
        unfilteredResult = [x for x in unfilteredResult if x.find("True") > -1]
        with open (args.temp + "_DetectedVariants.tsv", 'w') as f:
            f.write("ARO\tVariantClass\tVariantType\tResistantVariant\tSNP\tDepth\tAbsSupport\tRelativeSupport\tINFO\n")
            f.write("\n".join(unfilteredResult))

    FilterProteinVariants(unfilteredResult, filterDict)

if __name__ == "__main__":
    main()



