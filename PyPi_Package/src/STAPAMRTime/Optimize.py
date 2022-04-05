import ParseCardData
import GenerateReqFiles
import argparse
import GenerateSyntheticData


def Optimize():
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
    parser.add_argument(
        "-t",
        "--threads",
        default=16,
        type=int,
        help="output path for temp files",
    )
    parser.add_argument("--overwrite", action="store_true", help="overwrite existing files.")
    args = parser.parse_args()

    
    cardJsonPath = args.card_json #"/mnt/f/OneDrive/ProjectAMRSAGE/AMR_Metagenome_Simulator-master/full/card/card.json"
    cardSnpPath = args.card_snp #/mnt/f/OneDrive/ProjectAMRSAGE/AMR_Metagenome_Simulator-master/full/card/snps.txt"
    forward = args.forward
    reverse = args.reverse
    
    forward = "/mnt/h/AMRVariant_Genomes/Synthetic/synthetic_1.fq.gz"
    reverse = "/mnt/h/AMRVariant_Genomes/Synthetic/synthetic_2.fq.gz"

    CARD = ParseCardData.ParseJson(cardJsonPath) 
    proteinVariants, rnaVariants = ParseCardData.ParseSNP(cardSnpPath, CARD)
    variantsCollection = (rnaVariants) + proteinVariants
    #bowtie2 step.
    outputProtein = []
    outputRNA = []
    for mapq in range (32):
        generator = GenerateReqFiles.GenerateReqFiles(forward, reverse, outputDir=args.temp + str(mapq))
        generator.GenerateReferenceFasta(variantsCollection)
        generator.Bowtie2Align(args.threads)
        generator.SAMtoSortedBAM()
        generator.GenerateVariantFiles(quality=mapq)
        generator.GenerateUnfilteredVariantSummary(variantsCollection)
        
        #shutil.copy2(args.temp + str(mapq) + "/results.tsv", args.temp + str(mapq) + "_DetectedVariants.tsv")
        #optimize for protein
        for a in range(10):
            for r in range(0, 500, 2):
                with open (args.temp + str(mapq) + "/results.tsv", 'r') as f:
                    unfilteredResult = f.readlines()

                filterAbsolute = a
                filterRelative = float(r/10)
                for i in range(len(unfilteredResult)):
                    line = unfilteredResult[i].strip().split("\t")
                    if line[3] == "Resistant Variant":
                        #apply logic here
                        abs = min([float(x) for x in line[6].strip().split(";")])
                        rel = min([float(x) for x in line[7].strip().replace("%","").split(";")])

                        if (abs < filterAbsolute or rel < filterRelative):
                            line[3] = "False"
                        else:
                            line[3] = "True"
                    else:
                        line[3] = "False"
                    unfilteredResult[i] = "\t".join(line)

                matrix, totalRes, total = GenerateSyntheticData.CheckAccuracy(unfilteredResult, "./sim/SyntheticProteinVariants.fasta")
                metrics, header = matrix.GetAllMetrics()
                metrics = [str(round(x,4)) for x in metrics]
                
                outputProtein.append("{}\t{}\t{}\t{}\n".format(str(mapq), str(filterAbsolute), str(filterRelative), "\t".join(metrics)))
                print("MAPQ: {}; param: {},{}%".format(str(mapq), str(filterAbsolute), str(filterRelative)))
                print("Accuracy: {}; TPR: {}; FPR: {}".format(str(matrix.GetAccuracy()), str(matrix.GetTPR()), str(matrix.GetFPR())))


        #optimize for rna
        for a in range(10):
            for r in range(0, 500, 2):
                with open (args.temp + str(mapq) + "/results.tsv", 'r') as f:
                    unfilteredResult = f.readlines()

                filterAbsolute = a
                filterRelative = float(r/10)
                for i in range(len(unfilteredResult)):
                    line = unfilteredResult[i].strip().split("\t")
                    if line[3] == "Resistant Variant":
                        #apply logic here
                        abs = min([float(x) for x in line[6].strip().split(";")])
                        rel = min([float(x) for x in line[7].strip().replace("%","").split(";")])

                        if (abs < filterAbsolute or rel < filterRelative):
                            line[3] = "False"
                        else:
                            line[3] = "True"
                    else:
                        line[3] = "False"
                    unfilteredResult[i] = "\t".join(line)
                matrix, totalRes, total = GenerateSyntheticData.CheckAccuracy(unfilteredResult, "./sim/SyntheticRRnaVariants.fasta")
                metrics, header = matrix.GetAllMetrics()
                metrics = [str(round(x,4)) for x in metrics]

                outputRNA.append("{}\t{}\t{}\t{}\n".format(str(mapq), str(filterAbsolute), str(filterRelative),  "\t".join(metrics)))
                print("MAPQ: {}; param: {},{}%".format(str(mapq), str(filterAbsolute), str(filterRelative)))
                print("Accuracy: {}; TPR: {}; FPR: {}".format(str(matrix.GetAccuracy()), str(matrix.GetTPR()), str(matrix.GetFPR())))
    with open ("roc_protein.tsv", 'w') as f:
        f.write("MAPQ\tFilter_Abs\tFilter_Rel\t" + "\t".join(header) + "\n")
        f.writelines(outputProtein)
    with open ("roc_rna.tsv", 'w') as f:
        f.write("MAPQ\tFilter_Abs\tFilter_Rel\t" + "\t".join(header) + "\n")
        f.writelines(outputRNA)
Optimize()


