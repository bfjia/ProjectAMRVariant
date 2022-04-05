from urllib.request import urlretrieve
from . import DataTypes
import json
import tarfile
import os

def GetCardJsonFromWeb(version = "latest", downloadURL = None):
    print("Downloading CARD data...")
    if version == "latest":
        url = "https://card.mcmaster.ca/latest/data"
    else:
        url = "https://card.mcmaster.ca/download/0/broadstreet-v{}.tar.bz2".format(version)
    
    if (downloadURL != None):
        url = downloadURL

    try:
        urlretrieve(url, "card.tar.bz2")
    except Exception as e:
        print ("Error download card.json from {}: {}".format(url, str(e)))
        raise

    tar = tarfile.open("card.tar.bz2", "r:bz2")  
    tar.extractall("localDB")
    tar.close()

    os.remove("card.tar.bz2")
    print("Done.")
    


def ParseJson(jsonPath, version = "latest"): #discards metamodels
    #parseAROJson
    print("Parsing card.json...")
    cardJsonPath = jsonPath

    if not os.path.isfile(cardJsonPath):
        raise Exception("No card.json found at path specified")

    with open(cardJsonPath, 'r') as f:
        cardJson = json.loads(f.read())
    
    CARD = {}
    for model_id in (cardJson):
        try:
            accession = int(cardJson[model_id]["ARO_accession"])
            name = cardJson[model_id]["ARO_name"]
            cvterm = int(cardJson[model_id]["ARO_id"])
            seq = list(cardJson[model_id]["model_sequences"]["sequence"].values())[0]
            species = seq["NCBI_taxonomy"]["NCBI_taxonomy_name"]
            dna = seq["dna_sequence"]["sequence"]
            protein = seq["protein_sequence"]["sequence"]


            CARD[accession] = DataTypes.ARO(accession, name, cvterm, species, dna, protein)
        except Exception as e:
            try:
                print("Missing Sequences in CARD.json: ARO" + str(int(cardJson[model_id]["ARO_accession"])) )#+ str(e))
            except Exception:
                print(str(model_id) + str(e))
    return CARD

def ParseSNP(snpPath, CARD):
    print("Parsing available SNPs...")
    #lets parse a list of snps
    #snpPath = "/mnt/f/OneDrive/ProjectAMRSAGE/AMR_Metagenome_Simulator-master/full/card/snps.txt"
    with open(snpPath, 'r') as f:
        snpFile = f.readlines()

    proteinVariants = []
    rnaVariants = []
    for line in snpFile[1:]:
        try: 
            l = line.strip().split("\t")
            if (int(l[0]) in CARD):
                aro = CARD[int(l[0])]
            else:
                raise Exception("ARO{}: SNP exist in snps.txt but Data does not exist in CARD.json".format(str(l[0])))

            mt = l[3]
            if (mt == "single resistance variant"):
                mutationType = DataTypes.SingleMutationType()
            elif(mt == "multiple resistance variants"):
                mutationType = DataTypes.MultipleMutationType()
            elif(mt == "nonsense mutation"):
                mutationType = DataTypes.NonSenseMutationType()
            elif(mt == "frameshift mutation"):
                mutationType = DataTypes.FrameshiftMutationType()
            elif(mt == "co-dependent single resistance variant" or mutationType == "co-dependent insertion/deletion"):
                mutationType = DataTypes.CodependentMutationType()
            elif(mt.find("TB") > -1):
                if (l[4][0]== "+"):
                    mutationType = DataTypes.IndelMutationType()
                elif(l[4][0] == "-"):
                    mutationType = DataTypes.IndelMutationType()
                elif (l[4].find("\-") > -1):
                    mutationType = DataTypes.MutationType()
                elif (l[4][-2:].strip().lower() == "fs"):
                    mutationType = DataTypes.FrameshiftMutationType()
                elif (l[4][-4:].strip().lower() == "stop"):
                    mutationType = DataTypes.NonSenseMutationType()
                elif(l[4].find(",") > -1):
                    mutationType = DataTypes.MultipleMutationType()
                else:
                    mutationType = DataTypes.SingleMutationType()

            success,snp = mutationType.ParseSNP(str(aro.aro),l[4])
            if (success):
                if (l[2].strip().lower() == "protein variant model"):
                    variant = DataTypes.ProteinVariant(aro, snp, mutationType)
                    proteinVariants.append(variant)
                elif (l[2].strip().lower() == "protein overexpression model"):
                    variant = DataTypes.ProteinOverexpressionVariant(aro, snp, mutationType)
                    proteinVariants.append(variant)
                elif (l[2].strip().lower() == "rrna gene variant model"):
                    variant = DataTypes.RNAVariant(aro, snp, mutationType)
                    rnaVariants.append(variant)
                else:
                    raise Exception("Unsupported Model Type, got {}".format(l[2]))
            else:
                raise Exception("Could not parse the SNP, {}".format(snp))
            
        except Exception as e:
            #print(str(e))
            with open("SNPParsingErrors.txt", "a") as f:
                f.write(line.strip() + "\t" + str(e) + "\n")
            continue
    return proteinVariants, rnaVariants
