from itertools import count
import os
import DataTypes
import ParseCardData

def ParseVcfForVariants(variantCollection, vcfFile):
    vcfFile = [x for x in vcfFile if not x.startswith('#')]

    
    __variantCollection = variantCollection#variantCollection
    cur = 0
    tot = len(__variantCollection)

    output=[]
    output.append("ARO\tVariantType\tMutationType\tClassification\tSNP\tDepth\tAbs_Support\t%_Support\tINFO")
    for variant in __variantCollection:
        #scan the VCF file for this variant and checks if it exists. 
        print(str(cur) + "/" + str(tot))
        #variant = variant.FindIndelInVCF(vcfFile)
        variant = variant.FindVariantInPileup(vcfFile)
        classification, detectedSnp = variant.mutType.ClassifySNP(variant)
        
        if (detectedSnp == None):
            detectedSnp = [] 

        if (len(detectedSnp) > 0):
            if (any(isinstance(i, list) for i in detectedSnp)):
                countOut = ""
                percentOut = ""
                for group in detectedSnp:
                    counts = {"Wildtype":[], "Resistant Variant":[], "Other Variant" :[]}
                    percent = {"Wildtype":[], "Resistant Variant":[], "Other Variant" :[]}         
                    totalDepth = 0
                    if (group):
                        for snp in group:
                            counts["Wildtype"] = [0] * len(detectedSnp)
                            counts["Resistant Variant"] = [0] * len(detectedSnp)
                            counts["Other Variant"] = [0] * len(detectedSnp)


                            totalDepth = int(snp.totalDepth)
                            for i in range(len(variant.snp)):
                                if (snp.mut == variant.snp[i].mut):
                                    counts["Resistant Variant"][i] += snp.depth
                                elif (snp.mut == variant.snp[i].wt):
                                    counts["Wildtype"][i] += snp.depth
                            for i in range(len(variant.snp)):
                                counts["Other Variant"][i] = totalDepth - counts["Resistant Variant"][i] - counts["Wildtype"][i] 

                        for key in counts.keys():
                            percent[key] = [ str(round((int(x)*100)/totalDepth,2)) for x in counts[key]]
                            counts[key] = [str(x) for x in counts[key]]
                    
                if (classification == "Partial"):
                    countOut = countOut + ";".join(counts["Resistant Variant"])
                    percentOut = percentOut + ";".join(percent["Resistant Variant"]) + "%"
                    #countOut = countOut + "[wt:" + ";".join(counts["Wildtype"]) + ",res:" + ";".join(counts["Resistant Variant"]) + ",other:" + ";".join(counts["Other Variant"]) + "]"
                    #percentOut = percentOut +  "[wt:" + ";".join(counts["Wildtype"]) + ",res:" + ";".join(counts["Resistant Variant"]) + ",other:" + ";".join(counts["Other Variant"]) + "]"
                else:
                    countOut = countOut + ";".join(counts[classification]) 
                    percentOut = percentOut + ";".join(counts[classification]) + "%"
                debug = ""

            else:
                counts = {"Wildtype":0, "Resistant Variant":0, "Other Variant" :0}
                percent = {"Wildtype":0, "Resistant Variant":0, "Other Variant" :0}      
                totalDepth = 0
                for snp in detectedSnp:
                    totalDepth = int(snp.totalDepth)
                    if (snp.mut == variant.snp[0].mut):
                        counts["Resistant Variant"] += snp.depth
                    elif (snp.mut == variant.snp[0].wt):
                        counts["Wildtype"] += snp.depth
                counts["Other Variant"] = totalDepth - counts["Resistant Variant"] - counts["Wildtype"]
                for key in counts.keys():
                    percent[key] = str(round(((counts[key] / totalDepth) * 100),2)) + "%"
                    counts[key] = str(counts[key])
                countOut = counts[classification]
                percentOut = percent[classification]
                debug = ";".join([(str(x.mut) + ":" + str(x.depth)) for x in detectedSnp])
        else:
            totalDepth = 0
            countOut = ""
            percentOut = ""
            debug = ""

        out = "\t".join([str(variant.aro.aro), 
                        str(variant.GetType()),
                        str(variant.mutType), 
                        classification, 
                        ";".join([str(x) for x in (variant.snp)]), 
                        str(totalDepth),
                        countOut,#, z
                        percentOut,
                        debug
                        ])
        output.append(out)
        cur+=1

    return output
