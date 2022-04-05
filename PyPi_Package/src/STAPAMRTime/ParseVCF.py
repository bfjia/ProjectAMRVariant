def ParseVcfForVariants(variantCollection, vcfFile):
    vcfFile = [x for x in vcfFile if not x.startswith('#')]

    
    __variantCollection = variantCollection#variantCollection
    cur = 0
    tot = len(__variantCollection)

    output=[]
    output.append("ARO\tVariantType\tMutationType\tResistantVariant\tSNP\tDepth\tAbs_Support\t%_Support\tINFO")
    for variant in __variantCollection:
        #scan the VCF file for this variant and checks if it exists. 
        if (cur%250==0):
            print(str(cur) + "/" + str(tot))
        
        #variant = variant.FindIndelInVCF(vcfFile)
        variant = variant.FindVariantInPileup(vcfFile)
        classification, detectedSnp = variant.mutType.ClassifySNP(variant)
        debug = ""

        if (detectedSnp == None):
            detectedSnp = [] 
#        if (variant.aro.aro == 3003394 and str(variant.mutType) == "Frameshift" and classification == "Resistant Variant"):
#            print("this")
        if (len(detectedSnp) > 0):
            if (any(isinstance(i, list) for i in detectedSnp)):
                countOut = ""
                percentOut = ""
                counts = {"Wildtype":[], "Resistant Variant":[], "Other Variant" :[]}
                percent = {"Wildtype":[], "Resistant Variant":[], "Other Variant" :[]}      
                counts["Wildtype"] = [0] * len(detectedSnp)
                counts["Resistant Variant"] = [0] * len(detectedSnp)
                counts["Other Variant"] = [0] * len(detectedSnp)
                idx = 0
                totalDepth = 1
                for group in detectedSnp: #group represents 1 of the variants in a multi variant situation
                    if (group):
                        for snp in group: #codons in a variant
                            totalDepth = int(snp.totalDepth)
                            if (totalDepth == 0):
                                totalDepth = 1
                            for i in range(len(variant.snp)):
                                if (snp.position == variant.snp[i].position):
                                    if (snp.mut == "Indel"):
                                        counts["Resistant Variant"][idx] += snp.depth
                                        debug = debug + snp.type + ":" + ";".join(list(set(snp.sequences))) + ";"
                                    else:
                                        if (snp.mut == variant.snp[i].mut):
                                            counts["Resistant Variant"][idx] += snp.depth
                                        elif (snp.mut == variant.snp[i].wt):
                                            counts["Wildtype"][idx] += snp.depth
                                        else:
                                            counts["Other Variant"][idx] += snp.depth
                                        debug = debug + snp.mut + ":" + str(snp.depth) + ";"
                        #for i in range(len(variant.snp)):
                         #   counts["Other Variant"][i] = totalDepth - counts["Resistant Variant"][i] - counts["Wildtype"][i] 
                    idx += 1

                debug = debug[:-1] + ","
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
                    percentOut = percentOut + ";".join(percent[classification]) + "%"

            else:
                counts = {"Wildtype":0, "Resistant Variant":0, "Other Variant" :0}
                percent = {"Wildtype":0, "Resistant Variant":0, "Other Variant" :0}      
                totalDepth = 0
                for snp in detectedSnp:
                    totalDepth = int(snp.totalDepth)
                    if (totalDepth == 0):
                        totalDepth= 1

                        
                    if (snp.position == variant.snp[0].position):
                        if (snp.mut == "Indel"):
                            counts["Resistant Variant"] += snp.depth
                            debug = debug + snp.type + ":" + ";".join(list(set(snp.sequences))) + ";"
                        else:
                            if (snp.mut == variant.snp[0].mut):
                                counts["Resistant Variant"] += snp.depth #if there are multiple instances of a codon, the sum might be bigger than the depth
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
