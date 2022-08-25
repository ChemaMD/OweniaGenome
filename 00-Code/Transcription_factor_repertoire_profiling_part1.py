if __name__ == "__main__":

    import re
    specieslist = ["Owenia_fusiformis","Capitella_teleta","Dimorphilus_gyrociliatus"]

    # First generate the PFAM universe files for all species, might come in handy later on for background comparisons

    for species in specieslist:

        i = open("02-%s_annotation.txt"%(species), "r")
        o = open("03-%s_PFAM_universe.txt"%(species), "w")

        regex = re.compile(r'PF\d+')

        for line in i:
            pfam_matches = regex.findall(line)
            if species == 'Dimorphilus_gyrociliatus':
                Gene_ID = line.split(" ",1)[0]
            elif species == 'Capitella_teleta':
                Gene_ID = line.split("\t")[1]
            else:
                Gene_ID = line.split("\t",1)[0]
            if not pfam_matches == []:
                o.write(Gene_ID+"\t")
                for i, match in enumerate(pfam_matches):
                    if i+1 == len(pfam_matches):
                        o.write(match.strip("'")+"\n")
                    else:
                        o.write(match.strip("'")+", ")
        
        o.close()

    
    # Then subset the genes with a PFAM domain corresponding to a TF DNA-binding domain and generate stats files

    tf_input = open("01-TF_PFAM_domains.txt", "r")
    pfam_tf_list = []
    pfam_tf_dictionary = {}

    for i, line in enumerate(tf_input):
        if not i == 0:
            pfam_tf_list.append(line.split('\t')[4])
            pfam_tf_dictionary[line.split('\t')[4]] = line.split('\t')[3]
        

    for species in specieslist:

        i1 = open("02-%s_annotation.txt"%(species), "r")

        o1 = open("04-%s_TF_list_annotated.txt"%(species),"w")
        o2 = open("05-%s_TF_list_gene_ID_only.txt"%(species),"w")
        o3 = open("06-%s_TF_list_classified.txt"%(species),"w")

        # Subset TFs from the annotation file (o1)
        # Subset TFs from the annotation file, keeping the geneID only (o3)
        # Subset TFs from the annotation file, keeping the geneID, the PFAM entries that match a TF PFAM_ID, and the classification (o2)

        for line in i1:
            pfam_matches = regex.findall(line)
            if any(pfam_match in pfam_tf_list for pfam_match in pfam_matches):
                o1.write(line)
                copy_dict = pfam_tf_dictionary.copy()
                for key in copy_dict:
                    copy_dict[key] = 0
                for pfam_match in pfam_matches:
                    if pfam_match in pfam_tf_list:
                        copy_dict[pfam_match] = 1
                tf_class_list_for_gene = []
                pfam_id_list_for_gene = []
                for pfam_id in copy_dict:
                    if copy_dict[pfam_id] == 1:
                        tf_class_list_for_gene.append(pfam_tf_dictionary[pfam_id])
                        pfam_id_list_for_gene.append(pfam_id)
                if not tf_class_list_for_gene == []:
                    if not pfam_id_list_for_gene == []:
                        if species == 'Dimorphilus_gyrociliatus':
                            o2.write(line.split(' ')[0]+'\n')
                            o3.write(line.split(' ')[0]+'\t'+','.join(pfam_id_list_for_gene)+'\t'+','.join(tf_class_list_for_gene)+'\n')
                        elif species == 'Capitella_teleta':
                            o2.write(line.split('\t')[1]+'\n')
                            o3.write(line.split('\t')[1]+'\t'+','.join(pfam_id_list_for_gene)+'\t'+','.join(tf_class_list_for_gene)+'\n')
                        elif species == 'Owenia_fusiformis':
                            o2.write(line.split('\t')[0]+'\n')
                            o3.write(line.split('\t')[0]+'\t'+','.join(pfam_id_list_for_gene)+'\t'+','.join(tf_class_list_for_gene)+'\n')
                    else:
                        raise NameError("A gene is trying to be annotated as a TF but does not have a PFAM domain belonging to the TF classes defined in the input file")
                else:
                    raise NameError("A gene is trying to be annotated as a TF but does not have a PFAM domain belonging to the TF classes defined in the input file")

        o1.close()
        o2.close()
        o3.close()

        # Generate stats files with # of TFs by class (o4)

        i2 = open("06-%s_TF_list_classified.txt"%(species),"r")
        o4 = open("07-%s_TF_stats.txt"%(species),"w")

        count_dict = pfam_tf_dictionary.copy()
        genes_dict = pfam_tf_dictionary.copy()

        for key in count_dict.keys():
            count_dict[key] = 0
            genes_dict[key] = []

        for line in i2:
            pfam_matches = regex.findall(line)
            for pfam_match in pfam_matches:
                count_dict[pfam_match] += 1
                genes_dict[pfam_match].append(line.split('\t')[0])
        
        for key,value in count_dict.items():
            o4.write(key+'\t'+pfam_tf_dictionary[key]+'\t'+str(value)+'\t'+','.join(genes_dict[key])+'\n')

        o4.close()
        

