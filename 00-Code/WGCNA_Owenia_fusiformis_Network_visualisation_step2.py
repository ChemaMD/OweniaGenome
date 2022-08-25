i1 = open("WGCNA_30%_nodes.txt", "r")
i2 = open("CytoscapeInput-edges-black-red-brown-turquoise-blue-yellow-green-greenyellow-cyan-pink-midnightblue-magenta-grey-purple-tan-salmon.txt", "r")
o = open("WGCNA_30%_edges.txt","w")

genes = []
genes.append('fromNode')
genes.append('toNode')

for line in i1:
    if line.split('\t')[0] not in genes:
        genes.append(line.split('\t')[0])


for line in i2:
    if line.split('\t')[0] in genes:
        if line.split('\t')[1] in genes:
            o.write(line)

i1.close()
i2.close()
o.close()