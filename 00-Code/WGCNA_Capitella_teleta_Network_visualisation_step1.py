from random import sample

i = open("CytoscapeInput-nodes-yellow-turquoise-green-red-blue-tan-lightcyan-black-brown-greenyellow-lightyellow-magenta-salmon-grey60-midnightblue-lightgreen-pink-purple-cyan-grey.txt", "r")
o = open("WGCNA_30%_nodes.txt","w")

genes = []

for line in i:
    if line not in genes:
        genes.append(line)

desired_length = round(0.3*len(genes))
print(desired_length)

random_set = sample(genes,desired_length)

for line in random_set:
    o.write(line)

i.close()
o.close()