def gtf_parser():
    with open("Homo_sapiens.GRCh38.100.gtf", "r") as file:
        lines = file.readlines()

    genes = []
    for i in range(5, len(lines)):
        line = lines[i].split(" ")
        flag = line[0].split("\t")
        if flag[2] == "gene":
            ind = line.index("gene_name")
            line = line[ind+1][1:-2]
            genes.append(line)

    genes = list(set(genes))

    return genes
