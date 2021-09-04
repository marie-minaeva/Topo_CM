import pandas as pd


def signature_extractor(file):
    """
    Function is used to calculate signatures out of differentially expressed genes;
    :param file: a name of file with DE genes;
    :return: two lists: one with up-regulated  and the other with down-regulated genes;
    """

    # preprocessing data before
    # calculating up/down regulated genes

    stat = pd.read_table(file, sep=" ")
    stat = stat.apply(pd.to_numeric)

    up_genes = []
    down_genes = []
    names = stat.columns
    stat = stat.sort_values(by=["logFC", "PValue"], ascending=[False, True])

    # calculating up/down regulated genes

    if names[0] == "logFC" and names[2] == "PValue":
        data_up = stat[(stat["logFC"] >= 1.5) & (stat["PValue"] <= 1e-3)]
        data_down = stat[(stat["logFC"] <= -1.5) & (stat["PValue"] <= 1e-3)]
        up_genes = data_up.index
        down_genes = data_down.index
        logFC_up = data_up["logFC"]
        logFC_down = data_down["logFC"]

    if names[1] == "log2FoldChange" and names[5] == "padj":
        up_genes = stat[(stat["log2FoldChange"] >= 1.5) & (stat["padj"] <= 1e-3)].index
        down_genes = stat[(stat["log2FoldChange"] <= -1.5) & (stat["padj"] <= 1e-3)].index

    return up_genes, down_genes, logFC_up, logFC_down
