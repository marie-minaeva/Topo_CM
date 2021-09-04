
import pandas as pd


def signature_extractor(file):
    """
    Function is used to calculate signatures out of differentially expressed genes

    Args:
        file (str): path to the file with DE genes

    Returns:
        up_genes (:obj:`list` of :obj:`str`): List of up-regulated genes
        down_genes (:obj:`list` of :obj:`str`): List of down-regulated genes
        logFC_up (:obj:`list` of :obj:`float`): List of logFC for up-regulated genes
        logFC_down (:obj:`list` of :obj:`float`): List of logFC for down-regulated genes
    """

    # preprocessing data before
    # calculating up/down regulated genes

    stat = pd.read_table(file, sep=" ")
    stat = stat.apply(pd.to_numeric)

    up_genes = []
    down_genes = []
    logFC_up = []
    logFC_down = []
    names = stat.columns
    stat = stat.sort_values(by=["logFC", "PValue"], ascending=[False, True])

    # calculating up/down regulated genes
    if names[0] == "logFC" and names[2] == "PValue":
        data_up = stat[(stat["logFC"] >= 4.5) & (stat["PValue"] <= 1e-3)]
        data_down = stat[(stat["logFC"] <= -2.5) & (stat["PValue"] <= 1e-3)]
        if len(list(data_down.index)) > 5000:
            data_down = stat[(stat["logFC"] <= -2.5) & (stat["PValue"] <= 1e-3)]
        up_genes = data_up.index
        down_genes = data_down.index
        logFC_up = data_up["logFC"]
        logFC_down = data_down["logFC"]

    if names[1] == "log2FoldChange" and names[5] == "padj":
        up_genes = stat[(stat["log2FoldChange"] >= 1.5) & (stat["padj"] <= 1e-3)].index
        down_genes = stat[(stat["log2FoldChange"] <= -1.5) & (stat["padj"] <= 1e-3)].index

    with open("up_genes_ips.txt", "w") as file:
        file.write('\t'.join(up_genes))
    with open("down_genes_ips.txt", "w") as file:
        file.write('\t'.join(down_genes))
    return up_genes, down_genes, logFC_up, logFC_down


# Code for PWx project
"""
import pandas as pd

def signature_extractor(file):
    '''
    Function is used to calculate signatures out of differentially expressed genes;
   :param file: a name of file with DE genes;
   :return: two lists: one with up-regulated  and the other with down-regulated genes;
    '''

    # preprocessing data before
    # calculating up/down regulated genes

    stat = pd.read_table(file, sep="\t")
    stat = stat.apply(pd.to_numeric)
    print(stat)
    up_genes = []
    down_genes = []
    names = stat.columns
    stat = stat.sort_values(by=["FC", "P.adj"], ascending=[False, True])

    # calculating up/down regulated genes

    if names[1] == "FC" and names[0] == "P.adj":
        data_up = stat[(stat["FC"] >= 0.2) & (stat["P.adj"] <= 0.05)]
        data_down = stat[(stat["FC"] <= -0.2) & (stat["P.adj"] <= 0.05)]
        up_genes = data_up.index
        down_genes = data_down.index
        logFC_up = data_up["FC"]
        logFC_down = data_down["FC"]

    if names[1] == "log2FoldChange" and names[5] == "padj":
        up_genes = stat[(stat["log2FoldChange"] >= 0.2) & (stat["padj"] <= 5e-2)].index
        down_genes = stat[(stat["log2FoldChange"] <= -0.2) & (stat["padj"] <= 5e-2)].index
    print(len(up_genes), len(down_genes))
    up_string = '\t'.join(up_genes)
    down_string = '\t'.join(down_genes)
    up_out = ["PWx", "signature for PWX", up_string]
    down_out = ["PWx", "signature for PWX", down_string]
    with open("up_genes.gmt", "w") as file:
        file.write('\t'.join(up_out))
    with open("down_genes.gmt", "w") as file:
        file.write('\t'.join(down_out))
    return up_genes, down_genes, logFC_up, logFC_down
"""

# signature_extractor("~/Downloads/DE_IPS_fb_deseq2_edger.txt")
