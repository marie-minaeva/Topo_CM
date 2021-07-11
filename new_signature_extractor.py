import numpy as np
import pandas as pd

def new_signature_extractor(file1, file2):
    """
    This function is an upgrade of function signature_extractor to work with signatures of Lamb et al.

    Args:
        file1 (str): Path to the file with up genes
        file2 (str): Path to the file with up genes

    Returns:
        up_genes (:obj:`list` of :obj:`str`): List of up-regulated genes
        down_genes (:obj:`list` of :obj:`str`): List of down-regulated genes
        logfc_up (:obj:`list` of :obj:`float`): List of logFC for up-regulated genes
        logfc_down (:obj:`list` of :obj:`float`): List of logFC for down-regulated genes

    """
    up_genes = pd.read_table(file1, sep="\n", header=None)
    up_genes = up_genes.drop_duplicates()
    down_genes = pd.read_table(file2, sep="\n", header=None)
    down_genes = down_genes.drop_duplicates()
    logfc_up = np.repeat(np.nan, len(list(up_genes)))
    logfc_down = np.repeat(np.nan, len(list(down_genes)))
    return list(up_genes[0]), list(down_genes[0]), logfc_up, logfc_down