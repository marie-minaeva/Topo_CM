import pandas as pd
import numpy as np
from standard_chemicals_extraction import stand_chems


def in_GS(file, metadata, source, target):
    """
    Extracts positions of molecules of 'Golden Standard' for validation on CFM

    Args:
        file (str): Path to the file with TopoCMap output table
        metadata (str): Path to the file with drugs metadata
        source (str): Source cell type
        target (str): Target cell type

    """
    data = pd.read_csv(file)
    data = data.sort_values(by=["cosine_dist"], ascending=[True])
    data = data.reset_index(drop=True)
    chems = stand_chems(source, target)
    print(pd.DataFrame(chems).drop_duplicates())

    dictionary = dict.fromkeys(chems, [])
    for chem in chems:
        dictionary[chem] = np.array(data[data["pubchem_id"] == chem].index)
    print(dictionary)
    print(dictionary.values())
    less_50 = 0
    temp = []
    for chem in dictionary.values():
        temp.extend(chem)
    for val in temp:
        if val <= 2140:
            less_50 += 1
    print(less_50)
    temp.sort()
    print(temp)


# in_GS("../Data for thesis/small_mol_neuron_mes_general_coeffs.csv", None, "Mesenchymal Stem Cells", "Induced Neurons")
