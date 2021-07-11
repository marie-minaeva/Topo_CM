import pandas as pd
import numpy as np
from standard_chemicals_extraction import stand_chems

def in_GS(file, metadata, source, target):
    data = pd.read_csv(file)
    data = data.sort_values(by=["cosine_dist"], ascending=[True])
    data = data.reset_index(drop=True)
    chems = stand_chems(source, target)
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
        if val <= 50:
            less_50 += 1
    print(less_50)
    temp.sort()
    print(temp)
in_GS("~/Downloads/small_mol_neuron_correct_de_validated.csv", None, "Fibroblasts", "Induced Neurons")