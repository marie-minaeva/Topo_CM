import Cmap
import kolm_test
import BD_signature_parser
from signature_extractor import *
import pandas as pd
import numpy as np

bd_signatures, sig_names = BD_signature_parser.BD_signature_parser("/Users/littlequeen/Downloads/CD_signatures_binary_42809.gmt")
trial = Cmap.TopoCMap(bd_signatures[:200], "reverse", "~/Downloads/DE_neuron_fb_deseq2_edger.txt")
def validate(coeffs):
    print(len(coeffs))
    small_mol = trial.small_molec(sig_names[:100], "~/Downloads/CD_signature_metadata.csv", "~/Downloads/Drugs_metadata.csv", weights=coeffs)
    small_mol = pd.DataFrame(small_mol)
    small_mol = small_mol.T
    small_mol.columns = ["sign_id", "cosine_dist", "pert_id", "pubchem_id"]
    small_mol = small_mol.drop_duplicates()
    small_mol.to_csv("small_mol_neuron_" + '.csv', index=False)
    all_stat = kolm_test.kolm_test("small_mol_neuron_"  + '.csv', "test_1_neuron_" )
    mean = 0.0
    for j in range(len(all_stat)):
        print(all_stat[j][0])
        mean += all_stat[j][0]
    mean /= len(all_stat)

    all_stat = [" ". join(str(x)) for x in all_stat]
    coeffs = [" ". join(str(x)) for x in coeffs]
    with open("all_stat.txt", "a") as file:
        file.write('\n'.join(all_stat))
        file.write('\n')
        file.write('\n')


    with open("coeffs.txt", "a") as file:
        file.write('\n'.join(coeffs))
        file.write('\n')
        file.write('\n')

    return - mean
