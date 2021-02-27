import Cmap
import cholm_test
import BD_signature_parser
from signature_extractor import *
import pandas as pd
import numpy as np

up_request, down_request, FC_up, FC_down = signature_extractor("~/Downloads/DE_neuron_fb_deseq2_edger.txt")
bd_signatures, sig_names = BD_signature_parser.BD_signature_parser("/Users/littlequeen/Downloads"
                                                                   "/CD_signatures_binary_42809.gmt")
trial = Cmap.TopoCMap(bd_signatures[:200], "reverse", "~/Downloads/DE_neuron_fb_deseq2_edger.txt")
all_stat = []
coeffs = [[1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0]]

for i in range(1):
    #coeffs.append(list(np.random.uniform(1.0, 5.0, 6)))
    small_mol = trial.small_molec(sig_names[:100], "~/Downloads/CD_signature_metadata.csv", "~/Downloads"
                                                                                            "/Drugs_metadata.csv",
                                  weights=coeffs[i])
    small_mol = pd.DataFrame(small_mol)
    small_mol = small_mol.T
    small_mol.columns = ["sign_id", "cosine_dist", "pert_id", "pubchem_id"]
    all_stat.append(cholm_test.kolm_test(small_mol, "test_1_neuron_" + str(i)))
    small_mol.to_csv("small_mol_neuron_" + str(i) + '.csv', index=False)

all_stat = [" ". join(str(x)) for x in all_stat]
coeffs = [" ". join(str(x)) for x in coeffs]
with open("all_stat.txt", "w") as file:
    file.write('\n'.join(all_stat))


with open("coeffs.txt", "w") as file:
    file.write('\n'.join(coeffs))
