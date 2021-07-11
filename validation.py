import Cmap
import kolm_test
import BD_signature_parser
from signature_extractor import *
import pandas as pd
from gsea_metrics import *

bd_signatures, sig_names = BD_signature_parser.BD_signature_parser("/Users/littlequeen/Downloads/CD_signatures_binary_"
                                                                   "42809.gmt")
trial = Cmap.TopoCMap(bd_signatures, "mimic", "~/Downloads/DE_decitabine_counts_edger.txt")


def validate(coeffs):
    """
    Validates TopoCMap in different ways depending on validation dataset

    Args:
        coeffs (:obj:`list` of :obj:`float`): list of coefficients for inf_score expression

    Returns:
        float: value of function to be optimized in the optimization process (by default Normalized enrichment score)
    """
    source, target = "Fibroblasts", "Induced Pluripotent Stem Cells"
    small_mol = trial.small_molec(sig_names, "~/Downloads/CD_signature_metadata.csv", "~/Downloads/Drugs_metadata.csv",
                                  weights=coeffs)
    small_mol = pd.DataFrame(small_mol)
    small_mol = small_mol.T
    small_mol.columns = ["sign_id", "cosine_dist", "pert_id", "pubchem_id", "pert_iname"]
    small_mol = small_mol.drop_duplicates()
    small_mol = small_mol.sort_values(by=["cosine_dist"], ascending=[True])
    small_mol.to_csv("small_mol_decitabine_counts" + '.csv', index=False)
    # Code for validation on CFM below

    all_stat, mean_gold, mean_all = kolm_test.kolm_test("small_mol_IPS_mean_opt" + '.csv', "small_mol_IPS_mean_opt",
                                                        source, target, "~/Downloads/Drugs_metadata.csv")
    mean = 0.0
    for j in range(len(all_stat)):
        print(all_stat[j][0])
        mean += all_stat[j][0]
    mean /= len(all_stat)

    all_stat = [" ". join(str(x)) for x in all_stat]
    coeffs = [" ". join(str(x)) for x in coeffs]
    with open("../Data for thesis/all_stat_small_mol_IPS_mean_opt.txt", "a") as file:
        file.write('\n'.join(all_stat))
        file.write('\n')
        file.write('\n')

    with open("../Data for thesis/coeffs_small_mol_IPS_mean_opt.txt", "a") as file:
        file.write('\n'.join(coeffs))
        file.write('\n')
        file.write('\n')

    print(mean_gold - mean_all)
    return gsea_metrics(source, target, "small_mol_IPS_mean_opt" + '.csv',  "~/Downloads/Drugs_metadata.csv")


print(validate([5.223, 5.019, 6.973, 6.982, 6.466, 4.311, 5.772, 6.909, 0.715]))
