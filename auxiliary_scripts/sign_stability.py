import pandas as pd
from tqdm import tqdm
from BD_signature_parser import  BD_signature_parser
import time

def tanimoto (list1, list2):
  intersec = [common_item for common_item in list1 if common_item in list2]
  return intersec, float(len(intersec))/(len(list1) + len(list2) - len(intersec))

signatures, sig_names = BD_signature_parser("/Users/littlequeen/Downloads/CD_signatures_binary_42809.gmt")
metadata = pd.read_csv("~/Downloads/CD_signature_metadata.csv")
drug_meta = pd.read_csv("~/Downloads/Drugs_metadata.csv")
sig_for_drugs = []
cur = []

drug_ids = pd.unique(metadata["cell_id"])
drug_names = []
smiles = []
cell_line = []
for drug in drug_ids:
    data = metadata[metadata["pert_id"] == drug]
    drug_name = drug_meta[drug_meta["pert_id"] == drug]
    sig_for_drugs.append(list(data["sig_id"]))
    drug_names.append(drug_name["pert_iname"].iloc[0].lower())
    smiles.append(drug_name["canonical_smiles"])
    cell_line.append(drug_name["cell_id"])
print(len(sig_for_drugs))
#print(len(drug_names))
#print(len(smiles))
outdata = []

for k in tqdm(range(len(sig_for_drugs))):
    print(len(sig_for_drugs[k]))
    for i in range(len(sig_for_drugs[k])):
        start = time.time()
        for j in range(i, len(sig_for_drugs[k])):
            set1_up = signatures[sig_names.index(sig_for_drugs[k][i])]
            set1_down = signatures[sig_names.index(sig_for_drugs[k][i])+1]
            set2_up = signatures[sig_names.index(sig_for_drugs[k][j])]
            set2_down = signatures[sig_names.index(sig_for_drugs[k][j]) + 1]
            intersec_up, tc_up = tanimoto(set1_up, set2_up)
            intersec_down, tc_down = tanimoto(set1_down, set2_down)
            #sig_for_drug
            outdata.append([drug_ids[k], cell_line[k], drug_names[k],  sig_for_drugs[k][i], sig_for_drugs[k][j], list(smiles[k])[0], intersec_up, tc_up,  intersec_down, tc_down])
            #sig_for_cell
            #outdata.append([drug_ids[k], sig_for_drugs[k][i], sig_for_drugs[k][j], intersec_up, tc_up,  intersec_down, tc_down])
        end = time.time()
        print(end - start)

out_table = pd.DataFrame(outdata)
# sig_for_drug
out_table.columns = ["pert_id", "cell_id", "pert_name", "sig_id_1", "sig_id_2", "canonical_smiles", "List of up intersecting genes", "Tc Up", "List of down intersecting genes", "Tc Down"]
# sig_for_cell
#out_table.columns = ["cell_id", "sig_id_1", "sig_id_2", "List of up intersecting genes", "Tc Up", "List of down intersecting genes", "Tc Down"]
drugs = pd.unique(data["pert_id"])
out = []
for drug in drugs:
    data = out_table[out_table["pert_id"] == drug]
    print(data)
    mean = np.mean(list(data["Tc Up"])) + np.mean(list(data["Tc Down"]))
    print(mean)
    out.append(str(drug) + " " + str(mean * 0.5))

with open("mean_tanimoto.txt", "w") as file:
    file.write("\n".join(out))

out_table.to_csv("~/Downloads/sig_stability.csv", index=False)
