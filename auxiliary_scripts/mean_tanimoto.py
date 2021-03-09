import pandas as pd
import numpy as np
import tqdm

with open("sig_stability.csv", "r") as file:
    lines = file.readlines()

drugs = []
tanimotos = []
for i in tqdm.tqdm(range(1, len(lines))):
    cur = lines[i].split(",")
    cur[-1] = cur[-1][:-2]
    score = []
    for symb in cur:
        try:
            score.append(float(symb))
        except ValueError:
            continue
    drugs.append(cur[0])
    tanimotos.append((score[0] + score[1])*0.5)
out = []
data = pd.DataFrame([drugs, tanimotos])
data = data.T
data.columns = ["pert_id", "TC"]
drugs = list(set(drugs))

for drug in drugs:
    data_cur = data[data["pert_id"] == drug]
    mean = np.mean(data_cur["TC"])
    print(mean)
    out.append(str(drug) + " " + str(mean))

with open("mean_tanimoto.txt", "w") as file:
    file.write("\n".join(out))
