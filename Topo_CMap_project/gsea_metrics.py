import pandas as pd
import numpy as np

def gsea_metrics(small_mol, golden_stand):
    ess = []
    for i in range(0, small_mol.shape[0]):
        p_hit = 0.0
        p_miss = 0.0
        for ind in range(0, i):
            N_r = small_mol[small_mol["pert_id"].isin(golden_stand)]
            N_r = np.sum(list(N_r["cosine_dist"]))
            if list(small_mol["pert_id"])[ind] in golden_stand:
                p_hit += list(small_mol["cosine_dist"])[ind] / N_r
            else:
                p_miss += 1 / (len(small_mol["pert_id"]) - len(golden_stand))
        ess.append(p_hit - p_miss)

    ess = np.abs(ess)
    print(np.argmax(ess))
    return np.max(ess)