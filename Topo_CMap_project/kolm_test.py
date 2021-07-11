from scipy import stats
import pandas as pd
import numpy as np
import seaborn as sns
from standard_chemicals_extraction import *

def kolm_test(file, pref):
    dist = []
    cids = []
    cmap_db = pd.read_csv(file)
    drug_meta = pd.read_csv("~/Downloads/Drugs_metadata.csv")
    cids_cur = stand_chems("Fibroblasts", "Induced Neurons")
    cids_cur = [float(cid) for cid in cids_cur]
    pert_cur = []
    for ind, chem in enumerate(drug_meta['pubchem_cid']):
        for chem_1 in cids_cur:
            try:
                if int(chem) == int(chem_1):
                    pert_cur.append(drug_meta['pert_id'].loc[ind])
            except ValueError:
                continue


    for ind, chem in enumerate(cmap_db['pert_id']):
        for chem_1 in pert_cur:
            if chem == chem_1:
                cids.append(chem)
                dist.append(cmap_db['cosine_dist'].loc[ind])

    print(dist)
    statistics = []
    mean_1 = np.mean(dist)
    means = []
    for i in range(10):
        dist_rand = np.random.choice(list(cmap_db['cosine_dist']), len(dist))
        means.append(np.mean(dist_rand))
        stat, pval = stats.ks_2samp(dist, dist_rand)
        statistics.append((stat, pval))
        sns_plot = sns.displot([dist, dist_rand], kde=True)
        sns_plot.savefig(pref+str(i)+".png")
    mean_2 = np.mean(means)
    return statistics, mean_1, mean_2
