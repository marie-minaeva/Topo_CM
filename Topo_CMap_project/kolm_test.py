from scipy import stats
import pandas as pd
import numpy as np
import seaborn as sns
from standard_chemicals_extraction import *

def kolm_test(file, pref):
    dist = []
    cids = []
    cmap_db = pd.read_csv(file)
    cids_cur = stand_chems("Fibroblasts", "Induced Neurons")
    cids_cur = [float(cid) for cid in cids_cur]
    for ind, chem in enumerate(cmap_db['pubchem_id']):
        for chem_1 in cids_cur:
            print(chem, chem_1)
            try:
                if chem == chem_1:
                    print(chem, chem_1)
                    cids.append(chem)
                    dist.append(cmap_db['cosine_dist'].loc[ind])
            except ValueError:
                print(chem, type(chem))

    print(dist)
    statistics = []
    for i in range(10):
        dist_rand = np.random.choice(list(cmap_db['cosine_dist']), len(dist))
        stat, pval = stats.ks_2samp(dist, dist_rand)
        statistics.append((stat, pval))
        sns_plot = sns.displot([dist, dist_rand], kde=True)
        sns_plot.savefig(pref+str(i)+".png")

    return statistics
