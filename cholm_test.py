from scipy import stats
import pandas as pd
import numpy as np
import seaborn as sns

def kolm_test(cmap_db, pref):
    cfm_db = pd.read_csv("~/Topo_Camp/Data/intersection_of_CFM_and_L1000FWD.csv")
    dist = []
    cids = []
    cfm_db["CID"] = [str(cid) for cid in cfm_db["CID"]]
    for ind, chem in enumerate(cmap_db['pubchem_id']):
        if chem in list(cfm_db["CID"]):
            cids.append(chem)
            dist.append(cmap_db['cosine_dist'].loc[ind])
    statistics = []
    for i in range(10):
        dist_rand = np.random.choice(list(cmap_db['cosine_dist']), len(dist))
        stat, pval = stats.ks_2samp(dist, dist_rand)
        statistics.append((stat, pval))
        sns_plot = sns.displot([dist, dist_rand], kde=True)
        sns_plot.savefig(pref+str(i)+".png")

    return statistics