from scipy import stats
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


cfm_db = pd.read_csv("~/Topo_Camp/Data/intersection_of_CFM_and_L1000FWD.csv")
cmap_db = pd.read_csv("~/PycharmProjects/Topo_CM/small_mol.csv")
dist = []
cids = []
cfm_db["CID"] = [str(cid) for cid in cfm_db["CID"]]
for ind, chem in enumerate(cmap_db['pubchem_id']):
    if chem in list(cfm_db["CID"]):
        cids.append(chem)
cids = list(set(cids))
for ind, chem in enumerate(cmap_db['pubchem_id']):
    if chem in cids:
        dist.append(cmap_db['cosine_dist'].loc[ind])
        cids.remove(chem)
print(len(dist))
for i in range(10):
    dist_rand = np.random.choice(list(cmap_db['cosine_dist']), len(dist))
    stat, pval = stats.ks_2samp(dist, dist_rand)
    print(stats.ks_2samp(dist, dist_rand))
    sns.histplot(dist, color='pink', kde=True)
    sns.histplot(dist_rand, kde=True)
    plt.show()