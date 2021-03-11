from protein_network import *
from graph_tool.centrality import *
# from sklearn.preprocessing import *
import pandas as pd
from signature_extractor import *
from pathlib import Path
import argparse
import numpy as np


class centrality_metrics(protein_network) :

    def __init__(self, score, genes, FC) -> object :
        super().__init__(score, genes, FC)

    def metric_calculator(self):
        metrics = []
        self.creating_network()
        cur = []
        metrics.append(self.genes + self.not_in_STRING)
        for gene in metrics[0]:
            cur.append(np.abs(self.FC[gene]))
        metrics.append(cur)
        metrics.append([gene for gene in pagerank(self.graph)])
        metrics.append([gene for gene in betweenness(self.graph)[0]])  ##одинаковые
        metrics.append([gene for gene in eigenvector(self.graph, max_iter=1e4)[1]])  ##одинаковые
        metrics.append([gene for gene in closeness(self.graph)])
        metrics.append([gene for gene in katz(self.graph)])
        try:
            metrics.append([gene for gene in hits(self.graph, max_iter=1e4)[1]])  # одинаковые
            #metrics.append([gene for gene in hits(self.graph, max_iter=1e6)[2]])  # одинаковые
        except ZeroDivisionError:
            print("Cannot calculate Hits metrics because of a float division by zero")
        metrics.append([gene for gene in eigentrust(self.graph, self.graph.edge_properties["scores"], max_iter=1e4)])
        # metrics.append([gene for gene in trust_transitivity(self.graph, self.graph.edge_properties["scores"])])
        metrics_1 = pd.DataFrame(metrics[1:]).T
        #metrics_1 = metrics_1.fillna(1.0) ### убрать и сделать в подсчете инф_скора if inf_score==Nan inf_score = 1.0
        # scaler = StandardScaler()
        # metrics_1 = scaler.fit_transform(metrics_1)
        # ВРЕМЕННО УБРАТЬ НОРМАЛИЗАЦИЮ
        n_nzeros = metrics_1.ne(0).sum(axis=0)
        print(n_nzeros)
        # print(pd.DataFrame(metrics_1.T).describe())
        for i in range(1, 8):
            metrics[i] = metrics_1[i - 1]
        metrics = pd.DataFrame(metrics)
        print(metrics)

        return metrics
