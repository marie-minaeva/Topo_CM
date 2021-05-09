from protein_network import *
from graph_tool.centrality import *
from sklearn.preprocessing import *
import pandas as pd
from signature_extractor import *
from pathlib import Path
import argparse
import numpy as np
from scipy.stats import zscore


class centrality_metrics(protein_network) :

    def __init__(self, score, genes, FC) -> object :
        super().__init__(score, genes, FC)

    def metric_calculator(self):
        metrics = []
        self.creating_network()
        graph = self.visualising_graph("graph.png")
        cur = []
        metrics.append(self.genes + self.not_in_STRING)
        """
        for gene in metrics[0]:
            cur.append(np.abs(self.FC[gene]))
        cur = self.FC
        metrics.append(cur)
        """
        metrics.append([gene for gene in pagerank(self.graph)])
        for gene in self.not_in_STRING:
            metrics[1].append(0.0)
        metrics.append([gene for gene in betweenness(self.graph)[0]])  ##одинаковые
        for gene in self.not_in_STRING:
            metrics[2].append(0.0)
        metrics.append([gene for gene in eigenvector(self.graph, max_iter=1e4)[1]])  ##одинаковые
        for gene in self.not_in_STRING:
            metrics[3].append(0.0)
        metrics.append([gene for gene in closeness(self.graph)])
        for gene in self.not_in_STRING:
            metrics[4].append(0.0)
        metrics.append([gene for gene in katz(self.graph)])
        for gene in self.not_in_STRING:
            metrics[5].append(0.0)
        try:
            metrics.append([gene for gene in hits(self.graph, max_iter=1e4)[1]])  # одинаковые
        except ZeroDivisionError:
            metrics.append([np.nan for gene in range(len(self.genes))])  # одинаковые
            print("Cannot calculate Hits metrics because of a float division by zero")
        for gene in self.not_in_STRING:
            metrics[6].append(0.0)
        metrics.append([gene for gene in eigentrust(self.graph, self.graph.edge_properties["scores"], max_iter=1e4)])
        mean = np.mean(metrics[7][:len(self.genes)-len(self.in_biogrid)])
        print(mean)
        for gene in range(len(self.genes)-len(self.in_biogrid), len(self.genes)):
            metrics[7][gene] = mean

        for gene in self.not_in_STRING:
            metrics[7].append(0.0)
        metrics_1 = pd.DataFrame(metrics[1:]).T
        scaler = StandardScaler()
        cur = scaler.fit_transform(metrics_1)
        metrics_1[0] = np.abs(cur.T[0])

        #metrics_1 = metrics_1.fillna(1.0) ### убрать и сделать в подсчете инф_скора if inf_score==Nan inf_score = 1.0
        # ВРЕМЕННО УБРАТЬ НОРМАЛИЗАЦИЮ
        for i in range(1, 7):
            metrics[i] = metrics_1[i - 1]
        metrics = pd.DataFrame(metrics)
        print(metrics)

        return metrics
