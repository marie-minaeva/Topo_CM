from protein_network import *
from graph_tool.centrality import *
from sklearn.preprocessing import *
from pathlib import Path
import argparse


class centrality_metrics(protein_network):

    def __init__(self, file, mode, score) -> object:
        super().__init__(file, mode, score)

    @property
    def metric_calculator(self):
        metrics = []
        self.creating_network()
        cur = []
        if self.mode == "up_genes":
            metrics.append(self.up_genes + self.up_not_in_STRING)
            for gene in metrics[0]:
                cur.append(self.FC_up[gene])
            metrics.append(cur)
        else:
            metrics.append(self.down_genes + self.down_not_in_STRING)
            for gene in metrics[0]:
                cur.append(self.FC_down[gene])
            metrics.append(cur)
        metrics.append([gene for gene in pagerank(self.graph)])
        metrics.append([gene for gene in betweenness(self.graph)[0]])
        metrics.append([gene for gene in eigenvector(self.graph, max_iter=1e6)[1]])
        #metrics.append(central_point_dominance(self.graph, metrics[1][0]))
        metrics.append([gene for gene in closeness(self.graph)])
        metrics.append([gene for gene in katz(self.graph)])
        try:
            metrics.append([gene for gene in hits(self.graph, max_iter=1e6)[2]])
        except ZeroDivisionError:
            print("Cannot calculate Hits metrics because of a float division by zero")
        metrics.append([gene for gene in eigentrust(self.graph, self.trust_map)])
        metrics.append([gene for gene in trust_transitivity(self.graph, self.trust_map)])
        metrics_1 = pd.DataFrame(metrics[1:])
        metrics_1 = metrics_1.fillna(1.0)
        print(metrics_1.describe())
        scaler = StandardScaler()
        metrics_1 = scaler.fit_transform(metrics_1)
        print(pd.DataFrame(metrics_1).describe())
        print(metrics_1)
        for i in range(1,6):
            metrics[i] = metrics_1[i-1]
        metrics = pd.DataFrame(metrics)
        print(metrics.head())
        return metrics

