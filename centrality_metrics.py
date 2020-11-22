from protein_network import *
from graph_tool.centrality import *


class centrality_metrics(protein_network):

    def __init__(self, file, mode, score):
        super().__init__(file, mode, score)

    def metric_calculator(self):
        metrics = []
        self.creating_network()
        metrics.append(pagerank(self.graph))
        print("Pagerank:", metrics[0][0])
        metrics.append(betweenness(self.graph))
        print("Betweenness:", metrics[1][0][0])
        metrics.append(central_point_dominance(self.graph, metrics[1][0]))
        print("Central point dominance:", metrics[2])
        metrics.append(closeness(self.graph))
        print("Closeness:", metrics[3][0])
        metrics.append(eigenvector(self.graph))
        print("Eigenvector:", metrics[4][1][0])
        metrics.append(katz(self.graph))
        print("Katz:", metrics[5][0])
        try:
            metrics.append(hits(self.graph))
            print("Hits:", metrics[6])
        except ZeroDivisionError:
            print("Cannot calculate Hits metrics because of a float division by zero")

        return metrics


graph_up = centrality_metrics("~/Downloads/DE_heart_fb_deseq2_edger.txt", mode='up_genes', score=0.1)
metrics_up = graph_up.metric_calculator()

graph_down = centrality_metrics("~/Downloads/DE_heart_fb_deseq2_edger.txt", mode='down_genes', score=0.1)
metrics_down = graph_down.metric_calculator()
