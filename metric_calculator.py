from protein_network import *
from graph_tool.centrality import *
from pathlib import Path
import argparse


class centrality_metrics(protein_network):

    def __init__(self, file, mode, score):
        super().__init__(file, mode, score)

    @property
    def metric_calculator(self):
        metrics = []
        self.creating_network()
        if self.mode == "up_genes":
            metrics.append(self.up_genes)
        else:
            metrics.append(self.down_genes)
        metrics.append([gene for gene in pagerank(self.graph)])
        metrics.append([gene for gene in betweenness(self.graph)[0]])
        #metrics.append([gene for gene in eigenvector(self.graph)[1]])
        #metrics.append(central_point_dominance(self.graph, metrics[1][0]))
        metrics.append([gene for gene in closeness(self.graph)])
        metrics.append([gene for gene in katz(self.graph)])
        try:
            metrics.append([gene for gene in hits(self.graph, max_iter=1e6)[2]])
        except ZeroDivisionError:
            print("Cannot calculate Hits metrics because of a float division by zero")
        metrics = pd.DataFrame(metrics)
        metrics = metrics.fillna(1.0)
        return metrics


parser = argparse.ArgumentParser()
parser.add_argument('mode', nargs="?", help="Either up_genes or down_genes", default=["down_genes"])
parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=Path("~/Downloads"
                                                                                   "/DE_heart_fb_deseq2_edger.txt"))
parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'))
args = parser.parse_args()
