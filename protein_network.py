import requests
import tqdm
from graph_tool.all import *

from signature_extractor import *


class protein_network:
    """
    class protein_network could be used to represent PPI (protein-protein interaction) networks both in graph and\
        adjacency list forms;

    Attributes
    ----------
    up_genes: a list of overexpressed genes
    down_genes: a list of underexpressed genes
    interactions_up: a list of interactions for upregulated genes
    interactions_down: a list of interactions for downregulated genes
    adjac_list_up: a list of tuples where each tuple represents interaction of every single upregulated gene
    adjac_list_down: a list of tuples where each tuple represents interaction of every single downregulated gene


    Methods
    -------
    data_preprocessing(mode, species=9606)
        preprocesses a list of genes by checking their presence in STRING database;
        two modes are available:
            - "up_genes" for upregulated genes
            - "down_genes" for downregulated genes

    API_request(mode, species=9606)
        make a request to STRING DB to get interaction network

    creating_adj_list(mode, score)
        creates an ajacency list out list with interactions; possibly could filter only protein metting score\
         requirements

    creating_network(mode)
        creates a graph out of the adjacency matrix made below

    writing_adj_lists(mode)
        stores adjacency list into a file called mode.adjlist

    visualising_graph(mode)
        visualises PPI network

    """

    def __init__(self, file, mode, score):

        self.up_genes, self.down_genes = signature_extractor(file)
        self.interactions_up = []
        self.interactions_down = []
        self.adjac_list_up = []
        self.adjac_list_down = []
        self.mode = mode
        self.score = score
        self.graph = []
        self.up_not_in_STRING = None
        self.down_not_in_STRING = None

    def data_preprocessing(self, species=9606):

        string_api_url = "https://string-db.org/api"
        output_format = "tsv-no-header"
        method = "network"

        # Construct URL

        request_url = "/".join([string_api_url, output_format, method])

        # Set parameters

        my_genes = {
            'up_genes': self.up_genes[:100],
            'down_genes': self.down_genes[:100]
        }
        genes_in_string = {
            'up_genes': [],
            'down_genes': []
        }

        for gene in tqdm.tqdm(my_genes[self.mode]):

            params = {

                "identifiers": "%0d".join(gene),  # your protein
                "species": species,  # species NCBI identifier
                "caller_identity": "www.awesome_app.org"  # your app name

            }

            # Call STRING

            response = requests.post(request_url, data=params)

            flag = True
            for line in response.text.strip().split("\n"):
                l = line.strip().split("\t")
                try:
                    l[2]

                except IndexError:
                    flag = False
                    continue

            if flag:
                genes_in_string[self.mode].append(gene)

        if self.mode == 'up_genes':
            self.up_not_in_STRING = [x for x in self.up_genes if x not in genes_in_string]
            self.up_genes = genes_in_string['up_genes']
        else:
            self.down_not_in_STRING = [x for x in self.down_genes if x not in genes_in_string]
            self.down_genes = genes_in_string['down_genes']

        print(len(genes_in_string[self.mode]))

    def API_request(self):

        # Preprocessing genes
        species = 9606
        self.data_preprocessing()

        string_api_url = "https://string-db.org/api"
        output_format = "tsv-no-header"
        method = "network"

        request_url = "/".join([string_api_url, output_format, method])

        my_genes = {
            'up_genes': self.up_genes,
            'down_genes': self.down_genes
        }

        # Set parameters

        params = {

            "identifiers": "%0d".join(my_genes[self.mode]),  # your protein
            "species": species,  # species NCBI identifier
            "caller_identity": "www.awesome_app.org"  # your app name

        }

        # Call STRING

        response = requests.post(request_url, data=params)

        interactions = []

        for line in response.text.strip().split("\n"):
            l = line.strip().split("\t")
            try:
                p1, p2 = l[2], l[3]
                # filter the interaction according to experimental score
                experimental_score = float(l[10])
                interactions.append((p1, p2, experimental_score))
            except IndexError:
                print("Getting an error")
                print(response.text)

        if self.mode == 'up_genes':
            self.interactions_up = interactions
        else:
            self.interactions_down = interactions


    def creating_adj_list(self):

        self.API_request()
        adja_list = []
        genes = {
            'up_genes': self.up_genes[:1000],
            'down_genes': self.down_genes[:1000]
        }
        if self.mode == "up_genes":
            interactions = self.interactions_up[:1000]
        else:
            interactions = self.interactions_down[:1000]

        for ind_i, line in enumerate(genes[self.mode]):
            adja_list.append([])
            for ind_j, col in enumerate(genes[self.mode]):
                for ind_k, inter in enumerate(interactions):
                    if line == inter[0] and col == inter[1] and float(inter[2]) >= self.score:
                        if adja_list[ind_i] != [] and adja_list[ind_i][-1] != inter[1]:
                            adja_list[ind_i].append(inter[1])
                        elif adja_list[ind_i] == [] and inter[1] != "" and adja_list[ind_i - 1] != inter[1]:
                            adja_list[ind_i].append(inter[1])

            adja_list[ind_i] = tuple([genes[self.mode][ind_i]] + adja_list[ind_i])

        if self.mode == 'up_genes':
            for gene in self.up_not_in_STRING:
                adja_list.append((gene,))
            self.adjac_list_up = adja_list
        else:
            for gene in self.down_not_in_STRING:
                adja_list.append((gene,))
            self.adjac_list_down = adja_list


    def creating_network(self):

        self.creating_adj_list()
        g = Graph(directed=False)
        g.add_vertex()

        if self.mode == "up_genes":
            adj_list = self.adjac_list_up
            gene_set = self.up_genes + self.up_not_in_STRING
        else:
            adj_list = self.adjac_list_down
            gene_set = self.down_genes + self.down_not_in_STRING

        adj_list_int = []
        for ind_i, line in enumerate(adj_list):
            temp_tuple = []
            for ind_j, word in enumerate(line):
                for ind_k, gene in enumerate(gene_set):
                    if word == gene:
                        temp_tuple.append(ind_k)
            adj_list_int.append(tuple(temp_tuple))
        g.add_vertex(len(gene_set) - 1)
        g.add_edge_list(adj_list_int)
        self.graph = g

        return g

    def writing_adj_lists(self):

        if self.graph:
            graph = self.graph
            print("Came here")
        else:
            graph = self.creating_network()
        ## добавить директорию для записи
        if self.mode == "up_genes":
            graph.save("~/Topo_Camp/up_genes.gt")
        if self.mode == "down_genes":
            graph.save("~/Topo_Camp/down_genes.gt")

    def visualising_graph(self, file):

        if self.graph:
            graph = self.graph
        else:
            graph = self.creating_network()

        return graph_draw(graph, output=file)

