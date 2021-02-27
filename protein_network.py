import requests
import tqdm
from graph_tool.all import *
from collections import defaultdict
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

    def __init__(self, score, genes, FC):

        self.genes = genes
        self.FC = dict(zip(self.genes, FC))
        self.interactions = None
        self.adjac_list = None
        self.trust_array = None
        self.trust_map = None
        self.score = score
        self.graph = None
        self.not_in_STRING = None

    def data_preprocessing(self, species=9606):

        string_api_url = "https://string-db.org/api"
        output_format = "tsv-no-header"
        method = "network"

        # Construct URL

        request_url = "/".join([string_api_url, output_format, method])

        # Set parameters

        my_genes = self.genes[:100]
        genes_in_string = []

        for gene in tqdm.tqdm(my_genes):

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
                genes_in_string.append(gene)
        print(genes_in_string)
        print(self.genes)
        self.not_in_STRING = [x for x in self.genes if x not in genes_in_string]
        self.genes = genes_in_string

        print(len(genes_in_string))

    def API_request(self):

        # Preprocessing genes
        species = 9606
        self.data_preprocessing()

        string_api_url = "https://string-db.org/api"
        output_format = "tsv-no-header"
        method = "network"

        request_url = "/".join([string_api_url, output_format, method])

        my_genes = self.genes

        # Set parameters

        params = {

            "identifiers": "%0d".join(my_genes),  # your protein
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
                experimental_score = float(l[5])
                interactions.append((p1, p2, experimental_score))
            except IndexError:
                print("Getting an error")
                print(response.text)

        self.interactions = interactions


    def creating_adj_list(self):

        self.API_request()
        adja_list = defaultdict(list)
        scores = defaultdict(lambda: defaultdict(list))

        for ind_i, line in enumerate(self.genes):
            #adja_list[self.genes[ind_i]].append([])
            #scores.append([])
            for ind_j, col in enumerate(self.genes):
                for ind_k, inter in enumerate(self.interactions):
                    if line == inter[0] and col == inter[1] and float(inter[2]) >= self.score:
                        if adja_list[self.genes[ind_i]] != [] and adja_list[self.genes[ind_i]] [-1] != inter[1]:
                            adja_list[self.genes[ind_i]] .append(inter[1])
                            scores[self.genes[ind_i]][inter[1]].append(float(inter[2]))
                        elif adja_list[self.genes[ind_i]] == [] and inter[1] != "":
                            adja_list[self.genes[ind_i]].append(inter[1])
                            scores[self.genes[ind_i]][inter[1]].append(float(inter[2]))

            #adja_list[self.genes[ind_i]] = list(set(adja_list[self.genes[ind_i]]))

        for gene in self.not_in_STRING:
            for gene_1 in self.genes + self.not_in_STRING:
                adja_list[gene] = []
                scores[gene][gene_1] = []
        self.adjac_list = adja_list
        self.trust_array = scores



    def creating_network(self):

        self.creating_adj_list()
        gene_set = self.genes + self.not_in_STRING
        g = Graph(directed=False)
        g.add_vertex(n=len(gene_set))
        print(len(gene_set))
        vprop_proteins = g.new_vp("string")
        for i in range(len(gene_set)):
            vprop_proteins[i] = gene_set[i]
        g.vertex_properties['name_proteins'] = vprop_proteins
        """
        adj_list_int = []
        for ind_i, line in enumerate(self.adj_list.items()):
            temp_tuple = []
            for ind_j, word in enumerate(line):
                for ind_k, gene in enumerate(gene_set):
                    if word == gene:
                        temp_tuple.append(ind_k)

            adj_list_int.append(tuple(temp_tuple))
        """
        #g.add_edge_list(self.adjac_list)
        print(g.vertices())
        self.graph = g
        print(g.vertex(10))
        """
        self.trust_map = g.new_edge_property("double")
        if self.mode == "up_genes":
            for gene in zip(self.trust_array, g.edges()):
                for edge in gene:
                    #self.trust_map[g.vertex(ind)]
                    pass
        else:
            self.trust_map = g.new_edge_property("vector<double>")
            for ind, vert in enumerate(self.graph.vertices()):
                pass
        """
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

up_request, down_request, up_FC, down_FC = signature_extractor("~/Downloads/DE_neuron_fb_deseq2_edger.txt")
trial = protein_network(0.0, up_request, up_FC)
graph = trial.creating_network()
trial.visualising_graph("ex_1.png")