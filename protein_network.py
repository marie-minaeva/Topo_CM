from collections import defaultdict
import requests
import tqdm
from graph_tool.all import *


class ProteinNetwork:
    """
    A class is used to construct PPI network

    Attributes:
        genes (:obj:`list` of:obj:`str`): list of either up or down regulated genes of request signature
        FC (dict): list of logFC of either up or down regulated genes of request signature
        interactions (:obj:`list` of:obj:`tuple` of:obj:`str` or `float`): list of interactions with confidence
            scores
        adjac_list (:obj:`collections.defaultdict`): adjacency lists for each gene
        trust_array (:obj:`collections.defaultdict` of:obj:`collections.defaultdict`): a matrix of confident scores
            for interactions between each protein
        score (float): the lowest confident score taken into consideration
        graph (:obj: `graph_tool.Graph`): a graph of PPI network
        not_in_STRING (:obj:`list` of:obj:`str`): genes that were not found in STRING
        in_biogrid (:obj:`list` of:obj:`str`): genes that were found in BioGRID

    """

    def __init__(self, score, genes, FC):
        """
        Declares class attributes

        Args:
            score (float): the lowest confident score taken into consideration
            genes (:obj:`list` of:obj:`str`): list of either up or down regulated genes of request signature
            FC (dict): list of logFC of either up or down regulated genes of request signature

        """
        self.genes = genes
        self.FC = dict(zip(self.genes, FC))
        self.interactions = []
        self.adjac_list = []
        self.trust_array = []
        self.score = score
        self.graph = []
        self.not_in_STRING = []
        self.in_biogrid = None

    def data_preprocessing(self, species=9606):
        """
        Preprocesses genes, checks whether gene is in STRING or not

        Args:
            species (int, optional): specie identifier

        """

        string_api_url = "https://string-db.org/api"
        output_format = "tsv-no-header"
        method = "network"

        # Construct URL

        request_url = "/".join([string_api_url, output_format, method])

        # Set parameters
        my_genes = self.genes
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

        self.not_in_STRING = [x for x in self.genes if x not in genes_in_string]
        self.genes = genes_in_string

        print(len(genes_in_string))

    def API_request(self):
        """
        Conducts API request to STRING for preprocessed genes

        """
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

    def BIOGRID_request(self):
        """

        Conducts API request to BioGRID

        """
        request_url = "https://webservice.thebiogrid.org/" + "/interactions"

        # List of genes to search for
        geneList = self.genes  # Yeast Genes STE11 and NMD4
        evidenceList = ["POSITIVE GENETIC"]

        # These parameters can be modified to match any search criteria following
        # the rules outlined in the Wiki: https://wiki.thebiogrid.org/doku.php/biogridrest
        params = {
            "accesskey": "86ef3d2746f82a8c3644424ae4b5538c",
            "format": "tab2",  # Return results in TAB2 format
            "geneList": "|".join(geneList),  # Must be | separated
            "searchNames": "true",  # Search against official names
            "includeInteractors": "true",
            # Set to true to get any interaction involving EITHER gene, set to false to get interactions between genes
            "taxId": 559292,  # Limit to Homo Sapiens
            "evidenceList": "|".join(evidenceList),  # Exclude these two evidence types
            "includeEvidence": "true",
            # If false "evidenceList" is evidence to exclude, if true "evidenceList" is evidence to show
            "includeHeader": "true"
        }

        r = requests.post(request_url, params=params)
        interaction = r.text
        interactions = []
        genes = []
        for line in interaction.strip().split("\n"):
            l = line.strip().split("\t")
            try:
                p1, p2 = l[7], l[8]
                # filter the interaction according to experimental score
                experimental_score = 1.0
                genes.append(p1)
                interactions.append((p1, p2, experimental_score))
            except IndexError:
                print("Getting an error")
                print(interaction)
        print(interactions)
        self.interactions.extend(interactions[1:])
        print(genes)
        self.genes.extend(genes[1:])
        try:
            self.not_in_STRING.remove(genes[1:])
        except ValueError:
            pass
        self.in_biogrid = genes[1:]

    def creating_adj_list(self):
        """

        Parses self.interactions to create adjacency list and matrix of confidence scores

        """

        self.API_request()
        self.BIOGRID_request()
        adja_list = defaultdict(list)
        scores = defaultdict(lambda: defaultdict(list))
        genes = self.genes

        for ind_i, line in enumerate(genes):
            for ind_j, col in enumerate(genes):
                for ind_k, inter in enumerate(self.interactions):
                    if line == inter[0] and col == inter[1] and float(inter[2]) >= self.score:
                        if adja_list[genes[ind_i]] != [] and adja_list[genes[ind_i]][-1] != inter[1]:
                            adja_list[genes[ind_i]] .append(inter[1])
                            scores[genes[ind_i]][inter[1]].append(float(inter[2]))
                        elif adja_list[genes[ind_i]] == [] and inter[1] != "":
                            adja_list[genes[ind_i]].append(inter[1])
                            scores[genes[ind_i]][inter[1]].append(float(inter[2]))
        """
        for gene in self.not_in_STRING:
            for gene_1 in self.genes + self.not_in_STRING:
                adja_list[gene] = []
                scores[gene][gene_1] = []
        """
        self.adjac_list = adja_list
        self.trust_array = scores

    def creating_network(self):
        """
        Constructs a graph of PPI network based on adjacency list

        Returns:
            g (:obj:`graph_tool.Graph`): graph of reconstructed PPI network

        """
        self.creating_adj_list()
        # gene_set = self.genes + self.not_in_STRING
        gene_set = self.genes
        g = Graph(directed=False)
        g.add_vertex(n=len(gene_set))
        print(len(gene_set))
        vprop_proteins = g.new_vp("string")
        for i in range(len(gene_set)):
            vprop_proteins[i] = gene_set[i]
        g.vertex_properties['proteins'] = vprop_proteins
        eprop_scores = g.new_edge_property("double")
        for gene in self.adjac_list.items():
            for inter in gene[1]:
                edge = g.add_edge(self.genes.index(gene[0]), self.genes.index(inter))
                eprop_scores[edge] = self.trust_array[gene[0]][inter][0]

        g.edge_properties["scores"] = eprop_scores
        self.graph = g
        return g

    def writing_adj_lists(self):
        """

        Writes graph into file

        """

        if self.graph:
            graph = self.graph
        else:
            graph = self.creating_network()
        # добавить директорию для записи
        graph.save("~/Topo_Camp/genes.gt")

    def visualising_graph(self, file):
        """
        Visualises the graph

        Args:
            file (str): path to file where to save visualised graph

        Returns:
            function: function to draw the graph

        """

        if self.graph:
            graph = self.graph
        else:
            graph = self.creating_network()

        return graph_draw(graph, output=file)
