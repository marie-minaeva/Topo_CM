import time
from collections import *
from scipy.spatial import distance
import pandas as pd
import tqdm
from BD_signature_parser import *
from metric_calculator import *
import numpy as np
from signature_extractor import *
from new_signature_extractor import *


def decomposition(data):
    """
    Does decomposition of gene vectors on the bases of gene spaces

    Args:
        data (:obj:`list` of :obj:`list` of :obj:`str`): list of spaces (up space and down space)

    Returns:
        vector (:obj:`list` of :obj:`list` of :obj:`float`): list of 4 decomposed vectors

    """
    vector = []
    for j in range(len(data)):
        out_vector_1 = []
        out_vector_2 = []
        for gene, val in data[j].items():
            if len(val) == 2:
                out_vector_1.append(1.0)
                out_vector_2.append(1.0)
            if len(val) == 1 and val[0] == 0:
                out_vector_1.append(1.0)
                out_vector_2.append(0.0)
            if len(val) == 1 and val[0] == 1:
                out_vector_1.append(0.0)
                out_vector_2.append(1.0)

        vector.append(out_vector_1)
        vector.append(out_vector_2)
    return vector


def inf_score_expression(matrix, weights):
    """
    Calculates influence scores out of matrix of centrality metrics

    Args:
        matrix (:obj:`pandas.DataFrame`): matrix of centrality metrics
        weights (:obj:`list` of :obj:`float`): list of weights for each row of matrix

    Returns:
        vector (:obj:`list` of :obj:`float`): list with calculated influence scores

    """
    vector = (np.array(list(matrix.loc[1])) * weights[0] + 1) * (
            np.array(list(matrix.loc[2])) * weights[1] + 1) * (
                     np.array(list(matrix.loc[3])) * weights[2] + 1) * (
                     np.array(list(matrix.loc[4])) * weights[3] + 1) * (
                     np.array(list(matrix.loc[5])) * weights[4] + 1) * (
                     np.array(list(matrix.loc[6])) * weights[5] + 1) * (
                     np.array(list(matrix.loc[7])) * weights[6] + 1) * (
                     np.array(list(matrix.loc[8])) * weights[7] + 1) + weights[8]
    """
    Code for working with CMap signatures
    vector = (np.array(list(matrix.loc[1])) * weights[1] + 1) * (
                     np.array(list(matrix.loc[2])) * weights[2] + 1) * (
                     np.array(list(matrix.loc[3])) * weights[3] + 1) * (
                     np.array(list(matrix.loc[4])) * weights[4] + 1) * (
                     np.array(list(matrix.loc[5])) * weights[5] + 1) * (
                     np.array(list(matrix.loc[6])) * weights[6] + 1) * (
                     np.array(list(matrix.loc[7])) * weights[7] + 1) + weights[8]
    """
    return vector


class TopoCMap:
    """
    A class is used to perform Connectivity Map

    Attributes:
        up_request (:obj:`list` of :obj:`str`): list of up regulated genes of request signature
        down_request (:obj:`list` of :obj:`str`): list of down regulated genes of request signature
        up_FC (:obj:`list` of :obj:`float`): list of logFC of up regulated genes of request signature
        down_FC (:obj:`list` of :obj:`str`): list of logFC of down regulated genes of request signature
        db (:obj:`list` of :obj:`list` of :obj:`str`): list of list of genes from database
        mode (str): search mode; either "reverse" or "mimic"
        metrics_up (:obj: `pandas.DataFrame`): matrix of centrality metrics for up genes
        metrics_down (:obj: `pandas.DataFrame`): matrix of centrality metrics for down genes
        inf_score_up (:obj:`list` of :obj:`float`): list of influence scores for up genes
        inf_score_down (:obj:`list` of :obj:`float`): list of influence scores for down genes

    """

    def __init__(self, db, mode, file):
        """
        Declares class attributes

        Args:
            db (:obj:`list` of :obj:`list` of :obj:`str`): list of list of genes from database
            mode (str): search mode; either "reverse" or "mimic"
            file (str): path to the file with signatures from database

        """
        self.up_request, self.down_request, self.up_FC, self.down_FC = signature_extractor(file)
        self.db = db
        self.mode = mode
        self.metrics_up = centrality_metrics(0.1, self.up_request, self.up_FC).metric_calculator()
        self.metrics_down = centrality_metrics(0.1, self.down_request, self.down_FC).metric_calculator()
        self.inf_score_up = None
        self.inf_score_down = None

    def space_finding(self, db_up, db_down):
        """
        Finds gene spaces based on a signature from database and request signature

        Args:
            db_up (:obj:`list` of :obj:`str`): list of up genes from database signature
            db_down (:obj:`list` of :obj:`str`): list of down genes from database signature

        Returns: 
            ord_dict_1 (dict): Space constructed from DOWN genes from request signature and down or up genes,
                    in case of mimic or reverse mode, from database signature respectively 
            ord_dict_2 (dict): Space constructed from UP genes from request signature and up or down genes,
                    in case of mimic or reverse mode, from database signature respectively

        """
        set_1 = list(self.down_request)
        set_2 = list(self.up_request)
        if self.mode == "reverse":
            set_1.extend(list(set(db_up).difference(set(self.down_request))))
            set_2.extend(list(set(db_down).difference(set(self.up_request))))
            dict_1 = defaultdict(list)
            dict_2 = defaultdict(list)
            for gene in self.down_request:
                dict_1[gene].append(0)

            for gene in self.up_request:
                dict_2[gene].append(0)

            for gene in db_up:
                dict_1[gene].append(1)

            for gene in db_down:
                dict_2[gene].append(1)
            ord_dict_1 = {key: dict_1[key] for key in set_1}
            ord_dict_2 = {key: dict_2[key] for key in set_2}
        else:
            set_1.extend(list(set(db_down).difference(set(self.down_request))))
            set_2.extend(list(set(db_up).difference(set(self.up_request))))
            dict_1 = defaultdict(list)
            dict_2 = defaultdict(list)
            for gene in self.down_request:
                dict_1[gene].append(0)

            for gene in self.up_request:
                dict_2[gene].append(0)

            for gene in db_up:
                dict_2[gene].append(1)

            for gene in db_down:
                dict_1[gene].append(1)
            ord_dict_1 = {key: dict_1[key] for key in set_1}
            ord_dict_2 = {key: dict_2[key] for key in set_2}

        return ord_dict_1, ord_dict_2

    def influence_score(self, weights):
        """
        Calls functions for influence score calculations

        Args:
            weights (:obj:`list` of :obj:`float`): List of weights for the expression of influence scores

        """
        inf_scores_up = []
        inf_scores_down = []
        inf_scores_up.append(self.metrics_up.loc[0])
        inf_scores_up.append(list(inf_score_expression(self.metrics_up, weights)))
        inf_scores_down.append(self.metrics_down.loc[0])
        inf_scores_down.append(list(inf_score_expression(self.metrics_down, weights)))
        self.inf_score_up = inf_scores_up
        self.inf_score_down = inf_scores_down

    @staticmethod
    def current_scores(space, inf_score):
        """
        Rearranges list of influence scores for each space

        Args:
            space (dict): Current space based on current signature from database with either up or down genes
            inf_score (:obj:`list` of :obj:`float`): Influence scores list for either up or down genes

        Returns:
            output (:obj:`list` of :obj:`float`): Ordered list of influence scores according to the current space

        """
        inf_score[0] = list(inf_score[0])
        output = inf_score[1].copy()
        output.extend(np.ones(len(space) - len(inf_score[0])))
        return output

    def cmap(self, db_up, db_down):
        """
        Calculates cosine distances for the request signature with a signature from database

        Args:
            db_up (:obj:`list` of :obj:`str`): Up genes for the current database signature
            db_down (:obj:`list` of :obj:`str`): Down genes for the current database signature

        Returns:
            cosine_dist (float): Calculated cosine distance between signatures

        """
        set_1, set_2 = self.space_finding(db_up, db_down)
        # data = [(self.up_request, set_2), (db_down, set_2), (self.down_request, set_1), (db_up, set_1)]
        # data = [(self.up_request, space), (db_down, space), (self.down_request, space), (db_up, space)]
        spaces = [set_2, set_1]
        results = decomposition(spaces)
        scores = []
        data = [self.inf_score_up, self.inf_score_down]
        for ind in range(len(spaces)):
            scores.append(self.current_scores(spaces[ind], data[ind]))
        cos1 = distance.cosine(results[2], results[3], np.nan_to_num(scores[1], nan=1.0))
        cos2 = distance.cosine(results[0], results[1], np.nan_to_num(scores[0], nan=1.0))
        cosine_dist = 0.5 * (cos1 + cos2)

        return cosine_dist

    def small_molec(self, sig_names, file_meta, file_drug, weights):
        """
        Calculates cosine distances across all database signatures, extracts metadata for each molecule and ranges
                molecules by their cosine distance

        Args:
            sig_names (:obj:`list` of :obj:`str`): names of signatures
            file_meta (str): Path to the file with metadata
            file_drug (str): Path to the file with info about small molecules
            weights (:obj:`list` of :obj:`float`): List of weights for the expression of influence scores

        Returns:
            output (:obj:`pandas.DataFrame`): Table of molecules ranged by their cosine distance with some additional
                metadata such as signature name, pert_id, pert_name

        """
        self.influence_score(weights)
        cos_dist = []
        start = time.time()
        for i in tqdm.tqdm(range(0, len(self.db), 2)):
            cos_dist.append(self.cmap(self.db[i], self.db[i + 1]))
        end = time.time()
        print(end - start)
        output_data = pd.DataFrame([sig_names, cos_dist])
        metadata = pd.read_csv(file_meta)
        drugs = pd.read_csv(file_drug)
        pert_ids = []
        for sig in output_data.loc[0]:
            for ind, pert in enumerate(metadata['sig_id']):
                if sig == pert:
                    pert_ids.append(metadata["pert_id"].loc[ind])
        output_data = output_data.append(pd.Series(pert_ids), ignore_index=True)
        chems = []
        for pert in pert_ids:
            for ind, chem in enumerate(drugs['pert_id']):
                if pert == chem:
                    chems.append(drugs['pubchem_cid'].loc[ind])
        output_data = output_data.append(pd.Series(chems), ignore_index=True)
        chem_name = []
        for pert in pert_ids:
            for ind, chem in enumerate(drugs['pert_id']):
                if pert == chem:
                    chem_name.append(drugs['pert_iname'].loc[ind])
        output_data = output_data.append(pd.Series(chem_name), ignore_index=True)
        return output_data


# Lines for validation on CMap signatures
"""
coeffs = [7., 4.15797256, 2.12928013, 3.45806033, 8.71725902, 5.61280318, 1.21182083, 8.50291405, 0.67760011]
bd_signatures, sig_names = BD_signature_parser("/Users/littlequeen/Downloads/CD_signatures_binary_42809.gmt")
trial = TopoCMap(bd_signatures, "mimic", ["1132939s_sigS11_up.txt", "1132939s_sigS11_down.txt"])
small_mol = trial.small_molec(sig_names, "~/Downloads/CD_signature_metadata.csv", "~/Downloads/Drugs_metadata.csv",
                              weights=coeffs)
small_mol = pd.DataFrame(small_mol)
small_mol = small_mol.T
small_mol.columns = ["sign_id", "cosine_dist", "pert_id", "pubchem_id", "pert_iname"]
small_mol = small_mol.drop_duplicates()
small_mol = small_mol.sort_values(by=["cosine_dist"], ascending=[True])
small_mol.to_csv("1132939s_sigS11_mimic" + '.csv', index=False)
"""
