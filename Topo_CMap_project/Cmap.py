import time
from collections import *
from joblib import Parallel, delayed
from scipy.spatial import distance

from BD_signature_parser import *
from metric_calculator import *
#from cmap import decomposition
import numpy as np
#from gtf_parser import *
from signature_extractor import *
from new_signature_extractor import *


def decomposition(data):
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
    vector = (np.array(list(matrix.loc[1])) * weights[0] + 1) * (
                np.array(list(matrix.loc[2])) * weights[1] + 1) * (
                         np.array(list(matrix.loc[3])) * weights[2] + 1) * (
                         np.array(list(matrix.loc[4])) * weights[3] + 1) * (
                         np.array(list(matrix.loc[5])) * weights[4] + 1) * (
                         np.array(list(matrix.loc[6])) * weights[5] + 1) * (
                         np.array(list(matrix.loc[7])) * weights[6] + 1) * (
                         np.array(list(matrix.loc[8])) * weights[7] + 1) + weights[8]
    """
    vector = (np.array(list(matrix.loc[1])) * weights[1] + 1) * (
                     np.array(list(matrix.loc[2])) * weights[2] + 1) * (
                     np.array(list(matrix.loc[3])) * weights[3] + 1) * (
                     np.array(list(matrix.loc[4])) * weights[4] + 1) * (
                     np.array(list(matrix.loc[5])) * weights[5] + 1) * (
                     np.array(list(matrix.loc[6])) * weights[6] + 1) * (
                     np.array(list(matrix.loc[7])) * weights[7] + 1) + weights[8]

    return vector


class TopoCMap:

    def __init__(self, db, mode, file):
        self.up_request, self.down_request, self.up_FC, self.down_FC = new_signature_extractor(file[0], file[1])
        self.db = db
        self.mode = mode
        self.metrics_up = centrality_metrics(0.1, self.up_request, self.up_FC).metric_calculator()
        self.metrics_down = centrality_metrics(0.1, self.down_request, self.down_FC).metric_calculator()
        self.inf_score_up = None
        self.inf_score_down = None

    def space_finding(self, db_up, db_down):

        #if self.mode == "reverse":

        #set_1 = list(set(self.down_request + db_up))
        #set_2 = list(set(self.up_request + db_down))
            ## проверять что set1 не пересекается с set2
        #else:
            #set_1 = list(set(self.down_request + db_down))
            #set_2 = list(set(self.up_request + db_up))
        set_1 = list(self.down_request)
        set_2 = list(self.up_request)
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
        return ord_dict_1, ord_dict_2

    def influence_score(self, weights):
        inf_scores_up = []
        inf_scores_down = []
        inf_scores_up.append(self.metrics_up.loc[0])
        """
        for i in range(1, 7):
            print(np.std(list(self.metrics_up.loc[i])))
        """
        inf_scores_up.append(list(inf_score_expression(self.metrics_up, weights)))
        inf_scores_down.append(self.metrics_down.loc[0])
        inf_scores_down.append(list(inf_score_expression(self.metrics_down, weights)))
        self.inf_score_up = inf_scores_up
        self.inf_score_down = inf_scores_down

    @staticmethod
    def current_scores(space, inf_score):
        inf_score[0] = list(inf_score[0])
        output = inf_score[1].copy()
        output.extend(np.ones(len(space) - len(inf_score[0])))
        return output

    def cmap(self, db_up, db_down):

        set_1, set_2 = self.space_finding(db_up, db_down)
        #data = [(self.up_request, set_2), (db_down, set_2), (self.down_request, set_1), (db_up, set_1)]
        #data = [(self.up_request, space), (db_down, space), (self.down_request, space), (db_up, space)]
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
        self.influence_score(weights)
        cos_dist = []
        start = time.time()
        #space = gtf_parser()
        #scores = self.current_scores(space)
        for i in tqdm.tqdm(range(0, len(self.db), 2)):
            cos_dist.append(self.cmap(self.db[i], self.db[i+1]))
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

## отнормировать коэффициенты и посмотреть распределения коэффициентов для метрик
## баесовская оптимизация

