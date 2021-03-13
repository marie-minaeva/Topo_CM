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
    vector = (np.array(list(matrix.loc[1])) * weights[0] + 1) * (
                np.array(list(matrix.loc[2])) * weights[1] + 1) * (
                         np.array(list(matrix.loc[3])) * weights[2] + 1) * (
                         np.array(list(matrix.loc[4])) * weights[3] + 1) * (
                         np.array(list(matrix.loc[5])) * weights[4] + 1) * (
                         np.array(list(matrix.loc[6])) * weights[5] + 1) * (
                         np.array(list(matrix.loc[7])) * weights[6] + 1) * (
                         np.array(list(matrix.loc[8])) * weights[7] + 1) + weights[8]
    return vector


class TopoCMap:

    def __init__(self, db, mode, file):
        self.up_request, self.down_request, self.up_FC, self.down_FC = signature_extractor(file)
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

        set_1 = defaultdict(list)
        set_2 = defaultdict(list)
        for gene in self.down_request:
            set_1[gene].append(0)

        for gene in self.up_request:
            set_2[gene].append(0)

        for gene in db_up:
            set_1[gene].append(1)

        for gene in db_down:
            set_2[gene].append(1)
        #print(len(set_1), len(set_2))

        return set_1, set_2

    def influence_score(self, weights):
        inf_scores_up = []
        inf_scores_down = []
        print(self.metrics_up)
        inf_scores_up.append(self.metrics_up.loc[0])
        """
        for i in range(1, 7):
            print(np.std(list(self.metrics_up.loc[i])))
        """
        inf_scores_up.append(list(inf_score_expression(self.metrics_up, weights)))
        inf_scores_down.append(self.metrics_down.loc[0])
        inf_scores_down.append(list(inf_score_expression(self.metrics_down, weights)))
        print(inf_scores_up)
        self.inf_score_up = inf_scores_up
        self.inf_score_down = inf_scores_down


    def current_scores(self, spaces, inf_score):
        ## говорящие названия и описания классов
        #inf_scores = self.inf_score_up + self.inf_score_down
        output = []
        inf_score[0] = list(inf_score[0])
        for sp_gene in spaces:
            if sp_gene in inf_score[0]:
                output.append(inf_score[1][inf_score[0].index(sp_gene)])
            else:
                output.append(1.0)
        for ind, val in enumerate(output):
            if np.isnan(val):
                output[ind] = 1.0
        return output

    def cmap(self, db_up, db_down):

        set_1, set_2 = self.space_finding(db_up, db_down)
        #data = [(self.up_request, set_2), (db_down, set_2), (self.down_request, set_1), (db_up, set_1)]
        #data = [(self.up_request, space), (db_down, space), (self.down_request, space), (db_up, space)]
        spaces = [set_2, set_1]
        results = decomposition(spaces)
        scores = []
        data = [self.inf_score_down, self.inf_score_up]
        for ind in range(len(spaces)):
            scores.append(self.current_scores(spaces[ind], data[ind]))
        cos1 = distance.cosine(results[0], results[1], scores[0])
        cos2 = distance.cosine(results[2], results[3], scores[1])
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
        ## добавить pert_name
        return output_data

## отнормировать коэффициенты и посмотреть распределения коэффициентов для метрик
## баесовская оптимизация
