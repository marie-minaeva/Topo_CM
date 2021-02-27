import time

from joblib import Parallel, delayed
from scipy.spatial import distance

from BD_signature_parser import *
from metric_calculator import *
from cmap import decomposition
import numpy as np
from gtf_parser import *

class TopoCMap:

    def __init__(self, up_request, down_request, db, mode, file):
        self.up_request = up_request
        self.down_request = down_request
        self.db = db
        self.mode = mode
        self.metrics_up = centrality_metrics(file, mode="up_genes", score=0.1).metric_calculator
        self.metrics_down = centrality_metrics(file, mode="down_genes", score=0.1).metric_calculator
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

        set_1 = dict.fromkeys(list(set(self.down_request + db_up)), [])
        set_2 = dict.fromkeys(list(set(self.up_request + db_down)), [])
        print(set_1)
        print(set_2)
        for gene in self.down_request:
            set_1[gene].append(0)

        for gene in self.up_request:
            set_2[gene].append(0)

        for gene in db_up:
            set_1[gene].append(1)

        for gene in db_down:
            set_2[gene].append(1)


        return set_1, set_2

    def influence_score(self, weights):
        inf_scores_up = []
        inf_scores_down = []
        inf_scores_up.append(self.metrics_up.loc[0])
        """
        for i in range(1, 7):
            print(np.std(list(self.metrics_up.loc[i])))
        """
        vector = (np.array(list(self.metrics_up.loc[1]))*weights[0] + 1) * (np.array(list(self.metrics_up.loc[2]))*weights[1] + 1) * (np.array(list(self.metrics_up.loc[3]))*weights[2] + 1) * (np.array(list(self.metrics_up.loc[4]))*weights[3] + 1) * (np.array(list(self.metrics_up.loc[5]))*weights[4] + 1)*(np.array(list(self.metrics_up.loc[6]))*weights[5] + 1) + weights[6]
        inf_scores_up.append(list(vector))
        inf_scores_down.append(self.metrics_down.loc[0])
        vector = (np.array(list(self.metrics_down.loc[1]))*weights[0] + 1) * (np.array(list(self.metrics_down.loc[2]))*weights[1] + 1) * (np.array(list(self.metrics_down.loc[3]))*weights[2] + 1) * (np.array(list(self.metrics_down.loc[4]))*weights[3]  + 1) * (np.array(list(self.metrics_down.loc[5]))*weights[4] + 1)*(np.array(list(self.metrics_down.loc[6]))*weights[5] + 1) + weights[6]
        inf_scores_down.append(list(vector))
        self.inf_score_up = inf_scores_up
        self.inf_score_down = inf_scores_down


    def current_scores(self, spaces):
        ## говорящие названия и описания классов
        inf_scores = [self.inf_score_up, self.inf_score_down]
        #inf_scores = self.inf_score_up + self.inf_score_down
        output = []

        for  ind, inf_score in enumerate(inf_scores):
            output.append([])
            for sp_gene in spaces[ind]:
                if sp_gene in inf_score[0]:
                    output[ind].append(inf_score[1][inf_score[0].index(sp_gene)])
                else:
                    output[ind].append(1.0)
        """
        for sp_gene in spaces:
            if sp_gene in inf_scores[0]:
                output.append(inf_scores[1][inf_scores[0].index(sp_gene)])
            else :
                output.append(1.0)
        """
        return output

    def cmap(self, db_up, db_down):

        set_1, set_2 = self.space_finding(db_up, db_down)
        #data = [(self.up_request, set_2), (db_down, set_2), (self.down_request, set_1), (db_up, set_1)]
        #data = [(self.up_request, space), (db_down, space), (self.down_request, space), (db_up, space)]
        data = [set_2, set_1]
        results = decomposition(data)
        spaces = [set_2, set_1]
        scores = self.current_scores(spaces)
        cos1 = distance.cosine(results[1], results[3], scores[1])
        cos2 = distance.cosine(results[5], results[7], scores[0])
        cosine_dist = 0.5 * (cos1 + cos2)
        print(cosine_dist)

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
