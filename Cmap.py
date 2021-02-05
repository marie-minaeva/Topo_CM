import time

from joblib import Parallel, delayed
from scipy.spatial import distance

from BD_signature_parser import *
from metric_calculator import *
from cmap import decomposition

class TopoCMap:

    def __init__(self, up_request, down_request, db, mode):
        self.up_request = up_request
        self.down_request = down_request
        self.db = db
        self.mode = mode
        self.inf_score_up = None
        self.inf_score_down = None

    def space_finding(self, db_up, db_down):

        if self.mode == "reverse":
            set_1 = list(set(self.down_request + db_up))
            set_2 = list(set(self.up_request + db_down))
            ## проверять что set1 не пересекается с set2
        else:
            set_1 = list(set(self.down_request + db_down))
            set_2 = list(set(self.up_request + db_up))

        return set_1, set_2

    def influence_score(self):
        inf_scores_up = []
        inf_scores_down = []
        metrics_up = centrality_metrics(args.infile, mode="up_genes", score=0.1).metric_calculator
        metrics_down = centrality_metrics(args.infile, mode="down_genes", score=0.1).metric_calculator
        inf_scores_up.append(metrics_up.loc[0])
        inf_scores_up.append([(x[0] * x[1]+1) for x in zip(metrics_up.loc[2], metrics_up.loc[3])])
        inf_scores_down.append(metrics_down.loc[0])
        ## добавить сюда FC
        inf_scores_down.append([(x[0] * x[1]+1) for x in zip(metrics_down.loc[2], metrics_down.loc[3])])

        self.inf_score_up = inf_scores_up
        self.inf_score_down = inf_scores_down

    ## добавить в граф гены не из STRING и тогда не надо будет делать current_scores

    def current_scores(self, space):
        ## говорящие названия и описания классов
        inf_scores = self.inf_score_up + self.inf_score_down
        output = []
        for sp_gene in space:
            if sp_gene in inf_scores[0]:
                output.append(inf_scores[1][inf_scores[0].index(sp_gene)])
            else:
                output.append(1.0)

        return output

    """@staticmethod
    def decomposition(data):
        vector = []
        out_genes = []
        out_vector = []
        for sp_gene in data[1]:
            if sp_gene in data[0]:
                out_genes.append(sp_gene)
                out_vector.append(1.0)
            else:
                out_genes.append(sp_gene)
                out_vector.append(0.0)
        vector.append(out_genes)
        vector.append(out_vector)

        return vector
    """
    def cmap(self, db_up, db_down):
        start1 = time.time()
        set_1, set_2 = self.space_finding(db_up, db_down)
        data = [(self.up_request, set_2), (db_down, set_2), (self.down_request, set_1), (db_up, set_1)]
        results = decomposition(data)
        point1 = time.time()
        print(point1 - start1)
        data1 = [set_1, set_2]
        pool1 = Parallel(n_jobs=1, pre_dispatch='all')
        results1 = pool1(delayed(self.current_scores)(dd) for dd in data1)
        point2 = time.time()
        print(point2 - point1)
        print(len(results[1]), len(results1[1]))
        cos1 = distance.cosine(results[1], results[3], results1[1])
        cos2 = distance.cosine(results[5], results[7], results1[0])
        end1 = time.time()
        print(end1 - point2)
        cosine_dist = 0.5 * (cos1 + cos2)
        print(0.5 * (cos1 + cos2), ' ', end1 - start1)

        return cosine_dist

    def small_molec(self, sig_names, file_meta, file_drug):
        self.influence_score()
        cos_dist = []
        start = time.time()
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


start = time.time()
bd_signatures, sig_names = BD_signature_parser("/Users/littlequeen/Downloads/CD_signatures_binary_42809.gmt")
up_request, down_request = signature_extractor("~/Downloads/DE_heart_fb_deseq2_edger.txt")
trial = TopoCMap(list(up_request), list(down_request), bd_signatures[:200], 'reverse')
small_mol = trial.small_molec(sig_names, "~/Downloads/CD_signature_metadata.csv", "~/Downloads/Drugs_metadata.csv")
end = time.time()
print(end - start)
small_mol = pd.DataFrame(small_mol)
small_mol = small_mol.T
small_mol.columns = ["sign_id", "cosine_dist", "pert_id", "pubchem_id"]
small_mol.to_csv('~/Downloads/small_mol.csv', index=False)
