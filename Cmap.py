from metric_calculator import *
from BD_signature_parser import *
from scipy.spatial import distance
import time

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
        else:
            set_1 = list(set(self.down_request + db_down))
            set_2 = list(set(self.up_request + db_up))

        return set_1, set_2

    def influence_score(self):

        metrics_up = centrality_metrics(args.infile, mode="up_genes", score=0.1).metric_calculator
        metrics_down = centrality_metrics(args.infile, mode="down_genes", score=0.1).metric_calculator
        inf_scores_up = []
        inf_scores_down = []
        inf_scores_up.append(metrics_up.loc[0])
        inf_scores_up.append([x[0] * x[1] for x in zip(metrics_up.loc[2], metrics_up.loc[3])])
        inf_scores_down.append(metrics_down.loc[0])
        inf_scores_down.append([x[0] * x[1] for x in zip(metrics_down.loc[2], metrics_down.loc[3])])

        self.inf_score_up = inf_scores_up
        self.inf_score_down = inf_scores_down

    @staticmethod
    def decomposition(genes, space, scores):

        vector = []
        for sp_gene in space:
            flag = 0
            for ind, gene in enumerate(genes):
                flag1 = 0
                if sp_gene == gene:
                    try:
                        for ind_i, score in enumerate(scores[0]):
                            if score == gene:
                                if scores[1][ind_i] != 0:
                                    vector.append(scores[1][ind_i])
                                else:
                                    vector.append(1.0)
                            else:
                                flag1 += 1
                    except LookupError:
                        vector.append(1.0)
                else:
                    flag += 1
                try:
                    if flag1 == len(scores[0]):
                        vector.append(1.0)
                except LookupError:
                    continue
            if flag == len(genes):
                vector.append(0.0)

        return vector

    def cmap(self):
        cosine_dist = []
        self.influence_score()
        for i in range(0, len(self.db), 2):
            start1 = time.time()
            set_1, set_2 = self.space_finding(self.db[i], self.db[i+1])
            point1 = time.time()
            print(point1 - start1)
            vector_up_req = self.decomposition(self.up_request, set_1, self.inf_score_up)
            vector_down_db = self.decomposition(self.db[i+1], set_1, [])
            vector_down_req = self.decomposition(self.down_request, set_2, self.inf_score_down)
            vector_up_db = self.decomposition(self.db[i], set_2, [])
            point2 = time.time()
            print(point2 - point1)
            cos1 = distance.cosine(vector_up_req, vector_down_db)
            cos2 = distance.cosine(vector_down_req, vector_up_db)
            end1 = time.time()
            cosine_dist.append(0.5*(cos1+cos2))
            print(i, ' ', 0.5 * (cos1 + cos2), ' ', end1 - start1)

        return cosine_dist

    def small_molec(self, sig_names, file_meta, file_drug):
        cos_dist = self.cmap()
        print(sig_names[50:52])
        output_data = pd.DataFrame([sig_names[50:52], cos_dist])
        metadata = pd.read_csv(file_meta)
        drugs = pd.read_csv(file_drug)
        output_data.columns = ['sig_name', 'cos_dist']
        pert_ids = []
        for sig in output_data["sig_name"]:
            for ind, pert in enumerate(metadata['sig_id']):
                if sig == pert:
                    pert_ids.append(metadata["sig_id"].loc[ind])
        output_data["pert_id"] = pert_ids
        chems = []
        for pert in output_data["pert_id"]:
            for ind, chem in enumerate(drugs['pert_id']):
                if pert == chem:
                    chems.append(metadata['pubchem_cid'].loc[ind])
        output_data["pubchem_cid"] = chems

        return output_data

start = time.time()
bd_signatures, sig_names = BD_signature_parser("/Users/littlequeen/Downloads/CD_signatures_binary_42809.gmt")
print(len(bd_signatures))
up_request, down_request = signature_extractor("~/Downloads/DE_heart_fb_deseq2_edger.txt")
trial = TopoCMap(list(up_request), list(down_request), bd_signatures[100:104], 'reverse')
small_mol= trial.small_molec(sig_names, "~/Downloads/CD_signature_metadata.csv", "~/Downloads/Drugs_metadata.csv")
end = time.time()
print(end - start)
print(small_mol)

