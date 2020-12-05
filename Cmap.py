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
        # inf_scores_up.append(metrics_up.loc[0])
        inf_scores_up.append([x[0] * x[1] for x in zip(metrics_up.loc[2], metrics_up.loc[3])])
        # inf_scores_down.append(metrics_down.loc[0])
        inf_scores_down.append([x[0] * x[1] for x in zip(metrics_down.loc[2], metrics_down.loc[3])])

        self.inf_score_up = inf_scores_up
        self.inf_score_down = inf_scores_down

    @staticmethod
    def current_scores(genes, inf_score, vector):

        output = []
        for ind_i, cur_gene in enumerate(vector[0]):
            flag = 0
            for ind_j, gene in enumerate(genes):
                if gene == cur_gene:
                    try:
                        output.append(inf_score[ind_i])
                        break
                    except LookupError:
                        output.append(1.0)
                        break
                else:
                    flag += 1
            if flag == len(genes):
                output.append(1.0)

        return output

    @staticmethod
    def decomposition_1(genes, space):
        vector = []
        out_genes = []
        out_vector = []
        for sp_gene in space:
            flag = 0
            for ind, gene in enumerate(genes):
                if sp_gene == gene:
                    out_genes.append(gene)
                    out_vector.append(1.0)
                else:
                    flag += 1
            if flag == len(genes):
                out_genes.append(sp_gene)
                out_vector.append(0.0)
        vector.append(out_genes)
        vector.append(out_vector)

        return vector

    def cmap(self):
        cosine_dist = []
        self.influence_score()
        for i in range(0, len(self.db), 2):
            start1 = time.time()
            set_1, set_2 = self.space_finding(self.db[i], self.db[i + 1])
            vector_up_req = self.decomposition_1(self.up_request, set_1)
            vector_down_db = self.decomposition_1(self.db[i + 1], set_1)
            vector_down_req = self.decomposition_1(self.down_request, set_2)
            vector_up_db = self.decomposition_1(self.db[i], set_2)
            current_up_inf_score = self.current_scores(self.up_request, self.inf_score_up, vector_up_req)
            current_down_inf_score = self.current_scores(self.down_request, self.inf_score_down, vector_down_req)
            cos1 = distance.cosine(vector_up_req[1], vector_down_db[1], current_up_inf_score)
            cos2 = distance.cosine(vector_down_req[1], vector_up_db[1], current_down_inf_score)
            end1 = time.time()
            cosine_dist.append(0.5 * (cos1 + cos2))
            print(i, ' ', 0.5 * (cos1 + cos2), ' ', end1 - start1)

        return cosine_dist

    def small_molec(self, sig_names, file_meta, file_drug):
        cos_dist = self.cmap()
        output_data = pd.DataFrame([sig_names[:500], cos_dist])
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

        return output_data


start = time.time()
bd_signatures, sig_names = BD_signature_parser("/Users/littlequeen/Downloads/CD_signatures_binary_42809.gmt")
up_request, down_request = signature_extractor("~/Downloads/DE_heart_fb_deseq2_edger.txt")
trial = TopoCMap(list(up_request), list(down_request), bd_signatures[0:1000], 'reverse')
small_mol = trial.small_molec(sig_names, "~/Downloads/CD_signature_metadata.csv", "~/Downloads/Drugs_metadata.csv")
end = time.time()
print(end - start)
small_mol = pd.DataFrame(small_mol)
small_mol = small_mol.T
small_mol.columns = ["sign_id", "cosine_dist", "pert_id", "pubchem_id"]
small_mol.to_csv('~/Downloads/small_mol.csv', index=False)
