from rpy2.robjects import pandas2ri
import pandas as pd
from standard_chemicals_extraction import stand_chems
from rpy2 import robjects
pandas2ri.activate()


def gsea_metrics(source, target, file, meta_file):
    """
    Calculates GSEA score for validation on CFM

    Args:
        source (str): Source cell type
        target (str): Target cell type
        file (str): Path to the file with TopoCMap results table
        meta_file (str): Path to the file with drugs metadata

    Returns:
        float :Normalized enrichment score

    """
    drug_meta = pd.read_csv(meta_file)
    df = pd.read_csv(file)
    df = df.drop_duplicates()
    print(type(df))
    cids_cur = stand_chems(source, target)
    pert_cur = []
    for ind, chem in enumerate(drug_meta['pubchem_cid']):
        for chem_1 in cids_cur:
            try:
                if int(chem) == int(chem_1):
                    pert_cur.append(drug_meta['pert_id'].loc[ind])
            except ValueError:
                continue
    # Defining the R script and loading the instance in Python
    r = robjects.r
    r['source']('~/Downloads/fgsea-tutorial.R')
    # Loading the function we have defined in R.
    filter_country_function_r = robjects.globalenv['fgsea_analysis']
    # converting it into r object for passring into r function
    df_r = pandas2ri.py2rpy(df)
    pert_cur_r = robjects.vectors.FactorVector(pert_cur)
    # Invoking the R function and getting the result
    df_result_r = filter_country_function_r(df_r, pert_cur_r)
    print(df_result_r)
    # Converting it back to a pandas dataframe.
    return df_result_r["NES"]
