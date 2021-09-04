import pandas as pd


def stand_chems(source, target):
    """
    Extracts compounds of 'Golden Standard' for a certain cellular conversion
    Args:
        source (str): Source cell type
        target (str): Target cell type

    Returns:
        output (:obj:`list` of `float`): List of CIDs (PubChem identifiers) for compounds of 'Golden Standard' for a certain cellular conversion

    """
    data = pd.read_csv("~/Downloads/direct_reprogramming_non-genetics-structure.csv")
    # data = data[data["Species"] == "Homo sapiens"]
    data = data[data['Source Cell Type'] == source]
    data = data[data['Target Cell Type'] == target]
    cids = list(data['name of chemical 1,CID 1;name of chemical 2,CID 2'])
    output = []
    for cid in cids:
        line = cid.split(';')
        for word in line:
            word = word.split(',')
            output.append(word[1])
    if output[-1] == "_":
        output = output[:-1]
    meta = pd.read_csv("~/Downloads/intersection_of_CFM_and_L1000FWD.csv")
    cids_cfm = list(meta["CID"])
    cids_cfm = [str(i) for i in cids_cfm]
    output = list(set(output).intersection(set(cids_cfm)))
    print(len(output))
    return output
