import pandas as pd
def stand_chems(source, target):
    data = pd.read_csv("~/Downloads/direct_reprogramming_non-genetics_structure.csv")
    data = data[data["Species"] == "Homo sapiens"]
    data = data[data['Source Cell Type'] == source]
    data = data[data['Target Cell Type'] == target]
    cids = list(data['name of chemical 1,CID 1;name of chemical 2,CID 2'])
    output = []
    for cid in cids:
        line = cid.split(';')
        for word in line:
            word = word.split(',')
            output.append(word[1])

    print(output)
    return output
