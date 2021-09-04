def BD_signature_parser(file):
    with open(file) as gmt:
        lines = gmt.readlines()
    signatures = []
    sig_names = []
    for line in lines:
        signatures.append(line.split()[2:])
        sig_names.append(line.split()[0])
    print(len(signatures))
    sig_names = list(dict.fromkeys(sig_names))
    return signatures, sig_names
