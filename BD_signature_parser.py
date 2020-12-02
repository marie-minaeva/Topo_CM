def BD_signature_parser(file):
    with open(file) as gmt:
        lines = gmt.readlines()
    signatures = []
    sig_names = []
    for line in lines:
        signatures.append(line.split()[2:])
        sig_names.append(line.split()[0])

    sig_names = list(set(sig_names))
    return signatures, sig_names
