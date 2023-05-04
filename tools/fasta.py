def get_oligos(fasta_file):
    """Gets the oligo data from a fasta file"""
    headers, data = [], []
    with open(fasta_file, "r", encoding='utf-8') as f:
        tmp_data = f.readlines()
    for line in tmp_data:
        if line[0] in ["A", "T", "C", "G"]:
            if line[-1] == "\n":
                data.append(line[:-1])
            else:
                data.append(line)
        elif line[0] == ">":
            headers.append(line[:-1])
    return data
