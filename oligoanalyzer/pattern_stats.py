import argparse
from tools.fasta import get_oligos

PATTERN_SIZE_RANGE = (3,5)

def compute_data(data):
    """Compute simple pattern statistics on the oligo pool"""
    pattern_counts = []
    pattern_infos = []
    for i in range(len(data)):
        oligo = data[i]
        res = detect_patterns(oligo)
        pattern_counts.append(res[0])
        pattern_infos.append(res[1])
    return sum(pattern_counts)

def detect_patterns(oligo):
    """Detect the patterns in a specific oligo"""
    pattern_counts, pattern_info = 0, {}
    mask = [0] * len(oligo)
    for length in range(PATTERN_SIZE_RANGE[1], PATTERN_SIZE_RANGE[0]-1, -1):
        j = 0
        while j < len(oligo)-3*length:
            j_bis, j_ter = j+length, j+2*length
            if oligo[j:j+length] == oligo[j_bis:j_bis+length] == oligo[j_ter:j_ter+length] and mask[j:j+3*length] != [1]*3*length:
                pattern = oligo[j:j+length]
                new_info = {'position': {(j,j+length), (j_bis, j_bis+length), (j_ter, j_ter+length)}, 'population': 3}
                if pattern in pattern_info.keys():
                    update_pattern_info(pattern_info[pattern], new_info)
                else:
                    pattern_info[pattern] = new_info
                pattern_counts += 3
                mask[j:j+3*length] = [1]*3*length
                j += 3*length
                while (j+length < len(oligo) and oligo[j:j+length] == pattern and mask[j:j+length]):
                    new_info = {'position': {(j, j+length)}, 'population': 1}
                    update_pattern_info(pattern_info[pattern], new_info)
                    pattern_counts += 1
                    mask[j:j+length] = [1]*length
                    j += length
            else:
                j+=1
    return pattern_counts, pattern_info

def update_pattern_info(info, new_info):
    """Update the info object"""
    info['position'].update(new_info['position'])
    info['population'] += new_info['population']
    return info

def main():
    """Plots statistics on the patterns over all the oligos"""
    parser = argparse.ArgumentParser()
    parser.add_argument('FNAME',
                        type=str,
                        help='Filename of the input fasta file')
    args = parser.parse_args()
    fname = args.FNAME
    data = get_oligos(fname)
    n_patterns = compute_data(data)
    print(f'Number of patterns: {n_patterns}')

if __name__ == '__main__':
    main()
