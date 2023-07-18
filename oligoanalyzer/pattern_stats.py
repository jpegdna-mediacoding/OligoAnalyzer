import argparse
import numpy as np
import matplotlib.pyplot as plt
from tools.fasta import get_oligos

PATTERN_SIZE_RANGE = (3,7)
LEN_PATTERN_SIZE_RANGE = PATTERN_SIZE_RANGE[1] - PATTERN_SIZE_RANGE[0] + 1

def compute_data(data):
    """Compute simple pattern statistics on the oligo pool"""
    pattern_counts = []
    pattern_infos = []
    pattern_counts_2 = []
    pattern_infos_2 = []
    for i in range(len(data)):
        oligo = data[i]
        res = detect_patterns_old(oligo)
        pattern_counts.append(res[0])
        pattern_infos.append(res[1])
        res2 = detect_patterns(oligo)
        pattern_counts_2.append(res2[0])
        pattern_infos_2.append(res2[1])
    return pattern_counts, pattern_counts_2

def detect_patterns_old(oligo):
    """Detect the patterns in a specific oligo"""
    pattern_counts, pattern_info = [0, 0, 0, 0, 0], {}
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
                pattern_counts[length-PATTERN_SIZE_RANGE[0]] += 3
                mask[j:j+3*length] = [1]*3*length
                j += 3*length
                while (j+length < len(oligo) and oligo[j:j+length] == pattern and mask[j:j+length] != [1]*length):
                    new_info = {'position': {(j, j+length)}, 'population': 1}
                    update_pattern_info(pattern_info[pattern], new_info)
                    pattern_counts[length-PATTERN_SIZE_RANGE[0]] += 1
                    mask[j:j+length] = [1]*length
                    j += length
            else:
                j+=1
    return pattern_counts, pattern_info

def detect_patterns(oligo):
    """Detect the patterns in a specific oligo"""
    pattern_counts, pattern_info = [0, 0, 0, 0, 0], {}
    mask = [0] * len(oligo)
    for length in range(PATTERN_SIZE_RANGE[1], PATTERN_SIZE_RANGE[0]-1, -1):
        j = 0
        rep = 4 if (length >= 6) else 3
        while j < len(oligo)-rep*length:
            j_bis, j_ter = j+length, j+2*length
            if (rep==2 and
                oligo[j:j+length] == oligo[j_bis:j_bis+length] and
                mask[j:j+rep*length] != [1]*rep*length):
                pattern = oligo[j:j+length]
                new_info = {'position': {(j,j+length), (j_bis, j_bis+length)}, 'population': rep}
                if pattern in pattern_info.keys():
                    update_pattern_info(pattern_info[pattern], new_info)
                else:
                    pattern_info[pattern] = new_info
                pattern_counts[length-PATTERN_SIZE_RANGE[0]] += rep
                mask[j:j+rep*length] = [1]*rep*length
                j += rep*length
                while (j+length < len(oligo) and oligo[j:j+length] == pattern and mask[j:j+length] != [1]*length):
                    new_info = {'position': {(j, j+length)}, 'population': 1}
                    update_pattern_info(pattern_info[pattern], new_info)
                    pattern_counts[length-PATTERN_SIZE_RANGE[0]] += 1
                    mask[j:j+length] = [1]*length
                    j += length
            elif (rep==3 and
                oligo[j:j+length] == oligo[j_bis:j_bis+length] == oligo[j_ter:j_ter+length] and
                mask[j:j+rep*length] != [1]*rep*length):
                pattern = oligo[j:j+length]
                new_info = {'position': {(j,j+length), (j_bis, j_bis+length), (j_ter, j_ter+length)}, 'population': rep}
                if pattern in pattern_info.keys():
                    update_pattern_info(pattern_info[pattern], new_info)
                else:
                    pattern_info[pattern] = new_info
                pattern_counts[length-PATTERN_SIZE_RANGE[0]] += rep
                mask[j:j+rep*length] = [1]*rep*length
                j += rep*length
                while (j+length < len(oligo) and oligo[j:j+length] == pattern and mask[j:j+length] != [1]*length):
                    new_info = {'position': {(j, j+length)}, 'population': 1}
                    update_pattern_info(pattern_info[pattern], new_info)
                    pattern_counts[length-PATTERN_SIZE_RANGE[0]] += 1
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

def plot(stats):
    """Plot the stats on the data"""
    plt.figure()
    data = [[sum(el) for el in stats[0] if sum(el)!=0]]
    data2 = []
    legend = ["General"]
    for i in range(LEN_PATTERN_SIZE_RANGE):
        data.append([el[i] for el in stats[0] if el[i]!=0])
        data2.append(sum(data[-1]))
        legend.append(f"length={PATTERN_SIZE_RANGE[0]+i}")
    plt.subplot(221)
    plt.hist(data, bins=np.arange(3, max(data[0]), 1), label=legend)
    plt.legend(loc="upper right")
    plt.title("Nb of patterns per oligo")
    plt.xlabel("Nb of patterns")
    plt.ylabel("Nb of oligos")
    plt.subplot(222)
    plt.bar([str(el2) for el2 in range(PATTERN_SIZE_RANGE[0], PATTERN_SIZE_RANGE[1]+1)], data2)
    plt.title("Nb of patterns per size")
    plt.xlabel("Pattern Size")
    plt.ylabel("Nb of oligos")
    data = [[sum(el) for el in stats[1] if sum(el)!=0]]
    data2 = []
    legend = ["General"]
    for i in range(LEN_PATTERN_SIZE_RANGE):
        data.append([el[i] for el in stats[1] if el[i]!=0])
        data2.append(sum(data[-1]))
        legend.append(f"length={PATTERN_SIZE_RANGE[0]+i}")
    plt.subplot(223)
    plt.hist(data, bins=np.arange(3, max(data[0]), 1), label=legend)
    plt.legend(loc="upper right")
    plt.title("Nb of patterns per oligo")
    plt.xlabel("Nb of patterns")
    plt.ylabel("Nb of oligos")
    plt.subplot(224)
    plt.bar([str(el2) for el2 in range(PATTERN_SIZE_RANGE[0], PATTERN_SIZE_RANGE[1]+1)], data2)
    plt.title("Nb of patterns per size")
    plt.xlabel("Pattern Size")
    plt.ylabel("Nb of oligos")
    plt.show()

def main():
    """Plots statistics on the patterns over all the oligos"""
    parser = argparse.ArgumentParser()
    parser.add_argument('FNAME',
                        type=str,
                        help='Filename of the input fasta file')
    args = parser.parse_args()
    fname = args.FNAME
    data = get_oligos(fname)
    stats = compute_data(data)
    print(f'Number of old patterns: {sum([sum(el) for el in stats[0]])}')
    print(f'Number of patterns: {sum([sum(el) for el in stats[1]])}')
    plot(stats)

if __name__ == '__main__':
    main()
