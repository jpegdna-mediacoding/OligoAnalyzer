import argparse
import matplotlib.pyplot as plt


def compute_data(data):
    """Computes the data"""
    hom = []
    nuc_hom = []
    max_len = 4
    for line in data:
        if line[0] != ">":
            tmp_hom = 0
            tmp_nuc = 0
            i = 0
            while i < len(line)-3:
                if line[i] == line[i+1] == line[i+2] == line[i+3]:
                    tmp_size = 0
                    tmp_hom += 1
                    j = i
                    while  i < len(line) and line[j] == line[i]:
                        tmp_size += 1
                        tmp_nuc += 1
                        i += 1
                    max_len = max(max_len, tmp_size)
                else:
                    i+=1
            hom.append(tmp_hom)
            nuc_hom.append(tmp_nuc)
    avg_hom= []
    for i in range(len(hom)):
        if hom[i] != 0:
            avg_hom.append(nuc_hom[i]/hom[i])
        else:
            avg_hom.append(0)
    print(f"Total number of homopolymers: {sum(hom)}")
    print(f"Overall average size of the homopolymers: {sum(nuc_hom)/sum(hom):.3f}")
    print(f"Maximum length of the homopolymers: {max_len}")
    return (hom, nuc_hom, avg_hom)


def plot(stats):
    """Plot the stats on the data"""
    hom_stats, nuc_hom_stats, avg_hom_stats = stats
    fig, axs = plt.subplots(1, 3, sharey=True)
    axs[0].hist([el for el in hom_stats if el != 0], 200)
    axs[1].hist([el for el in nuc_hom_stats if el != 0], 200)
    axs[2].hist([el for el in avg_hom_stats if el != 0], 200)
    axs[0].set_title("Nb of homopolymers per oligo")
    axs[0].set_xlabel("Nb of homopolymers")
    axs[0].set_ylabel("Nb of oligos")
    axs[1].set_title("Nb of nucleotides belgonging to a homopolymer per oligo")
    axs[1].set_xlabel("Nb of nucleotides")
    axs[1].set_ylabel("Nb of oligos")
    axs[2].set_title("Avg homopolymer size per oligo")
    axs[2].set_xlabel("Avg homopolymer size")
    axs[2].set_ylabel("Nb of oligos")
    plt.show()



def main():
    """Plots the distribution of the GC content over all the oligos"""
    parser = argparse.ArgumentParser()
    parser.add_argument('FNAME',
                        type=str,
                        help='Filename of the input fasta file')
    args = parser.parse_args()
    fname = args.FNAME
    with open(fname, "r", encoding='utf-8') as f:
        data = f.readlines()
    stats = compute_data(data)
    plot(stats)


if __name__ == '__main__':
    main()
