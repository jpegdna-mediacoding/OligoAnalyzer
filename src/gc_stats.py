import argparse
import matplotlib.pyplot as plt


def compute_gc_content(dic):
    """Computes the GC content of a sequence"""
    gc_content = ((dic["G"] + dic["C"])/(dic["A"] + dic["T"] + dic["G"] + dic["C"]))
    return gc_content
    # return round(round(gc/0.05)*0.05,2)

def compute_data(data):
    """Computes the gc content of all the oligos in the fasta file"""
    n_oligos, n_bad_oligos = 0, 0
    stats = {"A":0, "T":0, "C":0, "G":0, "oligo_stats":[]}
    gc_content = []
    for line in data:
        if line[0] != ">":
            n_oligos += 1
            tmp_dic = {"A":0, "T":0, "C":0, "G":0}
            for el in line:
                if el in ["A", "T", "C", "G"]:
                    stats[el] += 1
                    tmp_dic[el] += 1
            stats["oligo_stats"].append(tmp_dic)
            tmp_gc = compute_gc_content(tmp_dic)
            if tmp_gc<=0.3 or tmp_gc>=0.6:
                n_bad_oligos += 1
            gc_content.append(tmp_gc)
    pb_rate = 100*n_bad_oligos/n_oligos
    print(f"Oligos with problematic GC content: {pb_rate:.2f}%")
    return stats, gc_content


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
    stats, gc_content = compute_data(data)
    fig = plt.figure()
    plt.hist(gc_content)
    plt.title("Distribution")
    plt.xlabel("GC content")
    plt.ylabel("Nb of oligos")
    plt.show()

if __name__ == '__main__':
    main()
