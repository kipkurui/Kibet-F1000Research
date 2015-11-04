
"""
This module contaisn functions to run a motif enrichment analysis using CentriMo, 
summarize and plot the data. 

Requires:
    meme 4.10.0 obtainable from http://meme-suite.org/doc/download.html?man_type=web

Usage:
    python run_centrimo.py <Tf_name> <chip-seq_list> <test_meme_file> <test_meme_file> >results_path>
"""

import os
import random

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def run_centrimo(tf, chip_seq_list, test_meme_input, files_path):
    
    # If chip_seq list is greater than 10, randomly sample 10.
    
    if len(chip_seq_list) > 10:
        random.seed(10)
        chip_seq_list = random.sample(chip_seq_list, 10)
        
    test_list = [] #get list of motifs in meme file
    
    test_list.append(["Motif"])
    with open(test_meme_input) as motifs:
        for motif_line in motifs:
            if motif_line.startswith("MOTIF"):
                test_list.append([motif_line.split()[1]])

    def get_centrimo_list(chip_name):
        '''
        Extracts important details form a CentriMo run output
        '''
        with open(chip_name) as cent:
            temp_dict = {}
            for line in cent:
                if line.startswith("#"):
                    continue
                else:
                    temp_dict[line.split()[1]] = float(line.split()[5])*-1
        
        for mot in range(len(test_list)):
            if test_list[mot][0] in temp_dict:
                test_list[mot].append(temp_dict[test_list[mot][0]])
            elif test_list[mot][0] == "Motif":
                continue
            else:
                test_list[mot].append(0)

    # Run an enrichment analsysis for each of the given sequences
    
    for chip_seq in chip_seq_list:
        file_name = chip_seq.split('/')[-1].split('.')[0]

        test_list[0].append(file_name)
        tab2fasta(chip_seq, '/tmp/'+file_name+'.fa', '/tmp/'+file_name+'.bg')

        os.system("fasta-get-markov  /tmp/%s.fa /tmp/%s.fa.bg" % (file_name, file_name))
        os.system("centrimo --oc /tmp/%s --verbosity 1 --local --optimize_score --score 5.0 --ethresh 100000.0 --neg /tmp/%s.bg --bgfile /tmp/%s.fa.bg /tmp/%s.fa %s" %
                  (file_name, file_name, file_name, file_name, test_meme_input))
        get_centrimo_list("/tmp/%s/centrimo.txt" % file_name)

    test_list[0].append("Average")

    for i in range(1, len(test_list)):
        test_list[i].append(np.mean(test_list[i][1:]))
    test_list.sort(key=lambda x: x[-1], reverse=True)
    with open('%s/%s_centrimo.txt' % (files_path, tf), 'w') as cent_out:
        for i in test_list:
            cent_out.writelines('\t'.join(map(str, i)) + '\n')
    plot_centrimo('%s/%s_centrimo.txt' % (files_path, tf), '%s/%s_centrimo.png' % (files_path, tf))
    plot_centrimo('%s/%s_centrimo.txt' % (files_path, tf), '%s/%s_centrimo.eps' % (files_path, tf))

    
def tab2fasta(posneg, fasta, background):
    '''
    Since CentriMo takes input in two separate fasta files, this function
    extracts the positive and the background sequences in two separate files
    '''
    i = 0
    with open(posneg) as tab:
        with open(fasta, 'w') as fa:
            with open(background, 'w') as bg:
                for line in tab:
                    details = line.split()
                    if len(details) == 2:
                        pos = 1
                    else:
                        pos = 2
                    if i < 500:
                        fa.write(">"+line.split()[0]+'\n'+line.split()[pos]+"\n")
                    else:
                        bg.write(">"+line.split()[0]+'\n'+line.split()[pos]+"\n")
                    i += 1


def plot_centrimo(centrimo_in, figure_output):
    centrimo_table = pd.read_table(centrimo_in, index_col=0)
    centrimo_table.sort(columns="Average", axis=0, ascending=False, inplace=True)
    sns.clustermap(centrimo_table, method='single', metric="euclidean",
                   z_score=None, row_cluster=False, col_cluster=True)
    f = plt.gcf()
    f.savefig(figure_output, bbox_inches='tight')


if __name__ == '__main__':
    if len(sys.argv) < 5:
        print __doc__
        sys.exit(1)
    tf = sys.argv[1]
    chip_seq_list = sys.argv[2]
    test_meme_input = sys.argv[3]
    results_path = sys.argv[4]
    run_centrimo(tf, chip_seq_list, test_meme_input, results_path)