from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import numpy as np
import pandas as pd

# define quantile normalized function
def quantile_normalize(dataframe):

    # copy dataframe and only use the columns with numerical values
    df = dataframe.copy()

    rank_mean = df.stack().groupby(df.rank(method='first').stack().astype(int)).mean()  

    result = df.rank(method='min').stack().astype(int).map(rank_mean).unstack()
    
    return result, rank_mean

# takes DNA sequence, outputs one-hot-encoded matrix with rows A, T, G, C
def one_hot_encoder(sequence):
    l = len(sequence)
    x = np.zeros((4,l),dtype = 'int8')
    for j, i in enumerate(sequence):
        if i == "A" or i == "a":
            x[0][j] = 1
        elif i == "T" or i == "t":
            x[1][j] = 1
        elif i == "G" or i == "g":
            x[2][j] = 1
        elif i == "C" or i == "c":
            x[3][j] = 1
        else:
            return "contains_N"
    return x

#read names and postions from bed file
def read_bed(filename):
    positions = defaultdict(list)
    with open(filename) as f:
        for line in f:
            name, chr, start, stop = line.split()
            positions[name].append((chr, int(start), int(stop)))

    return positions


# parse fasta file and turn into dictionary
def read_fasta(genome_dir):
    # chr_dict = dict()
    # for chr in range(1, num_chr):
    #     chr_file_path = genome_dir + "chr{}.fa".format(chr)
    #     chr_dict.update(SeqIO.to_dict(SeqIO.parse(open(chr_file_path), 'fasta')))
    chr_dict = dict()
    chr_file_path = genome_dir
    chr_dict.update(SeqIO.to_dict(SeqIO.parse(open(chr_file_path), 'fasta')))

    return chr_dict


#get sequences for peaks from reference genome
def get_sequences(positions, chr_dict):
    one_hot_seqs = []
    peak_seqs = []
    invalid_ids = []
    peak_names = []

    # target_chr = ['chr{}'.format(i) for i in range(1, num_chr)]

    for name in positions:
        for (chr, start, stop) in positions[name]:
            chr_seq = chr_dict[chr].seq
            peak_seq = str(chr_seq)[start - 1:stop].lower()
            one_hot_seq = one_hot_encoder(peak_seq)
            if isinstance(one_hot_seq, np.ndarray):
                peak_names.append(name)
                peak_seqs.append(peak_seq)
                one_hot_seqs.append(one_hot_seq)
            else:
                invalid_ids.append(name)

    one_hot_seqs = np.stack(one_hot_seqs)
    peak_seqs = np.stack(peak_seqs)
    peak_names = np.stack(peak_names)

    return one_hot_seqs, peak_seqs, invalid_ids, peak_names

def get_sequences_de_nove(de_nove_overlap, chr_dict):
    one_hot_seqs_normal = []
    one_hot_seqs_SNV = []
    peak_seqs_normal = []
    peak_seqs_SNV = []
    invalid_ids = []
    peak_names = []

    for de_nove_index in de_nove_overlap.shape[0]:
        chr = de_nove_overlap.loc[de_nove_index,'peak_chr']
        start = de_nove_overlap.loc[de_nove_index,'peak_start']
        stop = de_nove_overlap.loc[de_nove_index,'peak_end']
        Pos_SNV = de_nove_overlap.loc[de_nove_index,'Pos']
        
        chr_seq = chr_dict[chr].seq
        peak_seq_normal = str(chr_seq)[start - 1:stop].lower()
        peak_seq_SNV = peak_seq_normal[0:Pos_SNV-start] + de_nove_overlap.loc[de_nove_index,'Alt'].lower() + peak_seq_normal[Pos_SNV-start+1:]

        one_hot_seq_normal = one_hot_encoder(peak_seq_normal)
        one_hot_seq_SNV = one_hot_encoder(peak_seq_SNV)
        
        if isinstance(one_hot_seq_normal, np.ndarray):
            peak_names.append(de_nove_overlap.loc[de_nove_index,'peak_ID'])
            peak_seqs_normal.append(peak_seq_normal)
            peak_seqs_SNV.append(peak_seq_SNV)
            one_hot_seqs_normal.append(one_hot_seq_normal)
            one_hot_seqs_SNV.append(one_hot_seq_SNV)
        else:
            invalid_ids.append(de_nove_index)
        
        one_hot_seqs_normal = np.stack(one_hot_seqs_normal)
        one_hot_seqs_SNV = np.stack(one_hot_seqs_SNV)      
        peak_seqs_normal = np.stack(peak_seqs_normal)
        peak_seqs_SNV = np.stack(peak_seqs_SNV)
        peak_names = np.stack(peak_names)

    return one_hot_seqs_normal, one_hot_seqs_SNV, peak_seqs_normal, peak_seqs_SNV, invalid_ids, peak_names

def format_intensities(intensity_file, invalid_ids):
    
    multiomics_peak_height_norm = pd.read_table(intensity_file, sep=',', header=None)
    valid_ids = np.isin(multiomics_peak_height_norm[0],invalid_ids)
    cell_type_array = multiomics_peak_height_norm.loc[~valid_ids,1:]
    peak_names = multiomics_peak_height_norm.loc[~valid_ids,0]
    peak_names = np.stack(peak_names)

    return cell_type_array, peak_names



