import numpy as np
import json
import os
import _pickle as pickle
import preprocess_utils

################################ input files:

data_file = '../data/multiomics_peak_region_normalized.txt'
intensity_file = '../data/multiomics_peak_height_normalized.csv'
reference_genome = '../data/genome_hg38.fa'
output_directory = '../data/'

if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# read bed file with peak positions, and keep only entries with valid activity vectors
positions = preprocess_utils.read_bed(data_file)

# read reference genome fasta file into dictionary
if not os.path.exists('../data/chr_dict.pickle'):
    chr_dict = preprocess_utils.read_fasta(reference_genome)
    pickle.dump(chr_dict, open('../data/chr_dict.pickle', "wb"))
    
chr_dict = pickle.load(open('../data/chr_dict.pickle', "rb"))

one_hot_seqs, peak_seqs, invalid_ids, peak_names = preprocess_utils.get_sequences(positions, chr_dict)

# remove invalid ids from intensities file so sequence/intensity files match
cell_type_array, peak_names = preprocess_utils.format_intensities(intensity_file, invalid_ids)
cell_type_array = cell_type_array.astype(np.float32)

# write to file
np.save(output_directory + 'one_hot_seqs.npy', one_hot_seqs)
np.save(output_directory + 'peak_names.npy', peak_names)
np.save(output_directory + 'peak_seqs.npy', peak_seqs)
np.save(output_directory + 'cell_type_array.npy', cell_type_array)