# -*- coding: utf-8 -*-
"""
Created on Sun Jan  8 20:05:34 2023

@author: Sonric
"""
import numpy as np
from numpy import random
import pandas as pd
from collections import Counter

from sklearn.model_selection import KFold
import torch.utils.data
import matplotlib
import os
import sys
matplotlib.use('Agg')

import aitac as aitac
import plot_utils

filter_enrich_all = []
signif_threshold = 0.05
for model_index in (np.arange(10)+1):
    tomtom_model = pd.read_csv("../outputs/valid10x10/motifs/tomtom_compare_to_model_%s"%model_index + ".tsv",sep='\t')
    tomtom_model = tomtom_model[tomtom_model['q-value'] < signif_threshold]
    filter_enrich = list(set(list(tomtom_model['Query_ID'])))
    filter_enrich_all.extend(filter_enrich)

filter_enrich_counts = pd.DataFrame(pd.Series(Counter(filter_enrich_all)))
filter_enrich_counts = filter_enrich_counts.sort_values(by=0, ascending=False)

filter_enrich_counts['filter_name'] = filter_enrich_counts.index
filter_enrich_counts.columns = ['num','filter_name']
filter_enrich_counts = filter_enrich_counts[['filter_name','num']]
np.save('../outputs/test/motifs/filter_enrich_counts.npy',filter_enrich_counts)
