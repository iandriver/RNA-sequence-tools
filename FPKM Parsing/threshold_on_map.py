import cPickle as pickle
import numpy as np
import pandas as pd


def filter_by_mapping(path_to_align, cutoff_per_map = 50):
    c_to_del =[]
    with open(path_to_align, 'rb') as fp:
      a_data = pickle.load(fp)
      p_mapped = a_data['per_mapped']
      ind_list = p_mapped[p_mapped<cutoff_per_map]
      print ind_list.index.values


filter_by_mapping('/Volumes/Seq_data/results_pdgfra1_ctrl_pnx/results_pdgfra_align.p')
