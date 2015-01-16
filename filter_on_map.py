import cPickle as pickle
import numpy as np
import pandas as pd


def filter_by_mapping(path_to_align, cutoff_per_map = 50):
    c_to_del =[]
    with open(path_to_align, 'rb') as fp:
      a_data = pickle.load(fp)
      ind_list = a_data.index.tolist()
      cell_list = list(a_data.columns.values)
      for i,c in enumerate(cell_list):
          if a_data[c]['per_mapped'] < cutoff_per_map:
              c_to_del.append(c)
              print i, c
    print c_to_del

filter_by_mapping('/Volumes/Seq_data/pdgfra_align.p')
