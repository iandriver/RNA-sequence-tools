import pandas as pd
from pprint import pprint
import itertools
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt

qpcr_file ='/Users/iandriver/Documents/2016-04-06_results.csv'
qpcr_df = pd.read_csv(qpcr_file, sep=None, engine='python')

samples = list(set(qpcr_df['Sample Name']))
targets = list(set(qpcr_df['Target Name']))
sample_dict = {}
index_list = []
for samp in samples:
    target_dict = {}
    delta_ct_actb_dict = {}
    delta_ct_gapdh_dict = {}
    sample_df =  qpcr_df[qpcr_df['Sample Name'] == samp]
    for targ in targets:
        target_df = sample_df[(sample_df['Target Name'] == targ)]
        targ_mean = target_df['CT'].convert_objects(convert_numeric=True).mean()
        target_dict[targ] = targ_mean
        index_list.append(samp)
    for targ2 in targets:
        actb_mean = target_dict['ACTB']
        gapdh_mean = target_dict['GAPDH']
        if targ2 != 'ACTB':
            delta_ct_actb = actb_mean - target_dict[targ2]
        else:
            delta_ct_actb = 0
        if targ2 != 'GAPDH':
            delta_ct_gapdh = gapdh_mean - target_dict[targ2]
        else:
            delta_ct_gapdh = 0
        delta_ct_actb_dict[targ2] = delta_ct_actb
        delta_ct_gapdh_dict[targ2] = delta_ct_gapdh
    sample_dict[samp] = target_dict
    sample_dict['delta_ct_actb_'+samp] = delta_ct_actb_dict
    sample_dict['delta_ct_gapdh_'+samp] = delta_ct_gapdh_dict
delta_pairs = []
for samp1,samp2 in itertools.permutations(samples,2):
    if samp1 != samp2 and samp1[-2:] == samp2[-2:] and 'No' in samp2:
        delta_pairs.append((samp1,samp2))
results_df = pd.DataFrame.from_dict(sample_dict)
fc_results_df = 0
for p in delta_pairs:
    pow_dict = dict(zip(targets,[2 for t in targets]))
    ratio_dict = {'pos_dict_a':sample_dict['delta_ct_actb_'+p[0]],'neg_dict_a':sample_dict['delta_ct_actb_'+p[1]], 'pos_dict_g':sample_dict['delta_ct_gapdh_'+p[0]], 'neg_dict_g':sample_dict['delta_ct_gapdh_'+p[1]], 'pwer':pow_dict}
    pair_df = pd.DataFrame.from_dict(ratio_dict)
    pwer_df = pair_df['pwer']
    ratio_df_a = pd.DataFrame(pwer_df.pow(pair_df['pos_dict_a'])/pwer_df.pow(pair_df['neg_dict_a']), columns=[p[0]+'_to_'+p[1]+'_ratio_actb'])
    ratio_df_g = pd.DataFrame(pwer_df.pow(pair_df['pos_dict_g'])/pwer_df.pow(pair_df['neg_dict_g']), columns=[p[0]+'_to_'+p[1]+'_ratio_gapdh'])

    fc_all = pd.merge(ratio_df_a,ratio_df_g, right_index=True, left_index=True)
    all_results = pd.merge(results_df,fc_all, right_index=True, left_index=True)
    results_df = all_results.copy()
    if not isinstance(fc_results_df, pd.DataFrame):
        all_fc_results_df = fc_all.copy()
    else:
        all_fc_results_df = pd.merge(fc_results_df,fc_all, right_index=True, left_index=True)
    fc_results_df = all_fc_results_df.copy()
plot_df_dict = {}
target_list = all_fc_results_df.index.tolist()

plot_df_dict['Target'] = []
plot_df_dict['Ratio ITGA8+ to ITGA8-'] = []
plot_df_dict['Control'] = []
for  t in all_fc_results_df.columns.values:
    control_name = t.split('_')[-1]
    plot_df_dict['Target']= plot_df_dict['Target']+ target_list
    plot_df_dict['Ratio ITGA8+ to ITGA8-']= plot_df_dict['Ratio ITGA8+ to ITGA8-'] +all_fc_results_df[t].tolist()
    plot_df_dict['Control'] = plot_df_dict['Control']+[control_name]*len(target_list)
print(plot_df_dict)
plot_df = pd.DataFrame.from_dict(plot_df_dict)
print(plot_df)
sns.barplot(x='Target', y='Ratio ITGA8+ to ITGA8-', hue = 'Control', data=plot_df)
plt.show()
results_df.to_csv(os.path.join(os.path.dirname(qpcr_file), 'qpcr_results_'+qpcr_file.split('/')[-1]), sep='\t')
