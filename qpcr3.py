import pandas as pd
from pprint import pprint
import itertools
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

qpcr_file ='/Users/idriver/RockLab-files/Tcf_Belo-ctrl-090816_all.txt'
qpcr_df = pd.read_csv(qpcr_file, delimiter= '\t')

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
        targ_mean = pd.to_numeric(target_df['CT']).mean()
        target_dict[targ] = targ_mean
        index_list.append(samp)
    for targ2 in targets:
        actb_mean = target_dict['Actb']
        gapdh_mean = target_dict['Gapdh']
        if targ2 != 'Actb':
            delta_ct_actb = actb_mean - target_dict[targ2]
        else:
            delta_ct_actb = 0
        if targ2 != 'Gapdh':
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
    if samp1 != samp2 and 'Norm' in samp2 and 'Bleo' in samp1 and 'Socs3' not in samp1 and 'Norm2' not in samp2:
        delta_pairs.append((samp1,samp2))
results_df = pd.DataFrame.from_dict(sample_dict)
for p in delta_pairs:
    pow_dict = dict(zip(targets,[2 for t in targets]))
    ratio_dict = {'pos_dict_a':sample_dict['delta_ct_actb_'+p[0]],'neg_dict_a':sample_dict['delta_ct_actb_'+p[1]], 'pos_dict_g':sample_dict['delta_ct_gapdh_'+p[0]], 'neg_dict_g':sample_dict['delta_ct_gapdh_'+p[1]], 'pwer':pow_dict}
    pair_df = pd.DataFrame.from_dict(ratio_dict)
    pwer_df = pair_df['pwer']
    ratio_df_a = pd.DataFrame(pwer_df.pow(pair_df['pos_dict_a'])/pwer_df.pow(pair_df['neg_dict_a']), columns=[p[0]+'_to_'+p[1]+'_ratio_Actb'])
    ratio_df_g = pd.DataFrame(pwer_df.pow(pair_df['pos_dict_g'])/pwer_df.pow(pair_df['neg_dict_g']), columns=[p[0]+'_to_'+p[1]+'_ratio_Gapdh'])

    fc_all = pd.merge(ratio_df_a,ratio_df_g, right_index=True, left_index=True)

    all_results = pd.merge(results_df,fc_all, right_index=True, left_index=True)
    results_df = all_results.copy()
plot_df_dict = {}
target_list = results_df.index.tolist()

results_df.to_csv(os.path.join(os.path.dirname(qpcr_file), 'qpcr_results_'+qpcr_file.split('/')[-1]), sep='\t')

plot_df_dict['Target'] = []
plot_df_dict['Ratio Bleo to Saline Control'] = []
plot_df_dict['Control'] = []
selected_genes= ['Actb', 'Gapdh', 'Col3a1', 'Col1a2','Acta2','Il6', 'Pdgfra', 'G0s2', 'Tcf21']
ratio_df = results_df[[n for n in results_df.columns.values if 'ratio' in n]]
selected_ratio_df = ratio_df.loc[selected_genes]
selected_ratio_df.reindex(selected_genes)
for  t in [n for n in selected_ratio_df.columns.values if 'ratio' in n]:
    control_name = t.split('_')[-1]
    if control_name == 'Actb':
        plot_df_dict['Target']= plot_df_dict['Target']+ selected_genes
        plot_df_dict['Ratio Bleo to Saline Control']= plot_df_dict['Ratio Bleo to Saline Control'] +selected_ratio_df[t].tolist()
        plot_df_dict['Control'] = plot_df_dict['Control']+[control_name]*len(selected_genes)

plot_df = pd.DataFrame.from_dict(plot_df_dict)
plot_df.to_csv(os.path.join(os.path.dirname(qpcr_file), 'qpcr_plot_df_selected_'+qpcr_file.split('/')[-1]), sep='\t')

fig, ax = plt.subplots()
sns.boxplot(x='Target', y='Ratio Bleo to Saline Control', hue = 'Control', data=plot_df, ax = ax)
add_on =0
ymax = 30 #plot_df['Ratio Bleo to Saline Control'].max()
pos = np.arange(len(set(plot_df['Target'])))
ax.yaxis.set_ticks(np.arange(0,45,1))
ax.set_ylim([0,31])
for tick, label in zip(range(len(set(plot_df['Target']))), plot_df['Target']):

    df_1 = plot_df[(plot_df['Target']==label)&(plot_df['Control']=='Actb')]
    #df_2 = plot_df[(plot_df['Target']==label)&(plot_df['Control']=='Gapdh')]

    ratio_a = df_1['Ratio Bleo to Saline Control'].mean()
    #ratio_g = df_2['Ratio Bleo to Saline Control'].mean()
    ax.text(pos[tick], ymax + 0.1, "%.1f" % ratio_a,horizontalalignment='center', color='blue')
    #ax.text(pos[tick], ymax + 0.1, "%.1f" % ratio_g,horizontalalignment='center', color='green')
ax.legend_.remove()

plt.show()
