import pandas as pd
from pprint import pprint
import itertools
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
import scipy
import matplotlib

qpcr_file ='/Users/idriver/Downloads/2016-04-19_ITGA8_1_6_7_8.txt'
qpcr_df = pd.DataFrame.from_csv(qpcr_file, sep= '\t')

sns.set_context("notebook", font_scale=1.5, rc={"lines.linewidth": 2.5})
font = {'family' : 'normal',
        'size'   : 19}
matplotlib.rc('font', **font)

def stars(p):
    if p < 0.0001:
        return "****"
    elif (p < 0.001):
        return "***"
    elif (p < 0.01):
        return "**"
    elif (p < 0.05):
        return "*"
    else:
        return "ns"

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
    if samp1 != samp2 and samp1[-2:] == samp2[-2:] and 'neg' in samp2:
        delta_pairs.append((samp1,samp2))
results_df = pd.DataFrame.from_dict(sample_dict)
gene_df_list = []
for p in delta_pairs:
    pow_dict = dict(zip(targets,[2 for t in targets]))
    ratio_dict = {'pos_dict_a':sample_dict['delta_ct_actb_'+p[0]],'neg_dict_a':sample_dict['delta_ct_actb_'+p[1]], 'pos_dict_g':sample_dict['delta_ct_gapdh_'+p[0]], 'neg_dict_g':sample_dict['delta_ct_gapdh_'+p[1]], 'pwer':pow_dict}

    pair_df = pd.DataFrame.from_dict(ratio_dict)
    by_gene_df = pair_df.transpose()
    new_index_dict ={}
    for i in by_gene_df.index.tolist():
        new_index_dict[i]= i+'_'+p[0]+'_to_'+p[1]
    by_gene_df.rename(new_index_dict, inplace=True)
    gene_df_list.append(by_gene_df)
    pwer_df = pair_df['pwer']
    ratio_df_a = pd.DataFrame(pwer_df.pow(pair_df['pos_dict_a'])/pwer_df.pow(pair_df['neg_dict_a']), columns=[p[0]+'_to_'+p[1]+'_ratio_Actb'])
    ratio_df_g = pd.DataFrame(pwer_df.pow(pair_df['pos_dict_g'])/pwer_df.pow(pair_df['neg_dict_g']), columns=[p[0]+'_to_'+p[1]+'_ratio_Gapdh'])

    fc_all = pd.merge(ratio_df_a,ratio_df_g, right_index=True, left_index=True)

    all_results = pd.merge(results_df,fc_all, right_index=True, left_index=True)
    results_df = all_results.copy()
all_gene_df = pd.concat(gene_df_list)


plot_df_dict = {}
target_list = results_df.index.tolist()

results_df.to_csv(os.path.join(os.path.dirname(qpcr_file), 'qpcr_results_'+qpcr_file.split('/')[-1]), sep='\t')


plot_df_dict['Target'] = []
plot_df_dict['Ratio ITGA8+ to ITGA8-'] = []
plot_df_dict['Control'] = []


selected_genes= ['ACTB', 'GAPDH', 'ITGA8', 'PDGFRB','COL14A1','TCF21', 'G0S2', 'PLIN2']

pos_df_a = all_gene_df[all_gene_df.index.map(lambda x: 'pos_dict_a' in x)]
neg_df_a = all_gene_df[all_gene_df.index.map(lambda x: 'neg_dict_a' in x)]
pos_df_g = all_gene_df[all_gene_df.index.map(lambda x: 'pos_dict_g' in x)]
neg_df_g = all_gene_df[all_gene_df.index.map(lambda x: 'neg_dict_g' in x)]
stats_dict = {}
for gene in selected_genes:
    stats_dict[gene+'_a'] = scipy.stats.f_oneway(pos_df_a[gene],neg_df_a[gene])
    stats_dict[gene+'_g'] = scipy.stats.f_oneway(pos_df_g[gene],neg_df_g[gene])
stats_df = pd.DataFrame.from_dict(stats_dict)
print(stats_df)
ratio_df = results_df[[n for n in results_df.columns.values if 'ratio' in n]]
selected_ratio_df = ratio_df.loc[selected_genes]
selected_ratio_df.reindex(selected_genes)
plot_controls = 'Actb' #values can be 'Actb' 'Gapdh' or both

for  t in [n for n in selected_ratio_df.columns.values if 'ratio' in n]:
    control_name = t.split('_')[-1]
    if plot_controls == 'Actb':
        if control_name == 'Actb':
            plot_df_dict['Target']= plot_df_dict['Target']+ selected_genes
            plot_df_dict['Ratio ITGA8+ to ITGA8-']= plot_df_dict['Ratio ITGA8+ to ITGA8-'] +selected_ratio_df[t].tolist()
            plot_df_dict['Control'] = plot_df_dict['Control']+[control_name]*len(selected_genes)
    elif plot_controls == 'Gapdh':
        if control_name == 'Gapdh':
            plot_df_dict['Target']= plot_df_dict['Target']+ selected_genes
            plot_df_dict['Ratio ITGA8+ to ITGA8-']= plot_df_dict['Ratio ITGA8+ to ITGA8-'] +selected_ratio_df[t].tolist()
            plot_df_dict['Control'] = plot_df_dict['Control']+[control_name]*len(selected_genes)
    elif plot_controls == 'both':
        plot_df_dict['Target']= plot_df_dict['Target']+ selected_genes
        plot_df_dict['Ratio ITGA8+ to ITGA8-']= plot_df_dict['Ratio ITGA8+ to ITGA8-'] +selected_ratio_df[t].tolist()
        plot_df_dict['Control'] = plot_df_dict['Control']+[control_name]*len(selected_genes)

plot_df = pd.DataFrame.from_dict(plot_df_dict)
plot_df.to_csv(os.path.join(os.path.dirname(qpcr_file), 'qpcr_plot_df_selected_'+qpcr_file.split('/')[-1]), sep='\t')

fig, ax = plt.subplots()
sns.boxplot(x='Target', y='Ratio ITGA8+ to ITGA8-', hue = 'Control', data=plot_df, ax = ax)
add_on =0
ymax = plot_df['Ratio ITGA8+ to ITGA8-'].max()
y_min = plot_df['Ratio ITGA8+ to ITGA8-'].min()
pos = np.arange(len(set(plot_df['Target'])))
ax.yaxis.set_ticks(np.arange(0,ymax+5,5))
ax.set_ylim([0,ymax+1])
for tick, label in zip(range(len(set(plot_df['Target']))), plot_df['Target']):
    if plot_controls =='Actb':
        p_value_a = stats_df[label+'_a'][1]
        df_1 = plot_df[(plot_df['Target']==label)&(plot_df['Control']=='Actb')]
        ratio_a = df_1['Ratio ITGA8+ to ITGA8-'].mean()
        ax.text(pos[tick], ymax + 1.2, "%.1f" % ratio_a,horizontalalignment='center', color='blue')
        ax.text(pos[tick], ymax + abs(ymax - y_min)*0.1, stars(p_value_a),
                horizontalalignment='center',
                verticalalignment='center', color='black')
    elif plot_controls =='Gapdh':
        p_value_g = stats_df[label+'_g'][1]
        df_2 = plot_df[(plot_df['Target']==label)&(plot_df['Control']=='Gapdh')]
        ratio_g = df_2['Ratio ITGA8+ to ITGA8-'].mean()
        ax.text(pos[tick], ymax+1.2, "%.1f" % ratio_g,horizontalalignment='center', color='blue')
        ax.text(pos[tick], ymax + abs(ymax - y_min)*0.1, stars(p_value_g),
                horizontalalignment='center',
                verticalalignment='center', color='black')
    elif plot_controls =='both':
        df_1 = plot_df[(plot_df['Target']==label)&(plot_df['Control']=='Actb')]
        df_2 = plot_df[(plot_df['Target']==label)&(plot_df['Control']=='Gapdh')]
        p_value_a = stats_df[label+'_a'][1]
        p_value_g = stats_df[label+'_g'][1]
        ratio_a = df_1['Ratio ITGA8+ to ITGA8-'].mean()
        ratio_g = df_2['Ratio ITGA8+ to ITGA8-'].mean()
        ax.text(pos[tick]-0.15, ymax + 1.2, "%.1f" % ratio_a,
                horizontalalignment='center', color='blue')
        ax.text(pos[tick], ymax + -0.4 +abs(ymax - y_min)*0.1, stars(p_value_a),
                horizontalalignment='center',
                verticalalignment='center', color='blue')
        ax.text(pos[tick]+0.15, ymax + 1.2, "%.1f" % ratio_g,
                horizontalalignment='center', color='green')
        ax.text(pos[tick], ymax +0.2 + abs(ymax - y_min)*0.1, stars(p_value_g),
                horizontalalignment='center',
                verticalalignment='center', color='green')
if plot_controls == 'both':
    ax.text(pos[tick]+1, ymax + 0.2+ abs(ymax - y_min)*0.1, 'p-value',
            horizontalalignment='center',
            verticalalignment='center', color='green')
    ax.text(pos[tick]+1, ymax - 0.4 + abs(ymax - y_min)*0.1, 'p-value',
            horizontalalignment='center',
            verticalalignment='center', color='blue')
    ax.text(pos[tick]+1, ymax + 1.2, 'Ratio Mean',
            horizontalalignment='center',
            verticalalignment='center', color='blue')
else:
    ax.text(pos[tick]+1, ymax + abs(ymax - y_min)*0.1, 'p-value',
            horizontalalignment='center',
            verticalalignment='center', color='black')
    ax.text(pos[tick]+1, ymax + 1.2, 'Ratio Mean',
            horizontalalignment='center',
            verticalalignment='center', color='blue')

ax.legend_.remove()

plt.show()
