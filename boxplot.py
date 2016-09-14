import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats
import numpy as np
import matplotlib

path_to_file = '/Users/idriver/RockLab-files/tcf_pdgfra_edu_flow_data2.txt'
facs_df = pd.read_csv(path_to_file, sep='\t', header = 0)

sns.set_context("notebook", font_scale=1.5, rc={"lines.linewidth": 2.5})
font = {'family' : 'normal',
        'size'   : 22}

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

bp = sns.boxplot(x='Metric', y='Percent', hue='Condition', data=facs_df)
add_on = 0
for m in set(facs_df['Metric']):
    df_1 = facs_df[(facs_df['Metric']==m)&(facs_df['Condition']=='Saline')]
    df_2 = facs_df[(facs_df['Metric']==m)&(facs_df['Condition']=='Bleomycin')]
    z, p = scipy.stats.f_oneway(df_1['Percent'],df_2['Percent'])
    p_value = p * 2
    s = stars(p)
    print(m,p_value, s)
    bp.set_ylim(ymin=0)
    y_max = facs_df['Percent'].max()
    y_min = facs_df['Percent'].min()
    print(y_max)
    bp.annotate("", xy=(-.2+add_on, y_max+.15), xycoords='data',
                xytext=(.2+add_on, y_max+.15), textcoords='data',
                arrowprops=dict(arrowstyle="-", ec='black',
                                connectionstyle="bar,fraction=0.1", lw=1.3))
    bp.text(0+add_on, y_max + abs(y_max - y_min)*0.1, stars(p_value),
            horizontalalignment='center',
            verticalalignment='center', color='black')
    add_on+=1

plt.show()
