import pandas as pd
import numpy as np

gene_anova_df = pd.DataFrame.from_csv('/Users/idriver/RockLab-files/Pdgfra_Anova_w_gene_clusters.txt', sep='\t')
upin_2_1_df = gene_anova_df[(gene_anova_df['SC_2_1.Average'] > gene_anova_df['SC_3_3.Average']) & (gene_anova_df['SC_2_1.Average'] > gene_anova_df['SC_4_7.Average']) & (gene_anova_df['SC_2_1.Average'] > gene_anova_df['SC_4_8.Average'])]
downin_2_1_df = gene_anova_df[(gene_anova_df['SC_2_1.Average'] < gene_anova_df['SC_3_3.Average']) & (gene_anova_df['SC_2_1.Average'] < gene_anova_df['SC_4_7.Average']) & (gene_anova_df['SC_2_1.Average'] < gene_anova_df['SC_4_8.Average'])]
upin_2_1_df.to_csv('.txt', sep = '\t')
4
