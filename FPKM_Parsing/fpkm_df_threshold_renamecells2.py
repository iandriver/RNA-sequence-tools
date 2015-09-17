import os
import cPickle as pickle
import numpy as np
import pandas as pd
from collections import OrderedDict

def delete_cells(by_cell, cell_list, del_list):
    to_delete1 =[]
    for pos, cell_name in enumerate(cell_list):
        if cell_name in del_list:
            to_delete1.append(pos)
    to_delete = sorted(to_delete1, reverse=True)
    for pos in to_delete:
        print 'Deleted specific cell '+cell_list[pos]
        del cell_list[pos]
    n_by_cell = np.delete(by_cell, to_delete, axis=0)
    return cell_list, n_by_cell


def filter_by_mapping(path_to_align, cutoff_per_map = 150000):
    c_to_del =[]
    with open(path_to_align, 'rb') as fp:
        a_data = pickle.load(fp)
        p_mapped = a_data['mapped_L_num']
        ind_list = p_mapped[p_mapped<cutoff_per_map]
        c_to_del = ind_list.index.values
    return c_to_del

def filter_cells_sd(by_cell, cell_list, sd=2.8):
    average_gene_exp = []
    to_delete= []
    for cell_name, genes in zip(cell_list,by_cell):
        gen_exp = (genes >= 1).sum()
        #if the cell doesn't express at least 500 genes just delete it and exclude from average
        if gen_exp <=500:
            to_delete.append(cell_list.index(cell_name))
        else:
            average_gene_exp.append(gen_exp)
    np_av = np.array(average_gene_exp)
    averg = np.average(np_av)
    gene_sd = np.std(np_av)
    print averg, gene_sd
    #add cells that fall outside of stdev value sd
    for i1, exp_level in enumerate(np_av):
        if exp_level < (averg - (gene_sd*sd)) or exp_level > (averg + (gene_sd*sd)):
            if i1 not in to_delete:
                to_delete.append(i1)
    to_delete1 = sorted(to_delete, reverse = True)
    print to_delete1
    for pos in to_delete1:
        print 'Deleted outlier '+cell_list[pos]
        del cell_list[pos]
    n_by_cell = np.delete(by_cell, to_delete1, axis=0)
    print "Number of cells remaining: "+str(len(cell_list))
    naverage_gene_exp = []
    for ngenes in n_by_cell:
        ngen_exp = (ngenes >= 1).sum()
        naverage_gene_exp.append(ngen_exp)
    nnp_av = np.array(naverage_gene_exp)
    naverg = np.average(nnp_av)
    ngene_sd = np.std(nnp_av)
    print "New", naverg, ngene_sd
    return cell_list, n_by_cell


def threshold_genes(by_gene, gen_list, number_expressed=3):
    g_todelete = []
    for g1, gene in enumerate(by_gene):
        cells_exp = (gene >= 1.0).sum()
        if cells_exp < number_expressed:
            g_todelete.append(g1)
    g1_todelete = sorted(g_todelete, reverse = True)
    print by_gene.shape
    for pos in g1_todelete:
        if type(gen_list[pos]) != float:
            #print 'Gene '+gen_list[pos]+' not expressed in '+str(number_expressed)+' cells.'
            pass
        del gen_list[pos]
    n_by_gene = np.delete(by_gene, g1_todelete, axis=0)
    print n_by_gene.shape
    return gen_list, n_by_gene

#given a pandas dataframe of gene expression split out ERCC and return seperate dataframes of each
def sep_ERCC(pd_by_gene, gen_list):
    w_gene_list = list(gen_list)
    ERCC_list= []
    pop_list =[]
    to_del = []
    for i, gen in enumerate(w_gene_list):
        if 'ERCC-00' == gen[0:7] or 'RNASPIKE1-EC02' in gen:
            pop_list.append(i)
    to_del = sorted(to_del, reverse=True)
    for d in to_del:
        del w_gene_list[d]
    pop_list = sorted(pop_list, reverse=True)
    for pos in pop_list:
        ERCC_list.append(w_gene_list.pop(pos))
    pd_by_gene_no_ERCC = pd_by_gene[w_gene_list]
    pd_ERCC = pd_by_gene[ERCC_list]
    return pd_by_gene_no_ERCC.transpose(), pd_ERCC.transpose(), w_gene_list

<<<<<<< Updated upstream
path_to_file = '/Volumes/Seq_data/cuffnorm_spc_d0_4_7'
file_name = 'genes.fpkm_table'
base_name ='spc_d0_4_7'
data = pd.DataFrame.from_csv(os.path.join(path_to_file,file_name), sep='\t')
=======
path_to_file = '/Volumes/Seq_data/results_spc2_n2'
base_name ='spc2_count'
file_name = 'spc2_count_table.txt'
cuff_df = pd.DataFrame.from_csv(os.path.join(path_to_file,file_name), sep='\t')

with open(os.path.join(path_to_file,file_name.strip('.txt')+'.p'), 'wb') as fp:
  pickle.dump(cuff_df, fp)

with open(os.path.join(path_to_file,file_name.strip('.txt')+'.p'), 'rb') as fp:
  data = pickle.load(fp)
  gen_list = data.index.tolist()
  cell_list = [x.strip('_0') for x in list(data.columns.values)]
  path_to_align=os.path.join(path_to_file,'results_spc2_all_align.p')
  del_list=filter_by_mapping(path_to_align)
>>>>>>> Stashed changes

gen_list = data.index.tolist()
cell_list1 = [x[0:-2] for x in list(data.columns.values)]
cell_list = []
for ci in cell_list1:
    if ci[-1] == '_':
        cell_list.append(ci+'2')
    else:
        cell_list.append(ci)

data.columns=cell_list
print data
path_to_align=os.path.join(path_to_file,'results_spc_d0_4_7_align.p')
del_list=filter_by_mapping(path_to_align)
print del_list, 'del'
npdata = np.array(data.values, dtype='f')
by_cell1 = npdata.transpose()
rem_cell_list, rem_by_cell = delete_cells(by_cell1, cell_list, del_list)
npdata2 = rem_by_cell.transpose()
new_gene_list1, new_by_gene = threshold_genes(npdata2, gen_list)
by_cell = new_by_gene.transpose()
outlier_cell_list, outlier_by_cell = filter_cells_sd(by_cell, rem_cell_list)
final_by_gene = outlier_by_cell.transpose()
outlier_fpkm_dict = OrderedDict()
bulk_ctrl_dict = OrderedDict()
filter_on_lane = False
bulk = False
if filter_on_lane:
    for i, l in enumerate(outlier_by_cell):
        split_cell_list = outlier_cell_list[i].split('_')
        cell_name = outlier_cell_list[i]
        if 'neg' in cell_name or '+' in cell_name or '-' in cell_name:
            print cell_name, '1'
            if 'Ra' not in split_cell_list and 'pdgfra' not in split_cell_list:
                print cell_name, '2'
                bulk_ctrl_dict['_'.join(split_cell_list[1:])] = [float(lix) for lix in l]
        elif split_cell_list[0] == 'Lane1' and 'C' in split_cell_list[1] or split_cell_list[0] == 'Lane2' and 'C' in split_cell_list[1]:
            if 'Ra' not in split_cell_list and 'pdgfra' not in split_cell_list:
                if len(split_cell_list[1]) == 3:
                    cell_name = 'CTRL_'+split_cell_list[1]
                elif len(split_cell_list[1]) == 2:
                    cell_name = 'CTRL_'+split_cell_list[1][0]+'0'+split_cell_list[1][1]
                else:
                    print split_cell_list[1], 'wtf'
                outlier_fpkm_dict[cell_name] = [float(lx) for lx in l]
        elif split_cell_list[0] == 'Lane3' and 'C' in split_cell_list[1] or split_cell_list[0] == 'Lane4'and 'C' in split_cell_list[1]:
            if 'Ra' not in split_cell_list and 'pdgfra' not in split_cell_list:
                if len(split_cell_list[1]) == 3:
                    cell_name = 'PNX_'+split_cell_list[1]
                elif len(split_cell_list[1]) == 2:
                    cell_name = 'PNX_'+split_cell_list[1][0]+'0'+split_cell_list[1][1]
                else:
                    print split_cell_list[1], 'wtf'
                outlier_fpkm_dict[cell_name] = [float(lx) for lx in l]
else:
    for i, l in enumerate(outlier_by_cell):
        if bulk:
            cell_name = outlier_cell_list[i]
            if 'bulk' in cell_name:
                bulk_ctrl_dict[cell_name] = [float(lx) for lx in l]
            else:
                outlier_fpkm_dict[cell_name] = [float(lx) for lx in l]
        else:
            cell_name = outlier_cell_list[i]
            outlier_fpkm_dict[cell_name] = [float(lx) for lx in l]
if bulk:
    df_bulk = pd.DataFrame(bulk_ctrl_dict, index = new_gene_list1)
else:
    pass
fpkm_df_outlier1 = pd.DataFrame(outlier_fpkm_dict, index = new_gene_list1)
mod_gen_list = list(new_gene_list1)
fpkm_df_outlier, df_ERCC, new_gene_list  = sep_ERCC(fpkm_df_outlier1.transpose(), mod_gen_list)
mod_gen_list = list(new_gene_list1)
if bulk:
    bulk_ctrl_df, df_bulk_ERCC, bulk_gene_list = sep_ERCC(df_bulk.transpose(), mod_gen_list)
else:
    pass
outlier_cell_list = [x for x in list(fpkm_df_outlier.columns.values)]
df_ERCC.to_csv(os.path.join(path_to_file, base_name+'_ERCC.txt'), sep = '\t')
if bulk:
    df_bulk_ERCC.to_csv(os.path.join(path_to_file, base_name+'_bulkonly_ERCC.txt'), sep = '\t')
    bulk_ctrl_df.to_csv(os.path.join(path_to_file, base_name+'_bulk_ctrls.txt'), sep = '\t')
else:
    pass
fpkm_df_outlier.to_csv(os.path.join(path_to_file, base_name+'_outlier_filtered.txt'), sep = '\t')
with open(os.path.join(path_to_file,base_name+'_outlier_by_cell.p'), 'wb') as fp1:
    pickle.dump(fpkm_df_outlier.transpose(), fp1)
with open(os.path.join(path_to_file,base_name+'_outlier_cell_list.p'), 'wb') as fp2:
    pickle.dump(outlier_cell_list, fp2)
with open(os.path.join(path_to_file,base_name+'_outlier_by_gene.p'), 'wb') as fp3:
    pickle.dump(fpkm_df_outlier, fp3)
with open(os.path.join(path_to_file,base_name+'_outlier_gene_list.p'), 'wb') as fp4:
    pickle.dump(new_gene_list, fp4)
