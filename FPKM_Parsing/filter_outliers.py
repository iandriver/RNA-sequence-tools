import os
import fnmatch
import cPickle as pickle
import numpy as np
import pandas as pd
from collections import OrderedDict
import matplotlib.pyplot as plt

#This section will take fpkm matrix input and make pandas dataframe

#path to fpkm file (usually cuffnorm output)
path_to_file = '/Volumes/Seq_data/cuffnorm_pdgfra_1_and_2'
#default file name will use genes.fpkm_table from cuffnorm
file_name = 'genes.fpkm_table'
#provide base name for output files
base_name ='pdgfra2_all_n2'
#create pandas dataframe from fpkm files
data = pd.DataFrame.from_csv(os.path.join(path_to_file,file_name), sep='\t')

#This section is for making the alignment files


#the name of the file or files that contain the results from alignment, in a list
result_file_names = ['results_pdgfra_all_n2/results_Pdgfra-ctrl1', 'results_pdgfra_all_n2/results_Pdgfra-pnxd4', 'results_pdgfra_all_n2/results_Lane5_data', 'results_pdgfra_all_n2/results_Lane6_data', 'results_pdgfra_all_n2/results_Lane7_data', 'results_pdgfra_all_n2/results_Lane8_data']
#create full path to output file
path_to_align=os.path.join(path_to_file, 'results_'+base_name+'_align.p')
#change to True to force creation of new alignment even if file of same name already exists
run_even_if_align_exists = True


def make_align(result_file_names, base_name, path_to_file=path_to_file):
    #path to where the tophat result file or files are located, default is one level up from cuffnorm file
    path = os.path.dirname(path_to_file)
    cell_list =[]
    align_dict =OrderedDict()
    align_dict['input_L_num'] = []
    align_dict['mapped_L_num'] = []
    align_dict['input_R_num'] = []
    align_dict['mapped_R_num'] = []
    align_dict['per_mapped'] = []
    for rf in result_file_names:
        paths_from_list = os.path.join(path, rf)
        for root, dirnames, filenames in os.walk(paths_from_list):
            for filename in fnmatch.filter(filenames, 'align_summary.txt'):
                cell_name = (root.split('/')[-1])
                cell_list.append(cell_name)
                f = open(os.path.join(root,'align_summary.txt'), 'rw')
                for l in f:
                    if 'Left' in l:
                        side_s = 0
                    elif 'Right' in l:
                        side_s = 1
                    if "Input" in l and side_s == 0:
                        input_L_num = int(l.split(':')[-1])
                    if "Mapped" in l and side_s == 0:
                        mapped_L_1 = l.split(':')[-1]
                        mapped_L_num = int(mapped_L_1.split('(')[0].strip())
                    if "Input" in l and side_s == 1:
                        input_R_num = int(l.split(':')[-1])
                    if "Mapped" in l and side_s == 0:
                        mapped_R_1 = l.split(':')[-1]
                        mapped_R_num = int(mapped_R_1.split('(')[0].strip())
                    if "overall read mapping rate." in l:
                        per_mapped = float(l.split('%')[0])

                align_dict['input_L_num'].append(input_L_num)
                align_dict['mapped_L_num'].append(mapped_L_num)
                align_dict['input_R_num'].append(input_R_num)
                align_dict['mapped_R_num'].append(mapped_R_num)
                align_dict['per_mapped'].append(per_mapped)
                f.close()
    align_df = pd.DataFrame(align_dict, index = cell_list)
    align_df.to_csv(os.path.join(path_to_file,'results_'+base_name+'_align.txt'), sep = '\t')

    plt.hist(align_df['mapped_L_num'])
    plt.show()

    with open(os.path.join(path_to_file,'results_'+base_name+'_align.p'), 'wb') as fp:
      pickle.dump(align_df, fp)


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


def filter_by_mapping(path_to_align, cutoff_per_map = 100000, name_filter=False):
    c_to_del =[]
    if path_to_align[-2:] == '.p':
        with open(path_to_align, 'rb') as fp:
            a_data = pickle.load(fp)
    elif path_to_align[-4:] == '.txt':
        a_data = pd.DataFrame.from_csv(path_to_align, sep='\t')
    p_mapped = a_data['mapped_L_num']
    ind_list = p_mapped[p_mapped<cutoff_per_map]
    c_to_del = ind_list.index.values
    if name_filter:
        new_c_to_del = []
        for c in c_to_del:
            c2 = c.replace('-','_')
            nlist = c2.split('_')
            print nlist, 'nlist'
            if nlist[0] == 'pdgfra':
                if nlist[2] == 'ctrl':
                    new_c_to_del.append('Low_ctrl_'+nlist[3])
                elif nlist[2] == 'd4pnx':
                    new_c_to_del.append('Low_pnx_'+nlist[3])
            else:
                if len(nlist[0]) == 2:
                    new_c = 'C0'+nlist[0][-1]
                else:
                    new_c = nlist[0]
                if nlist[2] == 'ctrl1':
                    new_c_to_del.append('ctrl1_'+new_c)
                elif nlist[2] == 'pnxd4':
                    new_c_to_del.append('pnx1_'+new_c)
        return new_c_to_del
    else:
        return c_to_del

def filter_cells_sd(by_cell, cell_list, sd=3.8):
    average_gene_exp = []
    to_delete= []
    for cell_name, genes in zip(cell_list,by_cell):
        gen_exp = (genes >= 1).sum()
        print gen_exp, cell_list.index(cell_name), cell_name
        #if the cell doesn't express at least 500 genes just delete it and exclude from average
        if gen_exp <=1000:
            to_delete.append(cell_list.index(cell_name))
        else:
            average_gene_exp.append(gen_exp)
    print len(to_delete)
    np_av = np.array(average_gene_exp)
    averg = np.average(np_av)
    gene_sd = np.std(np_av)
    print averg, gene_sd
    #add cells that fall outside of stdev value sd
    for i1, exp_level in enumerate(np_av):
        if exp_level < (averg - (gene_sd*sd)) or exp_level > (averg + (gene_sd*sd)):
            if i1 not in to_delete:
                to_delete.append(i1)
    print len(to_delete)
    to_delete1 = sorted(to_delete, reverse = True)
    print len(to_delete1), 'td2'
    print to_delete1
    print len(by_cell)
    for pos in to_delete1:
        print len(cell_list)
        print 'Deleted outlier '+cell_list[pos]
        del cell_list[pos]
    to_delete2 = sorted(to_delete1, reverse = False)
    by_cell = np.delete(by_cell, to_delete2, axis=0)
    print len(by_cell), 'bc2'
    print "Number of cells remaining: "+str(len(cell_list))
    naverage_gene_exp = []
    for ngenes in by_cell:
        ngen_exp = (ngenes >= 1).sum()
        naverage_gene_exp.append(ngen_exp)
    nnp_av = np.array(naverage_gene_exp)
    naverg = np.average(nnp_av)
    ngene_sd = np.std(nnp_av)
    print "New", naverg, ngene_sd
    return cell_list, by_cell


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

def name_filtering(outlier_by_cell, outlier_cell_list):
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
    return outlier_fpkm_dict


#if alignment file doesn't already exist run make alignment file
if not os.path.isfile(path_to_align):
    make_align(result_file_names, base_name)
#if make new file is True
elif run_even_if_align_exists:
    make_align(result_file_names, base_name)
else:
    pass
#create cells to delete list. name_filter is an option for renaming if filenames are different between tophat and cuffnorm
del_list=filter_by_mapping(path_to_align, name_filter=False)
print del_list, 'del'
#create list of genes
gen_list = data.index.tolist()
#cell list
if file_name == 'genes.fpkm_table':
    cell_list = [x[:-2] for x in list(data.columns.values)]
else:
    cell_list = [x for x in list(data.columns.values)]

print cell_list , 'clist'
#make dataframe a numpy array for easy manipulation
npdata = np.array(data.values, dtype='f')
#transpose to delete whole cells
by_cell1 = npdata.transpose()
rem_cell_list, rem_by_cell = delete_cells(by_cell1, cell_list, del_list)
#transpose to delete by genes
npdata2 = rem_by_cell.transpose()
new_gene_list1, new_by_gene = threshold_genes(npdata2, gen_list)
#tranpose to do outlier analysis after deleting cells by mapping and genes without expression
by_cell = new_by_gene.transpose()
outlier_cell_list, outlier_by_cell = filter_cells_sd(by_cell, rem_cell_list)
#tranpose to create final gene oriented matrix
final_by_gene = outlier_by_cell.transpose()
outlier_fpkm_dict = OrderedDict()
bulk_ctrl_dict = OrderedDict()
#detailed funtion to tailor any name filtering. Change to True and edit name_filtering function
filter_on_lane = False
#function for removing cells by name, will be processed as a seperate matrix file
bulk = True
if filter_on_lane:
    outlier_fpkm_dict = name_filtering(outlier_by_cell, outlier_cell_list)
else:
    print len(outlier_by_cell)
    print len(outlier_cell_list)
    for i, l in enumerate(outlier_by_cell):
        print i
        if bulk:
            cell_name = outlier_cell_list[i]
            #create for exclusion from final matrix and inclusion into bulk matrix
            if 'bulk' in cell_name or '+' in cell_name or 'neg' in cell_name or '-' in cell_name:
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
