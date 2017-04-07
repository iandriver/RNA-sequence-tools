gene_svd_top = TruncatedSVD(n_components=5, n_iter=6)
    cell_svd_top = TruncatedSVD(n_components=5, n_iter=6)
    cell_svd_all = TruncatedSVD(n_components=5, n_iter=6)
    gene_svd_all = TruncatedSVD(n_components=5, n_iter=6)
    gene_svd_top.fit(np_top_gene)
    cell_svd_top.fit(np_top_cell)
    cell_svd_all.fit(np_by_cell)
    gene_svd_all.fit(np_by_gene)

    var_ratio_cell_all = cell_svd_all.explained_variance_ratio_ *100
    cum_var_exp_cell_all = np.cumsum(var_ratio_cell_all)
    var_ratio_cell_top = cell_svd_top.explained_variance_ratio_ *100
    cum_var_exp_cell_top = np.cumsum(var_ratio_cell_top)
    with plt.style.context('seaborn-whitegrid'):
        sns.set(context= 'poster', font_scale = 0.9)
        plt.figure(figsize=(6, 4))

        plt.bar(range(1,6), var_ratio_cell_all, alpha=0.5, align='center',
                label='all explained variance', color='b')
        plt.step(range(1,6), cum_var_exp_cell_all, where='mid',
                 label='cumulative explained variance all', color='b')
        plt.bar(range(1,6), var_ratio_cell_top, alpha=0.5, align='center',
                label='top explained variance', color='r')
        plt.step(range(1,6), cum_var_exp_cell_top, where='mid',
                 label='cumulative explained variance top', color='r')
        plt.ylabel('Explained variance ratio')
        plt.xlabel('Principal components')
        plt.legend(loc='best')
        plt.tight_layout()
        plt.savefig(os.path.join(matrix_data.new_filepath,'Cell_varience_explained.'+args.image_format), bbox_inches='tight')
        plt.close('all')
    var_ratio_gene_all = gene_svd_all.explained_variance_ratio_ *100
    cum_var_exp_gene_all = np.cumsum(var_ratio_gene_all)
    var_ratio_gene_top = gene_svd_top.explained_variance_ratio_ *100
    cum_var_exp_gene_top = np.cumsum(var_ratio_gene_top)
    with plt.style.context('seaborn-whitegrid'):
        sns.set(context= 'poster', font_scale = 0.9)
        plt.figure(figsize=(6, 4))

        plt.bar(range(1,6), var_ratio_gene_all, alpha=0.5, align='center',
                label='all explained variance ', color='b')
        plt.step(range(1,6), cum_var_exp_gene_all, where='mid',
                 label='cumulative explained variance all', color='b')
        plt.bar(range(1,6), var_ratio_gene_top, alpha=0.5, align='center',
             label='top explained variance', color='r')
        plt.step(range(1,6), cum_var_exp_gene_top, where='mid',
                          label='cumulative explained variance top', color='r')
        plt.ylabel('Explained variance ratio')
        plt.xlabel('Principal components')
        plt.legend(loc='best')
        plt.tight_layout()
        plt.savefig(os.path.join(matrix_data.new_filepath,'Gene_varience_explained.'+args.image_format), bbox_inches='tight')
        plt.close('all')
