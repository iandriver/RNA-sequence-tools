import os
import sys
import glob
import pandas as pd


def fpkm_filemerge(sort_on=False, sort_key= '', sort_cat = 'Sample_age'):
  #fpkm dictionary stored under cell names: Sample_age+'_'+repeat+'_'+sample_name
  path = input('Filepath to fpkm tracking files:')
  fpkm_matrix_dict = {}
  cell_names = []
  age_list =[]
  repeat_list = []
  sample_name_list = []
  n = 0
  for filename in glob.glob(os.path.join(path, '*.fpkm_tracking')):
    #example filename: GSM1271862_E18_1_C08_IL3230.sorted.genes.fpkm_tracking
    real_filename= filename.split('/')[-1]
    name_parse = real_filename.split('_')
    GSM_name = name_parse[0].strip('.').strip('/')
    Sample_age = name_parse[1]
    repeat = name_parse[2]
    sample_name = name_parse[3]
    if sort_on == True:
      if locals().get(sort_cat, '') == sort_key:
        print '1'
        cell_names.append(Sample_age+'_'+repeat+'_'+sample_name)
        age_list.append(Sample_age)
        repeat_list.append(repeat)
        sample_name_list.append(sample_name)

        file_header =['tracking_id', 'class_code', 'nearest_ref_id', 'gene_id', 'gene_short_name', 'tss_id', 'locus', 'length', 'coverage', 'q0_FPKM', 'q0_FPKM_lo', 'q0_FPKM_hi', 'q0_status']
        curr_file = open(filename, 'rw')
        curr_cell_fpkm = []
        if n != 1:
        	curr_cell_genes = []
        for line in curr_file:
        	curr_line = line.strip('\n').split('\t')
          #exclude RNA spike controls labeled as ERCC
        	if n == 0 and curr_line[3][:4] != 'ERCC':
        		curr_cell_genes.append(curr_line[3])
        	if curr_line[3][:4] != 'ERCC':
        		curr_cell_fpkm.append(curr_line[9])
        fpkm_matrix_dict[Sample_age+'_'+repeat+'_'+sample_name] = curr_cell_fpkm

        n=1
    elif sort_on == False:
      cell_names.append(Sample_age+'_'+repeat+'_'+sample_name)
      age_list.append(Sample_age)
      repeat_list.append(repeat)
      sample_name_list.append(sample_name)
      #not currently used but useful headers for fpkm tracking files
      file_header =['tracking_id', 'class_code', 'nearest_ref_id', 'gene_id', 'gene_short_name', 'tss_id', 'locus', 'length', 'coverage', 'q0_FPKM', 'q0_FPKM_lo', 'q0_FPKM_hi', 'q0_status']
      curr_file = open(filename, 'rw')
      curr_cell_fpkm = []
      if n != 1:
        curr_cell_genes = []
      for line in curr_file:
        curr_line = line.strip('\n').split('\t')
        #exclude RNA spike controls
        if n == 0 and curr_line[3][:4] != 'ERCC':
          curr_cell_genes.append(curr_line[3])
        if curr_line[3][:4] != 'ERCC':
          curr_cell_fpkm.append(curr_line[9])
      fpkm_matrix_dict[Sample_age+'_'+repeat+'_'+sample_name] = curr_cell_fpkm

      n=1
  annot_cells_dict = { 'Animal_age' : age_list, 'Repeat#' : repeat_list, 'cell_name' : sample_name_list}
  fpkm_df = pd.DataFrame(fpkm_matrix_dict, index = curr_cell_genes)
  annot_cells_df = pd.DataFrame(annot_cells_dict, index = cell_names)
  fpkm_df.to_csv('fpkm_matrix_'+sort_key+'.txt', sep = '\t')
  annot_cells_df.to_csv('cell_sample_sheet_'+sort_key+'.txt', sep = '\t')


fpkm_filemerge(sort_on=True, sort_key='E18')
