import os
import pandas as pd
import shutil
import sys
import hashlib

def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


pats = ['/Volumes/Seq_data/09142015_BU3_ips17_raw-data/BU3_fastq_renamed']
keep_file = '/Volumes/Drobo/Seq_data/rsem_results_h38_ips17_BU3/Samplelist_cells_new.txt'
insert_size_file = '/Volumes/Seq_data/rsem_results_h38_ips17_BU3_insert_metrics_picard_insert_metrics.txt'
new_filename = 'ips17_BU3_hg19_rsem_geo'
nextera_barcode_v = 1
py3 = sys.version_info[0] > 2
try:
    new_filename = os.path.join(os.path.dirname(pats[0]),new_filename)
    print("new_filename", new_filename)
    os.mkdir(new_filename)
except FileExistsError:
    if py3:
      response = input("New file already exists. Do you still want to continue (Y/N): ")
    else:
      response = raw_input("New file already exists. Do you still want to continue (Y/N): ")
    if response != 'Y':
        sys.exit("Canceled, please rerun with different new_filename.")

keep_df = pd.read_table(open(keep_file,'rU'), sep='\s+', engine='python')
insert_df = pd.read_table(open(insert_size_file,'rU'), sep='\t', engine='python', index_col=0)
samples_to_keep= keep_df["SampleID"].tolist()
#cells_to_keep = [x.split('_')[-2] for x in samples_to_keep]
sample_list_r1 = []
sample_list_r2 = []
md5_r1 = []
md5_r2 = []
sample_folder = []
index_list = []


for path_to_file in pats:
    for root, dirnames, filenames in os.walk(path_to_file):
        for f in filenames:
            if 'fastq.gz' in f:
                cell_num = os.path.basename(root)
                if cell_num in samples_to_keep:
                    new_name = os.path.basename(root)
                    try:
                        shutil.copy2(root+'/'+f,new_filename)
                    except FileExistsError:
                        print("alread moved "+new_name)
                    if "R1" in f:
                        sample_list_r1.append(f)
                        sample_folder.append(new_filename+'/'+new_name)
                        md5_r1.append(md5(os.path.join(new_filename,f)))
                    elif "R2" in f:
                        sample_list_r2.append(f)
                        index_list.append(cell_num)
                        md5_r2.append(md5(os.path.join(new_filename,f)))

sample_df = pd.DataFrame.from_dict({'index':index_list,'R1_md5':md5_r1, 'R2_md5':md5_r2, 'folder':sample_folder,'PE1':sample_list_r1, 'PE2':sample_list_r2})
sample_df.set_index('index', inplace=True)
filter_insert_df = insert_df.loc[index_list]
all_df = pd.concat([sample_df,filter_insert_df[['MEAN_INSERT_SIZE','STANDARD_DEVIATION']]], axis=1)

all_df.to_csv(os.path.join(path_to_file,os.path.basename(path_to_file)+'_samples.txt'), sep='\t')
