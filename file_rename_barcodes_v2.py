import subprocess
import os
import fnmatch
import pandas as pd
import sys
import shutil


def convert_bytes(num):
    """
    this function will convert bytes to MB.... GB... etc
    """
    for x in ['bytes', 'KB', 'MB', 'GB', 'TB']:
        if num < 1024.0:
            return "%3.1f %s" % (num, x)
        num /= 1024.0


def file_size(file_path):
    """
    this function will return the file size
    """
    if os.path.isfile(file_path):
        file_info = os.stat(file_path)
        return convert_bytes(file_info.st_size)

py3 = sys.version_info[0] > 2 #creates boolean value for test that Python major version > 2


result_file_names = ['Cfms_Pnx_d7_10242014']
path_to_file = '/Volumes/Drobo/Seq_data'
pats = []
new_filename = 'Cfms_Pnx_d7_10242014_renamed'
nextera_barcode_v = 1

try:
    new_filename = os.path.join(path_to_file,new_filename)
    os.mkdir(new_filename)
except FileExistsError:
    if py3:
      response = input("New file already exists. Do you still want to continue (Y/N): ")
    else:
      response = raw_input("New file already exists. Do you still want to continue (Y/N): ")
    if response != 'Y':
        sys.exit("Canceled, please rerun with different new_filename.")

if nextera_barcode_v == 1:
    barcode_map1 = pd.read_csv('/Volumes/Seq_data/v1-barcodes.txt', sep='\t', header = None, names = ['name', 'code'])
elif nextera_barcode_v == 2 :
    barcode_map1 = pd.read_csv('/Volumes/Seq_data/v2-barcodes.txt', sep='\t', header = None, names = ['name', 'code'])

map1_dict = {}
barcode_dict = {'TAAGGCGA':'N701','CGTACTAG':'N702','AGGCAGAA':'N703', 'TCCTGAGC':'N704',
                'GGACTCCT':'N705', 'TAGGCATG':'N706','CTCTCTAC':'N707', 'CAGAGAGG':'N708',
                'GCTACGCT':'N709', 'CGAGGCTG':'N710', 'AAGAGGCA':'N711', 'GTAGAGGA':'N712',
                'GCTCATGA':'N714', 'ATCTCAGG':'N715',
                'TAGATCGC':'S501', 'CTCTCTAT':'S502', 'TATCCTCT':'S503', 'AGAGTAGA':'S504',
                'GTAAGGAG':'S505', 'ACTGCATA':'S506', 'AAGGAGTA':'S507', 'CTAAGCCT':'S508',
                'CGTCTAAT':'S510', 'TCTCTCCG':'S511', 'TCGACTAG':'S513', 'GCGTAAGA':'S517'}
#nextera V1 index map
N_index1 = {'N701':'1', 'N702':'2', 'N703':'3', 'N704':'4', 'N705':'5', 'N706':'6', 'N707':'7', 'N708':'8',
            'N709':'9', 'N710':'10', 'N711':'11', 'N712':'12'}
S_index1 = {'S501':'A', 'S502':'B', 'S503':'C', 'S504':'D', 'S505':'E', 'S506':'F', 'S507':'G', 'S508':'H'}

#nextera V2 index map
N_index2 = {'N701':'1', 'N702':'2', 'N703':'3', 'N704':'4', 'N705':'5', 'N706':'6', 'N707':'7', 'N710':'8',
            'N711':'9', 'N712':'10', 'N714':'11', 'N715':'12'}
S_index2 = {'S502':'A', 'S503':'B', 'S505':'C', 'S506':'D', 'S507':'E', 'S508':'F', 'S510':'G', 'S511':'H'}

for n, c in zip(barcode_map1['name'].tolist(),barcode_map1['code'].tolist()):
    map1_dict[c]= n

for fname in result_file_names:
    pats.append(os.path.join(path_to_file, fname))
for p in pats:
    for root, dirnames, filenames in os.walk(p):
        x=0
        if dirnames == []:
            fastq_files = [files for files in filenames if '.fastq.gz' in files]

            name = os.path.basename(root)

            #subprocess.call('mkdir -p '+path+'/%s' % name, shell=True)
            if len(fastq_files) > 2:
                r1_files = [r1 for r1 in fastq_files if r1.strip('.fastq.gz').split('_')[-2] == 'R1']
                r1_files.sort()
                r2_files = [r2 for r2 in fastq_files if r2.strip('.fastq.gz').split('_')[-2] == 'R2']
                r2_files.sort()
                r1_cat_cmd = 'cat '
                r2_cat_cmd = 'cat '
                for r1,r2 in zip(r1_files,r2_files):
                    r1_cat_cmd+= os.path.join(root,r1)+' '
                    r2_cat_cmd+= os.path.join(root,r2)+' '
                r1_cat_cmd+= '> '+os.path.join(root,r1_files[0]+".temp")
                r2_cat_cmd+= '> '+os.path.join(root,r2_files[0]+".temp")
                print(r1_cat_cmd)
                if os.path.isfile(os.path.join(root,r1_files[0]+".temp")):
                    print(os.path.join(root,r1_files[0]+".temp"+' exists.'))
                    pass
                else:
                    process = subprocess.Popen(r1_cat_cmd, stdout=subprocess.PIPE, shell=True)
                    out1, err1 = process.communicate()
                    print(out1)
                    print(r2_cat_cmd)
                    process = subprocess.Popen(r2_cat_cmd, stdout=subprocess.PIPE, shell=True)
                    out2, err2 = process.communicate()
                    print(out2)
                pair_1 = os.path.join(root,r1_files[0]+".temp")
                pair_2 = os.path.join(root,r2_files[0]+".temp")
                shutil.move(os.path.join(root,r1_files[0]+".temp"), os.path.join(root,r1_files[0]))
                shutil.move(os.path.join(root,r2_files[0]+".temp"), os.path.join(root,r2_files[0]))
                for i in range(1,len(r1_files)):
                    os.remove(os.path.join(root,r1_files[i]))
                    os.remove(os.path.join(root,r2_files[i]))
            if len(fastq_files) == 2:
                for files in fastq_files:
                    num = files.strip('.fastq.gz').split('_')[-2]
                    if num == 'R1':
                        pair_1 = os.path.join(root,files)
                    if num == 'R2':
                        pair_2 = os.path.join(root,files)
            for filename in fastq_files:
                print(file_size(os.path.join(root,filename)),os.path.join(root,filename))
                if x==0:
                    g_cell_title = filename.split('_')[1]
                    barcode1 = g_cell_title[:8]
                    barcode2 = g_cell_title[8:]
                    code_lookup1 = barcode_dict[barcode1]
                    code_lookup2 = barcode_dict[barcode2]
                    code_lookup = code_lookup2+'/'+code_lookup1
                    skip= False
                    try:
                        cell_name = map1_dict[code_lookup]
                        if nextera_barcode_v == 1:
                            nindex = N_index1[code_lookup1]
                            sindex = S_index1[code_lookup2]
                        elif nextera_barcode_v == 2:
                            nindex = N_index2[code_lookup1]
                            sindex = S_index2[code_lookup2]
                    except KeyError:
                        print(code_lookup, 'file:',root, 'BarcodeError')
                        skip = True
                    if not skip:
                        print(code_lookup, 'present')
                        name = sindex+nindex+'_'+cell_name

                        #uncomment when you actually want to move
                        #os.mkdir(os.path.join(new_filename,name))
                        shutil.copytree(root+'/',os.path.join(new_filename,name))
                        print(root+'/',os.path.join(new_filename,name))
                        x+=1
