gtfToGenePred -genePredExt -geneNameAsName2 /Volumes/Seq_data/hg19_ERCC_bt2/Annotation/hg19_ERCC.gtf /Volumes/Seq_data/refFlat.tmp.txt
paste <(cut -f 12 /Volumes/Seq_data/refFlat.tmp.txt) <(cut -f 1-10 /Volumes/Seq_data/refFlat.tmp.txt) > /Volumes/Seq_data/refFlat_mm10ERS.txt
gzip /Volumes/Seq_data/refFlat_mm10ERS.txt
