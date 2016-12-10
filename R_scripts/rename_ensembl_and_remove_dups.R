rename_and_remove_dups <- function(df, keep_highest, species){

  library(biomaRt)
  if(species == 'mouse'){
    ensembl = useEnsembl(biomart="ensembl", dataset = 'mmusculus_gene_ensembl')
  }
  else{
    ensembl = useEnsembl(biomart="ensembl", dataset = 'hsapiens_gene_ensembl')
  }
  g_list <- getBM(filters= "ensembl_gene_id", attributes= c("external_gene_name", "ensembl_gene_id"),values=rownames(df),mart= ensembl)
  rename_df1 <- merge(as.data.frame(df),g_list,by.x="row.names",by.y="ensembl_gene_id")

  if(keep_highest==TRUE){
    dups <- rename_df1[duplicated(rename_df1[,'external_gene_name']),]
    dup_list <- unlist(as.list(dups$external_gene_name))
    if (length(dup_list) > 0){
      dup_df <- rename_df1[rename_df1$external_gene_name %in% dup_list,]
      i <- 1
      del_list <- c()
      for(g in dup_list){
        g_dup_df <- dup_df[dup_df$external_gene_name == g,]
        g_mean <- rowMeans(subset(g_dup_df, select = colnames(df), na.rm = TRUE))
        non_max_g <- as.integer(names(g_mean[names(g_mean) != names(which.max(g_mean))]))
        for(g_i in non_max_g){
          del_list[[i]] <- g_i
          i <- i + 1
        }

      }
      final_df <- rename_df1[-del_list,]
    }
    else{
      final_df <- rename_df1
    }
  }
  else{
    final_df <- rename_df1[!duplicated(rename_df1[,'external_gene_name']),]
  }
  row.names(final_df) <- final_df$external_gene_name
  final_df <- final_df[,!(names(final_df) %in% c('external_gene_name','Row.names'))]


return(final_df)
}
