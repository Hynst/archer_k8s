library(data.table)

args <- commandArgs(trailingOnly=TRUE)
run_ID <- args[1]
#data_dir = "/bibs/data/myelo/"
data_dir = args[2]

samples_list <- as.list(list.files(paste0(data_dir,run_ID,"/samples/")))

sample_var_list <-list()

sample_var_list <- lapply(samples_list, function(x) { 
  sample_var <- fread(paste0(data_dir,run_ID,"/samples/",x,"/variants/mutect2/",x,".mutect2.filt.norm.vep.csv"))
  sample_var <- cbind(sample = paste0(x), run_ID = paste0(run_ID), sample_var)
  return(sample_var) 
}
)

merge_var_df <- rbindlist(sample_var_list) 
print(head(merge_var_df))
merge_var_df[,Occurence_in_samples := .N, by = c("Start_Position","HGVSc")]

#write.table(merge_var_df, file=paste0("/bibs/data/myelo/", run_ID, "/", run_ID, ".merged_variant.table_NEW.tsv"), 
#            sep="\t", row.names = F, quote = F)

write.table(merge_var_df, file=paste0(data_dir, run_ID, "/", run_ID, ".merged_variant.table_NEW.tsv"),
            sep="\t", row.names = F, quote = F)
