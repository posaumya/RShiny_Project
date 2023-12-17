library('tidyverse')
library('fgsea')

id2_gene_file <- "data/human_id2gene.txt"
id2gene_df <- read.table(id2_gene_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE, col.names = c("gene_symbol","ensembl_id"))

gene_set_file <- "data/c3.tft.v2023.2.Hs.symbols.gmt"

deseq_file<-"data/GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.txt"
deseq_df <- read.csv(deseq_file,stringsAsFactors = FALSE, sep="\t", header = TRUE)


make_ranked_log2fc <- function(labeled_results, id2gene_df) {
  merged_data <- merge(labeled_results, id2gene_df, by.x = "symbol", by.y = "ensembl_id", all.x = TRUE)
  ranked_log2fc_df <- merged_data[order(merged_data$log2FoldChange, decreasing = TRUE), ]
  ranked_log2fc_named <- setNames(ranked_log2fc_df$log2FoldChange, ranked_log2fc_df$symbol)
  return(ranked_log2fc_named)
}

rnk_list <- make_ranked_log2fc(deseq_df,id2gene_df)

run_fgsea <- function(gmt_file_path, rnk_list, min_size, max_size) {
  rnk_list <- rnk_list[is.finite(rnk_list)]
  
  gene_sets <- gmtPathways(gmt_file_path)
  fgsea_results <- fgsea(pathways = gene_sets, 
                         stats = rnk_list, 
                         minSize = min_size, 
                         maxSize = max_size)
  
  return(as_tibble(fgsea_results))
}

fgsea_result <- run_fgsea(gene_set_file,rnk_list,15,1115)
# Convert list column to character
fgsea_result$leadingEdge <- sapply(fgsea_result$leadingEdge, paste, collapse = ", ")

print(head(fgsea_result))

write.csv(fgsea_result, "data/fgsea_results.csv", row.names = FALSE)
