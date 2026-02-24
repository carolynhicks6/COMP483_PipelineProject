library(sleuth)
library(dplyr)

# use samples and conditions to construct table
samples <- c("SRR5660030", "SRR5660033", "SRR5660044", "SRR5660045")
condition <- c("2dpi", "6dpi", "2dpi", "6dpi")

# get directories of files listed in snakemake input
kal_dirs <- unique(dirname(snakemake@input$abundance))

# make auxillary table
s2c <- data.frame(
    sample = samples,
    condition = condition,
    path = kal_dirs,
    stringsAsFactors = FALSE)

# sleuth prep
so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE)
so = sleuth_fit(so, ~condition, 'full')
so = sleuth_fit(so, ~1, 'reduced')
so = sleuth_lrt(so, 'reduced', 'full')

# make table of the sleuth results
sleuth_table = sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE) 

# filter out significant results
sleuth_significant = dplyr::filter(sleuth_table, qval <= 0.05) |> dplyr::arrange(pval) 

# write significant results to tsv
write.table(sleuth_significant, file=snakemake@output[[1]], sep="\t",quote = FALSE,row.names = FALSE)




