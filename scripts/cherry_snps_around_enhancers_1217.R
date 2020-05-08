library("ggplot2")

workingDir = "/data/atimms/cherry_snps_around_enhancers_1217";
setwd(workingDir);

tests = c("gnomad.combined_counts.0.001maf.counts", 
          "gnomad.combined_counts.0.01maf.counts",
          "gnomad.combined_counts.0.05maf.counts",
          "kaviar.combined_counts.0.001maf.counts",
          "kaviar.combined_counts.0.01maf.counts",
          "kaviar.combined_counts.0.05maf.counts",
          "1000g.combined_counts.0.001maf.counts",
          "1000g.combined_counts.0.01maf.counts",
          "1000g.combined_counts.0.05maf.counts")
for(i in 1:length(tests)){
  txt_file = paste(tests[i],"txt",sep=".")
  count_data = read.table(txt_file, header=T)
  ggplot(data=count_data, aes(x=window, y=count, group=variant_type, colour=variant_type))+ geom_line()
  pdf_file <- paste(tests[i],"pdf",sep=".")
  ggsave(file=pdf_file)
}

tests = c("gerp_window_averages",
          "hg38.phastCons7way.comb_results",
          "hg38.phastCons100way.comb_results",
          "hg38.phyloP100way.comb_results")
for(i in 1:length(tests)){
  txt_file = paste(tests[i],"txt",sep=".")
  count_data = read.table(txt_file, header=T)
  ggplot(data=count_data, aes(x=window, y=count))+ geom_line()
  pdf_file <- paste(tests[i],"pdf",sep=".")
  ggsave(file=pdf_file)
}