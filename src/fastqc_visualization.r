if(!require(optparse)){install.packages("optparse",repos = "http://cran.us.r-project.org")}
if(!require(tidyverse)){install.packages("tidyverse",repos = "http://cran.us.r-project.org")}


#load libraries
if (!requireNamespace("optparse", quietly = TRUE)) {
  write("R package 'optparse' needed for qc visualization. Please install it.\n",
       stderr())
  install.packages("optparse")
}  
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  stop("R package 'tidyverse' needed for qc visualization. Please install it.\n",
       stderr())
  install.packages("tidyverse")
}  


##optparser options
option_list <- list(
  make_option(c("-i", "--input_folder"), type="character", default=getwd(),
              help="The path to thefolder with the fastqc results"),
  make_option(c("-o", "--output_file"), type="character", default="QC_plots.pdf",
              help="Give the name of the pdf file where the plots are to be saved."),
  make_option(c("-p", "--print"), type="logical", default=TRUE,
              help="Print sample ids of samples that failed QC.")
  
)

opt <- parse_args(OptionParser(option_list=option_list))

#load fastQC summaries and create per test table
inp <- list.files(opt$input_folder, pattern = ".zip")


fastqc_results <- lapply(inp, function(k){
  unzip(paste(opt$input_folder, k, sep = "/"),exdir = opt$input_folder)
  inpu <- read_delim(paste(paste(gsub(".zip", "", paste(opt$input_folder,k, sep = "/"))), 
                           "summary.txt", sep = "/"), delim = "\t")
  out <- as_data_frame(t(inpu[, 1])) %>%
    mutate(sample.id = names(inpu)[3])
  names(out) <- c(gsub(" ", "_", unlist(inpu[,2])), "sample_id")
  unlink(x = paste(opt$input_folder, gsub(".zip", "", k), sep = "/"), 
         recursive = T, force = T)
  
  return(out)
})

outp <- do.call("rbind.data.frame", fastqc_results)%>%
  select(ID = sample_id,
         PBQ = Per_base_sequence_quality,
         PTQ = Per_tile_sequence_quality,
         PSQ = Per_sequence_quality_scores,
         PBC = Per_base_sequence_content,
         SGC = Per_sequence_GC_content,
         PBN = Per_base_N_content,
         SLD = Sequence_Length_Distribution,
         SDL = Sequence_Duplication_Levels,
         ORS = Overrepresented_sequences,
         AdC = Adapter_Content,
         KmC = Kmer_Content)

#change table format
ret <- outp %>% 
  group_by(ID) %>%
  gather(test, status, PBQ:KmC)

#plot how many samples failed the test
qc.fail <- ggplot()+
  geom_bar(data = ret, aes(x = test, fill = status), stat = 'count', position = 'dodge')+
  theme_bw()

#plot which sample failed which test
qc.samples <- ggplot()+
  geom_tile(data = ret, aes(y = ID, x = test, fill = as.factor(status)))+
  scale_fill_discrete(name = "status")+
  xlab("FastQC test")+
  ylab("Samples")+
  theme_bw()+
  theme(
    axis.text.y = element_blank()
  )

#plot pdf
pdf(opt$output_file)
print(qc.fail)
print(qc.samples)
dev.off()

png(gsub(".pdf", "1.png", opt$output_file))
print(qc.fail)
dev.off()

png(gsub(".pdf", "2.png", opt$output_file))
print(qc.samples)
dev.off()

#table with samples that faild a test
fail <- ret %>%
  filter(status == "FAIL")

#get the ID number of the failed samples
fail.samp <- fail %>%
  filter(!duplicated(ID)) %>%
  select(ID)%>%
  unlist() %>%
  parse_number()%>%
  unique() %>%
  sort()

if(opt$print){
  write(sprintf("The following sample failed at least one test: %s \n", fail.samp), stdout())
}
