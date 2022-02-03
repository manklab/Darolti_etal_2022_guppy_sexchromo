## This script normalizes snp density data based on the median snp density across the genome of each individual, in order to account for differences in overall coverage between individuals

rm(list=ls())
ls() 

#Load data
infiles <- Sys.glob("./reads/*.bygene_density")
for (infile in infiles) {
	sample <- read.csv(infile, sep="\t", header=F)
	#sort file by chromosome and position
	sample_sorted <- sample[order(sample$V1,sample$V7),]
	#normalize by median snp density
	sample_sorted$norm <- sample_sorted$V5/median(sample_sorted$V5)
	#write output file of normalized coverage
	output_file <- paste0(tools::file_path_sans_ext(infile),"_normalized.txt")
	write.table(sample_sorted, output_file, quote=F, sep="\t")
}