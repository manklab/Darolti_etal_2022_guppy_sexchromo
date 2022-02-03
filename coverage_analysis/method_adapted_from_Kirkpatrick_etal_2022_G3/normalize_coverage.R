## This script normalizes coverage data based on the median genomic coverage of each individual, in order to account for differences in overall coverage between individuals

rm(list=ls())
ls() 

#Load data
infiles <- Sys.glob("./reads/*multicov.txt")
for (infile in infiles) {
	sample <- read.csv(infile, sep="\t", header=F)
	#sort file by chromosome and start position
	sample_sorted <- sample[order(sample$V1,sample$V2),]
	#normalize by median coverage
	sample_sorted$norm <- sample_sorted$V4/median(sample_sorted$V4)
	#write output file of normalized coverage
	output_file <- paste0(tools::file_path_sans_ext(infile),"_normalized.txt")
	write.table(sample_sorted, output_file, quote=F, sep="\t")
}
