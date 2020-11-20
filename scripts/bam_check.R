bam <- snakemake@input[["bam"]]
outname <- snakemake@output[[1]]

log_vector <- c()
check_vector <- c()
for(i in bam){
	cmd <- paste0("samtools quickcheck -q ",i)
	if(system(cmd) == 0){
		log_vector <- append(log_vector, paste0("BAM valid - ",i))
		check_vector <- append(check_vector,FALSE)
	} else {
		log_vector <- append(log_vector,paste0("BAM invalid or missing - ",i))
       		check_vector <- append(check_vector,TRUE)
	}
}

if(any(check_vector)){
	outname <- gsub(pattern = "ok",replacement = "invalid",outname)
	writeLines(text = as.character(log_vector),con = outname)	
} else {
	writeLines(text = as.character(log_vector),con = outname)
}
