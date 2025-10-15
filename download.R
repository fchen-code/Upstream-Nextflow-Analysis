nums<-c(as.character(as.roman(1:16)),"Mito")
heading <- ("https://ftp.ensembl.org/pub/release-115/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.")
ending<-(".fa.gz")

for (num in nums) {
  url <- paste0(heading, num, ending)
  cmd <- paste("wget", url)
  system(cmd)
}


base = "Saccharomyces_cerevisiae.R64-1-1.dna.chromosome."
for (num in nums) {
  file <- paste0(base, num,".fa.gz")
  cmd <- paste0("zcat ", file, ">> Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.primary.assembly.fa")
  system(cmd)
}

folder <-  c("mkdir refs", "mv *primary.assembly.fa refs","mkdir data", 
	"fasterq-dump -p --split-files SRR13978643", "mv SRR13978643* data",
	"mkdir gtf", 
	"wget http://ftp.ensembl.org/pub/release-108/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.108.gtf.gz",
	"mv *gtf.gz gtf",
	"gunzip -k gtf/*.gz") 

for (i in 1:length(folder)) {
	system(folder[i])
}
system("rm *.fa.gz")
