setwd("H:/Documents/03_Syntenic_linkage/01_Data/OrthoFinder/input/")

## get species 
species <- list.files()[grep(list.files(), pattern = ".fasta")]
length(species) ## 126

## iterate through species
for(i in 1:length(species)) {
    
    ## truncate .fasta and 2 from species
    species_name <- strsplit(species[i], split = ".fasta")[[1]]
    
    ## load fasta
    fasta <- read.table(species[i], stringsAsFactors = FALSE)
    inds <- grep(fasta[,1], pattern = ">")
    fasta[inds, ] <- gsub(">", paste(">", species_name, "_", sep = ""), 
        fasta[inds, ])
    
    ## write fasta 
    write.table(fasta, file = paste("../", species[i], sep = ""), 
        col.names = FALSE, row.names = FALSE, quote = FALSE)
}
## no gff file for cycmi, equgi, pinpi, pinsy, taxba --> remove from analysis and remove manually

## in orysa rename orysa_ChrUn and orysa_ChrSy
j <- which(species == "orysa.fasta")

## truncate .fasta and 2 from species
species_name <- strsplit(species[j], split = ".fasta")[[1]]
species_name <- strsplit(species_name, split = "2")[[1]]

## load fasta 
fasta <- read.table(species[j], stringsAsFactors = FALSE)
inds <- grep(fasta[,1], pattern = ">")
fasta[inds, ] <- gsub(">", paste(">", species_name, "_", sep = ""), fasta[inds, ])
fasta[which(fasta[, 1] == ">orysa_ChrUn"),] <- paste(">orysa_ChrUn", 
    1:length(which(fasta[,1] == ">orysa_ChrUn")), sep = "_")
fasta[which(fasta[,1] == ">orysa_ChrSy"),] <- paste(">orysa_ChrSy", 
    1:length(which(fasta[,1] == ">orysa_ChrSy")), sep = "_")

## write fasta 
write.table(fasta, file = paste("../", species[j], sep = ""), 
    col.names = FALSE, row.names = FALSE, quote = FALSE)


## maldo contains some double sequences, remove the second duplicate
j <- which(species == "maldo.fasta")

## truncate .fasta and 2 from species
species_name <- strsplit(species[j], split = ".fasta")[[1]]
species_name <- strsplit(species_name, split = "2")[[1]]

## load fasta 
fasta <- read.table(species[j], stringsAsFactors = FALSE)
inds <- grep(fasta[, 1], pattern = ">")
duplicates <- fasta[inds[duplicated(fasta[inds, ])], ]

for (k in 1:length(duplicates)) {
    ## take 2nd duplicated
    ind_2 <- which(fasta[inds, ] == duplicates[k])[2] 
    
    ## replace by ""
    fasta[inds[ind_2]:c(inds[c(ind_2 + 1)]-1), ] <- "" 
}

fasta[inds, ] <- gsub(">", 
    paste(">", species_name, "_", sep = ""), fasta[inds, ])

## remove lines with ""
fasta <- fasta[-which(fasta[,1] == ""), ]

## write fasta
write.table(fasta, file = paste("../", species[j], sep = ""), col.names = FALSE, row.names = FALSE, quote = FALSE)
