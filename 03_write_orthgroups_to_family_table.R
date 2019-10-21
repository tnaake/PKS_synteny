setwd("/home/mpimp-golm.mpg.de/naake/winhome/Documents/03_Syntenic_linkage/01_Data/OrthoFinder/input/Results_Oct26/")

type <- "OrthoFinder" ## or MCL

if (type == "OrthoFinder") {
    og <- read.table("Orthogroups.txt", sep = " ", fill = TRUE, 
        stringsAsFactors = FALSE, header = FALSE)
    unassigned <- read.table("Orthogroups_UnassignedGenes.csv", sep = "\t", 
        dec = ".", fill = TRUE, stringsAsFactors = FALSE, header = TRUE)    
}

if (type == "MCL")  {
    ## start reading for mcl
    og <- read.table("mcl_families.processed.txt", sep = "\t", fill = TRUE, 
        stringsAsFactors = FALSE, header = FALSE)
    ## some genes are missing in mcl output, get these files from 
    ## Orthogroups_unassigned
    unassigned <- read.table("Results_Oct26/Orthogroups_UnassignedGenes.csv", 
        sep = "\t", dec = ".", fill = TRUE, stringsAsFactors = FALSE, 
        header = TRUE)
    ## write value to second column and switch column 1 and 2, that it is 
    ## compatible with the output of the script  
    for (i in 2:dim(unassigned)[2]) {
        unassigned[which(unassigned[,i] != ""), 2] <- unassigned[which(unassigned[,i] != ""), i]
    }
    unassigned <- unassigned[, 2:1]
    rownames(unassigned) <- NULL
    unassigned <- as.matrix(unassigned)

}

## create result file
family_file <- matrix(ncol = 2, nrow = 0)
family_all <- unlist(strsplit(as.character(og[, 1]), split = ":"))

dim_og <- dim(og)[1]


for (i in 1:dim_og) { 
    
    ## get the protein family
    if (type == "OrthoFinder") {
        family <- family_all[i]
    }
    if (type == "MCL") {
        family <- paste("group", i, sep = "_") 
    }
   
    ## set manually
    ## get the index of the last element, to speed up partition: 
    ## adjust the indices according to the data that it gets the filled 
    ## entries per line
    if (i <= 600) {
        og_i <- as.character(og[i, ])
        ind_last <- which(og_i == "")[1] - 1
    }
    if (i > 600 & i <= 2000) {
        og_i <- as.character(og[i, 1:800])
        ind_last <- which(og_i[1:800] == "")[1] - 1
    }
    if (i > 2000 & i <= 20000) {
        og_i <- as.character(og[i, 1:330]) 
        ind_last <- which(og_i[1:330] == "")[1] - 1
    }
    if (i > 20000 & i <= 60000) {
        og_i <- as.character(og[i, 1:10])
        ind_last <- which(og_i[1:10] == "")[1] - 1
    }
    if (i > 60000) {
        og_i <- as.character(og[i, 1:2])
        ind_last <- 1
    }
    
    if (is.na(ind_last)) {
        if (type == "OrthoFinder") {
            family_file <- rbind(family_file, 
                cbind(as.character(og_i[ 2:length(og_i) ]), 
                    rep(family, times = length(og_i)-1)))
        }
        if (type == "MCL") {
            family_file <- rbind(family_file, 
                cbind(as.character(og_i[ 1:length(og_i) ]), 
                    rep(family, times = length(og_i)))) 
        }
        
    } else { ## if !is.na(ind_last)
        if (type == "OrthoFinder") {
            family_file <- rbind(family_file, 
                cbind(as.character(og_i[ 2:ind_last ]), 
                    rep(family, times = length(2:ind_last)))) 
        }
        if (type == "MCL") {
            family_file <- rbind(family_file, 
                cbind(as.character(og_i[ 1:ind_last ]),
                    rep(family, times = length(1:ind_last)))) 
        }
    }
}

## for orthogroups.csv Orthofinder
if (type == "OrthoFinder") {
    family_file_tmp <- family_file[-grep(family_file[, 1], pattern = "OG0"), ]
    all(unassigned[, 1] %in% family_file_tmp[, 1])
    write.table(family_file_tmp, file = "family_orthogroups.txt", sep = "\t", 
                quote = FALSE, row.names = F, col.names = F)
}

## end for orthogroups.csv

## for orthogroups MCL

if (type == "MCL") {
    mcl <- family_file 
    which(duplicated(family_file[, 1]))
    all(unassigned[, 1] %in% mcl[, 1])
    ## add unassigned
    unassigned_add <- unassigned[which(!unassigned[, 1] %in% mcl[, 1]), ]
    group_end <- strsplit(as.character(mcl[dim(mcl)[1], 2]), split = "_")[[1]][2]
    group_end <- as.numeric(group_end)
    unassigned_add[, 2] <- paste("group", 
        c((group_end + 1):(group_end + dim(unassigned_add)[1])), sep="_")
    colnames(unassigned_add) <- c("V1", "V2")
    mcl <- rbind(mcl, unassigned_add)
    ## add sequences that are in OF but not in MCL
    family_file_tmp_add <- family_file_tmp[which(!family_file_tmp[, 1] %in% mcl[, 1]), ]
    group_end <- strsplit(as.character(mcl[dim(mcl)[1], 2]), split = "_")[[1]][2]
    group_end <- as.numeric(group_end)
    family_file_tmp_add[, 2] <- paste("group", 
        c((group_end + 1):(group_end + dim(family_file_tmp_add)[1])), sep = "_")
    colnames(family_file_tmp_add) <- c("V1", "V2")
    mcl <- rbind(mcl, family_file_tmp_add)
    write.table(mcl, file = "family_orthogroups_mcl.txt", sep="\t", 
        quote = FALSE, row.names = FALSE, col.names = FALSE)
}


## for i-ADHore
## read the tables
if (type == "OrthoFinder") {
    og <- read.table("family_orthogroups.txt")
    ## remove lines that contain retrotransposons (checked by BLAST)
    og <- og[!og[, 2] %in% c("OG0000001", "OG0000002"),]
    write.table(og, file = "family_orthogroups_iadhore.txt", sep = "\t", 
        quote = FALSE, row.names = FALSE, col.names = FALSE)
}

if (type == "MCL") {
    og_mcl <- read.table("family_orthogroups_mcl.txt")
    ## remove lines that contain retrotransposons (checked by BLAST)
    og_mcl <- og_mcl[!(og_mcl[, 2] %in% c("group_1", "group_3")),]
    write.table(og_mcl, file = "family_orthogroups_mcl_iadhore.txt", sep = "\t", 
        quote = FALSE, row.names = FALSE, col.names = FALSE)
}

##og <- read.table("family_orthogroups.txt") ## adjust for mcl or not
##og <- read.table("family_orthogroups_mcl.txt") ## adjust for mcl or not
## for family_orthogroups, remove OG0000001, OG0000002 since it contains retrotransposons
#og <- og[!og[,2] %in% c("OG0000001", "OG0000002"),]
## for family_orthogroups_mcl, remove group_1 and group_3 since it contains retrotransposons
#og <- og[!(og[,2] %in% c("group_1", "group_3")),]
##write.table(og, file="family_orthogroups_iadhore.txt", sep="\t", quote=FALSE, row.names=F, col.names=F)
##write.table(og, file="family_orthogroups_mcl_iadhore.txt", sep="\t", quote=FALSE, row.names=F, col.names=F)


## for MCScanX
## write orthogroups to homology file for each species condition
setwd("~/AG-Fernie/Thomas/Data/synteny/synteny_MCScanX")

if (type == "OrthoFinder") {
    og <- read.table("family_orthogroups.txt") 
    ## for family_orthogroups, remove OG0000001, OG0000002 since it 
    ## contains retrotransposons
    og <- og[!og[,2] %in% c("OG0000001", "OG0000002"),]
} 
if (type == "MCL") {
    og <- read.table("family_orthogroups_mcl.txt") 
    ## for family_orthogroups_mcl, remove group_1 and group_3 since it 
    ## contains retrotransposons
    og <- og[!(og[,2] %in% c("group_1", "group_3")),]
}

## get species combinations
species <- list.files()[grep(list.files(), pattern = "[.]gff")]
species <- unlist(lapply(strsplit(species, split = ".gff"), "[", 1))
species_combn <- cbind(combn(species, 2), rbind(species, species))
##species_combn <- matrix(c("papso", "zeama"), ncol=1)

og_1 <- as.character(og[, 1])
og_2 <- as.character(og[, 2])
for (i in 1:dim(species_combn)[2]) {
    
    species_combn_uniq <- unique(c(species_combn[1, i], species_combn[2, i]))
    
    species_combn_pattern <- if (length(species_combn_uniq) != 1) {
         paste(species_combn_uniq, collapse = "|") 
    } else {species_combn_uniq}
    
    og_ind <- grep(og_1, pattern = species_combn_pattern)
    og_gene <- as.character(og_1[og_ind])
    og_og <- as.character(og_2[og_ind])
    og_og_unique <- unique(og_og)

    df <- matrix(ncol = 2, nrow = 0)
    for (j in 1:length(og_og_unique)) {
        inds <- which(og_og == og_og_unique[j])
        if (length(inds) > 1) {
            og_gene_comb <- combn(
                og_gene[which(og_og == og_og_unique[j])], m = 2)
            og_gene_comb <- t(og_gene_comb)
            df <- rbind(df, og_gene_comb)
        }
    }
     
    ## adjust for OrthoFinder or MCL
    if (type == "OrthoFinder") {
        file_name <- paste0("MCScanX_gff/", species_combn[1, i], "_", 
            species_combn[2, i], ".homology")
    }
    if (type == "MCL") {
        file_name <- paste0("MCScanX_gff_mcl/", species_combn[1, i], "_", 
            species_combn[2, i], ".homology") 
    }
     
    write.table(df, file = file_name, quote = FALSE, row.names = FALSE, 
        col.names = FALSE, sep = "\t")
}