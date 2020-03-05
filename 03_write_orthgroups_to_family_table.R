setwd("~/AG-Fernie/Thomas/Data/synteny/Results_Nov05")

type <- "MCL" ## or MCL 

if (type == "OrthoFinder") { ## check Orthogroups.txt for lines
    og <- read.table("Orthogroups.txt", sep = " ", fill = TRUE, 
        stringsAsFactors = FALSE, header = FALSE, nrows = 63345)
	unassigned <- read.table("Orthogroups.txt", sep = " ", skip = 63345, 
		stringsAsFactors = FALSE, header = FALSE)
}

if (type == "MCL")  {
    ## start reading for mcl
    if (I == 1.5) {
        og <- read.table("out.seq.mci.I15", sep = "\t", fill = TRUE, 
                         stringsAsFactors = FALSE, header = FALSE, nrows = 61643)
        ## unassigned proteins
        unassigned <- read.table("out.seq.mci.I15", sep = "\t", dec = ".", 
            fill = TRUE, stringsAsFactors = FALSE, header = FALSE, skip = 61643)
    }
    if (I == 2.0) {
        og <- read.table("out.seq.mci.I20", sep = "\t", fill = TRUE, 
                         stringsAsFactors = FALSE, header = FALSE, nrows = 74059)
        ## unassigned proteins
        unassigned <- read.table("out.seq.mci.I20", sep = "\t", dec = ".", 
            fill = TRUE, stringsAsFactors = FALSE, header = FALSE, skip = 74059)
    }

}

## create result file
family_file <- matrix(ncol = 2, nrow = 0)

if (type == "OrthoFinder") 
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
        og_i <- as.character(og[i, 1:1100])
        ind_last <- which(og_i[1:1100] == "")[1] - 1
    }
    if (i > 2000 & i <= 20000) {
        og_i <- as.character(og[i, 1:390]) 
        ind_last <- which(og_i[1:390] == "")[1] - 1
    }
    if (i > 20000 & i <= 70000) {
        og_i <- as.character(og[i, 1:10])
        ind_last <- which(og_i[1:10] == "")[1] - 1
    }
    if (i > 70000) {
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

## for Orthofinder
if (type == "OrthoFinder") {
	
	unassigned[, 2] <- unlist(strsplit(unassigned[, 2], split = ":"))
    any(unassigned[, 1] %in% family_file[, 1])
	
    ## bind unassigned to family_file_tmp
	family_file_tmp <- rbind(family_file, unassigned)

    write.table(family_file_tmp, file = "family_orthogroups.txt", sep = "\t", 
                quote = FALSE, row.names = F, col.names = F)
}

## end for OrthoFinder

## for MCL
if (type == "MCL") {
    mcl <- family_file 
    which(duplicated(family_file[, 1]))
    all(unassigned[, 1] %in% mcl[, 1])
    ## add unassigned
    group_end <- strsplit(as.character(mcl[nrow(mcl), 2]), split = "_")[[1]][2]
    group_end <- as.numeric(group_end)
    unassigned_add <- cbind(unassigned, paste("group", 
        (group_end + 1):(group_end + nrow(unassigned)), sep="_"))
    colnames(unassigned_add) <- c("V1", "V2")
    mcl <- rbind(mcl, unassigned_add)
    
    ids <- read.table("../SequenceIDs.txt", sep = " ", stringsAsFactors = FALSE)
    ids[, 1] <- gsub(x = ids[, 1], pattern = ":", replacement = "")
    
    ## translate to proteins
    inds <- match(mcl[, 1], ids[, 1])
    mcl[, 1] <- ids[inds, 2]
    
    if (I == 1.5) {
        write.table(mcl, file = "family_orthogroups_mcl15.txt", sep="\t", 
                    quote = FALSE, row.names = FALSE, col.names = FALSE)    
    } 
    if (I == 2.0) {
        write.table(mcl, file = "family_orthogroups_mcl20.txt", sep="\t", 
                    quote = FALSE, row.names = FALSE, col.names = FALSE)    
    }
    
}


## for i-ADHore
## read the tables
if (type == "OrthoFinder") {
    og <- read.table("family_orthogroups.txt")
    ## remove lines that contain retrotransposons (checked by BLAST)
    og <- og[!og[, 2] %in% c("OG0000001", "OG0000002", "OG0000003"),] #######################checked 2019/11/10
    write.table(og, file = "family_orthogroups_iadhore.txt", sep = "\t", 
        quote = FALSE, row.names = FALSE, col.names = FALSE)
}

if (type == "MCL") {
    og_mcl <- read.table("family_orthogroups_mcl15.txt")
    ## remove lines that contain retrotransposons (checked by BLAST)
    og_mcl <- og_mcl[!(og_mcl[, 2] %in% c("group_1", "group_2")),] ##########################checked 2019/11/14
    write.table(og_mcl, file = "family_orthogroups_mcl_iadhore.txt", sep = "\t", 
        quote = FALSE, row.names = FALSE, col.names = FALSE)
}

## for MCScanX
## write orthogroups to homology file for each species condition
setwd("~/AG-Fernie/Thomas/Data/synteny/synteny_MCScanX")
type <- "MCL"

if (type == "OrthoFinder") {
    og <- read.table("family_orthogroups.txt") 
    ## for family_orthogroups, remove OG0000001, OG0000002 since it 
    ## contains retrotransposons
    og <- og[!og[, 2] %in% c("OG0000001", "OG0000002", "OG0000003"),] #######################checked 2019/11/10
} 
if (type == "MCL") {
    og <- read.table("family_orthogroups_mcl15.txt") 
    ## for family_orthogroups_mcl, remove group_1 and group_2 since it 
    ## contains retrotransposons
    og <- og[!(og[, 2] %in% c("group_1", "group_2")),] ######################checked 2019/11/14
}

## get species combinations
species <- list.files()[grep(list.files(), pattern = "[.]gff")]
species <- unlist(lapply(strsplit(species, split = ".gff"), "[", 1))
species_combn <- cbind(combn(species, 2), rbind(species, species))
##species_combn <- matrix(c("papso", "zeama"), ncol=1)

og_1 <- as.character(og[, 1])
og_2 <- as.character(og[, 2])

## 730:1000 ## 911
## 1001:2000 # 1083
## 2001:3000 # 2019
## 3001:4000 # 3080
## 4500:5000 # 4574
## 5001:6000 # 5072
## 6001:7000 # 6084
## 7001:8001 # 7029

for (i in 6562:7000) { ## 1:dim(species_combn)[2]
    
    species_combn_uniq <- unique(c(species_combn[1, i], species_combn[2, i]))
    species_combn_pattern <- paste(species_combn_uniq, collapse = "|") 
    
    
    og_ind <- grep(og_1, pattern = species_combn_pattern)
    og_gene <- og_1[og_ind]
    og_og <- og_2[og_ind]
    og_og_unique <- unique(og_og)

    df <- matrix(ncol = 2, nrow = 0)
    
    res <- lapply(og_og_unique, function(x) {
        inds <- which(og_og == x)
        if (length(inds) > 1) {
            og_gene_comb <- combn(og_gene[inds], m = 2)
            og_gene_comb <- t(og_gene_comb)
            return(og_gene_comb)
        } else {
            return(NULL)
        }
    })
    df <- do.call("rbind", res)
    
    ## adjust for OrthoFinder or MCL
    if (type == "OrthoFinder") {
        file_name <- paste0("MCScanX_gff/", species_combn[1, i], "_", 
            species_combn[2, i], ".homology")
    }
    if (type == "MCL") {
        file_name <- paste0("MCScanX_gff_mcl/00_qsub/00_check/", species_combn[1, i], "_", 
            species_combn[2, i], ".homology") 
    }
     
    write.table(df, file = file_name, quote = FALSE, row.names = FALSE, 
        col.names = FALSE, sep = "\t")
}
