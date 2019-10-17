## setwd
setwd("W:/Thomas/Data/synteny")

library(igraph)

## get pks genes from output of Orthofinder or MCL
genes_table <- read.table("./Results_Oct26/family_orthogroups.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)
pks_genes <- sort(genes_table[genes_table[, 2] == "OG0000260", 1]) ## for of 
genes_table_mcl <- read.table("./Results_Oct26/family_orthogroups_mcl.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
pks_genes_mcl <- sort(genes_table_mcl[genes_table_mcl[, 2] == "group_337", 1]) ## for mcl

## vector with all pks genes identified by orthofinder and mcl, use this to create the bin_mat_... matrices
pks_genes_all <- sort(unique(c(pks_genes, pks_genes_mcl)) )

## pinta_PITA_0064 and pinta_PITA_23328  are identical
pks_genes[which(pks_genes == "pinta_PITA_00664")] <- "pinta_PITA_00664/pinta_PITA_23328"
pks_genes[which(pks_genes == "pinta_PITA_23328")] <- "pinta_PITA_00664/pinta_PITA_23328"
pks_genes_mcl[which(pks_genes_mcl == "pinta_PITA_00664")] <- "pinta_PITA_00664/pinta_PITA_23328"
pks_genes_mcl[which(pks_genes_mcl == "pinta_PITA_23328")] <- "pinta_PITA_00664/pinta_PITA_23328"
pks_genes_all[which(pks_genes_all == "pinta_PITA_00664")] <- "pinta_PITA_00664/pinta_PITA_23328"
pks_genes_all[which(pks_genes_all == "pinta_PITA_23328")] <- "pinta_PITA_00664/pinta_PITA_23328"
pks_genes[which(pks_genes == "pinta_PITA_00679")] <- "pinta_PITA_00679/pinta_PITA_06018"
pks_genes[which(pks_genes == "pinta_PITA_06018")] <- "pinta_PITA_00679/pinta_PITA_06018"
pks_genes_mcl[which(pks_genes_mcl == "pinta_PITA_00679")] <- "pinta_PITA_00679/pinta_PITA_06018"
pks_genes_mcl[which(pks_genes_mcl == "pinta_PITA_06018")] <- "pinta_PITA_00679/pinta_PITA_06018"
pks_genes_all[which(pks_genes_all == "pinta_PITA_00679")] <- "pinta_PITA_00679/pinta_PITA_06018"
pks_genes_all[which(pks_genes_all == "pinta_PITA_06018")] <- "pinta_PITA_00679/pinta_PITA_06018"

## make the genes unique
pks_genes <- unique(pks_genes)
pks_genes_mcl <- unique(pks_genes_mcl)
pks_genes_all <- unique(pks_genes_all)

## load information on pks genes
supp <- xlsx::read.xlsx("~/Documents/03_Syntenic_linkage/01_Data/synteny_network_results/pks_genes_tree_id_type.xlsx", sheetIndex=1, stringsAsFactors=FALSE)
supp <- supp[,1:24]

## save the pks_genes
setwd("~/winhome/Documents/03_Syntenic_linkage/01_Data/synteny_network_results/")
save(pks_genes, pks_genes_mcl, pks_genes_all, file="./pks_genes.RData")

## i-ADHore
## function to find synteny by finding pks_genes in i-ADHore output
find_synteny <- function(files, bin_mat, pks_genes) {
    ## loop through files
    for (i in seq_along(files)) {
        files_i <- paste0(files[i], "/anchorpoints.txt")
        mult_pairs <- tryCatch(read.csv(files_i, header = FALSE, sep = "\t", 
            skip = 1, fill = TRUE, comment.char = "#", 
            stringsAsFactors = FALSE), error = function(e) NULL)
        if (!is.null(mult_pairs)) {
            gene_1_ind <- which(mult_pairs[, 4] %in% pks_genes)
            gene_2_ind <- which(mult_pairs[, 5] %in% pks_genes)
            gene_ind <- unique(c(gene_1_ind, gene_2_ind))
            
            links <- mult_pairs[gene_ind, 4:5]

            ## write 1 to collinear pair when there is synteny reported
            for (j in seq_len(nrow(links))) {
                bin_mat[links[j, 1], links[j, 2]] <- bin_mat[links[j, 2], links[j, 1]] <- 1
            }
            
        } else {print(files_i)}
    }
    return(bin_mat)
}

## function to find tandem by finding pks_genes in i-ADHoRe output
find_tandem <- function(files, bin_mat, pks_genes) {
    
    tandem_list <- vector("list", length(files))
    ## loop through files
    for (i in seq_along(files)) {
        tandem <- read.csv(paste(files[i], "genes.txt", sep = "/"), 
            header = TRUE, sep = "\t", fill = TRUE, comment.char = "#", 
            stringsAsFactors = FALSE)
        tandem <- tandem[tandem[, "id"] %in% pks_genes, c("id", "tandem_representative")]
        tandem_ind <- which(tandem[, "tandem_representative"] != "")
        tandem <- tandem[tandem_ind,]
        
        tandem_list[[i]] <- tandem
        ## write 1 to tandem when there is a tandem gene reported
        for (j in seq_len(nrow(tandem))) {
            bin_mat[tandem[j, 1], tandem[j, 2]] <- bin_mat[tandem[j, 2], tandem[j, 1]] <- 1
        }
    }
    return(list(bin_mat, tandem_list))
}

## load i-ADHore results from Orthofinder input
setwd("~/AG-Fernie/Thomas/Data/synteny/synteny_iadhore/output/")
output_files <- list.files()[grep(list.files(), pattern = "output_")]

## create binary matrix that will store if there exists a collinear or tandem
## relationship or not
bin_mat_tandem <- bin_mat <- matrix(data = 0, nrow = length(pks_genes_all), 
    ncol = length(pks_genes_all))
colnames(bin_mat) <- rownames(bin_mat) <- pks_genes_all
colnames(bin_mat_tandem) <- rownames(bin_mat_tandem) <- pks_genes_all
bin_mat <- find_synteny(output_files, bin_mat, pks_genes) ## use here pks_genes 
tandem <- find_tandem(output_files, bin_mat_tandem, pks_genes)

## load i-ADHore results from MCL input
setwd("../output_mcl")
output_files <- list.files()[grep(list.files(), pattern = "output_")]

## create binary matrix that will store if there exists a collinear or tandem
## relationship or not
bin_mat_tandem_mcl <- bin_mat_mcl <- matrix(data = 0, nrow = length(pks_genes_all), 
    ncol = length(pks_genes_all))
colnames(bin_mat_mcl) <- rownames(bin_mat_mcl) <- pks_genes_all
colnames(bin_mat_tandem_mcl) <- rownames(bin_mat_tandem_mcl) <- pks_genes_all
bin_mat_mcl <- find_synteny(output_files, bin_mat_mcl, pks_genes_mcl) ## use here pks_genes_mcl
tandem_mcl <- find_tandem(output_files, bin_mat_tandem_mcl, pks_genes_mcl)

## MCScanX
## function to find synteny by finding pks_genes in MSScanX output
find_synteny_mcscanx <- function(files, bin_mat, pks_genes) {
    ## loop through files
    for (i in 1:length(files)) {
        collinear_i <- NULL
        collinear_i <- tryCatch(read.csv(files[i], sep = "\t", fill = TRUE, 
            stringsAsFactors = FALSE, header = FALSE, comment.char = "#"), 
            error = function(x) NULL)
        if (!is.null(collinear_i)) {
            ## get index
            gene_1_ind <- which(collinear_i[, 2] %in% pks_genes)
            gene_2_ind <- which(collinear_i[, 3] %in% pks_genes)
            gene_ind <- unique(c(gene_1_ind, gene_2_ind))
            
            links <- collinear_i[gene_ind, 2:3]
            
            ## write 1 to collinear pair when there is synteny reported
            for (j in seq_len(nrow(links))) {
                bin_mat[links[j,1], links[j, 2]] <- bin_mat[links[j, 2], links[j, 1]] <- 1
            }    
        } else (print(files[i]))
    }
    return(bin_mat)
}

## function to find tandems by finding pks_genes in MCScanX output
find_tandem_mcscanx <- function(files, bin_mat, pks_genes) {
    ## loop through files
    tandem_list <- list()
    for (i in 1:length(files)) {
        file_i <- paste0(strsplit(files[i], split = "[.]collinearity")[[1]], ".tandem")
        tandem_i <- tryCatch(read.csv(file_i, sep = ",", 
            stringsAsFactors = FALSE, header = FALSE, comment.char = "#"), 
            error = function(x) NULL)
        if (!is.null(tandem_i)) {
            gene_1_ind <- which(tandem_i[,1] %in% pks_genes)
            gene_2_ind <- which(tandem_i[,2] %in% pks_genes)
            gene_ind <- unique(c(gene_1_ind, gene_2_ind)) 
            tandem <- tandem_i[gene_ind, ]
            
            tandem_list[[i]] <- tandem    
            ## write 1 to tandem pair when there is synteny reported
            for (j in seq_len(nrow(tandem))) {
                bin_mat[tandem[j, 1], tandem[j, 2]] <- bin_mat[tandem[j, 2], tandem[j, 1]] <- 1
            }
        } else (print(file_i))
    }
    return(list(bin_mat, tandem_list))
} ## replaces .collineary by .tandem

## load MCScanX results from Orthofinder input
setwd("~/AG-Fernie/Thomas/Data/synteny/synteny_MCScanX/MCScanX_gff")
collinearity_files <- list.files()[grep(list.files(), pattern = "[.]collinearity")]

## create binary matrix that will store if there exists a collinear or tandem
## relationship or not
bin_mat_tandem_mcscanx <- bin_mat_mcscanx <- matrix(0, 
    nrow = length(pks_genes_all), ncol = length(pks_genes_all))
colnames(bin_mat_mcscanx) <- rownames(bin_mat_mcscanx) <- pks_genes_all
colnames(bin_mat_tandem_mcscanx) <- rownames(bin_mat_tandem_mcscanx) <- pks_genes_all
bin_mat_mcscanx <- find_synteny_mcscanx(collinearity_files, bin_mat_mcscanx, pks_genes)
tandem_mcscanx <- find_tandem_mcscanx(collinearity_files, bin_mat_tandem_mcscanx, pks_genes)


## load MScanX results from MCL input
setwd("~/AG-Fernie/Thomas/Data/synteny/synteny_MCScanX/MCScanX_gff_mcl")
collinearity_files <- list.files()[grep(list.files(), pattern = "[.]collinearity")]

## create binary matrix that will store if there exists a collinear or tandem
## relationship or not
bin_mat_tandem_mcscanx_mcl <- bin_mat_mcscanx_mcl <- matrix(0, 
    nrow = length(pks_genes_all), ncol = length(pks_genes_all))
colnames(bin_mat_mcscanx_mcl) <- rownames(bin_mat_mcscanx_mcl) <- pks_genes_all
colnames(bin_mat_tandem_mcscanx_mcl) <- rownames(bin_mat_tandem_mcscanx_mcl) <- pks_genes_all
bin_mat_mcscanx_mcl  <- find_synteny_mcscanx(collinearity_files, bin_mat_mcscanx_mcl, pks_genes_mcl)
tandem_mcscanx_mcl <- find_tandem_mcscanx(collinearity_files, bin_mat_tandem_mcscanx_mcl, pks_genes_mcl)

## save and load tandem and bin_mat files
setwd("~/winhome/Documents/03_Syntenic_linkage/01_Data/synteny_network_results/")
save(bin_mat, bin_mat_mcl, bin_mat_mcscanx, bin_mat_mcscanx_mcl, file="./synteny_bin_mat.RData")
save(tandem, tandem_mcl, tandem_mcscanx, tandem_mcscanx_mcl, file="./synteny_tandem_mat.RData")
load("./synteny_bin_mat.RData")
load("./synteny_tandem_mat.RData")

##### TANDEM #####
## how to proceed with tandems? get tandems for all species from the four methods and use the lowest level union
## of tandem affiliation to concatenate/paste protein names, i.e. for
## method 1: x1 and x2 and x4 are a tandem
## method 2: x1 and x2 and x3 are a tandem 
## then x1 and x2 and x3 and x4 are a tandem
all(rownames(tandem[[1]]) == rownames(tandem_mcl[[1]]))
all(rownames(tandem[[1]]) == rownames(tandem_mcscanx[[1]]))
all(rownames(tandem[[1]]) == rownames(tandem_mcscanx_mcl[[1]]))
tandem_sum <- tandem[[1]] + tandem_mcl[[1]] + tandem_mcscanx[[1]] + tandem_mcscanx_mcl[[1]]
rownames(tandem_sum) <- rownames(tandem[[1]])
ind_keep <- apply(tandem_sum, 1, sum) > 0

tandem_sum <- tandem_sum[ind_keep, ind_keep]
##tandem_sum[which(tandem_sum>0)] <- 1
diag(tandem_sum) <- 0

## create a matrix that stores the source of tandem, i.e. if the information
## comes from i-ADHore + Orthofinder, i-ADHore + MCL, MCScanX + Orthofinder 
## and/or MCScanX + MCL
tandem_type <- matrix("", ncol=ncol(tandem[[1]]), nrow = nrow(tandem[[1]]))
rownames(tandem_type) <- rownames(tandem[[1]])
tandem_type[tandem[[1]] == 1] <- "iadhore/"
tandem_type[tandem_mcl[[1]] == 1] <- paste(tandem_type[tandem_mcl[[1]] == 1], "iadhore_mcl/", sep="")
tandem_type[tandem_mcscanx[[1]] == 1] <- paste(tandem_type[tandem_mcscanx[[1]] == 1], "mcscanx/", sep="")
tandem_type[tandem_mcscanx_mcl[[1]] == 1] <- paste(tandem_type[tandem_mcscanx_mcl[[1]] == 1], "mcscanx_mcl/", sep="")

## remove the final slash ("/")
tandem_type_vec <- lapply(as.vector(tandem_type), function(x) {
    paste(strsplit(x, split = "/")[[1]], collapse = "/")})
tandem_type_vec <- unlist(tandem_type_vec)

## write back to a matrix
tandem_type <- matrix(tandem_type_vec, ncol = 
    ncol(tandem_type), nrow = nrow(tandem_type))
rownames(tandem_type) <- colnames(tandem_type) <- rownames(tandem[[1]])

## plot tandem network
g <- igraph::graph_from_adjacency_matrix(tandem_sum, diag = FALSE)
plot(g, vertex.label.cex = 0.1, vertex.size = 0.1, edge.width = 1, 
     edge.arrow.size = 0.1)

## get components of complete graphs, i.e. members of each graphs, in this case
## all tandem genes per region
components_g <- igraph::components(g)$membership
components_g_unique <- unique(components_g)

## check reliability of connection by checking number of neighbouring genes
res <- vector("list", length(pks_genes_mcl))
setwd("~/AG-Fernie/Thomas/Data/synteny/synteny_iadhore")
gff_files <- list.files()
gff_files <- gff_files[grep(gff_files, pattern = "gff")]

## iterate through pks_genes_all
for(i in 1:length(pks_genes_all) ) {
    species <- unlist(lapply(strsplit(pks_genes_all[i], split = "_"), "[", 1))
    ## get chromosome of gene
    gene <- unlist(lapply(strsplit(
        pks_genes_all[i], split = paste(species, "_", sep = "")), "[", 2))
    gff <- read.table(gff_files[grep(gff_files, pattern = species)], 
        sep = "\t", stringsAsFactors = FALSE, quote = "")
    chr <- unique(gff[grep(gff[, 9], pattern = gene), 1])
    
    ## for species covsu and selmo strsplit the gene, since the genes were 
    ## pasted from different names (not directly from the gff files), get the 
    ## identifier that can be found in the gff
    if (species=="covsu") chr <- unique(gff[grep(gff[, 9], 
        pattern = unlist(lapply(strsplit(gene, split = "_"), "[", 2))), 1])
    if (species=="selmo") chr <- unique(gff[grep(gff[, 9], 
        pattern = unlist(lapply(strsplit(gene, split = "_"), "[", 2))), 1])
    ## load chromosome file lst and get position of gene in file
    res[[i]] <- list()
    if (!length(chr) == 0) {
        lst <- read.table(paste("lst_files/", species, "/", chr, ".lst", sep = ""), 
            stringsAsFactors = FALSE, quote="")
        lst <- unique(lst[, 1])
        genes_lst <- substring(lst, 1, nchar(lst)-1)
        position <- which(genes_lst == paste(species, gene, sep = "_"))
        end <- length(lst)
        ## write the following information to the list res: start end position, name of chromosome
        res[[i]][[1]] <- c(position-1, end-position)
        res[[i]][[2]] <- chr
    } else {
        res[[i]][[1]] <- NULL
        res[[i]][[2]] <- chr
    }
    print(c(pks_genes_all[i], res[[i]][[1]], res[[i]][[2]]))
}

names(res) <- pks_genes_all

## save res to neighbours_on_chromosomes.RData
setwd("~/winhome/Documents/03_Syntenic_linkage/01_Data/synteny_network_results")
save(res, file = "./neighbours_on_chromosomes.RData")
load("neighbours_on_chromosomes.RData")
## end check reliability

## some plots
## calculate length of chromosome/scaffold where PKS gene is located on
length_chr <- unlist(lapply(res, function(x) sum(x[[1]]) + 1)) 
hist(log(length_chr))
## calculate minimum to end to chromosome/scaffold
min_chr <- unlist(lapply(res, function(x) min(x[[1]]))) 
hist(log(min_chr))

## check tandem components if they are reported correctly --> tandems should 
## be on the same chromosome and within distance of 20 genes
for (i in 1:length(components_g_unique)) {
    name_i <- names(which(components_g == i))
    chrs <- unlist(lapply(res[name_i], "[", 2))

    if (length(unique(chrs)) != 1) {print(i)} else { 
        ## print i when genes are on different chr
        position <- lapply(res[name_i], "[", 1)
        
        ## get positions and calculate differences 
        position <- unlist(position)[c(TRUE, FALSE)] 
        position <- sort(position)
        position_diff <- numeric(length(position) -1)
        for (j in 1:length(position_diff)) position_diff[j] <- position[j + 1] - position[j]
        if (any(position_diff > 20)) print(i)
    }
}


## name comb will be the vector that stores the pasted names of all tandem genes
name_comb <- vector("numeric", igraph::components(g)$no)
for (i in seq_along(components_g_unique)) {
    names_i_ind <- components_g == components_g_unique[i]
    names_i <- names(components_g)[names_i_ind]
    name_comb[names_i_ind] <- paste(names_i, collapse = "/")
}

##### Syntenic links #####
## add all bin_mat* matrices

## some checks for bin_mat* matrices, they should all have the same rownames
all(rownames(bin_mat) == rownames(bin_mat_mcl))
all(rownames(bin_mat) == rownames(bin_mat_mcscanx))
all(rownames(bin_mat) == rownames(bin_mat_mcscanx_mcl))
bin_mat_complete <- 0.25 * bin_mat + 0.25 * bin_mat_mcl + 
    0.25 * bin_mat_mcscanx + 0.25 * bin_mat_mcscanx_mcl

## get source and write the source to mat_type
mat_type <- matrix("", ncol = ncol(bin_mat), nrow = nrow(bin_mat))
mat_type[bin_mat == 1] <- "iadhore/"
mat_type[which(bin_mat_mcl == 1)] <- paste(
    mat_type[which(bin_mat_mcl == 1)], "iadhore_mcl/", sep = "")
mat_type[which(bin_mat_mcscanx == 1)] <- paste(
    mat_type[which(bin_mat_mcscanx == 1)], "mcscanx/", sep = "")
mat_type[which(bin_mat_mcscanx_mcl == 1)] <- paste(
    mat_type[which(bin_mat_mcscanx_mcl == 1)], "mcscanx_mcl/", sep = "")

## remove the final slash ("/")
mat_type_vec <- lapply(as.vector(mat_type), 
    function(x) paste(strsplit(x, split = "/")[[1]], collapse="/"))
mat_type_vec <- unlist(mat_type_vec)

## write back a to matrix
mat_type <- matrix(mat_type_vec, ncol = ncol(mat_type), nrow = nrow(mat_type))

## assign rownames of bin_mat to col- and rownames of bin_mat_complete and 
## mat_type
colnames(bin_mat_complete) <- rownames(bin_mat_complete) <- rownames(bin_mat)
rownames(mat_type) <- colnames(mat_type) <- rownames(bin_mat)

## plot bin_mat_complete network
g <- graph_from_adjacency_matrix(bin_mat_complete, diag=FALSE)
plot(g, vertex.label.cex=0.1, vertex.size=0.1, edge.width=1, edge.arrow.size=0.1)


## rename tandem genes to pasted names
name <- names(components_g) ## name of tandem genes
name_comb_unique <- unique(name_comb)

## iterate through the unique name_comb
for (i in 1:length(name_comb_unique)) {
    
    name_i_ind <- name_comb_unique[[i]] == name_comb
    name_i <- name[name_i_ind]
    
    ## calculate the sum from all tandem genes 
    connection_sum <- apply(mat_type[name_i, ], 2, 
        function(x) paste(x, collapse = "/"))
    connection_sum <- paste(connection_sum, apply(mat_type[, name_i], 1, 
        function(x) paste(x, collapse = "/")), sep = "/")
    connection_sum <- lapply(strsplit(unlist(connection_sum), split = "/"), unique)                    
    connection_sum <- lapply(connection_sum, function(x) x[x != ""])
    connection_sum <- unlist(lapply(connection_sum, length))
    
    if (any(connection_sum>4)) stop(i)
    if (sum(connection_sum > 0) > 0) {
        ## combine all sources from all name_i and assign to mat_type
        if (sum(connection_sum>0) > 1) {
            mat_type_comb <- apply(mat_type[name_i, connection_sum > 0], 2, 
                function(x) paste(x, collapse = "/"))
            mat_type_comb <- paste(mat_type_comb, 
                apply(mat_type[connection_sum > 0, name_i], 1, 
                function(x) paste(x, collapse = "/")), sep = "/")
        } else {
            mat_type_comb <- paste(mat_type[name_i, connection_sum > 0], 
                collapse = "/")
            mat_type_comb <- paste(mat_type_comb, paste(
                mat_type[connection_sum > 0, name_i], collapse = "/"), sep = "/")
        }
        ## split by /
        mat_type_comb_l <- strsplit(mat_type_comb, split = "/") 
        ## make unique and sort
        mat_type_comb_l <- lapply(mat_type_comb_l, function(x) sort(unique(x))) 
        ## remove ""
        mat_type_comb_l <- lapply(mat_type_comb_l, function(x) x[x != ""])
        ## paste again
        mat_type_comb_l <- lapply(mat_type_comb_l, function(x) 
            paste(x, collapse = "/"))
        mat_type_comb <- as.vector(unlist(mat_type_comb_l))
        mat_type[name_i, connection_sum > 0] <- matrix(
            rep(mat_type_comb, times = length(name_i)), 
            ncol = sum(connection_sum > 0), byrow = T)
        mat_type[connection_sum > 0, name_i] <- matrix(
            rep(mat_type_comb, times = length(name_i)), 
            nrow = sum(connection_sum > 0), byrow = F)
        
        ## assign connection_sum to the first element in name i
        bin_mat_complete[name_i, ] <- matrix(
            rep(connection_sum / 4, times = length(name_i)), 
            ncol = length(connection_sum), byrow = T)
        bin_mat_complete[, name_i] <- matrix(
            rep(connection_sum / 4, times = length(name_i)), 
            nrow = length(connection_sum), byrow = F)
    }
    
    ## rename the first element to the combined name
    rownames(bin_mat_complete)[rownames(bin_mat_complete) %in% name_i] <- name_comb_unique[[i]] 
    colnames(bin_mat_complete)[colnames(bin_mat_complete) %in% name_i] <- name_comb_unique[[i]]
    rownames(mat_type)[rownames(mat_type) %in% name_i] <- name_comb_unique[[i]]
    colnames(mat_type)[colnames(mat_type) %in% name_i] <- name_comb_unique[[i]]
}

## remove all other elements from bin_mat_complete
ind_remove <- duplicated(rownames(bin_mat_complete)) 
## remove the duplicated rownames
bin_mat_complete <- bin_mat_complete[!ind_remove, !ind_remove]
mat_type <- mat_type[!ind_remove, !ind_remove]

## remove type of connection that are "iadhore_mcl/mcscanx", 
## "iadhore/mxscanx_mcl" (not compatible techniques and clustering)
table(mat_type)
bin_mat_complete[which(
    mat_type %in% c( "iadhore_mcl/mcscanx", "iadhore/mcscanx_mcl"))] <- 0
mat_type[which(
    mat_type %in% c("iadhore_mcl/mcscanx", "iadhore/mcscanx_mcl"))] <- ""

## for all nodes supported by only one technique 
## ("iadhore", "iadhore_mcl", "mcscanx", "mcscanx_mcl"): check in the 
## respective files and keep the connection if there are > n connections
## (do not remove connections in the end)
find_length_iadhore <- function(type = "iadhore_mcl", path = "output_mcl") {
    inds <- which(mat_type == type, arr.ind = TRUE)
    res <- vector("numeric", nrow(inds))
    
    ## iterate through the rows in inds
    for(i in 1:nrow(inds)) {
        inds_i <- inds[i, ]
        species <- lapply(
            strsplit(colnames(mat_type)[inds_i], split = "_"), "[", 1)
        species <- sort(unlist(species))
        species <- paste("output", species[1], species[2], sep = "_")
        
        ## paste the filename with path to adjust for the specific methods
        filepath <- paste("~/AG-Fernie/Thomas/Data/synteny/synteny_iadhore/", 
            path, sep = "")
        filepath <- paste(filepath, species, "anchorpoints.txt", sep = "/")
        file_iadhore <- read.table(filepath, sep = "\t", header = T)
        genes <- sort(colnames(bin_mat_complete)[inds_i])
        genes <- strsplit(genes, split = "/")
        ind_gene1 <- c(which(file_iadhore[, 4] %in% genes[[1]]), 
            which(file_iadhore[, 5] %in% genes[[1]]))
        ind_gene2 <- c(which(file_iadhore[, 4] %in% genes[[2]]), 
            which(file_iadhore[, 5] %in% genes[[2]]))
        
        ind_genes <- intersect(ind_gene1, ind_gene2)
        mp_i <- file_iadhore[ind_genes, "multiplicon"]
        
        len <- length(which(file_iadhore[, "multiplicon"] %in% mp_i))
        res[i] <- len
    }
    return(list(res, inds))
}

## apply the function
res_iadhore <- find_length_iadhore(type = "iadhore", path = "output")
res_iadhore_mcl <- find_length_iadhore(type = "iadhore_mcl", path = "output_mcl")

## write similarly a function for MSCanX files
find_length_mcscanx <- function(type = "mcscanx_mcl", path = "MCScanX_gff_mcl") {
    
    inds <- which(mat_type == type, arr.ind = TRUE)
    res <- vector("numeric", nrow(inds))
    
    ## iterate through the rows in inds
    for(i in 1:nrow(inds)) {
        inds_i <- inds[i, ]
        species <- lapply(strsplit(colnames(mat_type)[inds_i], split="_"), "[", 1)
        species <- sort(unlist(species))
        species <- paste(species[1], species[2], sep = "_")
        
        ## paste the filename with path to adjust for the specific methods
        filepath <- paste("~/AG-Fernie/Thomas/Data/synteny/synteny_MCScanX/", 
            path, "/", sep = "")
        filepath <- paste(filepath, species, ".collinearity", sep = "")
        file_mcscan <- read.table(filepath, sep = "\t", header = FALSE, 
            stringsAsFactors = FALSE)
        genes <- sort(colnames(bin_mat_complete)[inds_i])
        genes <- strsplit(genes, split = "/")
        ind_gene1 <- c(which(file_mcscan[, 2] %in% genes[[1]]), 
            which(file_mcscan[, 3] %in% genes[[1]]))
        ind_gene2 <- c(which(file_mcscan[, 2] %in% genes[[2]]), 
            which(file_mcscan[, 3] %in% genes[[2]]))
        
        ind_genes <- intersect(ind_gene1, ind_gene2)
        mp <- unlist(lapply(strsplit(
            as.character(file_mcscan[, 1]), split = "-"), "[", 1))
        mp <- gsub(mp, pattern = " ", replacement = "")
        
        mp_i <- unique(mp[ind_genes])
        
        len <- length(which(mp %in% mp_i))
        res[i] <- len
        if (len == 0) print(i)
    }
    return(list(res, inds))
}

## apply the function
res_mcscanx <- find_length_mcscanx(type = "mcscanx", path = "MCScanX_gff")
res_mcscanx_mcl <- find_length_mcscanx(type = "mcscanx_mcl", path = "MCScanX_gff_mcl")

## save
setwd("~/winhome/Documents/03_Syntenic_linkage/01_Data/synteny_network_results/")
save(res_iadhore, res_iadhore_mcl, res_mcscanx, res_mcscanx_mcl, 
    file = "filter_res_iadhore_mcscan.RData")
load("filter_res_iadhore_mcscan.RData")

## remove for iadhore
table(res_iadhore[[1]])
# for (i in 1:length(res_iadhore[[1]]))  {
#     if (res_iadhore[[1]][i] <= 5) {
#         inds_i <- res_iadhore[[2]][i, ]
#         bin_mat_complete[inds_i[1], inds_i[2]] <- 0
#         bin_mat_complete[inds_i[2], inds_i[1]] <- 0
#         mat_type[inds_i[1], inds_i[2]] <- ""
#         mat_type[inds_i[2], inds_i[1]] <- ""
#     }
# }
## remove for iadhore_mcl
table(res_iadhore_mcl[[1]])
# for (i in 1:length(res_iadhore_mcl[[1]]))  {
#     if (res_iadhore_mcl[[1]][i] <= 5) {
#         inds_i <- res_iadhore_mcl[[2]][i, ]
#         bin_mat_complete[inds_i[1], inds_i[2]] <- 0
#         bin_mat_complete[inds_i[2], inds_i[1]] <- 0
#         mat_type[inds_i[1], inds_i[2]] <- ""
#         mat_type[inds_i[2], inds_i[1]] <- ""
#     }
# }

## remove for mcscanx
table(res_mcscanx[[1]])
# for (i in 1:length(res_mcscanx[[1]]))  {
#     if (res_mcscanx[[1]][i] <= 7) {
#         inds_i <- res_mcscanx[[2]][i, ]
#         bin_mat_complete[inds_i[1], inds_i[2]] <- 0
#         bin_mat_complete[inds_i[2], inds_i[1]] <- 0
#         mat_type[inds_i[1], inds_i[2]] <- ""
#         mat_type[inds_i[2], inds_i[1]] <- ""
#     }
# }
## remove for mcscanx_mcl
table(res_mcscanx_mcl[[1]])
# for (i in 1:length(res_mcscanx_mcl[[1]]))  {
#     if (res_mcscanx_mcl[[1]][i] <= 7) {
#         inds_i <- res_mcscanx_mcl[[2]][i, ]
#         bin_mat_complete[inds_i[1], inds_i[2]] <- 0
#         bin_mat_complete[inds_i[2], inds_i[1]] <- 0
#         mat_type[inds_i[1], inds_i[2]] <- ""
#         mat_type[inds_i[2], inds_i[1]] <- ""
#     }
# }
## end of checking the length for a link that is only supported by one method

## plot the bin_mat_complete network
g <- graph_from_adjacency_matrix(bin_mat_complete, diag = FALSE)
plot(g, vertex.label.cex = 0.1, vertex.size = 0.1, edge.width = 1, 
     edge.arrow.size = 0.1)


## remove the proteins that do not link to others
inds_keep <- apply(bin_mat_complete, 1, sum) > 0
## keep also duplicated genes
inds_keep[rownames(bin_mat_complete) %in% name_comb_unique ] <- TRUE

## check how many genes are on the chromosome/scaffold for (not) linking ones 
## not linking ones
res_nlink <- res[rownames(bin_mat_complete)[!inds_keep]]
res_nlink <- lapply(1:length(res_nlink), function(x) res_nlink[[x]][[1]])
res_nlink_sum <- lapply(res_nlink, sum)
names(res_nlink_sum) <- rownames(bin_mat_complete)[!inds_keep]

## add 1 since the pks_genes is not counted
res_nlink_sum <- unlist(res_nlink_sum) + 1 

## create a vector that stores the colour depending on how many genes there 
## are on the chromosome/scaffold
reliability <- vector("character", dim(bin_mat_complete)[1])
names(reliability) <- rownames(bin_mat_complete)
tmp <- unlist(res_nlink_sum)
reliability[names(res_nlink_sum[tmp < 15])] <- "red"
reliability[names(res_nlink_sum[tmp >= 15 & tmp < 25])] <- "yellow"
reliability[names(res_nlink_sum[tmp >= 25])] <- "green"

## plot network that contains features that do not link to other
net <- graph_from_adjacency_matrix(
    bin_mat_complete[!inds_keep, !inds_keep], mode = "directed", diag = FALSE)
plot(net, vertex.label.cex = 0.3, vertex.size = 5, 
    vertex.color = reliability[!inds_keep],edge.arrow.size = 0.1)

## build network for linking features
res_link <- res[grep(names(res), 
    pattern = gsub(pattern = "/", replacement = "|", 
    paste(rownames(bin_mat_complete)[inds_keep], collapse = "/")))]
res_link <- lapply(1:length(res_link), function(x) res_link[[x]][[1]])
res_link_sum <- lapply(res_link, sum)
names(res_link_sum) <- names(res[grep(names(res), 
    pattern = gsub(pattern = "/", replacement = "|", 
    paste(rownames(bin_mat_complete)[inds_keep], collapse = "/")))])

## rename res_link_sum and truncate with names used in bin_mat_complete
for (i in 1:length(components_g_unique)) {
    new_name <- paste(names(components_g)[components_g == components_g_unique[i]], collapse="/")
    inds <- grep(names(res_link_sum), pattern = gsub("/", "|", new_name))
    names(res_link_sum)[inds[1]] <- new_name
    res_link_sum <- res_link_sum[ -inds[ 2:length(inds) ]]
}
## add 1 since the pks_genes is not counted
res_link_sum <- unlist(res_link_sum) + 1 


write.table(
    data.frame(names = names(unlist(res_link_sum)), value = unlist(res_link_sum)), 
    file = "./figure_S7_synteny_network_quality/sum_genes_scaffold_linking.txt", 
    quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(
    data.frame(names = names(unlist(res_nlink_sum)), value = unlist(res_nlink_sum)), 
    file = "./figure_S7_synteny_network_quality/sum_genes_scaffold_notlinking.txt", 
    quote = FALSE, row.names = FALSE, col.names = TRUE)

## plot for figure synteny_network_quality
df <- unlist(lapply(strsplit(names(res_nlink_sum), split = "_"), "[", 1))
df <- sort(table(df))
df <- data.frame(species = names(df), number = as.vector(df))
df$species <- factor(df$species, levels=df$species[order(df$number)])
g <- ggplot(df) + geom_bar(aes(x = species, y = number), stat = "identity") + 
    ylim(0,50) + theme_bw() + 
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(g, 
    file = "./figure_S7_synteny_network_quality/synteny_network_quality_notlinking_species_hist.pdf")

g <- ggplot(data.frame(value = unlist(res_link_sum)), aes(x = value)) + 
    ylim(0, 50) + xlim(0, 15000) + theme_bw() + geom_histogram(binwidth = 100) + 
    theme(panel.grid = element_blank())
ggsave(g, file = "./figure_S7_synteny_network_quality/sum_genes_scaffold_linking_hist.pdf")

g <- ggplot(data.frame(value = unlist(res_nlink_sum)), aes(x = value)) + 
    ylim(0,50) + xlim(0, 15000) + theme_bw() + geom_histogram(binwidth = 100) + 
    theme(panel.grid = element_blank())
ggsave(g, 
       file = "./figure_S7_synteny_network_quality/sum_genes_scaffold_notlinking_hist.pdf")

## set colour according to the number of neightbour genes
tmp <- unlist(res_link_sum)
reliability[names(tmp[tmp < 15])] <- "red"
reliability[names(tmp[tmp >= 15 & tmp < 25])] <- "yellow"
reliability[names(tmp[tmp >= 25])] <- "green"

## plot network that contains features that link to other
bin_mat_complete_cut <- bin_mat_complete[inds_keep, inds_keep]
mat_type_cut <- mat_type[inds_keep, inds_keep]
net <- graph_from_adjacency_matrix(bin_mat_complete_cut, weighted = TRUE, 
    mode = "undirected", diag = FALSE)
plot(net, vertex.label.cex = 0.1, vertex.size = 5, 
     vertex.color = reliability[inds_keep], edge.arrow.size = 0.1)

## validation with PLAZA
## 1 --> short protein of length 89 --> no colinear genes found for this input gene
## 2 --> --> no colinear genes found for this input gene
## 3 --> links to selmo?
## 4 --> vitvi_GSVIVT01032968001, --> deleted region Manes_09G185300-09G188900
## 5 duplicate, that does not link to other
i <- 18
colnames(bin_mat_complete_cut)[i]
mat_type_cut[which(bin_mat_complete_cut[, i]!=0), i]
df <- data.frame(names=names(unlist(sort(table(mat_type_cut)[-1]))), value=as.vector(sort(table(mat_type_cut)[-1])))
df$names <- factor(df$names, levels=df$names[order(df$value)])
g <- ggplot(df) + geom_bar(aes(x=names, y=value), stat="identity") + theme_bw() + theme(axis.text.x=element_text(angle = 90, hjust = 1))
ggsave(g, file="./figure_S7_synteny_network_quality/type_link_synteny.pdf")
##

setwd("H:/Documents/03_Syntenic_linkage/01_Data/synteny_network_results/")
setwd("~/winhome/Documents/03_Syntenic_linkage/01_Data/synteny_network_results/")
save(bin_mat_complete_cut, mat_type_cut, file = "bin_mat_complete_cut_tandemgenes.RData")
load("bin_mat_complete_cut_tandemgenes.RData")
attr <- openxlsx::read.xlsx(
    "~/winhome/Documents/03_Syntenic_linkage/01_Data/synteny_network_results/genes_synteny.xlsx", 
    sheet = "synteny_genes_attributes")

net <- graph_from_adjacency_matrix(bin_mat_complete_cut, 
    weighted = TRUE, mode = "undirected", diag = FALSE)

## write table with attributes
## weighted version of degree (summing up the edge weights of the adjacent 
## edges for each vertex)
net_strength <- strength(net, mode = "all", loops = FALSE)
## numer of adjacent edges
net_degree <- degree(net, mode = "all", loops = FALSE)
df <- data.frame(names = names(net_strength), strength = net_strength)
df$degree <- net_degree
df$cluster <- character(nrow(df))
df$R4C <- rep(0, nrow(df))
df$R2X <- rep(0, nrow(df))
df$R4A <- rep(0, nrow(df))
df$Other <- rep(0, nrow(df))
df$nd <- rep(0, nrow(df))
df$sum <- rep(0, nrow(df))

for (i in 1:nrow(df)) {
    ind <- which(attr[, 1] == rownames(df)[i] )
    df$cluster[i] <- attr[ind, "classification_community_apclust_dist"]
    type <- attr[ind, "type_pPAP_reviewed"]
    type <- strsplit(type, split = "/")[[1]]
    type_t <- table(type)
    r4c <- as.numeric(type_t["R-4-C"])
    r2x <- as.numeric(type_t["R-2-X"])
    r4a <- as.numeric(type_t["R-4-A"])
    other <- as.numeric(type_t["Other"])
    nd <- as.numeric(type_t["nd"])
    df$R4C[i] <- ifelse(is.na(r4c), 0, r4c)
    df$R2X[i] <- ifelse(is.na(r2x), 0, r2x)
    df$R4A[i] <- ifelse(is.na(r4a), 0, r4a)
    df$Other[i] <- ifelse(is.na(other), 0, other)
    df$nd[i] <- ifelse(is.na(nd), 0, nd)
    df$sum[i] <- sum(type_t)
}
write.table(df, file = "bin_mat_complete_cut_graphml_attributes.txt", 
    sep = "\t", dec = ".", quote = F, row.names = F)

## write graph to graphml format for use in Cytoscape
write_graph(net, file = "bin_mat_complete_cut_graphml.xml", format = "graphml")
net_nlink <- graph_from_adjacency_matrix(
    bin_mat_complete[!inds_keep, !inds_keep], mode = "undirected", diag = FALSE)
write_graph(net_nlink, 
    file = "./figure_S7_synteny_network_quality/bin_mat_complete_notlink_graphml.xml", 
    format = "graphml")

df <- data.frame(names = rownames(bin_mat_complete[!inds_keep, !inds_keep]))
df$type <- attr[match(df[, 1], attr[, 1]), 3]
write.table(df, 
    file = "./figure_S7_synteny_network_quality/bin_mat_complete_notlink_type.txt", 
    row.names = FALSE, col.names = TRUE, quote = FALSE)
g <- ggplot(df) + geom_bar(aes(x = type), stat = "count") + ylim(0, 200) + theme_bw()
ggsave(g, 
    file = "./figure_S7_synteny_network_quality/barplot_identity_notlinking.pdf")

## create again a network from the binary matrix 
net <- graph_from_adjacency_matrix(bin_mat_complete_cut, mode="undirected", 
    diag = FALSE, weighted = TRUE)
plot(net, vertex.label.cex = 0.1, vertex.size = 0.1, vertex.color = reliability,
     edge.arrow.size = 0.1, edge.width = E(net)$weight)

## only continue with nodes from the big component
names_net <- names(which(membership(components(net)) == 1))
net2 <- graph_from_adjacency_matrix(bin_mat_complete_cut[names_net, names_net], 
    mode = "undirected", diag = FALSE, weighted = TRUE)

## check for identity (type) of PKS sequences
gene_name_bin_mat <- strsplit(rownames(
    bin_mat_complete_cut[names_net, names_net]), split = "/")
supp <- xlsx::read.xlsx("../pks_genes_tree_id_type.xlsx", 
    sheetIndex = 1, stringsAsFactors = FALSE)
supp <- supp[, which(colnames(supp) == "ID"):which(colnames(supp) == "species")]
gene_name_bin_mat <- lapply(gene_name_bin_mat, 
    function(x) supp[match(x, supp[, "ID"]), "pPAP_reviewed"])
which(!unlist(lapply(1:length(gene_name_bin_mat), 
    function(z) { 
        tmp1 <- lapply(gene_name_bin_mat, length)[[z]] 
        tmp2 < lapply(strsplit(rownames(bin_mat_complete_cut), split = "/"), length)[[z]]
        tmp1 == tmp2
    })
))
colour_name_bin_mat <- gene_name_bin_mat

## colour according to pPAP classification/type
colour_name_bin_mat <- lapply(colour_name_bin_mat, function(x) ifelse(x=="R-4-A", "red", x))
colour_name_bin_mat <- lapply(colour_name_bin_mat, function(x) ifelse(x=="R-4-C", "blue", x))
colour_name_bin_mat <- lapply(colour_name_bin_mat, function(x) ifelse(x=="Other", "green", x))
colour_name_bin_mat <- lapply(colour_name_bin_mat, function(x) ifelse(x=="R-2-X", "yellow", x))
colour_name_bin_mat <- lapply(colour_name_bin_mat, function(x) ifelse(x=="nd", "black", x))

# ## colour according to BLAST classification/type
# colour_name_bin_mat <- lapply(colour_name_bin_mat, function(x) ifelse(x=="CHS", "red", x))
# colour_name_bin_mat <- lapply(colour_name_bin_mat, function(x) ifelse(x=="STS", "blue", x))
# colour_name_bin_mat <- lapply(colour_name_bin_mat, function(x) ifelse(x %in% c("PKS", "PKSA", "PKSB", "PKSC"), "red4", x))
# colour_name_bin_mat <- lapply(colour_name_bin_mat, function(x) ifelse(x=="BCS", "springgreen", x))
# colour_name_bin_mat <- lapply(colour_name_bin_mat, function(x) ifelse(x=="unknown", "black", x))
# colour_name_bin_mat <- lapply(colour_name_bin_mat, function(x) ifelse(x=="ADS", "sienna", x))
# colour_name_bin_mat <- lapply(colour_name_bin_mat, function(x) ifelse(x=="CCS", "orange", x))
# colour_name_bin_mat <- lapply(colour_name_bin_mat, function(x) ifelse(x=="PPS", "slateblue", x))
# colour_name_bin_mat <- lapply(colour_name_bin_mat, function(x) ifelse(x=="PPAS", "palegreen", x))
# colour_name_bin_mat <- lapply(colour_name_bin_mat, function(x) ifelse(x=="HCS", "violet", x))
# colour_name_bin_mat <- lapply(colour_name_bin_mat, function(x) ifelse(x=="BBS", "mediumorchid", x))
# colour_name_bin_mat <- lapply(colour_name_bin_mat, function(x) ifelse(x=="DPS", "palegreen", x))
# colour_name_bin_mat <- lapply(colour_name_bin_mat, function(x) ifelse(x=="PVS", "palegreen", x))

## cluster detection 
## https://stackoverflow.com/questions/9471906/what-are-the-differences-between-community-detection-algorithms-in-igraph
ec <- edge.betweenness.community(net2, modularity = TRUE)
fc <- fastgreedy.community(net2, modularity = TRUE)
wc <- walktrap.community(net2, modularity = TRUE, steps = 5)
le <- leading.eigenvector.community(net2, steps = 15)
sc <- spinglass.community(net2, spins = 100)
lp <- label.propagation.community(net2)

## combine all cluster detection results to a matrix 
members_comm <- cbind(ec = membership(ec), fc = membership(fc), 
    wc = membership(wc), le = membership(le), sc = membership(sc), 
    lp = membership(lp))

coords <- layout.auto(net2)
##coords <- layout.graphopt(net2)

par(mar=c(5.1, 4.1, 4.1, 2.1))
par(mar=c(0.1, 0.1, 0.9, 0.1))
pdf("~/winhome/Documents/graph_network%03d.pdf", paper="a4r", width=11, height=7)
par(mfrow=c(1, 3))
par(mar=c(0.1, 0.1, 0.8, 0.1))

## plot the network with the results from the different cluster algorithms
plot(net2, vertex.shape = "pie", 
     vertex.pie = lapply(colour_name_bin_mat, function(x) rep(1, length(x))), 
     vertex.pie.color = colour_name_bin_mat, vertex.label = "", 
     vertex.size = 1.5, edge.arrow.size = 0.1, edge.width = E(net2)$weight,
     vertex.frame.color = adjustcolor("black", alpha.f = 0.0), layout = coords)
plot(net2, vertex.color = membership(ec), vertex.size = 1.5, vertex.label = "", 
     layout = coords, main = "edge betweenness", edge.width = 0.2, 
     vertex.frame.color = adjustcolor("black", alpha.f = 0.0))
plot(net2, vertex.color = membership(fc), vertex.size = 1.5, vertex.label = "", 
     layout = coords, main = "fastgreedy", edge.width = 0.2, 
     vertex.frame.color = adjustcolor("black", alpha.f = 0.0))

plot(net2, vertex.shape = "pie", 
     vertex.pie = lapply(colour_name_bin_mat, function(x) rep(1, length(x))), 
     vertex.pie.color = colour_name_bin_mat, vertex.label = "", 
     vertex.size = 1.5, edge.arrow.size = 0.1, edge.width = E(net2)$weight,
     vertex.frame.color = adjustcolor("black", alpha.f = 0.0), layout = coords)
plot(net2, vertex.color = membership(wc), vertex.size = 1.5, vertex.label = "", 
     layout = coords, main = "walktrap", edge.width = 0.2, 
     vertex.frame.color = adjustcolor("black", alpha.f = 0.0))
plot(net2, vertex.color = membership(le), vertex.size = 1.5, vertex.label = "", 
     layout = coords, main = "leading eigenvector", edge.width = 0.2, 
     vertex.frame.color = adjustcolor("black", alpha.f = 0.0))
plot(net2, vertex.shape = "pie", 
     vertex.pie = lapply(colour_name_bin_mat, function(x) rep(1, length(x))), 
     vertex.pie.color = colour_name_bin_mat, vertex.label = "", 
     vertex.size = 1.5, edge.arrow.size = 0.1, edge.width = E(net2)$weight,
     vertex.frame.color = adjustcolor("black", alpha.f = 0.0), layout = coords)
plot(net2, vertex.color = membership(sc), vertex.size = 1.5, vertex.label = "", 
     layout = coords, main = "spinglass", edge.width = 0.2, 
     vertex.frame.color = adjustcolor("black", alpha.f = 0.0))
plot(net2, vertex.color = membership(lp), vertex.size = 1.5, vertex.label = "", 
     layout = coords, main = "label propagation", edge.width = 0.2, 
     vertex.frame.color = adjustcolor("black", alpha.f = 0.0))
dev.off()

## plot with gene identity
plot(net, vertex.shape = "pie", 
     vertex.pie = lapply(colour_name_bin_mat, function(x) rep(1, length(x))), 
     vertex.pie.color = colour_name_bin_mat, vertex.label = "", 
     vertex.size = 1.5, edge.arrow.size = 0.1, edge.width = E(net)$weight,
     vertex.frame.color = adjustcolor("black", alpha.f = 0.0), layout = coords)
plot(net, vertex.shape = "pie", 
     vertex.pie = lapply(colour_name_bin_mat, function(x) rep(1, length(x))), 
     vertex.pie.color = colour_name_bin_mat, vertex.label.cex = 0.05,
     vertex.size = 1.5, edge.arrow.size = 0.1, edge.width = 0.2,
     vertex.frame.color = adjustcolor("black", alpha.f = 0.0))

## create a plot that shows the number of identity of links (number of methods)
table_type <- sort(table(mat_type_cut)[-1])
names(table_type)[names(table_type) == "iadhore/iadhore_mcl/mcscanx/mcscanx_mcl"] <- "iadhore/iadhore_mcl/\nmcscanx/mcscanx_mcl"
names(table_type)[names(table_type) == "iadhore/iadhore_mcl/mcscanx_mcl"] <- "iadhore/iadhore_mcl/\nmcscanx_mcl"
names(table_type)[names(table_type) == "iadhore/mcscanx/mcscanx_mcl"] <- "iadhore/mcscanx/\nmcscanx_mcl"
names(table_type)[names(table_type) == "iadhore_mcl/mcscanx/mcscanx_mcl"] <- "iadhore_mcl/mcscanx/\nmcscanx_mcl"
names(table_type)[names(table_type) == "iadhore/iadhore_mcl/mcscanx"] <- "iadhore/iadhore_mcl/\nmcscanx"
barplot(table_type, las=2, cex.names=0.4)

## clustering: calculate distance between measures and make a union using apcluster
genes_attr <- openxlsx::read.xlsx("genes_synteny.xlsx", 
    sheet = "synteny_genes_attributes")

## function to calculate distance: if the cluster identity is identical, then
## distance is 0, if the cluster identity between two syntenic regions is not 
## identical, then distance is 1
calculate_dist_mat <- function(community_method = "edge_betweenness_community") {
    
    ## create matrix
    dist <- matrix(0, nrow = nrow(genes_attr), ncol = nrow(genes_attr))
    rownames(dist) <- colnames(dist) <- genes_attr[, "gene_tandem"]
    norm_factor <- length(unique(genes_attr[, community_method]))
    for (i in 1:nrow(dist)) {
        ## if the values match --> dist=0, else dist=1
        values <- genes_attr[, community_method]
        values_i <- genes_attr[i, community_method]
        dist[i, ] <- as.numeric(values == values_i) 
    }
    return(dist)
}

## apply the function
dist_ec <- calculate_dist_mat("edge_betweenness_community")
dist_fc <- calculate_dist_mat("fastgreedy_community")
dist_wc <-  calculate_dist_mat("walktrap_community")
dist_le <- calculate_dist_mat("leading_eigenvector_community")
dist_sc <- calculate_dist_mat("spinglass_community")
dist_lp <- calculate_dist_mat("label_propagation_community")

## get the sum of all 
dist_total <- dist_ec + dist_fc + dist_wc + dist_le + dist_sc + dist_lp 

## use apcluster to find consensus clusters using the correlation between 
## distances or between distances itself (distances itself will be used later)
dist_total_d <- cor(dist_total, method = "spearman")
ap_cor <- apcluster::apcluster(s = abs(dist_total_d), convits = 1000, maxits = 10000)
ap_dis <- apcluster::apcluster(s = dist_total/6, convits = 1000, maxits = 10000)

ap_cor_l <- lapply(ap_cor@clusters, function(x) names(x))
ap_cor_l <- lapply(1:length(ap_cor_l), function(x) {
    y <- rep(x, length(ap_cor_l[[x]])); names(y) <- ap_cor_l[[x]]; y})
ap_cor_l <- unlist(ap_cor_l)[sort(names(unlist(ap_cor_l)))]
ap_dis_l <- lapply(ap_dis@clusters, function(x) names(x))
ap_dis_l <- lapply(1:length(ap_dis_l), function(x) {
    y <- rep(x, length(ap_dis_l[[x]])); names(y) <- ap_dis_l[[x]]; y})
ap_dis_l <- unlist(ap_dis_l)[sort(names(unlist(ap_dis_l)))]
write.table(cbind(ap_cor_l, ap_dis_l), file = "new_partition_all_indices.csv", sep = ",")

## load the attribute xlsx file
attr <- openxlsx::read.xlsx("../synteny_network_results/genes_synteny.xlsx", 
    sheet="synteny_genes_attributes")
attr <- attr[order(attr[, 1]), ]
class_comm <- attr[, "classification_community_apclust_dist"]
attr <- attr[!class_comm == "-", ]
class_comm <- class_comm[!class_comm == "-"]

## colour the R-4-C-enriched clusters 2, 3 and 5
col_vertex <- rep("black", length(class_comm))
col_vertex[class_comm == 2] <- "red"
col_vertex[class_comm == 3] <- "blue"
col_vertex[class_comm == 5] <- "green" 

## plot
plot(net2, vertex.color = col_vertex, vertex.size = 1.5, vertex.label = "", 
     layout = coords, main = "classification_community", edge.width = 0.2, 
     vertex.frame.color = adjustcolor("black", alpha.f = 0.0))

## remove entries that do not have information on community apclust clusters
attr <- attr[attr[, "classification_community_apclust_dist"]  != "-",]
sort(table(attr[, "classification_community_apclust_dist"]))
cluster1 <- attr[attr[, "classification_community_apclust_dist"] == "1", ]
genes1 <- cluster1[, 1]
cluster2 <- attr[attr[, "classification_community_apclust_dist"] == "2", ]
genes2 <- cluster2[, 1] ## high
cluster3 <- attr[attr[, "classification_community_apclust_dist"] == "3", ]
genes3 <- cluster3[, 1] ## high
cluster4 <- attr[attr[, "classification_community_apclust_dist"] == "4", ]
genes4 <- cluster4[, 1] 
cluster5 <- attr[attr[, "classification_community_apclust_dist"] == "5", ]
genes5 <- cluster5[, 1] ## high
cluster6 <- attr[attr[, "classification_community_apclust_dist"] == "6", ]
genes6 <- cluster6[, 1]
cluster7 <- attr[attr[, "classification_community_apclust_dist"] == "7", ]
genes7 <- cluster7[, 1]
cluster8 <- attr[attr[, "classification_community_apclust_dist"] == "8", ]
genes8 <- cluster8[, 1]
cluster9 <- attr[attr[, "classification_community_apclust_dist"] == "9", ]
genes9 <- cluster9[, 1]
cluster10 <- attr[attr[, "classification_community_apclust_dist"] == "10", ]
genes10 <- cluster10[, 1]
cluster11 <- attr[attr[, "classification_community_apclust_dist"] == "11", ]
genes11 <- cluster11[, 1]
cluster12 <- attr[attr[, "classification_community_apclust_dist"] == "12", ]
genes12 <- cluster12[, 1]
cluster13 <- attr[attr[, "classification_community_apclust_dist"] == "13", ]
genes13 <- cluster13[, 1]
cluster14 <- attr[attr[, "classification_community_apclust_dist"] == "14", ]
genes14 <- cluster14[, 1]
cluster15 <- attr[attr[, "classification_community_apclust_dist"] == "15", ]
genes15 <- cluster15[, 1]
cluster16 <- attr[attr[, "classification_community_apclust_dist"] == "16", ]
genes16 <- cluster16[, 1]
cluster17 <- attr[attr[, "classification_community_apclust_dist"] == "17", ]
genes17 <- cluster17[, 1]
cluster18 <- attr[attr[, "classification_community_apclust_dist"] == "18", ]
genes18 <- cluster18[, 1]
cluster19 <- attr[attr[, "classification_community_apclust_dist"] == "19", ]
genes19 <- cluster19[, 1]
cluster20 <- attr[attr[, "classification_community_apclust_dist"] == "20", ]
genes20 <- cluster20[, 1]
cluster21 <- attr[attr[, "classification_community_apclust_dist"] == "21", ]
genes21 <- cluster21[, 1]
cluster22 <- attr[attr[, "classification_community_apclust_dist"] == "22", ]
genes22 <- cluster22[, 1]
cluster23 <- attr[attr[, "classification_community_apclust_dist"] == "23", ]
genes23 <- cluster23[, 1]
cluster24 <- attr[attr[, "classification_community_apclust_dist"] == "24", ]
genes24 <- cluster24[, 1]
cluster25 <- attr[attr[, "classification_community_apclust_dist"] == "25", ]
genes25 <- cluster25[, 1]
cluster26 <- attr[attr[, "classification_community_apclust_dist"] == "26", ]
genes26 <- cluster26[, 1]
cluster27 <- attr[attr[, "classification_community_apclust_dist"] == "27", ]
genes27 <- cluster27[, 1]
cluster28 <- attr[attr[, "classification_community_apclust_dist"] == "28", ]
genes28 <- cluster28[, 1]
cluster29 <- attr[attr[, "classification_community_apclust_dist"] == "29", ]
genes29 <- cluster29[, 1]
cluster30 <- attr[attr[, "classification_community_apclust_dist"] == "30", ]
genes30 <- cluster30[, 1]
cluster31 <- attr[attr[, "classification_community_apclust_dist"] == "31", ]
genes31 <- cluster31[, 1]
cluster32 <- attr[attr[, "classification_community_apclust_dist"] == "32", ]
genes32 <- cluster32[, 1]

genes_chs <- c(genes1, genes2, genes3, genes4, genes5, genes6, genes7,  genes8, genes9,#,  ##genes10,
                genes17, genes20, ## genes11, genes12, genes13, genes14, genes15, genes16, genes18, genes19,
               genes21, genes24, genes26, genes27, genes28, ##genes22, genes23,  genes25, genes29, genes30,
               genes31) #genes32)
type_chs <- rep("circle", length(genes_chs))
type_chs[attr[attr[,"gene_tandem"] %in% genes_chs, "conn_R-4-C_largenetwork"] == "yes" & attr[attr[, "gene_tandem"] %in% genes_chs, "gene_tandem"] %in% genes_chs] <- "square"

## plot the R-4-C-enriched clusters
cols <- rep("grey", length(genes_chs))
cols[genes_chs %in% genes2] <- "blue"
cols[genes_chs %in% genes3] <- "green"
cols[genes_chs %in% genes5] <- "yellow"
bin_12 <- bin_mat_complete_cut[genes_chs, genes_chs]
bin_12_g <- graph_from_adjacency_matrix(bin_12, mode = "undirected", 
    weighted = TRUE, diag = FALSE)
plot(bin_12_g, vertex.color = cols, vertex.label = "", vertex.size = 1, 
     vertex.shape = type_chs)
genes1 <- unlist(strsplit(cluster1[, 1], split = "/"))

## get information on when the nodes change in identity (e.g. by tandem 
## duplications) or when the type of PKS changes
bmat <- bin_mat_complete_cut[names_net, names_net]
res_change <- vector("list", nrow(bmat))
adj_mat_change <- matrix(0, nrow = nrow(bmat), ncol = ncol(bmat))
rownames(adj_mat_change) <- colnames(adj_mat_change) <- rownames(bmat)

## iterate through the large component of the network
for (i in 1:nrow(bmat)) {
    conn <- which(bmat[, i] >= 0.25)
    type <- gene_name_bin_mat
    type_notsame <- NULL
    if (length(conn) > 0) {
        notsame <- lapply(1:length(conn), function(x) {
            identical(type[[conn[x]]], type[[i]])
        })
        notsame <- unlist(notsame)
        notsame <- !notsame
        names_notsame <- names(conn[notsame])
        type_notsame <- type[conn[notsame]]
        names(type_notsame) <- names_notsame
        type_notsame <- lapply(type_notsame, function(x) paste(x, collapse="/"))
    }
    parent <- rbind(rownames(bmat)[i], paste(type[[i]], collapse = "/"), "")
    child <- rbind(names(type_notsame), unlist(type_notsame), 
                   bmat[rownames(bmat)[i], names(type_notsame)])
    
    adj_mat_change[rownames(bmat)[i], names(type_notsame)] <- bmat[rownames(bmat)[i], names(type_notsame)]
    
    res_change[[i]] <- cbind(parent, child)
}
mat_change <- matrix("", nrow = 3 * nrow(bmat), 
    ncol = max(unlist(lapply(res_change, ncol))))
for (i in 1:length(res_change)) {
    rc_i <- res_change[[i]]
    mat_change[c(i * 3 - 2):(i * 3), 1:ncol(rc_i)] <- rc_i
}
setwd("~/winhome/Documents/03_Syntenic_linkage/01_Data/synteny_network_results/")
write.table(mat_change, "syntenic_linkage_1_change_function.csv", sep="\t", quote=FALSE)

## plot the networks that contain information on the changes
adj_mat_change_075 <- adj_mat_change_05 <- adj_mat_change
adj_mat_change_05[adj_mat_change_05 < 0.5] <- 0
adj_mat_change_075[adj_mat_change_075 < 0.75] <- 0
net_change_025 <- graph_from_adjacency_matrix(adj_mat_change, mode="undirected", diag=FALSE, weighted=TRUE)
net_change_05 <- graph_from_adjacency_matrix(adj_mat_change_05, mode="undirected", diag=FALSE, weighted=TRUE)
net_change_075 <- graph_from_adjacency_matrix(adj_mat_change_075, mode="undirected", diag=FALSE, weighted=TRUE)

par(mfrow=c(2, 2))
par(mar=c(0.1, 0.1, 0.8, 0.1))
plot(net2, vertex.shape="pie", vertex.pie=lapply(colour_name_bin_mat, function(x) rep(1, length(x))), vertex.pie.color=colour_name_bin_mat, vertex.label="", 
     vertex.size=1.5, edge.arrow.size=0.1, edge.width=E(net2)$weight,
     vertex.frame.color=adjustcolor("black", alpha.f=0.0), layout=coords, rescale=T) 
plot(net_change_025, vertex.shape="pie", vertex.pie=lapply(colour_name_bin_mat, function(x) rep(1, length(x))), vertex.pie.color=colour_name_bin_mat, vertex.label="", 
     vertex.size=1.5, edge.arrow.size=0.1, edge.width=E(net_change)$weight,
     vertex.frame.color=adjustcolor("black", alpha.f=0.0), layout=coords, rescale=T, main=">=0.25")
plot(net_change_05, vertex.shape="pie", vertex.pie=lapply(colour_name_bin_mat, function(x) rep(1, length(x))), vertex.pie.color=colour_name_bin_mat, vertex.label="", 
     vertex.size=1.5, edge.arrow.size=0.1, edge.width=E(net_change)$weight,
     vertex.frame.color=adjustcolor("black", alpha.f=0.0), layout=coords, rescale=T, main=">=0.5")
plot(net_change_075, vertex.shape="pie", vertex.pie=lapply(colour_name_bin_mat, function(x) rep(1, length(x))), vertex.pie.color=colour_name_bin_mat, vertex.label="", 
     vertex.size=1.5, edge.arrow.size=0.1, edge.width=E(net_change)$weight,
     vertex.frame.color= "#00000000", layout=coords, rescale=T, main=">=0.75")
dev.off()
## end change network


## find the neighbours of BPS and STS from Vitis vinifera and Arachis genus
## BPS
neighbors(net, v="denof_XP_020690097.1/denof_XP_020690099.1", mode="all")

## STS vitvi
## find region names
gene_name_bin_mat[grep(rownames(bin_mat_complete_cut), 
    pattern = "vitvi")[c(2, 5)]]
name_bin_mat[grep(rownames(bin_mat_complete_cut), pattern = "vitvi")[c(2, 5)]]
tmp1 <- "vitvi_GSVIVT01026213001/vitvi_GSVIVT01026220001"
gene_name_bin_mat[neighbors(net, v = tmp, mode = "all")]
tmp2 <- "vitvi_GSVIVT01010554001/vitvi_GSVIVT01010556001/vitvi_GSVIVT01010557001/vitvi_GSVIVT01010561001/vitvi_GSVIVT01010563001/vitvi_GSVIVT01010565001/vitvi_GSVIVT01010568001/vitvi_GSVIVT01010570001/vitvi_GSVIVT01010572001/vitvi_GSVIVT01010574001/vitvi_GSVIVT01010578001/vitvi_GSVIVT01010579001/vitvi_GSVIVT01010580001/vitvi_GSVIVT01010581001/vitvi_GSVIVT01010582001/vitvi_GSVIVT01010583001/vitvi_GSVIVT01010584001/vitvi_GSVIVT01010585001/vitvi_GSVIVT01010589001/vitvi_GSVIVT01010590001"
gene_name_bin_mat[neighbors(net, v = tmp2, mode = "all")]
name_bin_mat[neighbors(net, v = tmp2, mode = "all")]

## STS aradu/araip
## find region names
gene_name_bin_mat[grep(rownames(bin_mat_complete_cut), 
    pattern = "aradu|araip")[c(2, 5, 6, 10, 13, 14)]]
name_bin_mat[grep(rownames(bin_mat_complete_cut), 
    pattern = "aradu|araip")[c(2, 5, 6, 10, 13, 14)]]

tmp3 <- "aradu_XP_015959399.1/aradu_XP_015959400.1/aradu_XP_015959401.1/aradu_XP_015959402.1/aradu_XP_015959403.1/aradu_XP_015959405.1/aradu_XP_015959406.1/aradu_XP_015959407.1/aradu_XP_015959435.1/aradu_XP_015959436.1/aradu_XP_015959437.1/aradu_XP_015959438.1/aradu_XP_015959439.2/aradu_XP_015959440.1/aradu_XP_015959441.1/aradu_XP_015959445.1/aradu_XP_015959446.1/aradu_XP_015959447.1/aradu_XP_015959448.1/aradu_XP_015959449.1/aradu_XP_015959453.1/aradu_XP_020995724.1/aradu_XP_020995733.1"
tmp4 <- "araip_XP_016197777.1/araip_XP_016197779.1/araip_XP_016197780.1/araip_XP_016197781.1/araip_XP_016197783.1/araip_XP_016197784.1/araip_XP_016197787.1/araip_XP_016197788.1/araip_XP_016197789.2/araip_XP_016197793.1/araip_XP_016197798.1/araip_XP_016197799.2/araip_XP_016197801.1/araip_XP_016197802.1/araip_XP_016197804.1/araip_XP_016197805.1/araip_XP_016197807.1/araip_XP_020976418.1/araip_XP_020976446.1/araip_XP_020978004.1"
names_STS_seed <- c(tmp1, tmp2, "aradu_XP_015940991.1", "aradu_XP_015958332.1",
    tmp3, "araip_XP_016180692.1", "araip_XP_016196906.1/araip_XP_016196907.1",
    tmp4)
names_STS <- c(names_STS_seed,
    rownames(bin_mat_complete_cut)[neighbors(net, v = tmp1, mode = "all")],
    rownames(bin_mat_complete_cut)[neighbors(net, v = tmp2, mode = "all")],
    rownames(bin_mat_complete_cut)[neighbors(net, v = "aradu_XP_015940991.1", mode = "all")],
    rownames(bin_mat_complete_cut)[neighbors(net, v = "aradu_XP_015958332.1", mode = "all")],
    rownames(bin_mat_complete_cut)[neighbors(net, v = tmp3, mode = "all")],
    rownames(bin_mat_complete_cut)[neighbors(net, v = "araip_XP_016180692.1", mode = "all")],
    rownames(bin_mat_complete_cut)[neighbors(net, v = "araip_XP_016196906.1/araip_XP_016196907.1", mode = "all")],
    rownames(bin_mat_complete_cut)[neighbors(net, v = tmp4, mode = "all")])
names_STS <- unique(names_STS)
bin_mat_sts <- bin_mat_complete_cut[names_STS, names_STS]
vertex_size <- rownames(bin_mat_sts) %in% names_STS_seed
net_STS <- graph_from_adjacency_matrix(bin_mat_sts, 
                        mode = "undirected", diag = FALSE, weighted = TRUE)
plot(net_STS, vertex.shape = "pie", 
     vertex.pie = lapply(colour_name_bin_mat[match(names_STS, rownames(bin_mat_complete_cut))], function(x) rep(1, length(x))), 
     vertex.pie.color = colour_name_bin_mat[match(names_STS, rownames(bin_mat_complete_cut))], vertex.label = "",
     vertex.size = ifelse(vertex_size, 2.5, 1.5), edge.arrow.size = 0.1, 
     edge.width = 0.2, vertex.frame.color = adjustcolor("black", alpha.f = 0.0))
rownames(bin_mat_complete_cut) %in% names_STS
