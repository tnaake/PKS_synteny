## setwd
setwd("~/AG-Fernie/Thomas/Data/synteny")

library(igraph)

## get pks genes from output of OrthoFinder or MCL
genes_table <- read.table("./Results_Oct26/family_orthogroups.txt", sep = "\t", 
    header = FALSE, stringsAsFactors = FALSE)
## for OrthoFinder
pks_genes <- sort(genes_table[genes_table[, 2] == "OG0000260", 1]) 
genes_table_mcl <- read.table("./Results_Oct26/family_orthogroups_mcl.txt", 
    sep = "\t", header = FALSE, stringsAsFactors = FALSE)
## forMCL
pks_genes_mcl <- sort(genes_table_mcl[genes_table_mcl[, 2] == "group_337", 1]) 

## vector with all pks genes identified by orthofinder and mcl, use this to 
## create the bin_mat_... matrices
pks_genes_all <- sort(unique(c(pks_genes, pks_genes_mcl)))


## i-ADHore
find_syntenic_regions_iadh <- function(files, pks_genes) {
    
    ## create vector to store syntenic genes (background)
    res <- vector("list", length(files))
    names(res) <- unlist(lapply(strsplit(files, split = "output_"), "[", 2))
    
    ## iterate through files
    for (i in seq_along(files)) {
        
        res[[i]] <- list(character(), character())
        names(res[[i]]) <- c("fg", "bg") ## stores back- and foreground
        
        ## create list to store syntenic genes surrounding pks genes
        res_pks <- vector("list", length(pks_genes_all))
        names(res_pks) <- pks_genes_all
        
        files_i <- paste0(files[i], "/anchorpoints.txt")
        mult_pairs <- NULL
        mult_pairs <- tryCatch(read.csv(files_i, header = FALSE, 
            sep = "\t", skip = 1, fill = TRUE, comment.char = "#", 
            stringsAsFactors = FALSE), error = function(e) NULL)
        
        if (!is.null(mult_pairs)) {
            region <- paste(mult_pairs[, 2], mult_pairs[, 3], sep = "_")
            gene_1_ind <- which(mult_pairs[, 4] %in% pks_genes)
            gene_2_ind <- which(mult_pairs[, 5] %in% pks_genes)

            if (length(gene_1_ind) != 0) {
            for (j in 1:length(gene_1_ind )) {
                ## get region of interest that contains syntenic genes
                roi <- region[gene_1_ind[j]]
                ## add genes in syntenic region to existing ones
                res_pks[[ mult_pairs[gene_1_ind[j], 4] ]] <- mult_pairs[region == roi, 4]
            }}
            
            if (length(gene_2_ind) != 0) {
            for (j in 1:length(gene_2_ind )) {
                ## get region of interest that contains syntenic genes
                roi <- region[gene_2_ind[j]]
                ## add genes in syntenic region to the pks gene slot
                res_pks[[ mult_pairs[gene_2_ind[j], 5] ]] <- mult_pairs[region == roi, 5]
            }}
            
            ## write all syntenic genes to res (background)
            res[[i]][["fg"]] <- res_pks
            res[[i]][["bg"]] <- c(mult_pairs[, 4], mult_pairs[, 5])
            
        } else {
            print(files_i)
            res[[i]][["fg"]] <- res_pks
            res[[i]][["bg"]] <- character()
        }
    }
    ## return list with proteins per pks gene (length == length(pks_genes_all))
    ## return vector with all syntenic genes 
    return(res)
}

## MCScanX
## functions 
find_syntenic_regions_mcsc <- function(files, pks_genes) {
    
    ## create vector to store syntenic genes (background)
    res <- vector("list", length(files)) 
    names(res) <- unlist(strsplit(files, split = "[.]collinearity"))
    
    ## iterate through files
    for (i in seq_along(files)) {
        
        res[[i]] <- list(character(), character())
        names(res[[i]]) <- c("fg", "bg") ## stores back- and foreground
        
        ## create list to store syntenic genes surrounding pks genes
        res_pks <- vector("list", length(pks_genes_all))
        names(res_pks) <- pks_genes_all
        
        coll_i <- NULL
        coll_i <- tryCatch(read.csv(files[i], sep = "\t", fill = TRUE, 
            stringsAsFactors = FALSE, header = FALSE, comment.char = "#"), 
            error = function(x) NULL)
        if (!is.null(coll_i)) {
            region <- unlist(lapply(strsplit(coll_i[, 1], split = "-"), "[", 1))
            ## get index
            gene_1_ind <- which(coll_i[, 2] %in% pks_genes)
            gene_2_ind <- which(coll_i[, 3] %in% pks_genes)
            
            if (length(gene_1_ind) != 0) {
                for (j in 1:length(gene_1_ind )) {
                    ## get region of interest that contains syntenic genes
                    roi <- region[ gene_1_ind[j] ]
                    ## add genes in syntenic region to the pks gene slot
                    res_pks[[ coll_i[gene_1_ind[j], 2] ]] <- coll_i[region == roi, 2]
            }}
            
            if (length(gene_2_ind) != 0) {
                for (j in 1:length(gene_2_ind )) {
                    ## get region of interest that contains syntenic genes
                    roi <- region[ gene_2_ind[j] ]
                    ## add genes in syntenic region to existing ones
                    res_pks[[ coll_i[gene_2_ind[j], 3] ]] <- coll_i[region == roi, 3]
            }}
            
            ## write all syntenic genes to res (background)
            res[[i]][["fg"]] <- res_pks
            res[[i]][["bg"]] <- c(coll_i[, 2], coll_i[, 3])
            
        } else {
            print(files[i])
            res[[i]][["fg"]] <- res_pks
            res[[i]][["bg"]] <- character()
        }
    }
    ## return list with proteins per pks gene (length == length(pks_genes_all))
    ## return vector with all syntenic genes 
    return(res)
}

## for OrthoFinder input
setwd("~/AG-Fernie/Thomas/Data/synteny/synteny_iadhore/output/")
output_files <- list.files()[grep(list.files(), pattern = "output_")]
## use here pks_genes 
reg_iadh_ofi <- find_syntenic_regions_iadh(output_files, pks_genes) 

## for MCL input
setwd("../output_mcl")
output_files <- list.files()[grep(list.files(), pattern = "output_")]
## use here pks_genes_mcl
reg_iadh_mcl <- find_syntenic_regions_iadh(output_files, pks_genes_mcl) 

## for orthoFinder input
setwd("~/AG-Fernie/Thomas/Data/synteny/synteny_MCScanX/MCScanX_gff")
collinearity_files <- list.files()[grep(list.files(), pattern = "[.]collinearity")]
## use here pks_genes 
reg_mcsc_ofi <- find_syntenic_regions_mcsc(collinearity_files, pks_genes)

## for MCL input 
setwd("~/AG-Fernie/Thomas/Data/synteny/synteny_MCScanX/MCScanX_gff_mcl")
collinearity_files <- list.files()[grep(list.files(), pattern = "[.]collinearity")]
reg_mcsc_mcl  <- find_syntenic_regions_mcsc(collinearity_files, pks_genes_mcl)

## make the union from all results --> one object that stores the syntenic 
## genes per species-pair and gene
setwd("~/winhome/Documents/03_Syntenic_linkage/01_Data/synteny_network_results/")
save(reg_iadh_ofi, reg_iadh_mcl, reg_mcsc_ofi, reg_mcsc_mcl, 
    file = "find_syntenic_regions_iadh_mcsc.RData")

## should not contain any NULL
fg_all <-   lapply(1:length(reg_iadh_ofi), function(x) {
    w <- lapply(1:length(reg_iadh_ofi[[x]][[1]]), function(y) {
        z <- unique(c(reg_iadh_ofi[[x]][[1]][[y]], reg_iadh_mcl[[x]][[1]][[y]], 
            reg_mcsc_ofi[[x]][[1]][[y]], reg_mcsc_mcl[[x]][[1]][[y]]))
        return(z)
    })
    names(w) <- names(reg_iadh_ofi[[x]][[1]])
    return(w)})
names(fg_all) <- names(reg_iadh_ofi)

bg_all <- lapply(1:length(reg_iadh_ofi), function(x) {unique(
    c(reg_iadh_ofi[[x]][[2]], reg_iadh_mcl[[x]][[2]], 
      reg_mcsc_ofi[[x]][[2]], reg_mcsc_mcl[[x]][[2]])) 
})
names(bg_all) <- names(reg_iadh_ofi)

## GO enrichment
## write GO terms to genes and create a list
setwd("~/winhome/Documents/03_Syntenic_linkage/01_Data/OrthoFinder/input/fasta_all")
go_files <- list.files()
go_files <- go_files[grep(go_files, pattern = "GO[.]out")]

go_l <- vector("list", length(go_files))
names(go_l) <- go_files

go_l_desc <- vector("list", length(go_files))
names(go_l_desc) <- go_files

for (j in seq_along(go_files)) {
    print(j)
    go_i <- read.csv(go_files[j], sep = "\t", stringsAsFactors = FALSE, 
        header = TRUE)
    sep0 <- lapply(nchar(go_i[, "goid"]), function(x) paste0(rep("0", 7-x), 
        collapse = ""))
    sep0 <- unlist(sep0)
    
    go_i[, "goid"] <- paste0(go_i[, "ontology"], ":", sep0, go_i[, "goid"], sep = "")
    
    go_id <- go_i[, "qpid"]
    go_id_u <- unique(go_id)
    
    go_l_species <- vector("list", length(go_id_u)) 
    names(go_l_species) <- go_id_u
    for (i in seq_along(go_id_u)) {
        inds <- which(go_id == go_id_u[i]) 
        go_l_species[[ go_id_u[i] ]] <- go_i[inds, "goid"]
        
    }
    
    go_l[[j]] <- go_l_species 
    go_l_desc[[j]] <- go_i[, c("qpid", "goid", "desc")]
}
save(go_l, go_l_desc, file = "go_terms_species_genes.RData")


## write function to perform Fisher test (test for enrichment)
fisher_test_all <- function(bg, fg) {
    
    universe <- unlist(bg)
    subset <- unlist(fg)
    
    ## get _terms
    subset_unique <- unique(subset)
    ## get table for the go terms in cluster
    table_subset <- table(subset)
    ## get table for the go terms in genome
    table_universe <- table(universe)
    
    p_value <- rep(NA, length(subset_unique))
    estimate <- rep(NA, length(subset_unique))
    
    res <- vector("list", length(subset_unique))
    
    ## do fisher_test_cluster for all clusters
    for (i in 1:length(subset_unique)) {
        ## get number of occurences of term in cluster and genome
        n_subset <- as.numeric(table_subset[subset_unique[i]])
        n_universe <- as.numeric(table_universe[subset_unique[i]])
            
        ## do fisher.exact test (only for known) to test for 
        ## overrepresentation (alternative="greater"):   
        ## number_of_go_term_in_cluster number_of_genes_in_cluster
        ## number_of_go_term_in_genome  number_of_genes_in_genome
        mm <- matrix(c(n_subset, length(subset) - n_subset, 
            n_universe - n_subset, 
            length(universe) - length(subset) - (n_universe - n_subset)), 
            nrow = 2, byrow = TRUE)
        res_i <- fisher.test(mm, alternative = "greater")
        p_value[i] <- res_i$p.value
        estimate[i] <- res_i$estimate
        names(p_value)[i] <- names(estimate)[i] <- subset_unique[i]
        res[[i]] <- res_i
    }
    names(res) <- subset_unique
    
    ## get number of total test conducted
    n_total <- length(unlist(lapply(res, "[", "p.value")))
    
    for (i in 1:length(subset_unique)) {
        ## apply multiple testing correction 
        res[[i]]$p_values_adj <- p.adjust(res[[i]]$p.value, method = "BH", 
        n = n_total)
    }
    return(res)
}

## open xlsx that contains information on syntenic clusters
attr <- openxlsx::read.xlsx(
    "~/winhome/Documents/03_Syntenic_linkage/01_Data/synteny_network_results/genes_synteny.xlsx", 
    sheet = "synteny_genes_attributes")
attr <- attr[attr[, "connecting"] == "yes", ]
attr <- attr[attr[, "classification_community_apclust_dist"] != "-", ]

## test for enrichment --> take clusters 2, 3, 5
clust_total <- table(unlist(strsplit(
    attr[!attr[, "classification_community_apclust_dist"] %in% "-", "type_pPAP_reviewed"], "/")))
clust_235 <- table(unlist(strsplit(
    attr[attr[, "classification_community_apclust_dist"] %in% c("2", "3", "5"), "type_pPAP_reviewed"], "/")))
freq_chs <- matrix(c(clust_235["R-4-C"], sum(clust_235) - clust_235["R-4-C"], 
    clust_total["R-4-C"] - clust_235["R-4-C"], 
    sum(clust_total) - sum(clust_235) - (clust_total["R-4-C"] - clust_235["R-4-C"])), 
    byrow = TRUE, ncol = 2)
colnames(freq_chs) <- c("R-4-C", "!R-4-C")
rownames(freq_chs) <- c("235", "non235")
fisher.test(freq_chs, alternative = "greater")
##        Fisher's Exact Test for Count Data
## data:  freq_chs
## p-value < 2.2e-16
## alternative hypothesis: true odds ratio is greater than 1
## 95 percent confidence interval:
##     9.605755      Inf
## sample estimates:
##    odds ratio
## 13.33392
## --> we can reject the null hypothesis, NH = true odds ratio is equal to or smaller than 1


## A) foreground #1 vs background (PKS containing regions=union of the four methods vs background)
## 1. get all genes in genes_fg
genes_fg1 <- attr[, "gene_tandem"]
genes_fg1 <- unlist(strsplit(genes_fg1, split = "/"))

## 2. get concatenate of all syntenic genes per gene for the species analysed: 
## 2.1 find first which species are included and truncate fg_all and bg_all
species_fg1 <- unique(unlist(lapply(strsplit(genes_fg1, split = "_"), "[", 1)))
inds_A <- unlist(lapply(1:length(fg_all), function(x) {
    all(strsplit(names(fg_all), split = "_")[[x]] %in% species_fg1)
}))
fg_all_A <- fg_all[inds_A]
bg_all_A <- bg_all[inds_A]

## 2.3 retain only for those genes that are in foreground 1
fg_all_A <- lapply(fg_all_A, function(x) x[genes_fg1])

## 2.4 concatenate over species combinations the genes 
fg_all_A_uniq <- vector("list", length(fg_all_A[[1]]))
names(fg_all_A_uniq) <- names(fg_all_A[[1]])

for (genes in seq_along(fg_all_A[[1]])) {
    for (i in seq_along(fg_all_A)) {
        if (is.null(fg_all_A[[i]][[genes]])) {
            fg_all_A[[i]][[genes]] <- character()
        }
        fg_all_A_uniq[[genes]] <- c(
            fg_all_A_uniq[[genes]], fg_all_A[[i]][[genes]])
    }
    
    fg_all_A_uniq[[genes]] <- sort(unique(fg_all_A_uniq[[genes]]))
}

## 2.5 combine tandem replicates
fg_all_A_uniq <- sort(unique(unlist(fg_all_A_uniq)))

## 2.6 convert genes to GO terms
## for foreground
fg_all_A_uniq_go <- lapply(1:length(fg_all_A_uniq), function(x) { 
    species <- lapply(strsplit(fg_all_A_uniq[x], split = "_"), "[", 1)
    species <- paste(unlist(species), "_GO.out", sep = "")
    unlist(go_l[[species]][ fg_all_A_uniq[[x]] ])
})

## for background
bg_all_A_uniq <- unique(unlist(bg_all_A))
bg_all_A_uniq_go <- lapply(1:length(bg_all_A_uniq), function(x) { 
    species <- lapply(strsplit(bg_all_A_uniq[x], split = "_"), "[", 1)
    species <- paste(unlist(species), "_GO.out", sep = "")
    unlist(go_l[[species]][ bg_all_A_uniq[[x]] ])
})

## 2.7 do enrichment test
A_fisher <- fisher_test_all(bg = bg_all_A_uniq_go, fg = fg_all_A_uniq_go)
A_fisher_p <- unlist(lapply(A_fisher, function(x) x$p_values_adj))
A_fisher_p[A_fisher_p < 0.05]

## B) foreground #2 vs background (CHS containing cluster regions=union of the 
## four methods vs background)
## 1. get all genes in genes_fg, clusters that contain high amount of CHS are 2, 3, 5
genes_fg2 <- attr[attr[, "classification_community_apclust_dist"] %in% c(2, 3, 5), "gene_tandem"]
genes_fg2 <- unlist(strsplit(genes_fg2, split = "/"))

## 2. get concatenate of all syntenic genes per gene for the species analysed: 
## 2.1 find first which species are included and truncate fg_all and bg_all
species_fg2 <- unique(unlist(lapply(strsplit(genes_fg2, split = "_"), "[", 1)))
inds_B <- unlist(lapply(1:length(fg_all), function(x) {
    all(strsplit(names(fg_all), split="_")[[x]] %in% species_fg2)
}))
fg_all_B <- fg_all[inds_B]

## 2.3 retain only for those genes that are in foreground 1
fg_all_B <- lapply(fg_all_B, function(x) x[genes_fg2])

## 2.4 concatenate over species combinations the genes 
fg_all_B_uniq <- vector("list", length(fg_all_B[[1]]))
names(fg_all_B_uniq) <- names(fg_all_B[[1]])
for (genes in seq_along(fg_all_B[[1]])) {
    for (i in seq_along(fg_all_B)) {
        if (is.null(fg_all_B[[i]][[genes]])) {
            fg_all_B[[i]][[genes]] <- character()}
        fg_all_B_uniq[[genes]] <- c(
            fg_all_B_uniq[[genes]], fg_all_B[[i]][[genes]])
    }
    fg_all_B_uniq[[genes]] <- sort(unique(fg_all_B_uniq[[genes]]))
}

## 2.5 combine tandem replicates
fg_all_B_uniq <- sort(unique(unlist(fg_all_B_uniq)))

## 2.6 convert genes to GO terms
fg_all_B_uniq_go <- lapply(1:length(fg_all_B_uniq), function(x) { 
    species <- lapply(strsplit(fg_all_B_uniq[x], split = "_"), "[", 1)
    species <- paste(unlist(species), "_GO.out", sep = "")
    unlist(go_l[[species]][ fg_all_B_uniq[[x]] ])
})

## 2.7 do enrichment test (bg is identical with bg_all_A_uniq_go)
B_fisher <- fisher_test_all(bg = bg_all_A_uniq_go, fg = fg_all_B_uniq_go)
B_fisher_p <- unlist(lapply(B_fisher, function(x) x$p_values_adj))
B_fisher_p[B_fisher_p < 0.05]

## C) foreground #2 vs foreground #1
C_fisher <- fisher_test_all(bg = fg_all_A_uniq_go, fg = fg_all_B_uniq_go)
C_fisher_p <- unlist(lapply(C_fisher, function(x) x$p_values_adj))
C_fisher_p[C_fisher_p < 0.05]

## create gmt file
## GO.id description
gmt <- do.call(rbind, go_l_desc)
gmt <- gmt[!duplicated(gmt[, "goid"]), ]
gmt <- gmt[, c("goid", "desc")]

## get relation: go terms --> gene list
names(fg_all_A_uniq_go) <- fg_all_A_uniq
names(fg_all_B_uniq_go) <- fg_all_B_uniq

## function to find from a go list corresponding genes
reverse_go <- function(go, term) {
    log <- lapply(go, function(x) unlist(any(x %in% term)))
    names(which(unlist(log)))
}

## create files per case
## A
go_A <- matrix("", nrow = length(A_fisher), ncol = 6)
colnames(go_A) <- c("GO.ID", "Description", "p.Val", "FDR", 
    "Phenotype", "gene.list")
go_A[, "GO.ID"] <- names(A_fisher_p)
go_A[, "Description"] <- gmt[match(go_A[, "GO.ID"], gmt[, "goid"]), "desc"]
go_A[, "p.Val"] <- unlist(lapply(A_fisher, function(x) x$p.value))
go_A[, "FDR"] <- unlist(lapply(A_fisher, function(x) x$p_values_adj))
go_A[, "Phenotype"] <- "+1" ## arbitrary
for (i in 1:nrow(go_A)) {
    go_A[i, "gene.list"] <- paste(
        reverse_go(fg_all_A_uniq_go, go_A[i, "GO.ID"]), collapse=",")
}

## B
go_B <- matrix("", nrow = length(B_fisher), ncol = 6)
colnames(go_B) <- c("GO.ID", "Description", "p.Val", "FDR", 
    "Phenotype", "gene.list")
go_B[, "GO.ID"] <- names(B_fisher_p)
go_B[, "Description"] <- gmt[match(go_B[, "GO.ID"], gmt[, "goid"]), "desc"]
go_B[, "p.Val"] <- unlist(lapply(B_fisher, function(x) x$p.value))
go_B[, "FDR"] <- unlist(lapply(B_fisher, function(x) x$p_values_adj))
go_B[, "Phenotype"] <- "+1" ## arbitrary
for (i in 1:nrow(go_B)) {
    go_B[i, "gene.list"] <- paste(
        reverse_go(fg_all_B_uniq_go, go_B[i, "GO.ID"]), collapse = ",")
}

## B
go_C <- matrix("", nrow = length(C_fisher), ncol = 6)
colnames(go_C) <- c("GO.ID", "Description", "p.Val", "FDR", 
    "Phenotype", "gene.list")
go_C[, "GO.ID"] <- names(C_fisher_p)
go_C[, "Description"] <- gmt[match(go_C[, "GO.ID"], gmt[, "goid"]), "desc"]
go_C[, "p.Val"] <- unlist(lapply(C_fisher, function(x) x$p.value))
go_C[, "FDR"] <- unlist(lapply(C_fisher, function(x) x$p_values_adj))
go_C[, "Phenotype"] <- "+1" ## arbitrary
for (i in 1:nrow(go_C)) {
    go_C[i, "gene.list"] <- paste(
        reverse_go(fg_all_B_uniq_go, go_C[i, "GO.ID"]), collapse = ",")
}

## retrieve for Molecular Function, Cellular compartment, Biological Process
## the terms
## A
go_A_MF <- go_A[grep(go_A[, "GO.ID"], pattern = "MF"), ]
go_A_CC <- go_A[grep(go_A[, "GO.ID"], pattern = "CC"), ]
go_A_BP <- go_A[grep(go_A[, "GO.ID"], pattern = "BP"), ]

## B
go_B_MF <- go_B[grep(go_B[, "GO.ID"], pattern = "MF"), ]
go_B_CC <- go_B[grep(go_B[, "GO.ID"], pattern = "CC"), ]
go_B_BP <- go_B[grep(go_B[, "GO.ID"], pattern = "BP"), ]

## C
go_C_MF <- go_C[grep(go_C[, "GO.ID"], pattern = "MF"), ]
go_C_CC <- go_C[grep(go_C[, "GO.ID"], pattern = "CC"), ]
go_C_BP <- go_C[grep(go_C[, "GO.ID"], pattern = "BP"), ]

## BP
df_BP <- data.frame(term = go_C_BP[, "Description"], 
    pvalue = as.numeric(go_C_BP[, "FDR"]))
df_BP <- df_BP[df_BP$pvalue < 0.05, ]
df_BP$pvalue <- -log10(df_BP$pvalue)
df_BP$type <- "BP"
df_BP <- within(df_BP, 
    term <- factor(term, levels = term[order(pvalue, decreasing = TRUE)]))

## MF
df_MF <- data.frame(term = go_C_MF[, "Description"], 
    pvalue = as.numeric(go_C_MF[, "FDR"]))
df_MF <- df_MF[df_MF$pvalue < 0.05, ]
df_MF$pvalue <- -log10(df_MF$pvalue)
df_MF$type <- "MF"
df_MF <- within(df_MF, 
    term <- factor(term, levels = term[order(pvalue, decreasing = TRUE)]))

## CC
df_CC <- data.frame(term = go_C_CC[, "Description"], 
    pvalue = as.numeric(go_C_CC[, "FDR"]))
df_CC <- df_CC[df_CC$pvalue < 0.05, ]
df_CC$pvalue <- -log10(df_CC$pvalue)
df_CC <- within(df_CC, 
    term <- factor(term, levels = term[order(pvalue, decreasing = TRUE)]))
df_CC$type <- "CC"

## dim(df_BP) 186 3
## dim(df_MF) 76 3 
## dim(df_CC) 34 3
## plot top 20 terms
df_all <- rbind(df_BP[order(df_BP[, 2], decreasing = TRUE)[1:20], ], 
    df_MF[order(df_MF[, 2], decreasing = TRUE)[1:20], ], 
    df_CC[order(df_CC[, 2], decreasing = TRUE)[1:20], ])
g <- ggplot(df_all) + 
    geom_bar(aes(x = term, y = pvalue, fill = type), stat = "identity", 
        width = 0.85) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0), 
        axis.text.y = element_text(angle = 90))  
ggsave(plot = g, device = "pdf", 
    filename = "figure_S2_GO_enrichment_terms/barplot_all_C.pdf")

## plot all terms for BP, MF and CC
## BP
g <- ggplot(df_BP) + 
    geom_bar(aes(x = term, y = pvalue), stat = "identity", width = 0.85) +
    ylim(0, 30) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0), 
        axis.text.y = element_text(angle = 90))  
ggsave(plot = g, device = "pdf", 
    filename = "figure_S2_GO_enrichment_terms/barplot_BP_C.pdf")

## MF
g <- ggplot(df_MF) + 
    geom_bar(aes(x = term, y = pvalue), stat = "identity", width = 0.85) + 
    ylim(0, 30) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0), 
        axis.text.y = element_text(angle = 90))
ggsave(plot = g, device = "pdf", 
    filename = "figure_S2_GO_enrichment_terms/barplot_MF_C.pdf")

## CC 
g <- ggplot(df_CC) + 
    geom_bar(aes(x = term, y = pvalue), stat = "identity", width = 0.85) + 
    ylim(0, 30) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0), 
        axis.text.y = element_text(angle = 90))
ggsave(plot = g, device = "pdf", 
    filename = "figure_S2_GO_enrichment_terms/barplot_CC_C.pdf")

## plot gene term frequencies per term
freq_terms_C <- unlist(lapply(strsplit(
    go_C[match(as.character(df_all[, 1]), go_C[, 2]), 6], split=","), length))
freq_terms_A <- unlist(lapply(strsplit(
    go_A[match(as.character(df_all[, 1]), go_A[, 2]), 6], split=","), length))

num_gene_C <- length(unique(unlist(strsplit(go_C[, 6], split = ","))))
num_gene_A <- length(unique(unlist(strsplit(go_A[, 6], split = ","))))

df_freq_C <- data.frame(name = as.character(df_all[, 1]), 
    freq = freq_terms_C / num_gene_C * 100, is = "C")
df_freq_A <- data.frame(name = as.character(df_all[, 1]), 
    freq = freq_terms_A / num_gene_A * 100, is = "A")
df_freq <- rbind(df_freq_A, df_freq_C)
df_freq <- within(df_freq, name <- factor(name, levels = as.character(df_freq_C$name)))

g <- ggplot(df_freq) + 
    geom_bar(aes(x = name, y = freq, fill = is), stat = "identity", 
        width = 0.85, position = position_dodge(width = 1)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0), 
        axis.text.y = element_text(angle = 90)) + ylim(0, 2)
ggsave(plot = g, device = "pdf", 
    filename = "figure_S2_GO_enrichment_terms/barplot_frequency_C_freq_greater2.pdf")

## write objects to files 
write.table(gmt, file = "go_gmt.gmt", sep = "\t", dec = ".", 
    row.names = FALSE, quote = FALSE)

## terms for A
write.table(go_A, file = "go_A.txt", sep = "\t", dec = ".", 
    row.names = FALSE, quote = FALSE)
write.table(go_A_MF, file = "go_A_MF.txt", sep = "\t", dec = ".", 
    row.names = FALSE, quote = FALSE)
write.table(go_A_CC, file = "go_A_CC.txt", sep = "\t", dec = ".", 
    row.names = FALSE, quote = FALSE)
write.table(go_A_BP, file = "go_A_BP.txt", sep = "\t", dec = ".", 
    row.names = FALSE, quote = FALSE) 

## terms for B
write.table(go_B, file = "go_B.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)
write.table(go_B_MF, file = "go_B_MF.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)
write.table(go_B_CC, file = "go_B_CC.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)
write.table(go_B_BP, file = "go_B_BP.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE) 

## terms for C
write.table(go_C, file = "go_C.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)
write.table(go_C_MF, file = "go_C_MF.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)
write.table(go_C_CC, file = "go_C_CC.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)
write.table(go_C_BP, file = "go_C_BP.txt", sep = "\t", dec = ".", row.names = FALSE, quote = FALSE)

###############################################################################
###############################################################################
## do the same as above but remove the pks genes from the gene set

## A) foreground #1 vs background
## foreground
fg_all_A_uniq_r <- fg_all_A_uniq[! fg_all_A_uniq %in% pks_genes_all]
fg_all_A_uniq_go_r <- lapply(1:length(fg_all_A_uniq_r), function(x) { 
    species <- lapply(strsplit(fg_all_A_uniq_r[x], split = "_"), "[", 1)
    species <- paste(unlist(species), "_GO.out", sep = "")
    unlist(go_l[[species]][ fg_all_A_uniq_r[[x]] ])
})

## background 
bg_all_A_uniq_r <- bg_all_A_uniq[! bg_all_A_uniq %in% pks_genes_all ]
bg_all_A_uniq_go_r <- lapply(1:length(bg_all_A_uniq_r), function(x) { 
    species <- lapply(strsplit(bg_all_A_uniq_r[x], split = "_"), "[", 1)
    species <- paste(unlist(species), "_GO.out", sep = "")
    unlist(go_l[[species]][ bg_all_A_uniq_r[[x]] ])
})
## 2.7 do enrichment test
A_fisher_r <- fisher_test_all(bg = bg_all_A_uniq_go_r, fg = fg_all_A_uniq_go_r)
A_fisher_p_r <- unlist(lapply(A_fisher_r, function(x) x$p_values_adj))
A_fisher_p_r[A_fisher_p_r < 0.05]

## B) foreground #2 vs background (CHS containing cluster regions=union of 
## the four methods vs background)
## foreground
fg_all_B_uniq_r <- fg_all_B_uniq[!fg_all_B_uniq %in% pks_genes_all]
fg_all_B_uniq_go_r <- lapply(1:length(fg_all_B_uniq_r), function(x) { 
    species <- lapply(strsplit(fg_all_B_uniq_r[x], split = "_"), "[", 1)
    species <- paste(unlist(species), "_GO.out", sep = "")
    unlist(go_l[[species]][ fg_all_B_uniq_r[[x]] ])
})
## background 
B_fisher_r <- fisher_test_all(bg = bg_all_A_uniq_go_r, fg = fg_all_B_uniq_go_r)
B_fisher_p_r <- unlist(lapply(B_fisher_r, function(x) x$p_values_adj))
B_fisher_p_r[B_fisher_p_r < 0.05]

## C) foreground #2 vs foreground #1
C_fisher_r <- fisher_test_all(bg = fg_all_A_uniq_go_r, fg = fg_all_B_uniq_go_r)
C_fisher_p_r <- unlist(lapply(C_fisher_r, function(x) x$p_values_adj))
C_fisher_p_r[C_fisher_p_r < 0.05]

## get relation: go terms --> gene list
names(fg_all_A_uniq_go_r) <- fg_all_A_uniq_r
names(fg_all_B_uniq_go_r) <- fg_all_B_uniq_r


## create files per case
## A
go_A_r <- matrix("", nrow = length(A_fisher_r), ncol = 6)
colnames(go_A_r) <- c("GO.ID", "Description", "p.Val", "FDR", 
    "Phenotype", "gene.list")
go_A_r[, "GO.ID"] <- names(A_fisher_p_r)
go_A_r[, "Description"] <- gmt[match(go_A_r[, "GO.ID"], gmt[, "goid"]), "desc"]
go_A_r[, "p.Val"] <- unlist(lapply(A_fisher_r, function(x) x$p.value))
go_A_r[, "FDR"] <- unlist(lapply(A_fisher_r, function(x) x$p_values_adj))
go_A_r[, "Phenotype"] <- "+1" ## arbitrary
for (i in 1:nrow(go_A_r)) {
    go_A_r[i, "gene.list"] <- paste(
        reverse_go(fg_all_A_uniq_go_r, go_A_r[i, "GO.ID"]), collapse = ",")
}

## B
go_B_r <- matrix("", nrow = length(B_fisher_r), ncol = 6)
colnames(go_B_r) <- c("GO.ID", "Description", "p.Val", "FDR", 
    "Phenotype", "gene.list")
go_B_r[, "GO.ID"] <- names(B_fisher_p_r)
go_B_r[, "Description"] <- gmt[match(go_B_r[, "GO.ID"], gmt[, "goid"]), "desc"]
go_B_r[, "p.Val"] <- unlist(lapply(B_fisher_r, function(x) x$p.value))
go_B_r[, "FDR"] <- unlist(lapply(B_fisher_r, function(x) x$p_values_adj))
go_B_r[, "Phenotype"] <- "+1"
for (i in 1:nrow(go_B_r)) {
    go_B_r[i, "gene.list"] <- paste(
        reverse_go(fg_all_B_uniq_go_r, go_B_r[i, "GO.ID"]), collapse = ",")
}

## C
go_C_r <- matrix("", nrow = length(C_fisher_r), ncol = 6)
colnames(go_C_r) <- c("GO.ID", "Description", "p.Val", "FDR", 
    "Phenotype", "gene.list")
go_C_r[, "GO.ID"] <- names(C_fisher_p_r)
go_C_r[, "Description"] <- gmt[match(go_C_r[, "GO.ID"], gmt[, "goid"]), "desc"]
go_C_r[, "p.Val"] <- unlist(lapply(C_fisher_r, function(x) x$p.value))
go_C_r[, "FDR"] <- unlist(lapply(C_fisher_r, function(x) x$p_values_adj))
go_C_r[, "Phenotype"] <- "+1"
for (i in 1:nrow(go_C_r)) {
    go_C_r[i, "gene.list"] <- paste(
        reverse_go(fg_all_B_uniq_go_r, go_C_r[i, "GO.ID"]), collapse = ",")
}

## retrieve for Molecular Function, Cellular compartment, Biological Process
## the terms
## A
go_A_MF_r <- go_A_r[grep(go_A_r[, "GO.ID"], pattern = "MF"), ]
go_A_CC_r <- go_A_r[grep(go_A_r[, "GO.ID"], pattern = "CC"), ]
go_A_BP_r <- go_A_r[grep(go_A_r[, "GO.ID"], pattern = "BP"), ]

## B
go_B_MF_r <- go_B_r[grep(go_B_r[, "GO.ID"], pattern = "MF"), ]
go_B_CC_r <- go_B_r[grep(go_B_r[, "GO.ID"], pattern = "CC"), ]
go_B_BP_r <- go_B_r[grep(go_B_r[, "GO.ID"], pattern = "BP"), ]

## C
go_C_MF_r <- go_C_r[grep(go_C_r[, "GO.ID"], pattern = "MF"), ]
go_C_CC_r <- go_C_r[grep(go_C_r[, "GO.ID"], pattern = "CC"), ]
go_C_BP_r <- go_C_r[grep(go_C_r[, "GO.ID"], pattern = "BP"), ]

## BP
df_BP_r <- data.frame(term = go_C_BP_r[, "Description"], 
    pvalue = as.numeric(go_C_BP_r[, "FDR"]))
df_BP_r <- df_BP_r[df_BP_r$pvalue < 0.05, ]
df_BP_r$pvalue <- -log10(df_BP_r$pvalue)
df_BP_r$type <- "BP"
df_BP_r <- within(df_BP_r, 
    term <- factor(term, levels = term[order(pvalue, decreasing = TRUE)]))

## MF
df_MF_r <- data.frame(term = go_C_MF_r[, "Description"], 
    pvalue = as.numeric(go_C_MF_r[, "FDR"]))
df_MF_r <- df_MF_r[df_MF_r$pvalue < 0.05, ]
df_MF_r$pvalue <- -log10(df_MF_r$pvalue)
df_MF_r$type <- "MF"
df_MF_r <- within(df_MF_r, 
    term <- factor(term, levels = term[order(pvalue, decreasing = TRUE)]))

## C
df_CC_r <- data.frame(term = go_C_CC_r[, "Description"], 
    pvalue = as.numeric(go_C_CC_r[, "FDR"]))
df_CC_r <- df_CC_r[df_CC_r$pvalue < 0.05, ]
df_CC_r$pvalue <- -log10(df_CC_r$pvalue)
df_CC_r <- within(df_CC_r, 
    term <- factor(term, levels = term[order(pvalue, decreasing = TRUE)]))
df_CC_r$type <- "CC"

## plot top 20 terms 
df_all_r <- rbind(df_BP_r[order(df_BP_r[,2], decreasing = TRUE)[1:20], ], 
    df_MF_r[order(df_MF_r[,2], decreasing = TRUE)[1:20], ], 
    df_CC_r[order(df_CC_r[,2], decreasing = TRUE)[1:20], ])
g <- ggplot(df_all_r) +
    geom_bar(aes(x = term, y = pvalue, fill = type), stat = "identity", 
        width = 0.85) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0), 
         axis.text.y = element_text(angle = 90))
ggsave(plot = g, filename = "barplot_all_C_removedPKS.pdf", device = "pdf")

## plot gene term frequencies per term
freq_terms_C_r <- unlist(lapply(strsplit(
    go_C_r[match(as.character(df_all_r[,1]), go_C_r[,2]), 6], split = ","), length))
freq_terms_A_r <- unlist(lapply(strsplit(
    go_A_r[match(as.character(df_all_r[,1]), go_A_r[,2]), 6], split = ","), length))

num_gene_C_r <- length(unique(unlist(strsplit(go_C_r[, 6], split = ","))))
num_gene_A_r <- length(unique(unlist(strsplit(go_A_r[, 6], split = ","))))

df_freq_C_r <- data.frame(name = as.character(df_all_r[, 1]), 
    freq = freq_terms_C_r/num_gene_C_r * 100, is = "C")
df_freq_A_r <- data.frame(name = as.character(df_all_r[, 1]), 
    freq = freq_terms_A_r/num_gene_A_r * 100, is = "A")
df_freq_r <- rbind(df_freq_A_r, df_freq_C_r)
df_freq_r <- within(df_freq_r, name <- factor(name, levels = as.character(df_freq_C_r$name)))

g <- ggplot(df_freq_r) + 
    geom_bar(aes(x = name, y = freq, fill = is), stat = "identity", 
        width = 0.85, position = position_dodge(width = 1)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0, size = 5), 
        axis.text.y = element_text(angle = 90)) + ylim(0, 2)
ggsave(plot = g, device = "pdf", 
    filename = "barplot_frequency_C_freq_removedPKS.pdf")


## plot all terms 
## BP
g <- ggplot(df_BP_r) + 
    geom_bar(aes(x = term, y = pvalue), stat = "identity", width = 0.85) + 
    ylim(0,30) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0), 
        axis.text.y = element_text(angle = 90))  
ggsave(plot = g, filename = "barplot_BP_C_removedPKS.pdf", device = "pdf")
## MF
g <- ggplot(df_MF_r) + 
    geom_bar(aes(x = term, y = pvalue), stat = "identity", width = 0.85) + 
    ylim(0,30) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0, size = 5), 
          axis.text.y = element_text(angle = 90))  
ggsave(plot = g, filename = "barplot_MF_C_removedPKS.pdf", device = "pdf")
## CC
g <- ggplot(df_CC_r) + 
    geom_bar(aes(x = term, y = pvalue), stat = "identity", width = 0.85) + 
    ylim(0,30) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0), 
        axis.text.y = element_text(angle = 90))  
ggsave(plot = g, filename = "barplot_CC_C_removedPKS.pdf", device = "pdf")

## plot gene term frequencies per term
freq_terms_C_r <- unlist(lapply(strsplit(
    go_C_r[match(levels(df_all_r[, 1]), go_C_r[, 2]), 6], split = ","), length))
freq_terms_A_r <- unlist(lapply(strsplit(
    go_A_r[match(levels(df_all_r[, 1]), go_A_r[, 2]), 6], split = ","), length))

num_gene_C_r <- length(unique(unlist(strsplit(go_C_r[, 6], split = ","))))
num_gene_A_r <- length(unique(unlist(strsplit(go_A_r[, 6], split = ","))))

df_freq_C_r <- data.frame(name = levels(df_all_r[, 1]), 
    freq = freq_terms_C_r/num_gene_C_r * 100, is = "C")
df_freq_A_r <- data.frame(name = levels(df_all_r[, 1]), 
    freq = freq_terms_A_r/num_gene_A_r * 100, is = "A")
df_freq_r <- rbind(df_freq_A_r, df_freq_C_r)
df_freq_r <- within(df_freq_r, 
    name <- factor(name, levels = as.character(df_freq_C_r$name)))

g <- ggplot(df_freq_r) + 
    geom_bar(aes(x = name, y = freq, fill = is), stat = "identity", 
        width = 0.85, position = position_dodge(width = 1)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0, size = 5), 
        axis.text.y = element_text(angle = 90)) + ylim(0, 14)
ggsave(plot = g, device = "pdf", 
    filename = "barplot_frequency_C_freq_smaller_removedPKS.pdf")

## write objects to files 
## for A
write.table(go_A_r, file = "go_A_PKSremoved.txt", sep = "\t", dec = ".", 
    row.names = FALSE, quote = FALSE)
write.table(go_A_MF_r, file = "go_A_MF_PKSremoved.txt", sep = "\t", dec = ".", 
    row.names = FALSE, quote = FALSE)
write.table(go_A_CC_r, file = "go_A_CC_PKSremoved.txt", sep = "\t", dec = ".", 
    row.names = FALSE, quote = FALSE)
write.table(go_A_BP_r, file = "go_A_BP_PKSremoved.txt", sep = "\t", dec = ".", 
    row.names = FALSE, quote = FALSE)

## for B
write.table(go_B_r, file = "go_B_PKSremoved.txt", sep = "\t", dec = ".", 
    row.names = FALSE, quote = FALSE)
write.table(go_B_MF_r, file = "go_B_MF_PKSremoved.txt", sep = "\t", dec = ".", 
    row.names = FALSE, quote = FALSE)
write.table(go_B_CC_r, file = "go_B_CC_PKSremoved.txt", sep = "\t", dec = ".", 
    row.names = FALSE, quote = FALSE)
write.table(go_B_BP_r, file = "go_B_BP_PKSremoved.txt", sep = "\t", dec = ".", 
    row.names = FALSE, quote = FALSE)

## for C
write.table(go_C_r, file = "go_C_PKSremoved.txt", sep = "\t", dec = ".", 
    row.names = FALSE, quote = FALSE)
write.table(go_C_MF_r, file = "go_C_MF_PKSremoved.txt", sep = "\t", dec = ".", 
    row.names = FALSE, quote = FALSE)
write.table(go_C_CC_r, file = "go_C_CC_PKSremoved.txt", sep = "\t", dec = ".", 
    row.names = FALSE, quote = FALSE)
write.table(go_C_BP_r, file = "go_C_BP_PKSremoved.txt", sep = "\t", dec = ".", 
    row.names = FALSE, quote = FALSE)


################################################################################
################################################################################
## use STRING to get co-expression within gene region 
## use arath, zeama, orysa and solyc for this kind of analysis
## write the gene list containing the regions to a file
#### for clusters 2, 3 and 5 ####
## arath_AT5G13930
at5g13930 <- fg_all_B_uniq[grep(fg_all_B_uniq, pattern = "arath_AT5G")]
## truncate the region
at5g13930 <- at5g13930[(grep(at5g13930, pattern = "AT5G13930")-399):(grep(at5g13930, pattern = "AT5G13930")+399)]
## zeama_Zm00001d007400/zeama_Zm00001d007403 
zm00001d007403 <- fg_all_B_uniq[grep(fg_all_B_uniq, 
    pattern = "zeama_Zm00001d007[2|3|4|5]")]
## zeama_Zm00001d052673/zeama_Zm00001d052675/zeama_Zm00001d052676
zm00001d052673 <- fg_all_B_uniq[grep(fg_all_B_uniq, 
    pattern = "zeama_Zm00001d052[5|6|7]")]
## orysa_LOC_Os05g12180/orysa_LOC_Os05g12190/orysa_LOC_Os05g12210/orysa_LOC_Os05g12240 (no CHS)
os05g12180 <- fg_all_B_uniq[grep(fg_all_B_uniq, pattern = "orysa_LOC_Os05g")]
## orysa_LOC_Os11g32540/orysa_LOC_Os11g32550/orysa_LOC_Os11g32580/orysa_LOC_Os11g32610/orysa_LOC_Os11g32620/orysa_LOC_Os11g32650
os11g32650 <- fg_all_B_uniq[grep(fg_all_B_uniq, pattern = "orysa_LOC_Os11g")]
## solyc_Solyc05g053550.3.1
solyc05g053550 <- fg_all_B_uniq[grep(fg_all_B_uniq, pattern = "solyc_Solyc05g")]
## solyc_Solyc09g091510.3.1
solyc09g091510 <- fg_all_B_uniq[grep(fg_all_B_uniq, pattern = "solyc_Solyc09g")]
## solyc_Solyc12g098100.2.1 (no CHS)
solyc12g098100 <- fg_all_B_uniq[grep(fg_all_B_uniq, pattern = "solyc_Solyc12g")]
## vitvi_GSVIVT01032968001
gsvivt01032968001 <- fg_all_B_uniq[grep(fg_all_B_uniq, 
    pattern = "vitvi_GSVIVT0103")]

#### for other clusters #####
##arath_AT1G02050
at1g02050 <- fg_all_A_uniq[grep(fg_all_A_uniq, pattern = "arath_AT1G")]
##arath_AT4G00040
at4g00040 <- fg_all_A_uniq[grep(fg_all_A_uniq, pattern = "arath_AT4G0")]
##arath_AT4G34850
at4g34850 <- fg_all_A_uniq[grep(fg_all_A_uniq, pattern = "arath_AT4G3")]
##solyc_Solyc01g090600.3.1
solyc01g090600 <- fg_all_A_uniq[grep(fg_all_A_uniq, pattern = "solyc_Solyc01g")]
solyc01g090600 <- solyc01g090600[(grep(solyc01g090600, pattern = "01g090600")-399):(grep(solyc01g090600, pattern = "01g090600")+399)]
##solyc_Solyc01g111070.3.1
solyc01g111070 <- fg_all_A_uniq[grep(fg_all_A_uniq, pattern = "solyc_Solyc01g")]
solyc01g111070 <- solyc01g111070[(grep(solyc01g111070, pattern = "01g111070")-399):length(solyc01g111070)]
##solyc_Solyc05g053170.3.1/solyc_Solyc05g053550.3.1
solyc05g053170 <- fg_all_A_uniq[grep(fg_all_A_uniq, pattern = "solyc_Solyc05g")]
##solyc_Solyc06g043120.2.1
solyc06g043120 <- fg_all_A_uniq[grep(fg_all_A_uniq, pattern = "solyc_Solyc06g")]
##vitvi_GSVIVT01000521001
gsvivt01000521001 <- fg_all_A_uniq[grep(fg_all_A_uniq, 
    pattern = "vitvi_GSVIVT0100")]
##vitvi_GSVIVT01010554001/vitvi_GSVIVT01010556001/vitvi_GSVIVT01010557001/vitvi_GSVIVT01010561001/vitvi_GSVIVT01010563001/vitvi_GSVIVT01010565001/
##vitvi_GSVIVT01010568001/vitvi_GSVIVT01010570001/vitvi_GSVIVT01010572001/vitvi_GSVIVT01010574001/vitvi_GSVIVT01010578001/vitvi_GSVIVT01010579001/
##vitvi_GSVIVT01010580001/vitvi_GSVIVT01010581001/vitvi_GSVIVT01010582001/vitvi_GSVIVT01010583001/vitvi_GSVIVT01010584001/vitvi_GSVIVT01010585001/
##vitvi_GSVIVT01010589001/vitvi_GSVIVT01010590001 (STS)
gsvivt01010554001 <- fg_all_A_uniq[grep(fg_all_A_uniq, 
    pattern = "vitvi_GSVIVT01010")]
##vitvi_GSVIVT01018219001
gsvivt01018219001 <- fg_all_A_uniq[grep(fg_all_A_uniq, 
    pattern = "vitvi_GSVIVT0101[7|8]")]
##vitvi_GSVIVT01024107001
gsvivt01024107001 <- fg_all_A_uniq[grep(fg_all_A_uniq, 
    pattern = "vitvi_GSVIVT0102[3|4]")]
##vitvi_GSVIVT01026213001/vitvi_GSVIVT01026220001
gsvivt01026213001 <- fg_all_A_uniq[grep(fg_all_A_uniq, 
    pattern = "vitvi_GSVIVT0102[6|7|8|9]")]
##zeama_Zm00001d006445
zm00001d000438 <- fg_all_A_uniq[grep(fg_all_A_uniq, 
    pattern = "zeama_Zm00001d00[5|6]")]
##zeama_Zm00001d007717
Zm00001d007717 <- fg_all_A_uniq[grep(fg_all_A_uniq, 
    pattern = "zeama_Zm00001d007[5|6|7|8]")]
##zeama_Zm00001d013991
Zm00001d013991 <- fg_all_A_uniq[grep(fg_all_A_uniq, 
    pattern = "zeama_Zm00001d013")]
##zeama_Zm00001d019478
Zm00001d019478 <- fg_all_A_uniq[grep(fg_all_A_uniq, 
    pattern = "zeama_Zm00001d019")]
##zeama_Zm00001d021562
Zm00001d021562 <- fg_all_A_uniq[grep(fg_all_A_uniq, 
    pattern = "zeama_Zm00001d02")]
##zeama_Zm00001d032662
Zm00001d032662 <- fg_all_A_uniq[grep(fg_all_A_uniq, 
    pattern = "zeama_Zm00001d03")]
##zeama_Zm00001d040479/zeama_Zm00001d040483/zeama_Zm00001d040497
Zm00001d040479 <- fg_all_A_uniq[grep(fg_all_A_uniq, 
    pattern = "zeama_Zm00001d04")]
## zeama_Zm00001d052915/zeama_Zm00001d052916
Zm00001d052915 <- fg_all_A_uniq[grep(fg_all_A_uniq, 
    pattern = "zeama_Zm00001d0529")]
##orysa_LOC_Os01g41834
os01g41834 <- fg_all_A_uniq[grep(fg_all_A_uniq, pattern = "Os01g")]
##orysa_LOC_Os03g47000
os03g47000 <- fg_all_A_uniq[grep(fg_all_A_uniq, pattern = "Os03g")]
##orysa_LOC_Os04g23940
os04g23940 <- fg_all_A_uniq[grep(fg_all_A_uniq, pattern = "Os04g")]
##orysa_LOC_Os05g41645
os05g41645 <- fg_all_A_uniq[grep(fg_all_A_uniq, pattern = "Os05g4")]
##orysa_LOC_Os07g11440 + orysa_LOC_Os07g17010 + orysa_LOC_Os07g22850 + orysa_LOC_Os07g31750/orysa_LOC_Os07g31770 + orysa_LOC_Os07g34140/orysa_LOC_Os07g34190/orysa_LOC_Os07g34260
os07g11440 <- fg_all_A_uniq[grep(fg_all_A_uniq, pattern = "Os07g")]
##orysa_LOC_Os10g07616 + orysa_LOC_Os10g08620/orysa_LOC_Os10g08670/orysa_LOC_Os10g08680/orysa_LOC_Os10g08710 + orysa_LOC_Os10g09860
os10g07616 <- fg_all_A_uniq[grep(fg_all_A_uniq, pattern = "Os10g[0|1]")]
##orysa_LOC_Os10g34360 + orysa_LOC_Os10g36972
os10g34360 <- fg_all_A_uniq[grep(fg_all_A_uniq, pattern = "Os10g[2|3]")]
##orysa_LOC_Os11g35930
os11g35930 <- fg_all_A_uniq[grep(fg_all_A_uniq, pattern = "Os11g3")]
##orysa_LOC_Os12g07690
os12g07690 <- fg_all_A_uniq[grep(fg_all_A_uniq, pattern = "Os12g")]
## end other clusters ## 


## strsplit (cluster 2, 3, 5)
## A. thaliana
at5g13930 <- unlist(lapply(strsplit(at5g13930, split = "arath_"), "[", 2))

## Z. mays
zm00001d007403 <- unlist(lapply(
    strsplit(zm00001d007403, split = "zeama_"), "[", 2))
zm00001d052673 <- unlist(lapply(
    strsplit(zm00001d052673, split = "zeama_"), "[", 2))

## O. sativa
os05g12180 <- unlist(lapply(
    strsplit(os05g12180, split = "orysa_"), "[", 2))
os11g32650 <- unlist(lapply(
    strsplit(os11g32650, split = "orysa_"), "[", 2))

## S. lycopersicum
solyc05g053550 <- unlist(lapply(
    strsplit(solyc05g053550, split = "solyc_"), "[", 2))
solyc09g091510 <- unlist(lapply(
    strsplit(solyc09g091510, split = "solyc_"), "[", 2))
solyc12g098100 <- unlist(lapply(
    strsplit(solyc12g098100, split = "solyc_"), "[", 2))
gsvivt01032968001 <- unlist(lapply(
    strsplit(gsvivt01032968001, split = "vitvi_"), "[", 2))

## strsplit (other clusters)
## A. thaliana
at1g02050 <- unlist(lapply(strsplit(at1g02050, split = "arath_"), "[", 2))
at4g00040 <- unlist(lapply(strsplit(at4g00040, split = "arath_"), "[", 2))
at4g34850 <- unlist(lapply(strsplit(at4g34850, split = "arath_"), "[", 2))

## S. lycopersicum
solyc01g090600 <- unlist(lapply(
    strsplit(solyc01g090600, split = "solyc_"), "[", 2))
solyc01g111070 <- unlist(lapply(
    strsplit(solyc01g111070, split = "solyc_"), "[", 2))
solyc05g053170 <- unlist(lapply(
    strsplit(solyc05g053170, split = "solyc_"), "[", 2))
solyc06g043120 <- unlist(lapply(
    strsplit(solyc06g043120, split = "solyc_"), "[", 2))

## V. vinifera
gsvivt01000521001 <- unlist(lapply(
    strsplit(gsvivt01000521001, split = "vitvi_"), "[", 2))
gsvivt01010554001 <- unlist(lapply(
    strsplit(gsvivt01010554001, split = "vitvi_"), "[", 2))
gsvivt01018219001 <- unlist(lapply(
    strsplit(gsvivt01018219001, split = "vitvi_"), "[", 2))
gsvivt01024107001 <- unlist(lapply(
    strsplit(gsvivt01024107001, split = "vitvi_"), "[", 2))
gsvivt01026213001 <- unlist(lapply(
    strsplit(gsvivt01026213001, split = "vitvi_"), "[", 2))

## Z. mays 
zm00001d000438 <- unlist(lapply(
    strsplit(zm00001d000438, split = "zeama_"), "[", 2))
Zm00001d007717 <- unlist(lapply(
    strsplit(Zm00001d007717, split = "zeama_"), "[", 2))
Zm00001d013991 <- unlist(lapply(
    strsplit(Zm00001d013991, split = "zeama_"), "[", 2))
Zm00001d019478 <- unlist(lapply(
    strsplit(Zm00001d019478, split = "zeama_"), "[", 2))
Zm00001d021562 <- unlist(lapply(
    strsplit(Zm00001d021562, split = "zeama_"), "[", 2))
Zm00001d032662 <- unlist(lapply(
    strsplit(Zm00001d032662, split = "zeama_"), "[", 2))
Zm00001d040479 <- unlist(lapply(
    strsplit(Zm00001d040479, split = "zeama_"), "[", 2))
Zm00001d052915 <- unlist(lapply(
    strsplit(Zm00001d052915, split = "zeama_"), "[", 2))
os01g41834 <- unlist(lapply(strsplit(os01g41834, split = "orysa_"), "[", 2))

## O. sativa
os03g47000 <- unlist(lapply(strsplit(os03g47000, split = "orysa_"), "[", 2))
os04g23940 <- unlist(lapply(strsplit(os04g23940, split = "orysa_"), "[", 2))
os05g41645 <- unlist(lapply(strsplit(os05g41645, split = "orysa_"), "[", 2))
os07g11440 <- unlist(lapply(strsplit(os07g11440, split = "orysa_"), "[", 2))
os10g07616 <- unlist(lapply(strsplit(os10g07616, split = "orysa_"), "[", 2))
os10g34360 <- unlist(lapply(strsplit(os10g34360, split = "orysa_"), "[", 2))
os11g35930 <- unlist(lapply(strsplit(os11g35930, split = "orysa_"), "[", 2))
os12g07690 <- unlist(lapply(strsplit(os12g07690, split = "orysa_"), "[", 2))

## function to retrieve fasta from list of gene names
get_fasta <- function(gene_list, species = "zeama") {
    
    p <- "~/winhome/Documents/03_Syntenic_linkage/01_Data/OrthoFinder/input/"
    fasta <- read.csv(paste(p, species, ".fasta", sep = ""), 
        header = FALSE, stringsAsFactors = FALSE)
    beg <- grep(fasta[,1], pattern = ">")
    end <- c(beg[2:length(beg)]-1, nrow(fasta))
    gene_list_p <- paste(">", species, "_", gene_list, sep = "")
    res <- character()
    
    for (i in 1:length(gene_list)) {
        ind <- which(fasta[,1] == gene_list_p[i])
        ind_beg <- which(ind == beg)
        res <- c(res, fasta[beg[ind_beg]:end[ind_beg],1])
    }
    return(res)
}

## apply get_fasta on the gene_list for genes surrounding PKS genes
## clusters 2, 3 and 5
## A. thaliana 
fasta_at5g13930 <- get_fasta(at5g13930, "arath")

## Z. mays
fasta_zm00001d007403 <- get_fasta(zm00001d007403, "zeama")
fasta_zm00001d052673 <- get_fasta(zm00001d052673, "zeama")

## O. sativa
fasta_os05g12180 <- get_fasta(os05g12180, "orysa")
fasta_os11g32650 <- get_fasta(os11g32650, "orysa")

## S. lycopersicum
fasta_solyc05g053550 <- get_fasta(solyc05g053550, "solyc")
fasta_solyc09g091510 <- get_fasta(solyc09g091510, "solyc")
fasta_solyc12g098100 <- get_fasta(solyc12g098100, "solyc")

## V. vinifera
fasta_gsvivt01032968001 <- get_fasta(gsvivt01032968001, "vitvi")

## other clusters 
## A. thaliana 
fasta_at1g02050 <- get_fasta(at1g02050, "arath")
fasta_at4g00040 <- get_fasta(at4g00040, "arath")
fasta_at4g34850 <- get_fasta(at4g34850, "arath")

## S. lycopersicum 
fasta_solyc01g090600 <- get_fasta(solyc01g090600, "solyc")
fasta_solyc01g111070 <- get_fasta(solyc01g111070, "solyc")
fasta_solyc05g053170 <- get_fasta(solyc05g053170, "solyc")
fasta_solyc06g043120 <- get_fasta(solyc06g043120, "solyc")

## V. vinifera
fasta_gsvivt01000521001 <- get_fasta(gsvivt01000521001, "vitvi")
fasta_gsvivt01010554001 <- get_fasta(gsvivt01010554001, "vitvi")
fasta_gsvivt01018219001 <- get_fasta(gsvivt01018219001, "vitvi")
fasta_gsvivt01024107001 <- get_fasta(gsvivt01024107001, "vitvi")
fasta_gsvivt01026213001 <- get_fasta(gsvivt01026213001, "vitvi")

## Z. mays 
fasta_zm00001d000438 <- get_fasta(zm00001d000438, "zeama")
fasta_Zm00001d007717 <- get_fasta(Zm00001d007717, "zeama")
fasta_Zm00001d013991 <- get_fasta(Zm00001d013991, "zeama")
fasta_Zm00001d019478 <- get_fasta(Zm00001d019478, "zeama")
fasta_Zm00001d021562 <- get_fasta(Zm00001d021562, "zeama")
fasta_Zm00001d032662 <- get_fasta(Zm00001d032662, "zeama")
fasta_Zm00001d040479 <- get_fasta(Zm00001d040479, "zeama")
fasta_Zm00001d052915 <- get_fasta(Zm00001d052915, "zeama")

## O. sativa
fasta_os01g41834 <- get_fasta(os01g41834, "orysa")
fasta_os03g47000 <- get_fasta(os03g47000, "orysa")
fasta_os04g23940 <- get_fasta(os04g23940, "orysa")
fasta_os05g41645 <- get_fasta(os05g41645, "orysa")
fasta_os07g11440 <- get_fasta(os07g11440, "orysa")
fasta_os10g07616 <- get_fasta(os10g07616, "orysa")
fasta_os10g34360 <- get_fasta(os10g34360, "orysa")
fasta_os11g35930 <- get_fasta(os11g35930, "orysa")
fasta_os12g07690 <- get_fasta(os12g07690, "orysa")

## write fasta to files 
## clusters 2, 3 and 5
## A. thalian 
write.table(fasta_at5g13930, quote = FALSE, row.names = FALSE, col.names = FALSE, 
    file = "figure_S5_co_expression_STRING/list_at5g13930_chs.txt")

## Z. mays
write.table(fasta_zm00001d007403, quote = FALSE, row.names = FALSE, col.names = FALSE,
    file = "figure_S5_co_expression_STRING/list_zm00001d007403_chs.txt")
write.table(fasta_zm00001d052673, quote = FALSE, row.names = FALSE, col.names = FALSE,
    file = "figure_S5_co_expression_STRING/list_zm00001d052673_chs.txt")

## O. sativa
write.table(fasta_os05g12180, quote = FALSE, row.names = FALSE, col.names = FALSE,
    file = "figure_S5_co_expression_STRING/list_os05g12180_no_chs.txt")
write.table(fasta_os11g32650, quote = FALSE, row.names = FALSE, col.names = FALSE,
    file = "figure_S5_co_expression_STRING/list_os11g32650_chs.txt")

## S. lycopersicum
write.table(fasta_solyc05g053550, quote = FALSE, row.names = FALSE, col.names = FALSE,
    file = "figure_S5_co_expression_STRING/list_solyc05g053550_chs.txt")
write.table(fasta_solyc09g091510, quote = FALSE, row.names = FALSE, col.names = FALSE,
    file = "figure_S5_co_expression_STRING/list_solyc09g091510_chs.txt")
write.table(fasta_solyc12g098100, quote = FALSE, row.names = FALSE, col.names = FALSE,
    file = "figure_S5_co_expression_STRING/list_solyc12g098100_no_chs.txt")

## V. vinifera
write.table(fasta_gsvivt01032968001, quote = FALSE, row.names = FALSE, col.names = FALSE,
    file = "figure_S5_co_expression_STRING/list_gsvivt01032968001_chs.txt")

## other clusters 
## A. thalian 
write.table(fasta_at1g02050, quote = FALSE, row.names = FALSE, col.names = FALSE,
    file = "figure_S5_co_expression_STRING/other_clusters/list_at1g02050.txt")
write.table(fasta_at4g00040, quote = FALSE, row.names = FALSE, col.names = FALSE,
    file = "figure_S5_co_expression_STRING/other_clusters/list_at4g00040.txt")
write.table(fasta_at4g34850, quote = FALSE, row.names = FALSE, col.names = FALSE,
    file = "figure_S5_co_expression_STRING/other_clusters/list_at4g34850.txt")

## S. lycopersicum
write.table(fasta_solyc01g090600, quote = FALSE, row.names = FALSE, col.names = FALSE, 
    file = "figure_S5_co_expression_STRING/other_clusters/list_solyc01g090600.txt")
write.table(fasta_solyc01g111070, quote = FALSE, row.names = FALSE, col.names = FALSE, 
    file = "figure_S5_co_expression_STRING/other_clusters/list_solyc01g111070.txt")
write.table(fasta_solyc05g053170, quote = FALSE, row.names = FALSE, col.names = FALSE, 
    file = "figure_S5_co_expression_STRING/other_clusters/list_solyc05g053170.txt")
write.table(fasta_solyc06g043120, quote = FALSE, row.names = FALSE, col.names = FALSE, 
    file = "figure_S5_co_expression_STRING/other_clusters/list_solyc06g043120.txt")

## V. vinifera
write.table(fasta_gsvivt01000521001, quote = FALSE, row.names = FALSE, col.names = FALSE, 
    file = "figure_S5_co_expression_STRING/other_clusters/list_gsvivt01000521001.txt")
write.table(fasta_gsvivt01010554001, quote = FALSE, row.names = FALSE, col.names = FALSE, 
    file = "figure_S5_co_expression_STRING/other_clusters/list_gsvivt01010554001.txt")
write.table(fasta_gsvivt01018219001, quote = FALSE, row.names = FALSE, col.names = FALSE,
    file = "figure_S5_co_expression_STRING/other_clusters/list_gsvivt01018219001.txt")
write.table(fasta_gsvivt01024107001, quote = FALSE, row.names = FALSE, col.names = FALSE,
    file = "figure_S5_co_expression_STRING/other_clusters/list_gsvivt01024107001.txt")
write.table(fasta_gsvivt01026213001, quote = FALSE, row.names = FALSE, col.names = FALSE, 
    file = "figure_S5_co_expression_STRING/other_clusters/list_gsvivt01026213001.txt")

## Z. mays 
write.table(fasta_zm00001d000438, quote = FALSE, row.names = FALSE, col.names = FALSE,
    file = "figure_S5_co_expression_STRING/other_clusters/list_zm00001d000438.txt")
write.table(fasta_Zm00001d007717, quote = FALSE, row.names = FALSE, col.names = FALSE,
    file = "figure_S5_co_expression_STRING/other_clusters/list_Zm00001d007717.txt")
write.table(fasta_Zm00001d013991, quote = FALSE, row.names = FALSE, col.names = FALSE,
    file = "figure_S5_co_expression_STRING/other_clusters/list_Zm00001d013991.txt")
write.table(fasta_Zm00001d019478, quote = FALSE, row.names = FALSE, col.names = FALSE,
    file = "figure_S5_co_expression_STRING/other_clusters/list_Zm00001d019478.txt")
write.table(fasta_Zm00001d021562, quote = FALSE, row.names = FALSE, col.names = FALSE,
    file = "figure_S5_co_expression_STRING/other_clusters/list_Zm00001d021562.txt")
write.table(fasta_Zm00001d032662, quote = FALSE, row.names = FALSE, col.names = FALSE,
    file = "figure_S5_co_expression_STRING/other_clusters/list_Zm00001d032662.txt")
write.table(fasta_Zm00001d040479, quote = FALSE, row.names = FALSE, col.names = FALSE,
    file = "figure_S5_co_expression_STRING/other_clusters/list_Zm00001d040479.txt")
write.table(fasta_Zm00001d052915, quote = FALSE, row.names = FALSE, col.names = FALSE,
    file = "figure_S5_co_expression_STRING/other_clusters/list_Zm00001d052915.txt")

## O. sativa
write.table(fasta_os01g41834, quote = FALSE, row.names = FALSE, col.names = FALSE,
    file = "figure_S5_co_expression_STRING/other_clusters/list_os01g41834.txt")
write.table(fasta_os03g47000, quote = FALSE, row.names = FALSE, col.names = FALSE,
    file = "figure_S5_co_expression_STRING/other_clusters/list_os03g47000.txt")
write.table(fasta_os04g23940, quote = FALSE, row.names = FALSE, col.names = FALSE,
    file = "figure_S5_co_expression_STRING/other_clusters/list_os04g23940.txt")
write.table(fasta_os05g41645, quote = FALSE, row.names = FALSE, col.names = FALSE,
    file = "figure_S5_co_expression_STRING/other_clusters/list_os05g41645.txt")
write.table(fasta_os07g11440, quote = FALSE, row.names = FALSE, col.names = FALSE,
    file = "figure_S5_co_expression_STRING/other_clusters/list_os07g11440.txt")
write.table(fasta_os10g07616, quote = FALSE, row.names = FALSE, col.names = FALSE,
    file = "figure_S5_co_expression_STRING/other_clusters/list_os10g07616.txt")
write.table(fasta_os10g34360, quote = FALSE, row.names = FALSE, col.names = FALSE,
    file = "figure_S5_co_expression_STRING/other_clusters/list_os10g34360.txt")
write.table(fasta_os11g35930, quote = FALSE, row.names = FALSE, col.names = FALSE,
    file = "figure_S5_co_expression_STRING/other_clusters/list_os11g35930.txt")
write.table(fasta_os12g07690, quote = FALSE, row.names = FALSE, col.names = FALSE,
    file = "figure_S5_co_expression_STRING/other_clusters/list_os12g07690_.txt")
