## barplot
setwd("/home/thomas/Documents/PhD/Polyketide_synthase/")

## load xlsx with synteny information
attr <- openxlsx::read.xlsx("./HMMer_pHMM/genes_synteny.xlsx", 
    sheet = "synteny_genes_attributes")
attr <- attr[attr[, "classification_community_apclust_dist"] != "-",]

cluster <- as.numeric(attr[, "classification_community_apclust_dist"])
cluster_u <- unique(cluster)

gene <- attr[, "gene_tandem"]
type <- attr[, "type_pPAP_reviewed"]
species <- unlist(lapply(strsplit(gene, "_"), "[", 1))
species_u <- unique(species)

## order species according to species tree
order_species <- c("selmo", "marpo", "sphfa", "phypa", "amtri", "spipo",
    "zosma", "phaeq", "denof", "musac", "anaco", "orysa", "brast", "bradi",
    "oroth", "zeama", "sorbi", "setvi", "setit", "panvi", "panha", "papso",
    "aquco", "kalla", "kalfe", "betvu", "amahy", "vitvi", "vacco", "camsi",
    "panno", "dauca", "helan", "oleeu", "salmi", "mumgu", "petin", "petax",
    "nicbe", "nicta", "nicsy", "capan", "soltu", "solyc", "solpe", "pungr",
    "eucgr", "araip", "aradu", "phavu", "glyma", "lotja", "glyur", "tripr",
    "medtr", "momch", "citla", "cucsa", "cucme", "quero", "jugre", "zizju",
    "ruboc", "frave", "prupe", "maldo", "linus", "salpu", "poptr", "ricco",
    "manes", "jatcu", "citsi", "citcl", "gosra", "theca", "corol", "corca",
    "carpa", "thepa", "eutsa", "brara", "braol", "lepme", "boest", "capru",
    "capgr", "arath", "araly", "araha")

## create binary matrix cols = cluster
mat <- matrix(0, nrow=length(cluster_u), ncol=length(species_u), 
    dimnames = list(cluster_u, species_u))

## iterate through colnames
for (i in colnames(mat)) {
    mat[unique(as.character(cluster[species == i])), i] <- 1
}

## sort according to order_species
mat <- mat[, order_species]

## sort according to number of regions
num_regions <- table(cluster)
mat <- mat[names(sort(num_regions, descrasing = TRUE)), ]

library(pheatmap)
pdf("cluster_species_presenceabsence.pdf")
pheatmap(t(mat), cluster_cols = FALSE, cluster_rows = FALSE, scale = "none", 
    border_color = alpha("darkgrey", 1))
dev.off()

## barplot of number of regions per cluster
df <- data.frame(num_regions)
df$cluster <- factor(x = df$cluster, levels = names(sort(num_regions)))
g <- ggplot(df) + geom_bar(aes(y = Freq, x = cluster), stat = "identity") + 
    ylim(0, 100) + theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank())
ggsave(g, file="cluster_numberregions.pdf")

## barplot of number of genes per cluster + identity
types <- strsplit(type, "/")
cluster_gene_l <- vector("list", length(unique(cluster)))
names(cluster_gene_l) <- sort(unique(cluster))
number_genes <- vector("list", length(cluster_gene_l))

## create data.frame to store results
df <- data.frame(cluster=sort(unique(cluster)), 
    r4c = 0, r4a = 0, r2x = 0, other = 0, nd = 0)

## iterate through clusters and store number of type 
for (i in seq_along(cluster_gene_l)) {
    number_genes[[i]] <- table(unlist(types[which(cluster == i)]))
    df$r4c[i] <- ifelse(
        is.na(number_genes[[i]]["R-4-C"]), 0, number_genes[[i]]["R-4-C"])
    df$r4a[i] <- ifelse(
        is.na(number_genes[[i]]["R-4-A"]), 0, number_genes[[i]]["R-4-A"])
    df$r2x[i] <- ifelse(
        is.na(number_genes[[i]]["R-2-X"]), 0, number_genes[[i]]["R-2-X"])
    df$other[i] <- ifelse(
        is.na(number_genes[[i]]["Other"]), 0, number_genes[[i]]["Other"])
    df$nd[i] <- ifelse(
        is.na(number_genes[[i]]["nd"]), 0, number_genes[[i]]["nd"])
}

## plot
df$cluster <- as.character(df$cluster)
df <- reshape::melt(df)
df$cluster <- factor(x = df$cluster, levels = names(sort(num_regions)))
g <- ggplot(df) + 
    geom_bar(aes(y = value, x = cluster, fill = variable), stat = "identity") +
    ylim(0, 200) + theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_blank())
ggsave(g, file = "cluster_numbergenes_identity.pdf")

## barplot of number genes per species + identity
region_species <- lapply(gene, function(x) {
    unique(unlist(lapply(strsplit(x, split = "_"), "[", 1)))
})
region_species_u <- sort(unique(unlist(region_species)))

## create data.frame to store results
df <- data.frame(species=region_species_u, 
    r4c = 0, r4a = 0, r2x = 0, other = 0, nd = 0)

## iterate through clusters and store number of type 
for (i in seq_along(region_species_u)) {
    genes_t <- table(
        unlist(types[which(unlist(region_species) == region_species_u[i])]))
    df$r4c[i] <- ifelse(
        is.na(genes_t["R-4-C"]), 0, genes_t["R-4-C"])
    df$r4a[i] <- ifelse(
        is.na(genes_t["R-4-A"]), 0, genes_t["R-4-A"])
    df$r2x[i] <- ifelse(
        is.na(genes_t["R-2-X"]), 0, genes_t["R-2-X"])
    df$other[i] <- ifelse(
        is.na(genes_t["Other"]), 0, genes_t["Other"])
    df$nd[i] <- ifelse(
        is.na(genes_t["nd"]), 0, genes_t["nd"])
}

## plot 
df <- reshape::melt(df)
df$species <- factor(df$species, levels = order_species)
g <- ggplot(df) + 
    geom_bar(aes(y = value, x = species, fill = variable), stat = "identity") +
    ylim(0, 40) + theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))
ggsave(g, file = "species_numbergenes_identity.pdf")

## barplot number  regions per species
df <- data.frame(table(unlist(region_species)))
df$Var1 <- factor(df$Var1, levels = order_species)
g <- ggplot(df) + geom_bar(aes(x = Var1, y = Freq), stat = "identity") + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0), 
          panel.border = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), axis.line = element_blank())
ggsave(g, file = "species_numberregions.pdf")