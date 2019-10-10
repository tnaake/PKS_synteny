## Figure expression 
setwd("~/Documents/03_Syntenic_linkage/01_Data/synteny_network_results/figure_S6_coexpression_five_species_cluster/")
library(openxlsx)
tkm <- read.xlsx("expression_values_summary.xlsx", sheet = "CoNekT_final_custom")
tkm <- cbind(tkm, 
    species = unlist(lapply(strsplit(tkm[, "gene"], split="_"), "[", 1)))

ind_values <- which(colnames(tkm) == "roots/rhizoids"):which(colnames(tkm) == "pollen")
tkm <- tkm[!apply(tkm[, ind_values], 1, function(x) all(x == 0)), ]

## log2 transform
tkm[, ind_values] <- log2(tkm[, ind_values] + 1)

## remove PKs of type "nd"
tkm <- tkm[tkm[, "type"] != "nd",]

tkm_m <- as.matrix(tkm[, ind_values])
rownames(tkm) <- tkm[, 1]
rownames(tkm_m) <- tkm[, 1]
heatmap(tkm_m)

## calculate Pearson correlation values (pairwise complete)
cor_tkm <- cor(t(tkm_m), use = "pairwise.complete.obs")

## rm col and row where there are NA values
which(is.na(cor_tkm), arr.ind = TRUE)
cor_tkm_rm <- cor_tkm##[-27,-27]
dist_cor_tkm <- dist(cor_tkm_rm)

## affinity propagation clustering on correlation values
ap_dis <- apcluster::apcluster(s = cor_tkm_rm, convits = 1000, 
    maxits = 10000, nonoise = TRUE, seed = 1000)

## plot 
tkm <- cbind(tkm, apclust="black")
tkm$apclust <- as.character(tkm$apclust)
tkm[tkm[, 1] %in% names(ap_dis@clusters[[1]]), "apclust"] <- "1"
tkm[tkm[, 1] %in% names(ap_dis@clusters[[2]]), "apclust"] <- "2"
tkm[tkm[, 1] %in% names(ap_dis@clusters[[3]]), "apclust"] <- "3"
tkm[tkm[, 1] %in% names(ap_dis@clusters[[4]]), "apclust"] <- "4"

## create dendrogram
dend <- as.dendrogram(hclust(dist_cor_tkm))
labels_colors(dend) <- tkm[order.dendrogram(dend), "apclust"]
plot(dend)
#dist_matrix <- as.matrix(dist_cor_tkm)

library(pheatmap)
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(-1, 0, length.out = ceiling(paletteLength/2) + 1), 
    seq(max(cor_tkm, na.rm = T) / paletteLength, 1, 
    length.out = floor(paletteLength / 2)))
p <- pheatmap(cor_tkm, border_color = NA, fontsize = 5, breaks = myBreaks, 
    color = myColor, annotation_col = tkm[, c("type", "apclust")], 
    annotation_row = tkm[, c("type", "apclust", "species")])
ggsave(p, file="pks_expression_pheatmap_pks.pdf")

## plot according to clusters
tkm_m_z <- t(apply(tkm_m, 1, scale))
colnames(tkm_m_z) <- colnames(tkm_m)

tkm_m_melt <- reshape::melt(tkm_m_z)  
tkm_m_melt <- cbind(tkm_m_melt, clusters = 0)
tkm_m_melt <- cbind(tkm_m_melt, lwid = 0.5)
tkm_m_melt <- cbind(tkm_m_melt, colour = scales::alpha("grey", .3))
tkm_m_melt$colour <- as.character(tkm_m_melt$colour)
tkm_m_melt[tkm_m_melt[, "X1"] %in% names(ap_dis@clusters[[1]]), "clusters"] <- 1
tkm_m_melt[tkm_m_melt[, "X1"] %in% names(ap_dis@clusters[[2]]), "clusters"] <- 2
tkm_m_melt[tkm_m_melt[, "X1"] %in% names(ap_dis@clusters[[3]]), "clusters"] <- 3
tkm_m_melt[tkm_m_melt[, "X1"] %in% names(ap_dis@clusters[[4]]), "clusters"] <- 4
tkm_m_melt[tkm_m_melt[, "X1"] %in% names(ap_dis@exemplars), "lwid"] <- 2
tkm_m_melt[tkm_m_melt[, "X1"] %in% names(ap_dis@exemplars), "colour"] <- "black"
levels(tkm_m_melt$clusters) <- c("1", "2", "3", "4")
tkm_m_melt$X2 <- as.character(tkm_m_melt$X2)

library(ggplot2)
ggplot() + 
    geom_point(data = tkm_m_melt[tkm_m_melt[, "colour"] == alpha("grey", .3), ], 
        aes(x = X2, y = value, group = X1, size = lwid, colour = colour)) + 
    geom_point(data = tkm_m_melt[tkm_m_melt[, "colour"] == "black",], 
        aes(x = X2, y = value, group = X1, size = lwid, colour = colour)) + 
    geom_line(data = tkm_m_melt[tkm_m_melt[, "colour"] == alpha("grey", .3), ], 
        aes(x = X2, y = value, group = X1, colour = colour)) + 
    geom_line(data = tkm_m_melt[tkm_m_melt[, "colour"] == "black",], 
        aes(x = X2, y = value, group = X1, colour = colour)) + 
    #geom_ribbon(data=predframe,aes(ymin=lwr,ymax=upr),alpha=0.3))
    ##geom_point(aes(x=X2, y=value, group=X1, colour=colour, size=lwid)) + #, colour=colour)) + 
    facet_grid(. ~ clusters) + ylim(-3, 3) + 
    theme(axis.text.x = element_text(angle = -90)) + 
    scale_colour_manual(name = "grp",values = unique(tkm_m_melt$colour))


## check expression for LAP5/6 orthologs

## create vector of gene names that are LAP5/6 orthologs
lap <- c("arath_AT1G02050", "arath_AT4G00040", "arath_AT4G34850",
    "orysa_LOC_Os07g22850", "orysa_LOC_Os10g34360", 
    "selmo_Smo122361|PACid_15419808", "selmo_Smo231846|PACid_15419824",
    "solyc_Solyc01g090600.3.1", "solyc_Solyc01g111070.3.1",
    "vitvi_GSVIVT01018219001", "vitvi_GSVIVT01024107001",
    "zeama_Zm00001d013991", "zeama_Zm00001d019478", "zeama_Zm00001d032662"
)

## plot 
ggplot(tkm_m_melt[tkm_m_melt[, "X1"] %in% lap,]) + 
    geom_line(aes(x = X2, y = value, group = X1)) + 
    theme(axis.text.x = element_text(angle = -90)) 

## use geom_ribbon --> calculate 25% lower and 75% upper quantile for each i in cluster
tkm_rib <- data.frame(cluster = numeric(), condition = character(), 
    qu_l = numeric(), qu_u = numeric(), examplar = numeric())

## iterate through clusters
for (i in 1:4) {
    tkm_i <- tkm_m_melt[tkm_m_melt[, "clusters"] == i, ]
    
    for (j in unique(tkm_i[, "X2"])) {
        qu_l <- quantile(tkm_i[tkm_i[, "X2"] == j, "value"], 0.25, na.rm = TRUE)
        qu_u <- quantile(tkm_i[tkm_i[, "X2"] == j, "value"], 0.75, na.rm = TRUE)
        value_examplar <- 
            tkm_i[as.character(tkm_i[, "X1"]) == names(ap_dis@exemplars[i]) & tkm_i[, "X2"] == j, "value"]
        
        ## add values to data.frame
        tkm_rib <- data.frame(
            cluster = c(tkm_rib$cluster, i),
            condition = as.character(c(as.character(tkm_rib$condition), j)),
            qu_l = c(tkm_rib$qu_l, qu_l),
            qu_u = c(tkm_rib$qu_u, qu_u),
            examplar = c(tkm_rib$examplar, value_examplar)
        )
    }
}

## reorder factors
tkm_rib$condition <- factor(tkm_rib$condition, 
    levels = colnames(tkm_m_z)[amap::hcluster(t(tkm_m_z))$order])

## plot 
g <- ggplot(tkm_rib) + 
    geom_ribbon(aes(x = condition, ymin = qu_l, ymax = qu_u, group = cluster), 
        fill = "lightgrey") + 
    geom_point(aes(x = condition, y = examplar, group = cluster), size = 0.5) + 
    geom_line(aes(x = condition, y = examplar, group = cluster), size = 0.9) + 
    facet_grid(. ~ cluster) + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) 
ggsave(g, file = "pks_expression_cluster_ribbon.pdf")

## cluster 1
cluster_level <- c("2", "3", "4", "7", "9", "11", "13", "14", "15", "18", 
    "19", "20", "21", "25", "26", "28", "29", "30", "-")

tkm_clust1 <- tkm[tkm[, "gene"] %in% names(ap_dis@clusters[[1]]), ][, 1:3]
tkm_clust1 <- rbind(tkm_clust1, c(NA, "2", "0"))
tkm_clust1 <- rbind(tkm_clust1, c(NA, "4", "0"))
tkm_clust1 <- rbind(tkm_clust1, c(NA, "7", "0"))
tkm_clust1 <- rbind(tkm_clust1, c(NA, "11", "0"))
tkm_clust1 <- rbind(tkm_clust1, c(NA, "21", "0"))
tkm_clust1 <- rbind(tkm_clust1, c(NA, "25", "0"))
tkm_clust1 <- rbind(tkm_clust1, c(NA, "26", "0"))
tkm_clust1 <- rbind(tkm_clust1, c(NA, "28", "0"))
tkm_clust1 <- rbind(tkm_clust1, c(NA, "29", "0"))
tkm_clust1 <- rbind(tkm_clust1, c(NA, "30", "0"))
tkm_clust1 <- rbind(tkm_clust1, c(NA, "-", "0"))

## reorder factors 
tkm_clust1$cluster <- factor(tkm_clust1$cluster, levels=cluster_level)

## plot
g <- ggplot(tkm_clust1) + 
    geom_bar(aes(x = cluster, fill = type, colour = type)) + ylim(0, 12)
ggsave(g, file = "pks_expression_heatmap_clust1_all.pdf")

## plot for CHS
tkm_clust1_chs <- tkm_clust1
tkm_clust1_chs[which(tkm_clust1_chs[, "type"] != "R-4-C"), "type"] <- NA
g <- ggplot(tkm_clust1_chs) + 
    geom_bar(aes(x = cluster, fill = type, colour = type)) + ylim(0, 6)
ggsave(g, file = "pks_expression_heatmap_clust1_chs.pdf")

## cluster2
tkm_clust2 <- tkm[tkm[, "gene"] %in% names(ap_dis@clusters[[2]]), ][, 1:3]
tkm_clust2 <- rbind(tkm_clust2, c(NA, "9", "0"))
tkm_clust2 <- rbind(tkm_clust2, c(NA, "13", "0"))
tkm_clust2 <- rbind(tkm_clust2, c(NA, "14", "0"))
tkm_clust2 <- rbind(tkm_clust2, c(NA, "15", "0"))
tkm_clust2 <- rbind(tkm_clust2, c(NA, "18", "0"))
tkm_clust2 <- rbind(tkm_clust2, c(NA, "19", "0"))
tkm_clust2 <- rbind(tkm_clust2, c(NA, "28", "0"))

## reorder factors 
tkm_clust2$cluster <- factor(tkm_clust2$cluster, levels = cluster_level)

## plot
g <- ggplot(tkm_clust2) + 
    geom_bar(aes(x = cluster, fill = type, colour = type)) + ylim(0, 12)
ggsave(g, file = "pks_expression_heatmap_clust2_all.pdf")

## plot for CHS
tkm_clust2_chs <- tkm_clust2
tkm_clust2_chs[which(tkm_clust2_chs[, "type"] != "R-4-C"), "type"] <- "0"
g <- ggplot(tkm_clust2_chs) + 
    geom_bar(aes(x = cluster, fill = type, colour = type)) + ylim(0, 6)
ggsave(g, file = "pks_expression_heatmap_clust2_chs.pdf")

## cluster 3
tkm_clust3 <- tkm[tkm[, "gene"] %in% names(ap_dis@clusters[[3]]), ][, 1:3]
tkm_clust3 <- rbind(tkm_clust3, c(NA, "7", "0"))
tkm_clust3 <- rbind(tkm_clust3, c(NA, "11", "0"))
tkm_clust3 <- rbind(tkm_clust3, c(NA, "14", "0"))
tkm_clust3 <- rbind(tkm_clust3, c(NA, "15", "0"))
tkm_clust3 <- rbind(tkm_clust3, c(NA, "21", "0"))
tkm_clust3 <- rbind(tkm_clust3, c(NA, "25", "0"))
tkm_clust3 <- rbind(tkm_clust3, c(NA, "26", "0"))
tkm_clust3 <- rbind(tkm_clust3, c(NA, "28", "0"))
tkm_clust3 <- rbind(tkm_clust3, c(NA, "29", "0"))
tkm_clust3 <- rbind(tkm_clust3, c(NA, "30", "0"))
tkm_clust3 <- rbind(tkm_clust3, c(NA, "26", "-"))

## reorder factors 
tkm_clust3$cluster <- factor(tkm_clust3$cluster, levels = cluster_level)

## plot 
g <- ggplot(tkm_clust3) + 
    geom_bar(aes(x = cluster, fill = type, colour = type)) + ylim(0, 12)
ggsave(g, file = "pks_expression_heatmap_clust3_all.pdf")

## plot for cHS
tkm_clust3_chs <- tkm_clust3
tkm_clust3_chs[which(tkm_clust3_chs[, "type"] != "R-4-C"), "type"] <- NA
g <- ggplot(tkm_clust3_chs) + 
    geom_bar(aes(x = cluster, fill = type, colour = type)) + ylim(0, 6)
ggsave(g, file = "pks_expression_heatmap_clust3_pks.pdf")

## cluster 4
tkm_clust4 <- tkm[tkm[, "gene"] %in% names(ap_dis@clusters[[4]]), ][, 1:3]
tkm_clust4 <- rbind(tkm_clust4, c(NA, "2", "0"))
tkm_clust4 <- rbind(tkm_clust4, c(NA, "7", "0"))
tkm_clust4 <- rbind(tkm_clust4, c(NA, "9", "0"))
tkm_clust4 <- rbind(tkm_clust4, c(NA, "13", "0"))
tkm_clust4 <- rbind(tkm_clust4, c(NA, "14", "0"))
tkm_clust4 <- rbind(tkm_clust4, c(NA, "15", "0"))
tkm_clust4 <- rbind(tkm_clust4, c(NA, "18", "0"))
tkm_clust4 <- rbind(tkm_clust4, c(NA, "19", "0"))
tkm_clust4 <- rbind(tkm_clust4, c(NA, "21", "0"))
tkm_clust4 <- rbind(tkm_clust4, c(NA, "26", "0"))
tkm_clust4 <- rbind(tkm_clust4, c(NA, "30", "0"))

## reorder factors 
tkm_clust4$cluster <- factor(tkm_clust4$cluster, levels = cluster_level)

## plot 
g <- ggplot(tkm_clust4) + 
    geom_bar(aes(x = cluster, fill = type, colour = type)) + ylim(0, 12)
ggsave(g, file = "pks_expression_heatmap_clust4_all.pdf")

## plot for CHS
tkm_clust4_chs <- tkm_clust4
tkm_clust4_chs[which(tkm_clust4_chs[, "type"] != "R-4-C"), "type"] <- "0"
g <- ggplot(tkm_clust4_chs) + 
    geom_bar(aes(x=cluster, fill = type, colour = type)) + ylim(0, 6)
ggsave(g, file = "pks_expression_heatmap_clust4_pks.pdf")

table(tkm[tkm[, "gene"] %in% names(ap_dis@clusters[[1]]), "cluster"])
table(tkm[tkm[, "gene"] %in% names(ap_dis@clusters[[2]]), "cluster"])
table(tkm[tkm[, "gene"] %in% names(ap_dis@clusters[[3]]), "cluster"])
table(tkm[tkm[, "gene"] %in% names(ap_dis@clusters[[4]]), "cluster"])

##M.pca <- prcomp(cor_tkm_rm, scale = FALSE, center = FALSE)

