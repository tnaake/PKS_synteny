    ## visualization
## 1) load R and package ggtree
library(ggtree)
library(treeio)

setwd("/home/thomas/Documents/PhD/Polyketide_synthase/")
tree <- read.newick("HMMer_pHMM/bootstrap_n1000_outgroup_addedMCL_penium/booster_addedMCL.nw", node.label = "support")	

## load supplementary information
pks_genes <- openxlsx::read.xlsx("./HMMer_pHMM/pks_genes_tree_id_type.xlsx", sheet = 1)
attr <- openxlsx::read.xlsx("./HMMer_pHMM/genes_synteny.xlsx", sheet = "synteny_genes_attributes")

tipLabel <- tree@phylo$tip.label
tipLabel <- gsub(tipLabel, pattern = "[.]", replacement = "_")
pks_genes_ID <- pks_genes[, "ID"]
pks_genes_ID <- gsub(pks_genes_ID, pattern = "[.]", replacement = "_")
attr_ID <- gsub(attr[, "gene_tandem"], pattern = "[.]", replacement = "_")

## find indices in attr that match to tipLabel (gives row indices in attr where the entries are for tipLabel)
inds_attr <- lapply(tipLabel, function(x) 
    grep(attr_ID, pattern = paste(strsplit(x, split = "/")[[1]], collapse = "|")))
## set to NA if the sequence is not found in the network
inds_attr <- lapply(inds_attr, function(x) ifelse(length(x) == 0, NaN, x)) 
names(inds_attr) <- tipLabel


## get cluster from "classification_community_apclust_dist"
cluster <- lapply(inds_attr, function(x) 
    paste(unique(attr[x, "classification_community_apclust_dist"]), collapse = "/"))
## set cluster to "not" if not in 2,4,5,14 (but not if it is NA)
cluster24514 <- lapply(cluster, function(x) ifelse(x == "NA", x, 
    ifelse(x %in% c("2", "4", "5", "14"), x, "not")))
cluster <- unlist(cluster)
cluster24514 <- unlist(cluster24514)


## find indices in pks_genes that match to tipLabel (gives row indices in 
## pks_genes where the entries are for tipLabel)
inds <- lapply(tipLabel, function(x) 
    grep(pks_genes_ID, pattern=paste(strsplit(x, split="/")[[1]], collapse = "|"))) 
names(inds) <- tipLabel

## get pPAP_reviewed
pPAP <- unlist(lapply(inds, function(x) paste(unique(pks_genes[x, "pPAP_reviewed"]), collapse = "/")))
## get class
class <- unlist(lapply(inds, function(x) 
    paste(unique(pks_genes[x, "Class"]), collapse = "/")))
## get unranked
unranked <- unlist(lapply(inds, function(x) 
    paste(unique(pks_genes[x, which(colnames(pks_genes) == "unranked")[2]]), collapse = "/")))
## get Order
Order <- unlist(lapply(inds, function(x) 
    paste(unique(pks_genes[x, "Order"]), collapse = "/")))
## get Family
family <- unlist(lapply(inds, function(x) 
    paste(unique(pks_genes[x, "Family"]), collapse = "/")))
## get number of exon
numberexon <- unlist(lapply(inds, function(x) 
    paste(unique(pks_genes[x, "number_exon_PLAZA"]), collapse = "/")))

p <- ggplot(tree, mapping = aes(x,y), layout = "circular", ladderize = TRUE, 
        right = FALSE, branch.length = "branch.length", ndigits = NULL) + 
    geom_tree(size = 0.1) + theme_tree() + geom_treescale() + 
    theme(legend.position = "right") #+# geom_nodepoint(aes(color=support), size=0.02, pch=20) 
p <- p + ggtree:::layout_circular()
p <- p + geom_tiplab(aes(angle = angle), size = 0.2) + 
    #theme(legend.position="right") + 
    scale_color_continuous(low="red", high="green") 

phenotype <- data.frame(cluster = cluster, cluster24514 = cluster24514, 
    is = pPAP, exon = numberexon, class = class, unranked = unranked, 
    order = Order, family = family)
rownames(phenotype) <- tree@phylo$tip.label

library(tidyr)

## modify function gheatmap from ggtree
gheatmap_mod <- function (p, data, offset = 0, width = 1, low = "green", 
    high = "red", color = "white", colnames = TRUE, colnames_position = "bottom", 
    colnames_angle = 0, colnames_level = NULL, colnames_offset_x = 0, 
    colnames_offset_y = 0, font.size = 4, hjust = 0.5, legend=TRUE) {
    
    #colnames_position %<>% match.arg(c("bottom", "top"))
    variable <- value <- lab <- y <- NULL
    width <- width * (p$data$x %>% range(na.rm = TRUE) %>% diff) / ncol(data)
    isTip <- x <- y <- variable <- value <- from <- to <- NULL
    df <- p$data
    df <- df[df$isTip, ]
    start <- max(df$x, na.rm = TRUE) + offset
    dd <- as.data.frame(data)
    i <- order(df$y)
    i <- i[!is.na(df$y[i])]
    lab <- df$label[i]
    dd <- dd[match(lab, rownames(dd)), , drop = FALSE]
    dd$y <- sort(df$y)
    dd$lab <- lab
    dd <- gather(dd, variable, value, -c(lab, y))
    i <- which(dd$value == "")
    if (length(i) > 0) {
        dd$value[i] <- NA
    }
    if (is.null(colnames_level)) {
        dd$variable <- factor(dd$variable, levels = colnames(data))
    }
    else {
        dd$variable <- factor(dd$variable, levels = colnames_level)
    }
    V2 <- start + as.numeric(dd$variable) * width
    mapping <- data.frame(from = dd$variable, to = V2)
    mapping <- unique(mapping)
    dd$x <- V2
    dd$width <- width
    dd[[".panel"]] <- factor("Tree")
    if (is.null(color)) {
        p2 <- p + geom_tile(data = dd, aes(x, y, fill = value),
            width = width, height=1, inherit.aes = FALSE)
    }
    else {
        p2 <- p + geom_tile(data = dd, aes(x, y, fill = value), 
            width = width, height=1, color = color, inherit.aes = FALSE)
    }
    if (is(dd$value, "numeric")) {
        p2 <- p2 + scale_fill_gradient(low = low, high = high, na.value = NA)
    }
    else {
        p2 <- p2 + scale_fill_discrete(na.value = NA)
    }
    if (colnames) {
        if (colnames_position == "bottom") {
            y <- 0
        }
        else {
            y <- max(p$data$y) + 1
        }
        mapping$y <- y
        mapping[[".panel"]] <- factor("Tree")
        p2 <- p2 + geom_text(data = mapping, aes(x = to, y = y, 
            label = from), size = font.size, inherit.aes = FALSE, 
            angle = colnames_angle, nudge_x = colnames_offset_x, 
            nudge_y = colnames_offset_y, hjust = hjust)
    }
    if (!legend) {
        p2 <- p2 + theme(legend.position = "none", legend.title = element_blank())
        if (!colnames) {
            p2 <- p2 + scale_y_continuous(expand = c(0, 0))
        }
        attr(p2, "mapping") <- mapping
        return(p2)
    } else {
        ## Function to extract legend
        g_legend <- function(a.gplot){ 
            tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
            leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
            legend <- tmp$grobs[[leg]] 
            legend
        } 
        p2 <- p2 + theme(legend.text=(element_text(size=2.5)), legend.title=element_blank())
        
        legend <- g_legend(p2) 
        return(legend)
        #grid.newpage()
        # grid.draw(legend) 
    }
    
}

p2 <- gheatmap_mod(p, phenotype, offset = 0.1, width = 0.1, font.size = 1, 
    colnames_angle = -45, hjust = 0, legend = F)
ggsave(plot = p2, filename = "HMMer_pHMM/phylogenetic_annotation.pdf", 
    scale = 2.5, useDingbats = FALSE)

p2 <- gheatmap_mod(p, phenotype, offset = 0.1, width = 0.1, font.size = 1, 
    colnames_angle = -45, hjust = 0, legend = T)
ggsave(plot = p2, filename = "HMMer_pHMM/phylogenetic_annotation_legend.pdf", 
    scale = 2.5, useDingbats = FALSE)

## plot tree with bootstrap values
p <- ggplot(tree, mapping = aes(x,y), layout = "circular", ladderize = TRUE, 
    right = FALSE, branch.length = "branch.length", ndigits = NULL) + 
    geom_tree(size = 0.1) + theme_tree() + geom_treescale() + 
    theme(legend.position = "right") + 
    geom_nodepoint(aes(color = support), size = 0.02, pch = 20) 
p <- p + ggtree:::layout_circular()
p <- p + geom_tiplab(aes(angle = angle), size = 0.2) + 
    scale_color_continuous(low = "red", high = "green") 
ggsave(plot = p, filename = "supp_fig_bootstrap_values/HMMerphylogenetic_bootstrapsupport.pdf", 
    scale = 2.5, useDingbats = FALSE)

q <- qplot(tree@data$support, geom = "histogram", binwidth = 0.025, 
    main = "histogram for node support", xlab = "support") + 
    theme_bw() +  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.border = element_blank())
ggsave(plot = q, 
    filename = "supp_fig_bootstrap_values/HMMerphylogenetic_bootstrapsupport_hist.pdf")
