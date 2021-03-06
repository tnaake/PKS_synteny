## setwd
setwd("W:/Thomas/Data/synteny")

## get pks genes from output of Orthofinder or MCL
## for Orthofinder
genes_table <- read.table("./Results_Nov05/family_orthogroups.txt", sep = "\t", 
    header = FALSE, stringsAsFactors = FALSE)
pks_genes <- sort(genes_table[genes_table[, 2] == "OG0000256", 1])
## for MCL
genes_table_mcl <- read.table("./Results_Nov05/family_orthogroups_mcl15.txt", 
    sep = "\t", header = FALSE, stringsAsFactors = FALSE)
pks_genes_mcl <- sort(genes_table_mcl[genes_table_mcl[, 2] == "group_365", 1]) 

## vector with all pks genes identified by orthofinder and mcl
pks_genes_all <- sort(unique(c(pks_genes, pks_genes_mcl)) )

length(pks_genes)
length(pks_genes_all)
which(!pks_genes %in% pks_genes_mcl)
distribution_pks_all <- table(unlist(lapply(strsplit(pks_genes_all, split="_"), "[", 1)))


## create a vector that contains abbreviations and full names of species
species_abbr <- rbind(c("actch", "Actinidia chinensis"), 
    c("amahy", "Amaranthus hypochondriacus"), c("amtri", "Amborella trichopoda"), 
    c("anaco", "Ananas comosus"), c("aquco", "Aquilegia coerulea"), 
    c("aradu", "Arachis duranensis"), c("araip", "Arachis ipaensis"), 
    c("araly", "Arabidopsis lyrata"), c("arath", "Arabidopsis thaliana"), 
    c("artan", "Arthemisia annua"), c("aspof", "Asparagus officinalis"), 
    c("auran", "Aureococcus aneophagefferens"), c("azofi", "Azolla filiculoides"), 
    c("betpe", "Betula pendula"),c("betvu", "Beta vulgaris"), 
    c("bradi", "Brachipodium distachyon"), c("brana", "Brassica napus"), 
    c("braol", "Brassica oleracea"), c("brara", "Brassica rapa"), 
    c("cajca", "Cajanus cajan"), c("camsa", "Camelina sativa"), 
    c("camsi", "Camellia sinensis"), c("cansa", "Cannabis sativa"), 
    c("capan", "Capsicum annuum"), c("capgr", "Capsella grandiflora"), 
    c("capru", "Capsella rubella"), c("carpa", "Carica papaya"), 
    c("chequ", "Chenopodium quinoa"), c("cicar", "Cicer arietinum"), 
    c("citcl", "Citrus clementina"), c("citla", "Citrullus lanatus"), 
    c("citsi", "Citrus sinensis"), c("cofca", "Coffea arabica"), 
    c("corca", "Corchorus capsularis"), c("corol", "Corchorus olitorius"), 
    c("covsu", "Coccomyxa subellipsoidea C-169"), c("cucme", "Cucumis melo"), 
    c("cucsa", "Cucumis sativus"), c("dauca", "Daucus carota"), 
    c("denof", "Dendrobium officinale"), c("ectsi", "Ectocarpus siliculosus"), 
    c("elagu", "Elaeis guineensis"), c("eucgr", "Eucalyptus grandis"), 
    c("eutsa", "Eutrema salsugineum"), c("frave", "Fragaria vesca"), 
    c("genau", "Genlisea aurea"), c("ginbi", "Ginkgo biloba"), 
    c("glyma", "Glycine max"), c("glyur", "Glycyrrhiza uralensis"), 
    c("gnemo", "Gnetum montanum"), c("gosra", "Gossypium raimondii"), 
    c("helan", "Helianthus annuus"), c("horvu", "Hordeum vulgare"),
    c("humlu", "Humulus lupulus"), c("iponi", "Ipomoea nil"),
    c("jatcu", "Jatropha curcas"), c("jugre", "Juglans regia"), 
    c("kalfe", "Kalanchoe fedtschenkoi"), c("leepe", "Leersia perrieri"),
    c("lepme", "Lepidium meyenii"), c("linus", "Linum usitatissimum"), 
    c("lotja", "Lotus japonicus"), c("lupan", "Lupinus angustifolius"),
    c("maldo", "Malus domestica"), c("manes", "Manihot esculenta"), 
    c("marpo", "Marchantia polymorpha"), c("medtr", "Medicago truncatula"), 
    c("mimgu", "Mimulus guttatus"), c("momch", "Momordica charantia"), 
    c("morno", "Morus notabilis"), c("musac", "Musa acuminata"), 
    c("nelnu", "Nelumbo nucifera"), c("nicat", "Nicotiana attenuata"), 
    c("nicbe", "Nicotiana benthamiana"), c("nicsy", "Nicotiana sylvestris"), 
    c("nicta", "Nicotiana tabacum"), c("oleeu", "Olea europaea"), 
    c("oroth", "Oropetium thomaeum"), c("oryru", "Oryza rufipogon"), 
    c("orysa", "Oryza sativa"), c("ostlu", "Ostreococcus lucimarinus"), 
    c("papso", "Papaver somniferum"), c("petax", "Petunia axillaris"), 
    c("petin", "Petunia inflata"), c("phaeq", "Phalaenopsis equestris"), 
    c("phavu", "Phaseolus vulgaris"), c("phoda", "Phoenix dactylifera"), 
    c("phypa", "Physcomitrella patens"), c("picab", "Picea abies"), 
    c("pinta", "Pinus taeda"), c("poptr", "Populus trichocarpa"), 
    c("prupe", "Prunus persica"), c("pseme", "Pseudotsuga menziesii"), 
    c("pungr", "Punica granatum"), c("pyrbr", "Pyrus x bretschneideri"),
    c("quero", "Quercus robur"), c("ricco", "Ricinus communis"), 
    c("ruboc", "Rubus occidentalis"), c("salcu", "Salvinia cucullata"),
    c("salmi", "Salvia miltiorrhiza"), c("selmo", "Selaginella moellendorffii"), 
    c("setit", "Setaria italica"), c("solpe", "Solanum pennellii"), 
    c("soltu", "Solanum tuberosum"), c("solyc", "Solanum lycopersicum"), 
    c("sorbi", "Sorghum bicolor"), c("spipo", "Spirodela polyrhiza"), 
    c("tarha", "Tarenaya hassleriana"), c("theca", "Theobroma cacao"), 
    c("thepa", "Thelungiella parvula"), c("triae", "Triticum aestivum"),
    c("tripr", "Trifolium pratense"), c("vacco", "Vaccinium corymbosum"), 
    c("vitvi", "Vitis vinifera"), c("zeama", "Zea mays"),
    c("zizju", "Ziziphus jujuba"), c("zosma", "Zostera marina"))

## rename the abbreviations to full names
names(distribution_pks_all)[!names(distribution_pks_all) %in% species_abbr[, 1]]
species_abbr[, 1][!species_abbr[, 1] %in% names(distribution_pks_all)]
names(distribution_pks_all) <- species_abbr[match(species_abbr[,1], names(distribution_pks_all)), 2]

names_sort_all <- names(sort(distribution_pks_all))

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

pks_genes <- unique(pks_genes)
pks_genes_mcl <- unique(pks_genes_mcl)
pks_genes_all <- unique(pks_genes_all)

## read supplementary information
supp <- openxlsx::read.xlsx("H:/Documents/03_Syntenic_linkage/01_Data/synteny_network_results/pks_genes_tree_id_type.xlsx", 
    sheet = "pks_genes_tree_id_type")

## create data.frame with information on pks genes and species
df <- data.frame(pks_genes_all = pks_genes_all, 
    species = unlist(lapply(strsplit(pks_genes_all, split = "_"), "[", 1)), 
    ofi = pks_genes_all %in% pks_genes, mcl = pks_genes_all %in% pks_genes_mcl, 
    length = supp[match(pks_genes_all, supp[, "ID"]), "length"] >= 200,
    pPAP = supp[match(pks_genes_all, supp[, "ID"]), "pPAP_reviewed"])

## get df_200 that only contains seqs that have an AA length >= 200
df_200 <- df[df[, "length"] == TRUE,]

## get df_ofi, df_mcl that contain only seqs that were identified by OrthoFinder
## and MCL
df_ofi <- df[df[, "ofi"] == TRUE, ]
df_mcl <- df[df[, "mcl"] == TRUE, ]

## get data.frame that only contains seqs that have an AA length >= 200
df_200_ofi <- df_200[df_200[, "ofi"] == TRUE, ]
df_200_mcl <- df_200[df_200[, "mcl"] == TRUE, ]

## define function overlap to find Orthofinder-, MCL-specific and shared PKS
overlap <- function(df_ofi, df_mcl) {
    spec_ofi <- df_ofi[, "species"]
    spec_mcl <- df_mcl[, "species"]
    spec_u <- unique(c(as.character(spec_ofi), as.character(spec_mcl)))
    res <- matrix(ncol = 3, nrow = 117) ## nrow is equal to numbre of species
    rownames(res) <- as.character(unique(df[, "species"]))
    colnames(res) <- c("ofi_specific", "mcl_specific", "shared")
    for (i in 1:length(spec_u)) {
        gene_ofi <- df_ofi[spec_ofi == spec_u[i], "pks_genes_all"]
        gene_mcl <- df_mcl[spec_mcl == spec_u[i], "pks_genes_all"]
        ofi_specific <- sum(!gene_ofi %in% gene_mcl)
        mcl_specific <- sum(!gene_mcl %in% gene_ofi)
        shared <- length(intersect(gene_ofi, gene_mcl))
        res[as.character(spec_u[i]), 1:3] <- c(ofi_specific, mcl_specific, shared)
    }
    return(res)
}

## apply the function
overlap_all <- overlap(df_ofi, df_mcl)
overlap_200 <- overlap(df_200_ofi, df_200_mcl)

overlap_all <- melt(overlap_all)
overlap_200 <- melt(overlap_200)

## set ostlu in overlap_200 to 0, was NA
overlap_200[overlap_200[, 1] == "ostlu", "value"] <- 0

overlap_df <- rbind(overlap_200, overlap_all)

overlap_all[, 1] <- species_abbr[match(overlap_all[, 1], species_abbr[, 1]), 2]
overlap_200[, 1] <- species_abbr[match(overlap_200[, 1], species_abbr[, 1]), 2]

## order according to increasing values in overlap_200
df_200_mcl_spec <- melt(table(df_200_mcl[, 2]))
overlap_all$X1 <-factor(overlap_all$X1, 
    levels = overlap_all$X1[order(df_200_mcl_spec[, 2])])
overlap_200$X1 <-factor(overlap_200$X1, 
    levels = overlap_200$X1[order(df_200_mcl_spec[, 2])])

## do the plotting
setwd("~/Documents/03_Syntenic_linkage/01_Data/synteny_network_results/figure_1_distribution")
pdf("distribution_all.pdf")
ggplot(overlap_all) + geom_bar(aes(x = X1, y = value, fill = X2), 
    stat = "identity", width = 0.5) + ylim(c(0, 100)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0, size = 5)) 
dev.off()
pdf("distribution_200.pdf")
ggplot(overlap_200) + geom_bar(aes(x = X1, y = value, fill = X2), 
    stat = "identity", width = 0.5) + ylim(c(0, 100)) +  
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0, size = 5)) 
dev.off()

## identity plot 
identity_percent <- function(df = df_ofi) {
    spec <- as.character(df[, "species"])
    spec_u <- levels(df[, "species"])
    
    res <- matrix(0, nrow = length(spec_u), ncol = 5)
    colnames(res) <- c("nd", "Other", "R-2-X", "R-4-A", "R-4-C")
    rownames(res) <- unique(as.character(df_ofi[, "species"]))
    
    for (i in 1:length(spec_u)) {
        type <- df[which(spec == spec_u[i]), "pPAP"]
        if (sum(table(type)) > 0)
            res[as.character(spec_u[i]), ] <- table(type) / sum(table(type))
    }
    return(res) 
}

df_ofi_type <- identity_percent(df_ofi)
df_mcl_type <- identity_percent(df_mcl)
df_200_ofi_type <- identity_percent(df_200_ofi)
df_200_mcl_type <- identity_percent(df_200_mcl)

rownames(df_ofi_type) <- species_abbr[
    match(rownames(df_ofi_type), species_abbr[, 1]), 2]
rownames(df_mcl_type) <- species_abbr[
    match(rownames(df_mcl_type), species_abbr[, 1]), 2]
rownames(df_200_ofi_type) <- species_abbr[
    match(rownames(df_200_ofi_type), species_abbr[, 1]), 2]
rownames(df_200_mcl_type) <- species_abbr[
    match(rownames(df_200_mcl_type), species_abbr[, 1]), 2]

## melt the data.frames
df_ofi_type <- melt(df_ofi_type)
df_mcl_type <- melt(df_mcl_type)
df_200_ofi_type <- melt(df_200_ofi_type)
df_200_mcl_type <- melt(df_200_mcl_type)

## reorder the factors
df_ofi_type$X1 <- factor(df_ofi_type$X1, 
    levels = df_ofi_type$X1[order(df_200_mcl_spec[, 2])])
df_mcl_type$X1 <- factor(df_mcl_type$X1, 
    levels = df_mcl_type$X1[order(df_200_mcl_spec[, 2])])
df_200_ofi_type$X1 <-factor(df_200_ofi_type$X1, 
    levels = df_200_ofi_type$X1[order(df_200_mcl_spec[, 2])])
df_200_mcl_type$X1 <-factor(df_200_mcl_type$X1, 
    levels = df_200_mcl_type$X1[order(df_200_mcl_spec[, 2])])

## write to pdfs
pdf("distribution_type_ofi_all.pdf")
ggplot(df_ofi_type) + 
    geom_bar(aes(x = X1, y = value, fill = X2), stat = "identity", width = 0.5) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0, size = 5))
dev.off()
pdf("distribution_type_mcl_all.pdf")
ggplot(df_mcl_type) + 
    geom_bar(aes(x = X1, y = value, fill = X2), stat = "identity", width = 0.5) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0, size = 5))
dev.off()
pdf("distribution_type_ofi_200.pdf")
ggplot(df_200_ofi_type) + 
    geom_bar(aes(x = X1, y = value, fill = X2), stat = "identity", width = 0.5) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0, size = 5))
dev.off()
pdf("distribution_type_mcl_200.pdf")
ggplot(df_200_mcl_type) + 
    geom_bar(aes(x = X1, y = value, fill = X2), stat = "identity", width = 0.5) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0, size = 5))
dev.off()

