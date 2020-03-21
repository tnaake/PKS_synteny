setwd("AG-Fernie/Thomas/Data/synteny/synteny_MCScanX")
tmp <- c(789, 790, 1116, 2028, 2029, 3249, 4988:4992, 5876, 5877)

 "araip_araly" 

"actch_actch"


######################### 1 #######################

## load orthogroups

type <- "MCL" ## "MCL"

if (type == "OrthoFinder") {
    og <- read.table("family_orthogroups.txt")
} 
if (type == "MCL") {
    og <- read.table("family_orthogroups_mcl15.txt")
}


species <- list.files()[grep(list.files(), pattern = "[.]gff")]
species <- unlist(lapply(strsplit(species, split = ".gff"), "[", 1))
species_combn <- cbind(combn(species, 2), rbind(species, species))



######################### 2 #######################
## parse gff to MCScanX_gff compatible file
setwd("AG-Fernie/Thomas/Data/synteny/synteny_MCScanX")

## two functions for parsing adjusting for the different files and how the
## gene ids need to be constructed
write_MCScanX_gff <- function(species_name, is_type_col = 3, is_type, 
    first_sep, second_sep, starting_position_col = 4, ending_position_col = 5, 
    gene_col = 9, uniq_id_col = 9, uniq_id_first, 
    uniq_id_second = " ", append = "", front = FALSE, ...) {
    
    gff <- read.table(species_name, header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "", ...)    
    gff <- gff[which(as.character(gff[, is_type_col]) == is_type), ]   
    
    ## make unique
    uniq_id <- as.character(gff[, uniq_id_col])
    if (!is.null(uniq_id_first)) uniq_id <- unlist(lapply(strsplit(uniq_id, split = uniq_id_first), "[", 2))
    if (!is.null(uniq_id_second)) uniq_id <- unlist(lapply(strsplit(uniq_id, split = uniq_id_second), "[", 1))
    gff <- gff[!duplicated(uniq_id),]
    
    ## get species name and paste with chromosome 
    species_name_new <- strsplit(x = species_name, split = ".gff")[[1]][1]
    species_chr <- paste(species_name_new, gff[, 1], sep = "_")
    ## get gene 
    gene <- as.character(gff[, gene_col])
    if (!is.null(first_sep)) gene <- unlist(lapply(strsplit(gene, split = first_sep), "[", 2))
    if (!is.null(second_sep)) gene <- unlist(lapply(strsplit(gene, split = second_sep), "[", 1))
    if (front) gene <- paste0(append, gene)
    if (!front) gene <- paste0(gene, append)
    
    gene <- paste(species_name_new, gene, sep = "_")
    
    ## get start and end position
    starting_position <- gff[, starting_position_col]
    ending_position <- gff[, ending_position_col]
    
    
    df <- data.frame(species_chr, gene, starting_position, ending_position)
    ## order according to the first column and third column
    df <- df[with(df, order(species_chr, starting_position)), ]

    ##df <- df[!duplicated(df[, "gene"]),]
    
    write.table(df, file = paste0("MCScanX_gff/", species_name_new, ".gff"),
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")    
}

write_MCScanX_gff_combine <- function(species_name, is_type_col=3, is_type, first_sep_first, second_sep_first, first_sep_second, second_sep_second, 
    starting_position_col=4, ending_position_col=5, gene_col=9, uniq_id_col=9, uniq_id_first, uniq_id_second=" ", append="", ...) {
    
    gff <- read.table(species_name, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", ...)    
    gff <- gff[which(as.character(gff[, is_type_col]) == is_type), ]   
    
    ## make unique
    uniq_id <- unlist(lapply(strsplit(as.character(gff[, uniq_id_col]), split = uniq_id_first), "[", 2))
    uniq_id <- unlist(lapply(strsplit(uniq_id, split = uniq_id_second), "[", 1))
    gff <- gff[!duplicated(uniq_id),]
    
    ## get species name and paste with chromosome 
    species_name_new <- strsplit(x = species_name, split = ".gff")[[1]][1]
    species_chr <- paste(species_name_new, gff[, 1], sep = "_")
    ## get gene 
    gene <- as.character(gff[, gene_col])
    first <- unlist(lapply(strsplit(gene, split = first_sep_first), "[", 2))
    first <- unlist(lapply(strsplit(first, split = second_sep_first), "[", 1))
    second <- unlist(lapply(strsplit(gene, split = first_sep_second), "[", 2))
    second <- unlist(lapply(strsplit(second, split = second_sep_second), "[", 1))
    gene <- paste0(species_name_new, "_", append, first, "|", "PACid_", second)
    
    ## get start and end position
    starting_position <- gff[, starting_position_col]
    ending_position <- gff[, ending_position_col]
    
    
    df <- data.frame(species_chr, gene, starting_position, ending_position)
    ## order according to the first column and third column
    df <- df[with(df, order(species_chr, starting_position)), ]
    ##df <- df[!duplicated(df[, "gene"]),]
    
    write.table(df, file = paste0("MCScanX_gff/", species_name_new, ".gff"), 
        col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")    
    
    
}

## start here with the parsing
write_MCScanX_gff("actch.gff3", is_type = "mRNA", first_sep = "Parent=", 
    second_sep = ";", uniq_id_first = "Parent=", uniq_id_second = "[.]")
write_MCScanX_gff("amahy.gff3", is_type="mRNA", first_sep="ID=", 
    second_sep="[.]", uniq_id_first="ID=", uniq_id_second="[.]") 
write_MCScanX_gff("amtri.gff3", is_type="mRNA", first_sep="ID=", 
    second_sep=";", uniq_id_first="ID=", uniq_id_second=";")
write_MCScanX_gff("anaco.gff3", is_type="mRNA", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";") 
write_MCScanX_gff("aquco.gff3", is_type="mRNA", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";", append=".p")
write_MCScanX_gff("aradu.gff", is_type="CDS", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";", skip=5)
write_MCScanX_gff("araip.gff", is_type="CDS", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff("araly.gff3", is_type="mRNA", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff("arath.gff", is_type="gene", first_sep="ID=", 
    second_sep=";", uniq_id_first="ID=", uniq_id_second=";")
write_MCScanX_gff("artan.gff", is_type="CDS", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff("aspof.gff", is_type = "CDS", first_sep = "Name=", 
    second_sep = ";", uniq_id_first = "Name=", uniq_id_second=";")
write_MCScanX_gff("azofi.gff", is_type = "gene", first_sep = "Name=", 
    second_sep = ";", uniq_id_first = "Name=", uniq_id_second = ";")
write_MCScanX_gff("auran.gff", is_type="CDS", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff("betpe.gff", is_type = "mRNA", first_sep = "ID=", 
    second_sep = ";", uniq_id_first = "ID=", uniq_id_second = ";")
write_MCScanX_gff("betvu.gff3", is_type="mRNA", first_sep="ID=", 
    second_sep=";", uniq_id_first="ID=", uniq_id_second=";")
write_MCScanX_gff("bradi.gff3", is_type="mRNA", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";", append=".p")
write_MCScanX_gff("brana.gff3", is_type = "mRNA", first_sep = "Alias=", 
    second_sep = ";", uniq_id_first = "Alias=", uniq_id_second = ";")
write_MCScanX_gff("braol.gff3", is_type = "mRNA", first_sep = "transcript_id=", 
    second_sep = ";", uniq_id_first = "transcript_id=", uniq_id_second = ";")
write_MCScanX_gff("brara.gff3", is_type="mRNA", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";", append=".p")
write_MCScanX_gff("cajca.gff", is_type = "mRNA", first_sep = "ID=", 
    second_sep = ";", uniq_id_first = "ID=", uniq_id_second = ";")
write_MCScanX_gff("camsa.gff", is_type = "CDS", first_sep = "Name=", 
    second_sep = ";", uniq_id_first = "Name=", uniq_id_second = ";")
write_MCScanX_gff("camsi.gff3", is_type="mRNA", first_sep="ID=", 
    second_sep=";", uniq_id_first="ID=", uniq_id_second=";")
write_MCScanX_gff("cansa.gff", is_type = "CDS", first_sep = "Name=", 
    second_sep = ";", uniq_id_first = "Name=", uniq_id_second = ";")
write_MCScanX_gff("capan.gff", is_type="mRNA", first_sep="ID=", 
    second_sep=";", uniq_id_first="ID=", uniq_id_second=";")
write_MCScanX_gff("capgr.gff3", is_type="mRNA", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";", append=".p")
write_MCScanX_gff("capru.gff3", is_type="mRNA", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff("carpa.gff3", is_type="mRNA", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff("chabr.gff", is_type = "CDS", first_sep = "Name=", 
    second_sep = ";", uniq_id_first = "Name=", uniq_id_second = ";")
write_MCScanX_gff("chequ.gff3", is_type = "mRNA", first_sep = "Name=", 
    second_sep = ";", uniq_id_first = "Name=", uniq_id_second = ";")
write_MCScanX_gff("chlre.gff3", is_type="mRNA", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff("cicar.gff", is_type = "mRNA", first_sep = "ID=", 
    second_sep = ";", uniq_id_first = "ID=", uniq_id_second = ";")
write_MCScanX_gff("citcl.gff", is_type="CDS", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff("citla.gff", is_type="mRNA", first_sep="ID=", 
    second_sep=";", uniq_id_first="ID=", uniq_id_second=";")
write_MCScanX_gff("citsi.gff3", is_type="mRNA", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff("cofca.gff3", is_type = "mRNA", first_sep = "Name=", 
    second_sep = ";", uniq_id_first = "Name=", uniq_id_second = ";")
write_MCScanX_gff("corca.gff", is_type="CDS", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff("corol.gff", is_type="CDS", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff_combine("covsu.gff3", is_type="mRNA", first_sep_first="Name=", 
    second_sep_first=";", first_sep_second="pacid=", second_sep_second=";", 
    uniq_id_first="pacid=", uniq_id_second=";", append="Csu")
write_MCScanX_gff("cucme.gff3", is_type="transcript", first_sep="ID=", 
    second_sep="T1", uniq_id_first="ID=", uniq_id_second=";", append="P1")
write_MCScanX_gff("cucsa.gff", is_type="CDS", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff("cyame.gff3", is_type="mRNA", first_sep="ID=transcript:", 
    second_sep=";", uniq_id_first="ID=", uniq_id_second=";")
write_MCScanX_gff("cyapa.gff", is_type="mRNA", first_sep="ID=", 
    second_sep=";", uniq_id_first="ID=", uniq_id_second=";", append="Cpa|", front=TRUE) 
write_MCScanX_gff("dauca.gff3", is_type="mRNA", first_sep="ID=", 
    second_sep="[.]v1", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff("denof.gff", is_type="CDS", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff("dunsa.gff3", is_type="mRNA", first_sep="ID=", 
    second_sep="[.]v1[.]0", uniq_id_first="ID=", uniq_id_second=";")
write_MCScanX_gff("ectsi.gff3", is_type="gene", first_sep="ID=", 
    second_sep=";", uniq_id_first="ID=", uniq_id_second=";")
write_MCScanX_gff("elagu.gff", is_type = "CDS", first_sep = "Name=", 
    second_sep = ";", uniq_id_first = "Name=", uniq_id_second = ";")
write_MCScanX_gff("eucgr.gff3", is_type="mRNA", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";", append=".p")
write_MCScanX_gff("eutsa.gff3", is_type="mRNA", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff("frave.gff3", is_type="mRNA", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff("genau.gff", is_type="CDS", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff("ginbi.gff", is_type="mRNA", first_sep="ID=", 
    second_sep=";", uniq_id_first="ID=", uniq_id_second=";")
write_MCScanX_gff("glyma.gff3", is_type="mRNA", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";", append=".p")
write_MCScanX_gff("glyur.gff3", is_type="mRNA", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff("gnemo.gff", is_type="mRNA", first_sep="ID=", 
    second_sep=";", uniq_id_first="ID=", uniq_id_second=";")
write_MCScanX_gff("gosra.gff3", is_type="mRNA", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff("helan.gff3", is_type="mRNA", first_sep="ID=", 
    second_sep=";", uniq_id_first="ID=", uniq_id_second=";")
write_MCScanX_gff("horvu.gff3", is_type = "mRNA", first_sep = "ID=transcript:", 
    second_sep = ";", uniq_id_first = "ID=transcript:", uniq_id_second = ";")
write_MCScanX_gff("humlu.gff3", is_type="mRNA", first_sep="ID=", 
    second_sep=";", uniq_id_first="ID=", uniq_id_second=";")
write_MCScanX_gff("iponi.gff", is_type = "CDS", first_sep = "Name=", 
    second_sep = ";", uniq_id_first = "Name=", uniq_id_second = ";")
write_MCScanX_gff("jatcu.gff", is_type="CDS", first_sep="JCDB_ID=", 
    second_sep=";", uniq_id_first="JCDB_ID=", uniq_id_second=";")
write_MCScanX_gff("jugre.gff", is_type="gene", first_sep="J", 
    second_sep=";", uniq_id_first="J", uniq_id_second=";", append="J", front=TRUE)
write_MCScanX_gff("kalfe.gff3", is_type="mRNA", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";", append=".p")
write_MCScanX_gff("kleni.gff", is_type="mRNA", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";") ## remove manually Note=... from the build kleni.gff
write_MCScanX_gff("leepe.gff3", is_type = "mRNA", first_sep = "transcript_id=", 
    second_sep = ";", uniq_id_first = "transcript_id=", uniq_id_second = ";")
write_MCScanX_gff("lepme.gff3", is_type="mRNA", first_sep="ID=", 
    second_sep=";", uniq_id_first="ID=", uniq_id_second=";")
write_MCScanX_gff("linus.gff3", is_type="mRNA", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff("lotja.gff3", is_type="mRNA", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff("lupan.gff", is_type = "CDS", first_sep = "Name=", 
    second_sep = ";", uniq_id_first = "Name=", uniq_id_second = ";")
write_MCScanX_gff("maldo.gff", is_type="CDS", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";") ## check if duplicated genes
write_MCScanX_gff("manes.gff", is_type="CDS", first_sep="Phytozome:", 
    second_sep=",", uniq_id_first="Phytozome:", uniq_id_second=",")
write_MCScanX_gff("marpo.gff3", is_type="mRNA", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff("medtr.gff", is_type="CDS", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff("mimgu.gff3", is_type = "mRNA", first_sep = "Name=", 
    second_sep = ";", uniq_id_first = "Name=", uniq_id_second = ";", append = ".p")
write_MCScanX_gff("morno.gff", is_type = "CDS", first_sep = "Name=", 
    second_sep = ";", uniq_id_first = "Name=", uniq_id_second = ";")
write_MCScanX_gff("momch.gff", is_type="CDS", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff("musac.gff3", is_type="polypeptide", first_sep="ID=", 
    second_sep=";", uniq_id_first="ID=", uniq_id_second=";")
write_MCScanX_gff("nelnu.gff", is_type = "CDS", first_sep = "Name=", 
    second_sep = ";", uniq_id_first = "Name=", uniq_id_second = ";")
write_MCScanX_gff("nicat.gff", is_type = "CDS", first_sep = "Name=", 
    second_sep = ";", uniq_id_first = "Name=", uniq_id_second = ";")
write_MCScanX_gff("nicbe.gff", is_type="mRNA", first_sep="ID=", 
    second_sep=";", uniq_id_first="ID=", uniq_id_second=";")
write_MCScanX_gff("nicsy.gff", is_type="CDS", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff("nicta.gff", is_type="mRNA", first_sep="ID=", 
    second_sep=";", uniq_id_first="ID=", uniq_id_second=";")
write_MCScanX_gff("oleeu.gff", is_type="mRNA", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff("oroth.gff3", is_type="mRNA", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff("orysa.gff3", is_type="CDS", first_sep="Parent=", 
    second_sep=".1;", uniq_id_first="pacid=", uniq_id_second=" ") 
write_MCScanX_gff("oryru.gff3", is_type = "mRNA", first_sep = "transcript_id=", 
    second_sep = ";", uniq_id_first = "transcript_id=", uniq_id_second = ";")
write_MCScanX_gff_combine("ostlu.gff3", is_type="mRNA", first_sep_first="Name=", 
    second_sep_first=";", first_sep_second="pacid=", second_sep_second=";", 
    uniq_id_first="pacid=", uniq_id_second=";", append="Olu")
write_MCScanX_gff("papso.gff", is_type="CDS", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff("petax.gff", is_type="mRNA", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";") 
write_MCScanX_gff("petin.gff", is_type="mRNA", first_sep="ID=", 
    second_sep=";", uniq_id_first="ID=", uniq_id_second=";") 
write_MCScanX_gff("phaeq.gff", is_type="CDS", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";") 
write_MCScanX_gff("phavu.gff3", is_type = "mRNA", first_sep = "transcript_id=", 
    second_sep = ";", uniq_id_first = "transcript_id=", uniq_id_second = ";")
write_MCScanX_gff("phoda.gff", is_type = "CDS", first_sep = "Name=", 
    second_sep = ";", uniq_id_first = "Name=", uniq_id_second = ";")
write_MCScanX_gff("phypa.gff3", is_type="mRNA", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";", append=".p") 
write_MCScanX_gff("picab.gff", is_type="gene", first_sep="ID=", 
    second_sep=";", uniq_id_first="ID=", uniq_id_second=";")
write_MCScanX_gff("pinta.gff", is_type="gene", first_sep="P", 
    second_sep=";", uniq_id_first="P", uniq_id_second=";", append="P", front=TRUE)
write_MCScanX_gff("poptr.gff3", is_type="mRNA", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";") 
write_MCScanX_gff("porpu.gff3", is_type="mRNA", first_sep="ID=", 
    second_sep=";", append="Popu|", uniq_id_first="ID=", uniq_id_second=";", front=TRUE)  
write_MCScanX_gff("prupe.gff3", is_type="mRNA", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";", append=".p") 
write_MCScanX_gff("pseme.gff3", is_type="mRNA", first_sep="ID=", 
    second_sep=";", uniq_id_first="ID=", uniq_id_second=";")
write_MCScanX_gff("pungr.gff", is_type="CDS",  first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff("pyrbr.gff", is_type = "mRNA", first_sep = "ID=", 
    second_sep = ";", uniq_id_first = "ID=", uniq_id_second = ";")
write_MCScanX_gff("quero.gff", is_type="mRNA", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";") ## replace Qrob_T by Qrob_P
write_MCScanX_gff("ricco.gff", is_type="mRNA", first_sep="mRNA ", 
    second_sep=";", uniq_id_first="mRNA ", uniq_id_second=";") 
write_MCScanX_gff("ruboc.gff3", is_type="gene", first_sep="ID=", 
    second_sep=";", uniq_id_first="ID=", uniq_id_second=";") 
write_MCScanX_gff("salcu.gff", is_type = "gene", first_sep = "ID=", 
    second_sep = ";", uniq_id_first = "ID=", uniq_id_second = ";")
write_MCScanX_gff("salmi.gff3",  is_type="mRNA", first_sep="ID=", 
    second_sep=";", uniq_id_first="ID=", uniq_id_second=";") 
write_MCScanX_gff_combine("selmo.gff3", is_type="mRNA", first_sep_first="Name=", 
    second_sep_first=";", first_sep_second="pacid=", second_sep_second=";", 
    uniq_id_first="pacid=", uniq_id_second=";", append="Smo") 
write_MCScanX_gff("setit.gff", is_type="mRNA", first_sep="ID=", 
    second_sep=";", uniq_id_first="ID=", uniq_id_second=";") 
write_MCScanX_gff("solpe.gff3", is_type="mRNA", first_sep="ID=", 
    second_sep=";", uniq_id_first="ID=", uniq_id_second=";") 
write_MCScanX_gff("soltu.gff", is_type="mRNA", first_sep="Parent=", 
    second_sep=";", uniq_id_first="Parent=", uniq_id_second=";") 
write_MCScanX_gff("solyc.gff", is_type="mRNA", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff("sorbi.gff3", is_type="mRNA", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";", append=".p")
write_MCScanX_gff("spipo.gff3", is_type="mRNA", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff("synec.gff", is_type="CDS", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff("tarha.gff", is_type = "CDS", first_sep = "Name=", 
    second_sep = ";", uniq_id_first = "Name=", uniq_id_second = ";")
write_MCScanX_gff("theca.gff3", is_type="mRNA", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff("thepa.gff", is_type="CDS", first_sep=NULL, 
    second_sep=NULL, uniq_id_first=NULL, uniq_id_second=NULL)
write_MCScanX_gff("triae.gff3", is_type = "mRNA", first_sep = "Name=", 
    second_sep = ";", uniq_id_first = "Name=", uniq_id_second = ";")
write_MCScanX_gff("tripr.gff3", is_type="mRNA", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff("vacco.gff3", is_type="mRNA", first_sep="ID=", 
    second_sep=";", uniq_id_first="ID=", uniq_id_second=";")
write_MCScanX_gff("vitvi.gff", is_type="mRNA", first_sep="mRNA ", 
    second_sep=" ;", uniq_id_first="mRNA ", uniq_id_second=" ;")
write_MCScanX_gff("volca.gff3", is_type="mRNA", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff("zeama.gff3", is_type="CDS", first_sep="protein_id=", 
    second_sep="_T00", uniq_id_first="protein_id", uniq_id_second="_T00")
write_MCScanX_gff("zizju.gff", is_type="CDS", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")
write_MCScanX_gff("zosma.gff3", is_type="mRNA", first_sep="Name=", 
    second_sep=";", uniq_id_first="Name=", uniq_id_second=";")

######################## 3 ##########################
############## combine gff files ####################
setwd("~/AG-Fernie/Thomas/Data/synteny/synteny_MCScanX/")
species <- list.files()[grep(list.files(), pattern = "[.]gff")]
species <- unlist(lapply(strsplit(species, split = ".gff"), "[", 1)) ## 126
species_combn <- cbind(combn(species, 2), rbind(species, species)) ## 2 8001

# new_gff <- c("cyapa")
# species_combn <- species_combn[, unique(c(which(species_combn[1,] %in% new_gff), which(species_combn[2,] %in% new_gff)))]


for (i in 1:dim(species_combn)[2]) {
    species1 <- read.table(paste0("MCScanX_gff/00_gff_files/", species_combn[1,i], ".gff"), stringsAsFactors = FALSE, quote = "", header = FALSE)
    species2 <- read.table(paste0("MCScanX_gff/00_gff_files/", species_combn[2,i], ".gff"), stringsAsFactors = FALSE, quote = "", header = FALSE)
    file_name <- paste0("MCScanX_gff/00_gff_files/", species_combn[1,i], "_", species_combn[2,i], ".gff")
    write.table(rbind(species1, species2), file = file_name, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
}
## run the above in MCScanX_gff and copy all gff files to MCScanX_gff_mcl/

## copy homology files to the folder 00_rerun
combn_lst <- apply(species_combn, 2, function(x) paste(x[1], x[2], sep = "_"))
combn_lst <- matrix(combn_lst, ncol = 1)
combn_lst_homology <- paste0("../", combn_lst, ".homology")
file.copy(from = combn_lst_homology, to = ".")

## start MCScanX in console by
ls *.homology | cut -c -11 > MCScanX_files.lst ## create MCScanX_files.lst
while read p; do .pathtoMCScanX./MCScanX_h $p -b 0; done < MCScanX_files.lst



x <- read.table("MCScanX_gff/arath_orysa.gff")
x[,1] <- gsub("arath", "at", x[,1])
x[,1] <- gsub("orysa", "os", x[,1])
x[,1] <- gsub("_Chr", "", x[,1])
x[,2] <- gsub("_", "|", x[,2])
write.table(x, file = "MCScanX_gff/arath_orysa.gff", col.names=F, row.names=F, quote=F)

y <- read.table("MCScanX_gff/arath_orysa.homology")
y[,1] <- gsub("_", "|", y[,1])
y[,2] <- gsub("_", "|", y[,2])
y <- cbind(y, rep(1, dim(y)[1]))
write.table(y, file = "MCScanX_gff/arath_orysa.homology", col.names=F, row.names=F, quote=F)
