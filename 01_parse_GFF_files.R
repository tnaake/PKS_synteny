## a script to create a file with genes and their orientations compatible with i-ADHore

setwd("H:/Documents/03_Syntenic_linkage/01_Data/OrthoFinder/input/gff_files")
setwd("W:/Thomas/Data/synteny/synteny_iadhore")

## define the functions write_file,  write_file_gff and
## write_file_combine_different_entries
write_file <- function(file, column_name = 9,  column_scaffold = 1, 
    column_orientation = 7, first_pattern = "ID=", second_pattern = ";", 
    keep_pattern_cut = "aaaaaaa", add_pattern = "", filter_take = "mRNA", 
    replace_pattern = FALSE, pattern = "abc", replacement = "abc") {
    
    species <- unlist(lapply(strsplit(file, split = ".gff"), "[", 1))
    gff <- read.table(file, sep = "\t", stringsAsFactors = FALSE, quote = "")
    
    ## filter for genes/mRNA, etc.
    gff <- gff[gff[, 3] == filter_take, ]
    name_tr <- unlist(lapply(strsplit(as.character(gff[, column_name]), first_pattern), "[", 2))
    name_tr <- unlist(lapply(strsplit(name_tr, second_pattern), "[", 1))
    # ind_remove <- grep(name_tr, pattern = exclude_pattern)
    # gff <- gff[-ind_remove,]
    # name_tr <- unlist(lapply(strsplit(gff[, column_name], first_pattern), "[", 2))
    # name_tr <- unlist(lapply(strsplit(name_tr, ";"), "[", 1)
    #ind_keep <- grep(name_tr, pattern = keep_pattern)
    ## truncate the gff that it contains only the kept features
    #gff <- gff[ind_keep,]
    ## truncate name_tr that it contains only the kept features
    #name_tr <- name_tr[ind_keep]
    name_tr <- unlist(lapply(strsplit(
        name_tr, split = keep_pattern_cut), "[", 1))
    name_tr <- paste(name_tr, add_pattern, sep = "")
    ## paste species in front of gene name
    name_tr <- paste(species, "_", name_tr, sep = "")
    
    if (replace_pattern) {
        name_tr <- gsub(x = name_tr, replacement = replacement, pattern = pattern)
    }
    
    ## read fasta file
    fasta <- read.table(file = paste("../", species, ".fasta", sep = ""), 
        header = FALSE, stringsAsFactors = FALSE)
    ind_fasta <- grep(fasta[, 1], pattern = ">")
    fasta <- fasta[ind_fasta, ]
    fasta <- as.character(gsub(">", "", fasta))
    
    ## match values in fasta with the ones in name_tr and get indices of the matching ones 
    inds <- match(name_tr, fasta)
    inds_na <- which(!is.na(inds))
    
    ## index name_tr and gff
    gff <- gff[inds_na, ]
    name_tr <- name_tr[inds_na]
    
    value_scaffold <- as.character(gff[, column_scaffold])
    unique_scaffold <- unique(value_scaffold)
    
    for (i in 1:length(unique_scaffold)) {
        inds <- which(value_scaffold == unique_scaffold[i])
        name_tr_i <- name_tr[ inds ]
        orientation_i <- gff[inds, column_orientation]
        m_i <- matrix(paste(name_tr_i, orientation_i, sep = ""), ncol = 1)
        m_i <- m_i[order(gff[inds, 4]), ]
        ## return matrix with gene + orientation
        #if (!dir.exists(species)) dir.create(species)
        file_name <- gsub("[|]", "_", unique_scaffold[i])
        write.table(m_i, quote = FALSE, 
            file = paste("./lst_files/", species, "/", file_name, ".lst", sep = ""), 
            row.names = F, col.names = F)    
    }
}

write_file_gff <- function(file, fasta_file = file, column_name = 9, 
    column_orientation = 7, column_scaffold = 1, filter_take, first_pattern, 
    second_pattern) {
    
    species <- strsplit(file, split = ".gff")[[1]][1]
    species2 <- strsplit(fasta_file, split = ".gff|.fasta")[[1]][1]
    fasta <- read.table(paste("../", species2, ".fasta", sep = ""))
    inds <- grep(fasta[, 1], pattern = ">")
    fasta <- gsub(">", "", fasta[inds, ])
    gff <- read.csv(file, comment.char = "#", header = FALSE, 
        stringsAsFactors = FALSE, sep = "\t")
    gff <- gff[gff[, 3] == filter_take, ]
    name_tr <- unlist(lapply(strsplit(
        gff[,column_name], split = first_pattern), "[", 2))
    name_tr <- unlist(lapply(strsplit(
        name_tr, split = second_pattern), "[", 1))
    
    ## paste species in front of gene name
    name_tr <- paste(species, "_", name_tr, sep = "")
    
    inds <- match(name_tr, fasta)
    inds_na <- which(!is.na(inds))
    ## index name_tr and gff
    gff <- gff[inds_na, ]
    name_tr <- name_tr[inds_na]

    
    scaffold <- gff[, column_scaffold] 
    scaffold_unique <- unique(scaffold)
    
    for (i in 1:length(scaffold_unique)) {
        inds <- which(scaffold == scaffold_unique[i])
        
        ## paste species in front of gene name, add gene gene name and orientation
        m_i <- matrix(paste(
            name_tr[inds], gff[inds, column_orientation], sep = ""), ncol = 1)
        m_i <- m_i[order(gff[inds, 4]), ]
        write.table(m_i, 
            file = paste("lst_files/", species, "/", scaffold_unique[i], ".lst", sep=""),
            col.names = FALSE, row.names = FALSE, quote = FALSE)
    }
}

write_file_combine_different_entries <- function(file, fasta_file, first_pattern1 = "Name=", 
    first_pattern2 = ";", second_pattern1 = "pacid=", second_pattern2 = ";",
    connector1 = "Csu", connector2 = "|PACid_", column_scaffold = 1, 
    column_orientation = 7, column_name = 9, filter_take = "mRNA") {
    
    species <- unlist(strsplit(file, split = ".gff"))[1]
    
    ## load gff
    gff <- read.table(file, comment.char = "#", header = FALSE, 
        stringsAsFactors = FALSE, sep = "\t")
    gff <- gff[gff[, 3] == filter_take, ]
    gff_name <- as.character(gff[, column_name])
    
    ## strsplit the patterns and combine
    pattern1 <- unlist(lapply(
        strsplit(gff_name, split = first_pattern1), "[", 2))
    pattern1 <- unlist(lapply(
        strsplit(pattern1, split = first_pattern2), "[", 1))
    pattern2 <- unlist(lapply(
        strsplit(gff_name, split = second_pattern1), "[", 2))
    pattern2 <- unlist(lapply(
        strsplit(pattern2, split = second_pattern2), "[", 1))
    combinedName <- paste(connector1, pattern1, connector2, pattern2, sep="")
    
    ## paste species in front of gene name
    combinedName <- paste(species, "_", combinedName, sep = "")
    
    ## read fasta 
    fasta <- read.table(fasta_file)
    fasta <- fasta[grep(fasta[, 1], pattern = ">"), 1]
    fasta <- gsub(">", "", fasta)
    
    ## match fasta and combinedName
    inds <- match(combinedName, fasta)
    inds_na <- which(!is.na(inds))

    ## truncata gff and combinedName
    gff <- gff[inds_na, ]
    combinedName <- combinedName[inds_na]
    
    scaffold <- gff[, column_scaffold] 
    scaffold_unique <- unique(scaffold)
    
    ## iterate through scaffold_unique
    for (i in 1:length(scaffold_unique)) {
        inds <- scaffold == scaffold_unique[i]
        ## paste species in front of gene name, add gene gene name and orientation
        m_i <- matrix(
            unique(paste(combinedName[inds], gff[inds, column_orientation], sep="")), ncol=1)
        m_i <- m_i[order(gff[inds, 4]), ]
        write.table(m_i, 
            file = paste("lst_files/", species, "/", scaffold_unique[i], ".lst", sep=""),
            col.names=FALSE, row.names=FALSE, quote = FALSE)
    }
}

## apply the functions for all species 

## actch
write_file("actch.gff3", first_pattern = "Parent=", second_pattern = ";", 
    keep_pattern_cut = "aaa")

## amahy
write_file("amahy.gff3", first_pattern = "ID=", second_pattern = ";", 
    keep_pattern_cut = ".v2.1")

## amtri
write_file("amtri.gff3", first_pattern = "ID=", second_pattern = ";", 
    keep_pattern_cut = "aaa")

## anaba ## no gff file

## anaco
write_file("anaco.gff3", first_pattern = "ID=", second_pattern = ";", 
    keep_pattern_cut = ".v3")

## aquco 
write_file("aquco.gff3", first_pattern = "ID=", second_pattern = ";", 
    keep_pattern_cut = ".v3.1", add_pattern = ".p")

## aradu
write_file_gff("aradu.gff", filter_take = "CDS", first_pattern = "Genbank:", 
    second_pattern = ";")

## araip
write_file_gff("araip.gff", filter_take = "CDS", first_pattern = "Genbank:", 
    second_pattern = ";")

## araly
write_file("araly.gff3", first_pattern = "Name=", second_pattern = ";")

## arath 
write_file("arath.gff", first_pattern = "Parent=", second_pattern = ";")

## artan
write_file("artan.gff", filter_take = "CDS", first_pattern = "Name=", 
    second_pattern = ";")

## aspof
write_file_gff("aspof.gff", first_pattern = "Name=", second_pattern = ";", 
    filter_take = "CDS")

## azofi
write_file_gff("azofi.gff", first_pattern = "ID=", second_pattern = ";", 
    filter_take = "gene")

## auran
write_file("auran.gff", filter_take = "CDS", first_pattern = "Name=", 
    second_pattern = ";")

## betpe
write_file_gff("betpe.gff", first_pattern = "ID=", second_pattern = ";", 
    filter_take = "mRNA")

## betvu
write_file("betvu.gff3", first_pattern = "ID=", second_pattern = ";", 
    keep_pattern_cut = "aaaa")

## bradi
write_file("bradi.gff3", first_pattern = "ID=", second_pattern = ";", 
    keep_pattern_cut = ".v3.1", add_pattern = ".p")

## brana
write_file("brana.gff3", first_pattern = "Alias=", second_pattern = ";", 
    filter_take = "mRNA")

## braol
write_file("braol.gff3", first_pattern = "transcript_id=", second_pattern = ";", 
    filter_take = "mRNA")

## brara
write_file("brara.gff3", first_pattern = "ID=", second_pattern = ";", 
    keep_pattern_cut = ".v1.3", add_pattern = ".p")

## cajca
write_file_gff("cajca.gff", first_pattern = "ID=", second_pattern = ";", 
    filter_take = "mRNA")

## camsa
write_file_gff("camsa.gff", first_pattern = "Name=", second_pattern = ";", 
    filter_take = "CDS")

## camsi
write_file("camsi.gff3", first_pattern = "ID=", second_pattern = ";", 
    filter_take = "mRNA")

## cansa
write_file_gff("cansa.gff", first_pattern = "Name=", second_pattern = ";", 
    filter_take = "CDS")

## capan
write_file("capan.gff", first_pattern = "=", second_pattern = ";", 
    filter_take = "mRNA")

## capgr
write_file("capgr.gff3", first_pattern = "Name=", second_pattern = ";", 
    add_pattern = ".p")

## capru
write_file("capru.gff3", first_pattern = "Name=", second_pattern = ";")

## carpa
write_file("carpa.gff3", first_pattern = "Name=", second_pattern = ";")

## chabr
write_file_gff("chabr.gff", first_pattern = "Name=", second_pattern = ";", 
    filter_take = "CDS")

## chequ
write_file("chequ.gff3", first_pattern = "Name=", second_pattern = ";", 
    filter_take = "mRNA")

## chlre 
write_file_gff("chlre.gff3", filter_take = "mRNA", first_pattern = "Name=", 
    second_pattern = ";") 

## cicar
write_file_gff("cicar.gff", first_pattern = "ID=", second_pattern = ";", 
    filter_take = "mRNA")

## citcl
write_file_gff("citcl.gff", fasta_file = "citcl.fasta", filter_take = "CDS", 
    first_pattern = "Name=", second_pattern = ";")

## citla
write_file("citla.gff", first_pattern = "ID=", second_pattern = ";")

## citsi
write_file("citsi.gff3", first_pattern = "ID=", second_pattern = ";")

## cofca
write_file("cofca.gff3", first_pattern = "Name=", second_pattern = ";", 
    filter_take = "mRNA")

## corca 
write_file("corca.gff", filter_take = "CDS", first_pattern = "Name=", 
    second_pattern = ";")

## corol
write_file("corol.gff", filter_take = "CDS", first_pattern = "Name=", 
    second_pattern = ";")

## covsu 
write_file_combine_different_entries("covsu.gff3", "../covsu.fasta", 
    first_pattern1 = "Name=", first_pattern2 = ";", second_pattern1 = "pacid=", 
    second_pattern2 = ";", connector1 = "Csu", connector2 = "|PACid_", 
    filter_take = "mRNA")

## cucme
write_file("cucme.gff3", first_pattern = "Target=", second_pattern = " ", 
    filter_take = "CDS", replace_pattern = TRUE, pattern = "T", 
    replacement = "P")

## cucsa
write_file_gff(file = "cucsa.gff", fasta_file = "cucsa.fasta", 
    filter_take = "CDS", first_pattern = "Genbank:", second_pattern = ";")

## cyame
write_file("cyame.gff3", first_pattern = "ID=transcript:", second_pattern = ";")

## cyapa
write_file_combine_different_entries("cyapa.gff", "../cyapa.fasta", 
    first_pattern1 = "ID=", first_pattern2 = ";", second_pattern1 = "ID=", 
    second_pattern2 = "e", connector1 = "Cpa|", connector2 = "", 
    filter_take = "mRNA")

## cycmi no gff

## dauca
write_file("dauca.gff3", first_pattern = "ID=", second_pattern = ";", 
    keep_pattern_cut = ".v1.0.388")

## denof
write_file_gff("denof.gff", filter_take = "CDS", first_pattern = "Genbank:", 
    second_pattern = ";")

## dunsa 
write_file("dunsa.gff3", first_pattern = "ID=", second_pattern = ";", 
    keep_pattern_cut = ".v1.0")

## ectsi
write_file("ectsi.gff3", first_pattern = "ID=", second_pattern = ";", 
    filter_take = "gene")

## elagu
write_file_gff("elagu.gff", first_pattern = "Name=", second_pattern = ";", 
    filter_take = "CDS")

## equgi no gff

## eucgr
write_file("eucgr.gff3", first_pattern = "Name=", second_pattern = ";", 
    add_pattern = ".p")

## eutsa
write_file("eutsa.gff3", first_pattern = "Name=", second_pattern = ";")

## frave
write_file("frave.gff3", first_pattern = "Name=", second_pattern = ";")

## genau
write_file_gff("genau.gff", filter_take = "CDS", first_pattern = "Name=", 
    second_pattern = ";")

## ginbi 
write_file("ginbi.gff", first_pattern = "ID=", second_pattern = ";")

## glyma
write_file("glyma.gff3", first_pattern = "Name=", second_pattern = ";", 
    add_pattern = ".p")

## glyur
write_file("glyur.gff3", first_pattern = "ID=", second_pattern = ";")

## gnemo 
write_file("gnemo.gff", first_pattern = "ID=", second_pattern = ";")

## gosra
write_file("gosra.gff3", first_pattern = "Name=", second_pattern = ";")

## helan
write_file("helan.gff3", first_pattern = "ID=", second_pattern = ";")

## horvu
write_file("horvu.gff3", first_pattern = "ID=transcript:", second_pattern = ";", 
    filter_take = "mRNA")

## humlu
write_file("humlu.gff3", first_pattern = "ID=", second_pattern = ";")

## iponi
write_file_gff("iponi.gff", first_pattern = "Name=", second_pattern = ";", 
    filter_take = "CDS")

## jatcu
write_file("jatcu.gff", filter_take = "CDS", first_pattern = "JCDB_ID=",
    second_pattern = ";")

## jugre 
write_file("jugre.gff", filter_take = "gene", first_pattern = "J", 
    second_pattern = ";", replace_pattern = TRUE, pattern = "ure", 
    replacement = "Jure")

## kalfe
write_file("kalfe.gff3", first_pattern = "Name=", second_pattern = ";", 
    add_pattern = ".p")

## kleni
write_file("kleni.gff", first_pattern = "ID=", second_pattern = ";")

## leepe
write_file("leepe.gff3", first_pattern = "transcript_id=", second_pattern = ";", 
    filter_take = "mRNA")

## lepme
write_file("lepme.gff3", first_pattern = "D=", second_pattern = ";")

## linus
write_file("linus.gff3", first_pattern = "Name=", second_pattern = ";")

## lotja
write_file("lotja.gff3", first_pattern = "ID=", second_pattern = ";")

## lupan
write_file_gff("lupan.gff", first_pattern = "Name=", second_pattern = ";", 
    filter_take = "CDS")

## maldo
write_file_gff(file = "maldo.gff", filter_take = "CDS", 
    first_pattern = "Genbank:", second_pattern = ";")

## manes
write_file("manes.gff", first_pattern = "Dbxref=Phytozome:", 
    second_pattern = ",", filter_take = "CDS")

## marpo
write_file("marpo.gff3", first_pattern = "Name=", second_pattern = ";")

## medtr
write_file_gff(file = "medtr.gff", filter_take = "CDS", 
    first_pattern = "Genbank:", second_pattern = ";")

## mimgu
write_file("mimgu.gff3", first_pattern = "Name=", second_pattern = ";", 
    filter_take = "mRNA", add_pattern=".p")

## morno
write_file_gff("morno.gff", first_pattern = "Name=", second_pattern = ";", 
    filter_take = "CDS")

## musac
write_file("musac.gff3", first_pattern = "ID=", second_pattern = ";", 
    replace_pattern = TRUE, pattern = "_t", replacement = "_p")

## momch
write_file("momch.gff", filter_take = "CDS", first_pattern = "Name=", 
    second_pattern = ";")

## nelnu
write_file_gff("nelnu.gff", first_pattern = "Name=", second_pattern = ";", 
    filter_take = "CDS")

## nicat
write_file_gff("nicat.gff", first_pattern = "Name=", second_pattern = ";", 
    filter_take = "CDS")

## nicbe
write_file("nicbe.gff", first_pattern = "ID=", second_pattern = ";")

## nicsy 
write_file_gff("nicsy.gff", fasta_file = "nicsy.fasta", filter_take = "CDS", 
    first_pattern = "Name=", second_pattern = ";")

## nicta
write_file("nicta.gff", first_pattern = "ID=", second_pattern = ";")

## oleeu
write_file("oleeu.gff", first_pattern = "Name=", second_pattern = ";")

## oroth
write_file("oroth.gff3", first_pattern = "Name=", second_pattern = ";")

## orysa
write_file("orysa.gff3", first_pattern = "Parent=", second_pattern = ".1;", 
    filter_take = "CDS")

## oryru
write_file("oryru.gff3", first_pattern = "transcript_id=", second_pattern = ";", 
    filter_take = "mRNA")

## ostlu
write_file_combine_different_entries("ostlu.gff3", "../ostlu.fasta", 
    first_pattern1 = "Name=", first_pattern2 = ";", second_pattern1 = "pacid=", 
    second_pattern2 = ";", connector1 = "Olu", connector2 = "|PACid_")

## papso 
write_file_gff("papso.gff", filter_take = "CDS", first_pattern = "Name=", 
    second_pattern = ";")

## petax
write_file("petax.gff", first_pattern = "ID=", second_pattern = ";")

## petin
write_file("petin.gff", first_pattern = "ID=", second_pattern = ";")

## phaeq
write_file_gff("phaeq.gff", filter_take = "CDS", first_pattern = "Genbank:", 
    second_pattern = ";")

## phavu 
write_file("phavu.gff3", first_pattern = "transcript_id=", second_pattern = ";", 
    filter_take = "mRNA")

## phoda
write_file_gff("phoda.gff", first_pattern = "Name=", second_pattern = ";", 
    filter_take = "CDS")

## phypa 
write_file("phypa.gff3", first_pattern = "Name=", second_pattern = ";", 
    add_pattern = ".p")

## picab 
write_file("picab.gff",  filter_take = "gene", first_pattern = "ID=", 
    second_pattern = " ")

## pinpi no gff

## pinsy no gff

## pinta 
write_file("pinta.gff", filter_take = "gene", first_pattern = "P", 
    second_pattern = ";", replace_pattern = TRUE, pattern = "ITA", 
    replacement = "PITA")

## poptr
write_file("poptr.gff3", first_pattern = "Name=", second_pattern = ";")

## porpu 
write_file_combine_different_entries("porpu.gff3", "../porpu.fasta", 
    first_pattern1 = "ID=", first_pattern2 = ";", second_pattern1 = "ID=", 
    second_pattern2 = "e", connector1 = "Popu|", connector2 = "")

## prupe 
write_file("prupe.gff3", first_pattern = "Name=", second_pattern = ";", 
    add_pattern = ".p")

## pseme 
write_file("pseme.gff3", first_pattern = "ID=", second_pattern = ";")

## pungr
write_file("pungr.gff", filter_take = "CDS", first_pattern = "Name=", 
    second_pattern = ";")

## pyrbr
write_file_gff("pyrbr.gff", first_pattern = "ID=", second_pattern = ";", 
    filter_take = "mRNA")

## quero
write_file("quero.gff", first_pattern = "Name=", second_pattern = ";", 
    replace_pattern = TRUE, pattern = "T", replacement = "P")

## ricco
write_file("ricco.gff", first_pattern = "mRNA ", second_pattern = "; ")

## ruboc
write_file("ruboc.gff3", first_pattern = "ID=", second_pattern = ";", 
    filter_take = "gene")

## salcu
write_file_gff("salcu.gff", first_pattern = "ID=", second_pattern = ";", 
    filter_take = "gene")

## salmi
write_file("salmi.gff3", first_pattern = "ID=", second_pattern = ";", 
    add_pattern = "")

## selmo
write_file_combine_different_entries("selmo.gff3", "../selmo.fasta", 
    first_pattern1 = "Name=", first_pattern2 = ";", second_pattern1 = "pacid=", 
    second_pattern2 = ";", connector1 = "Smo", connector2 = "|PACid_")

## setit
write_file("setit.gff", first_pattern = "ID=", second_pattern = ";")

## solpe
write_file("solpe.gff3", first_pattern = "ID=", second_pattern = ";")

## soltu 
write_file("soltu.gff", first_pattern = "Parent=", second_pattern = ";")

## solyc
write_file("solyc.gff", first_pattern = "ID=mRNA:", second_pattern = ";")

## sorbi
write_file("sorbi.gff3", first_pattern = "Name=", second_pattern = ";", 
    add_pattern = ".p")

## spipo
write_file("spipo.gff3", first_pattern = "Name=", second_pattern = ";")

## synec
write_file("synec.gff", filter_take = "CDS", first_pattern = "Name=", 
    second_pattern = ";")

## tarha
write_file_gff("tarha.gff", first_pattern = "Name=", second_pattern = ";", 
    filter_take = "CDS")

## taxba no gff

## theca
write_file("theca.gff3", first_pattern = "Name=", second_pattern = ";")

## thepa
write_file_combine_different_entries("thepa.gff", "../thepa.fasta", 
    first_pattern1 = "Tp", first_pattern2 = ";", second_pattern1 = "T", 
    second_pattern2 = "p", connector1 = "Tp", connector2 = "")

## triae
write_file("triae.gff3", first_pattern = "Name=", second_pattern = ";", 
    filter_take = "mRNA")

## tripr
write_file("tripr.gff3", first_pattern = "Name=", second_pattern = ";")

## vacco
write_file_gff("vacco.gff3", first_pattern = "ID=", second_pattern = ";", 
    filter_take = "mRNA")

## vitvi
write_file("vitvi.gff", first_pattern = "mRNA ", second_pattern = " ;")

## volca
write_file("volca.gff3", first_pattern = "Name=", second_pattern = ";")

## zeama
write_file("zeama.gff3", first_pattern = "transcript:", second_pattern = "_", 
    filter_take = "CDS")

## zizju
write_file_gff("zizju.gff", filter_take = "CDS", first_pattern = "Genbank:",
    second_pattern = ";")

## zosma
write_file("zosma.gff3", first_pattern = "Name=", second_pattern = ";")
