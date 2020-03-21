setwd("Documents/03_Syntenic_linkage/01_Data/OrthoFinder/input/gff_files")
setwd("W:/Thomas/Data/synteny/synteny_iadhore/")

species <- unlist(lapply(strsplit(
    list.files()[grep(list.files(), pattern=".gff")], split=".gff"), "[", 1)) 
length(species) ## 126
##species <- c("arath", "zeama", "orysa", "theca")

## get all combinations
combinations <- combn(species, 2)
combinations <- cbind(combinations, rbind(species, species)) 
dim(combinations) ## 2 8001


setwd("W:/Thomas/Data/synteny/synteny_iadhore/")

type <- "MCL" ## or "MCL"

if (type == "OrthoFinder") {
    setwd("lst_files_orthofinder/") 
}
if (type == "MCL") {
    setwd("lst_files_mcl/")
}

## iterate through combinations and write ini file
for (i in 1:dim(combinations)[2]) {
    
    if (type == "OrthoFinder") {
        genome_path <- "./lst_files_orthofinder"
        blast_table  <- "blast_table= ./family_orthogroups_iadhore.txt"
        output_path <- "./output/output"
        final_path <- "../ini_files/"
    }
    if (type == "MCL") {
        genome_path <- "./lst_files_mcl"
        blast_table <- "blast_table= ./family_orthogroups_mcl_iadhore.txt"
        output_path <- "./output_mcl/output"
        final_path <- "../ini_files_mcl/"
    }
    
    
    final <- rbind(
        matrix(
            c(paste("genome=", combinations[1, i], sep = " "),
            paste(
                gsub(list.files(combinations[1, i]), pattern = ".lst", 
                    replacement = ""), 
                paste(genome_path, combinations[1, i], 
                    list.files(combinations[1, i]), sep = "/"), sep = " "), 
                ""
            )),
        matrix(
            c(paste("genome=", combinations[2, i], sep = " "),
            paste(
                gsub(list.files(combinations[2, i]), pattern = ".lst", 
                    replacement = ""), 
                paste(genome_path, combinations[2, i], 
                    list.files(combinations[2,i]), sep = "/"), sep = " "), 
                ""
            )),
        matrix(
            c(blast_table, 
            "",
            "table_type= family",
            "cluster_type= collinear",
            paste("output_path=", paste(output_path, combinations[1, i], 
                combinations[2, i], sep = "_"), sep = " "), 
            "", 
            "alignment_method=gg2",
            "gap_size=15", 
            "cluster_gap=20", 
            "max_gaps_in_alignment=20",
            "q_value=0.9",
            "prob_cutoff=0.001",
            "anchor_points=5",
            "level_2_only=true",
            "write_stats=true",
            "number_of_threads=4"
            ))
    )
    write.table(final, 
        file = paste(final_path, "iadhore_", 
            paste(combinations[1 ,i], combinations[2, i], sep = "_"), 
                ".ini", sep = ""), 
            quote = FALSE, col.names = FALSE, row.names = FALSE) 
}
