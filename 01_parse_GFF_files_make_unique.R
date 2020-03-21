## the .lst files do not have unique sequence, will cause errors in the 
## following analysis --> make them unique 
setwd("~/AG-Fernie/Thomas/Data/synteny/Results_Nov05/")
groups_of <- read.table("family_orthogroups_iadhore.txt") 
groups_mcl <- read.table("family_orthogroups_mcl_iadhore.txt") 

make_unique <- function(path, new_path, groups) {
    setwd(path)
    file_list <- list.files()[grep(list.files(), pattern = ".lst")]
    new_path <- gsub("lst_files", new_path, path)
    
    ## create directory
    dir.create(file.path(new_path), showWarnings = FALSE)
    
    #setwd(path)
    for (i in 1:length(file_list)) {
        
        x <- read.table(file_list[i], stringsAsFactors = FALSE)
        x <- unique(x[, 1])
        x_mod <- substr(x, 1, nchar(x) - 1) ## remove +/-
        x <- x[x_mod %in% groups]
        
        if (length(x) > 0) {
            write_file_path <- paste(new_path, file_list[i], sep = "/")
            write.table(x, write_file_path, row.names = FALSE, 
                col.names = FALSE, quote = FALSE)
        }
    }
}
setwd("~/AG-Fernie/Thomas/Data/synteny/synteny_iadhore/lst_files/")
species <- list.files("..")[grep(list.files(".."), pattern = "[.]gff")]
##species <- c("camsi", "chlre", "covsu", "cyame", "dunsa", "ginbi", 
##      "nicsy", "picab", "pinta", "porpu", "pseme", "salmi", "soltu", "volca")
species <- lapply(strsplit(species, split = ".gff"), "[", 1)
species <- unlist(species)

species_path <- paste0(getwd(), "/", species)

## apply the functions
for (i in 1:length(species_path)) {
    make_unique(species_path[i], "lst_files_orthofinder", groups_of[, 1])
}
for (i in 1:length(species_path)) {
    make_unique(species_path[i], "lst_files_mcl", groups_mcl[, 1])
}

