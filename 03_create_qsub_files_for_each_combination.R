setwd("/home/mpimp-golm.mpg.de/naake/winhome/Documents/03_Syntenic_linkage/01_Data/OrthoFinder/input/gff_files")
setwd("W:/Thomas/Data/synteny/synteny_iadhore")

species <- unlist(lapply(strsplit(
    list.files()[grep(list.files(), pattern = ".gff")], split = ".gff"), "[", 1)) ## 126 species
##species <- c("arath", "zeama", "orysa", "theca")

## gett all combinations
combinations <- combn(species, 2)
combinations <- cbind(combinations, rbind(species, species)) 
dim(combinations) ## 2 8001
combinations <- apply(combinations, 2, function(x) paste(x[1], x[2], sep = "_"))

write.table(combinations, file = "iadhore_list.txt", quote = FALSE, 
    col.names = FALSE, row.names = FALSE)

## iterate through combinations and write qsub file 
for (i in 1:dim(combinations)[2]) {
    final <- rbind(
        matrix(
           "#!/bin/sh
#
## --------------------------------------------------
## Please note:
## Lines beginning with '#$'
## are not commented out
## Lines beginning with '##'
## are...
## Name of the queue the job is supposed to run in:
## right now there are 'regular' and 'first'.
## If you leave it blank, any queue is used
#$ -q regular
## --------------------------------------------------
## Name of the host/node
## If you want to run the job on a certain node:
#$ -l h=rock11
## --------------------------------------------------
## Execute the job from the current working directory:
#$ -cwd
## Standard error  stream  of the job is merged into
## the standard output stream:
#$ -j y
## --------------------------------------------------
## Specifies the interpreting shell for the job:
#$ -S /bin/sh
## --------------------------------------------------
## Name our job as you wish: "),
    matrix(paste("#$ -N", paste("i-adhore", combinations[1, i], combinations[2, i], sep = "_"))),
    matrix("## --------------------------------------------------
## Put standard output of your job in this file:
#$ -o OUT_$JOB_NAME.$JOB_ID
## Put standard error output of your job in this file:
#$ -e ERR_$JOB_NAME.$JOB_ID
## --------------------------------------------------
## E-Mail options:
## Specify mail address to send job mails to
## 'b' job has begun
## 'e' job has ended
## 'a' job has got aborted or rescheduled
#$ -m a
## --------------------------------------------------
## Your E-Mail address
#$ -M naake@mpimp-golm.mpg.de
## --------------------------------------------------
## Here starts the 'real' stuff:
## --------------------------------------------------
## write the starting date and time
## into the output file
date
## write the execution host
## into the output file
hostname
## load the correct module environment - in this case
## for ncbi
module load biotools/i-adhore-3.0.01
## finally start the actual work"),
        matrix(
            paste("i-adhore ", paste("../ini_files/iadhore", combinations[1, i], combinations[2, i], sep = "_"), ".ini", sep = "")
        ), ## ini_files_mcl/iadhore
        matrix("## write the end date and time
## into the outputfile
date")
)
    write.table(final,
        file = paste("qsub/qsub_iadhore_", 
            paste(combinations[1, i], combinations[2, i], sep = "_"), ".sh", sep = ""),
        quote = FALSE, col.names = FALSE, row.names = FALSE) ## adjust qsub_mcl/qsub_mcl_iadhore
}

## switch to golem and copy the relevant files (*.lst, *.ini, qsub_iadhore_*.sh)
## run in the qsub directory
## qsub qsub_iadhore.sh
## qsub qsub_iadhore_mcl.sh

## for i in ls qsub_iadhore_*.sh
## do 
##   qsub $i
## done 