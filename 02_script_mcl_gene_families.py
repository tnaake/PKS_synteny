"""
script to get gene families based on MCL 
"""

## This part is taken and modified from orthology.py from LSTrAP/pipeline
"""
Runs MCL clustering on OrthoFinder output to obtain homologous families (without re-running blast)
"""
## run with python3


__template = """#!/bin/bash
#
#$ -N %s
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -o OUT_$JOB_NAME.$JOB_ID
#$ -e ERR_$JOB_NAME.$JOB_ID
#email
%s
#
%s
date
hostname
%s
date
"""


def build_template(name, module, cmd):
    """
    Generates submit script for a normal job.
    :param name: Name for the job
    :param email: Email address of the user, set to None to disable email
    :param module: Module to load, separate multiple modules using spaces in case more than one module is required
    :param cmd: The command to execute, separate multiple commands using newlines
    :return: The completed template
    """
    include_email = "" ## if email is None else "#$ -m bea\n#$ -M " + email
    load_module = "" if module is None else "module load " + module
    return __template % (name, include_email, load_module, cmd)




def write_submission_script(jobname, module, command, filename): ## (self, ...)
	"""
	Writes a job submission script that includes a timestamp, required to keep track if a job is running or not
	param jobname: name of the job include %d for the timestamp !
	param module: Module to load, separate multiple modules using spaces in case more than one module is required
	:param command: The command to execute, separate multiple commands using newlines
	:param filename: filename for the script include %d for the timestamp !
	:return: tuple with stamped_filename and stamped_jobname
	"""
    timestamp = int(time.time())
    stamped_filename = str(filename % timestamp)
    stamped_jobname = str(jobname % timestamp)
    template = build_template(stamped_jobname, module, command)
    with open(stamped_filename, "w") as f:
        print >> f, template #print(template, file=f)
    return stamped_filename, stamped_jobname




#orthofinder_dir = self.dp['GLOBAL']['orthofinder_output'] ## set path to orthofinder_output
#os.chdir("/scratch2/evorepro/OrthoFinder/input/")
import os
import subprocess
import time
import sys



os.chdir("/home/naake/03_synteny/new_fasta/")
orthofinder_dir = os.getcwd()
orthofinder_results_dir = list(filter(lambda x: 'Results_' in x, os.listdir(orthofinder_dir)))[0]

# Concatenate OrthoFinder blast files
working_dir = os.path.join(orthofinder_dir, orthofinder_results_dir, 'WorkingDirectory')
orthofinder_blast_files = list(filter(lambda x: x.startswith('Blast'), os.listdir(working_dir)))
full_blast = os.path.join(working_dir, 'full_blast.out')
full_blast_abc = os.path.join(working_dir, 'full_blast.abc')
mcl_families_out = os.path.join(orthofinder_dir, 'mcl_families.unprocessed.txt')

#with open(full_blast, 'w') as outfile:
#	for fname in orthofinder_blast_files:
#		with open(os.path.join(working_dir, fname), encoding="utf8", errors='ignore') as infile:
#			for line in infile:
#				outfile.write(line)
                
import gzip
with open(full_blast, 'w') as outfile:
	for fname in orthofinder_blast_files:
		with gzip.open(os.path.join(working_dir, fname), 'rt') as infile:
			for line in infile:
				outfile.write(line)




filename, jobname = write_submission_script("mcl_%d", 'biotools/mcl-14.137', ## mcl_module=biotools/mcl-14.137
              'mcxdeblast_cmd=perl /apps/biotools/mcl-14.137/bin/mcxdeblast --m9 --line-mode=abc ${blast_in} > ${abc_out}' + '\n' + ## mcxdeblast_cmd=perl /apps/biotools/mcl-14.137/bin/mcxdeblast --m9 --line-mode=abc ${blast_in} > ${abc_out}
              'mcl ${in} --abc -o ${out} -te 4', "mcl_%d.sh") ## mcl_cmd=mcl ${in} --abc -o ${out} -te 4

filename = os.getcwd() + "/" + filename
## run this in the console
# submit job
#command = ["qsub"] + [""] + ["-v " + "blast_in=" + full_blast + ",abc_out=" + full_blast_abc + ",in=" + full_blast_abc + ",out=" + mcl_families_out + " " + filename] ## ["qsub"] + self.qsub_mcxdeblast + "-v"\
blast_in=/home/naake/03_synteny/Results_Oct29/WorkingDirectory/full_blast.out
abc_out=/home/naake/03_synteny/Results_Oct29/WorkingDirectory/full_blast.abc
in=/home/naake/03_synteny/Results_Oct29/WorkingDirectory/full_blast.abc
out=/home/naake/03_synteny/mcl_families.unprocessed.txt

# run two jobs

##module load biotools/mcl-14.137
# copy mcxdeblat in current working directory (which mcxdeblast) and edit first line to the correct path of perl (find out by which perl) (/usr/bin/perl)
## old (next two lines
##./mcxdeblast /apps/biotools/mcl-14.137/bin/mcxdeblast --m9 --line-mode=abc ${blast_in} > ${abc_out}
## mcl ${in} --abc -o ${out} -te 4

cut -f 1,2,11 full_blast.out > seq.abc
## check for comments and filter using grep
mcxload -abc seq.abc --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o seq.mci -write-tab seq.tab
#mcl seq.mci -I 2
#mcxdump -icl out.seq.mci.I20 -tabr seq.tab -o dump.seq.mci.I20
mcl seq.mci -I 2 -use-tab seq.tab -te 10 ## -I 2 is also default
mcl seq.mci -I 1.5 -use-tab seq.tab -te 10 ## -I 2 is also default
cp out.seq.mci.I20 Results_Oct29/WorkingDirectory/mcl_families2.unprocessed.txt
cp out.seq.mci.I15 Results_Oct29/WorkingDirectory/mcl_families15.unprocessed.txt

cd Results_Oct26/WorkingDirectory/
cp out.seq.mci.I20 mcl_families.unprocessed.txt

# move the file to working_dir
## end: run this in the console 

## run in Python 2.7
## convert ids
id_conversion = {}
with open(os.path.join(working_dir, 'SequenceIDs.txt')) as infile:
    for line in infile:
        parts = line.strip().split()
        id = parts[0].strip(':')
        gene = parts[1]
        id_conversion[id] = gene


with open(mcl_families_out, 'r') as infile, open(os.path.join(orthofinder_dir, 'mcl_families.processed.txt'), 'w') as outfile:
    for l in infile:
        parts = [id_conversion[id] if id in id_conversion.keys() else '!error!' for id in l.strip().split()]
        print >> outfile, '\t'.join(parts)

# remove the submission script
os.remove(filename)

# remove OUT_ files
PipelineBase.clean_out_files(jobname)

print("Done\n\n")


## id conversion using R 
R:
setwd("~/03_synteny/Results_Oct29/WorkingDirectory/mcl")
ids <- read.table("../SequenceIDs.txt", sep = " ", stringsAsFactors = FALSE)
ids[, 1] <- gsub(x = ids[, 1], pattern = ":", replacement = "")
rownames(ids) <- ids[, 1]

I <- 2.0
if (I == 1.5) {
	n <- 61643 ## row number where there is last pair of orthologous genes
	mcl1 <- read.table("out.seq.mci.I15", fill = TRUE, sep = "\t", stringsAsFactors = FALSE, 
		comment.char="", colClasses = "character", nrows = 1)
	mcl2 <- read.table("out.seq.mci.I15", fill = TRUE, sep = "\t", stringsAsFactors = FALSE, 
		comment.char="", colClasses = "character", nrows = 2 - 1, skip = 1)
	mcl3 <- read.table("out.seq.mci.I15", fill = TRUE, sep = "\t", stringsAsFactors = FALSE, 
		comment.char="", colClasses = "character", nrows = 3 - 2, skip = 2)
	mcl4 <- read.table("out.seq.mci.I15", fill = TRUE, sep = "\t", stringsAsFactors = FALSE, 
		comment.char="", colClasses = "character", nrows = 4 - 3, skip = 3)
	mcl5 <- read.table("out.seq.mci.I15", fill = TRUE, sep = "\t", stringsAsFactors = FALSE, 
		comment.char="", colClasses = "character", nrows = 5 - 4, skip = 4)
	mcl6 <- read.table("out.seq.mci.I15", fill = TRUE, sep = "\t", stringsAsFactors = FALSE, 
		comment.char="", colClasses = "character", nrows = 6 - 5, skip = 5)
	mcl7 <- read.table("out.seq.mci.I15", fill = TRUE, sep = "\t", stringsAsFactors = FALSE, 
		comment.char="", colClasses = "character", nrows = 7 - 6, skip = 6)
	mcl8 <- read.table("out.seq.mci.I15", fill = TRUE, sep = "\t", stringsAsFactors = FALSE, 
		comment.char="", colClasses = "character", nrows = 8 - 7, skip = 7)
	mcl9 <- read.table("out.seq.mci.I15", fill = TRUE, sep = "\t", stringsAsFactors = FALSE, 
		comment.char="", colClasses = "character", nrows = 9 - 8, skip = 8)
	mcl10 <- read.table("out.seq.mci.I15", fill = TRUE, sep = "\t", stringsAsFactors = FALSE, 
		comment.char="", colClasses = "character", nrows = 10 - 9, skip = 9)	
	mcl <- read.table("out.seq.mci.I15", fill = TRUE, sep = "\t", stringsAsFactors = FALSE, 
		comment.char="", colClasses = "character", nrows = n - 10, skip = 10)
	unassigned <- read.table("out.seq.mci.I15", fill = TRUE, sep = "\t", stringsAsFactors = FALSE, 
		comment.char="", colClasses = "character", skip = n)
}

if (I == 2.0) {
	n <- 74059 ## row number where there is last pair of orthologous genes
	mcl1 <- read.table("out.seq.mci.I20", fill = TRUE, sep = "\t", stringsAsFactors = FALSE, 
		comment.char="", colClasses = "character", nrows = 1)
	mcl2 <- read.table("out.seq.mci.I20", fill = TRUE, sep = "\t", stringsAsFactors = FALSE, 
		comment.char="", colClasses = "character", nrows = 2 - 1, skip = 1)
	mcl3 <- read.table("out.seq.mci.I20", fill = TRUE, sep = "\t", stringsAsFactors = FALSE, 
		comment.char="", colClasses = "character", nrows = 3 - 2, skip = 2)
	mcl4 <- read.table("out.seq.mci.I20", fill = TRUE, sep = "\t", stringsAsFactors = FALSE, 
		comment.char="", colClasses = "character", nrows = 4 - 3, skip = 3)
	mcl5 <- read.table("out.seq.mci.I20", fill = TRUE, sep = "\t", stringsAsFactors = FALSE, 
		comment.char="", colClasses = "character", nrows = 5 - 4, skip = 4)
	mcl6 <- read.table("out.seq.mci.I20", fill = TRUE, sep = "\t", stringsAsFactors = FALSE, 
		comment.char="", colClasses = "character", nrows = 6 - 5, skip = 5)
	mcl7 <- read.table("out.seq.mci.I20", fill = TRUE, sep = "\t", stringsAsFactors = FALSE, 
		comment.char="", colClasses = "character", nrows = 7 - 6, skip = 6)
	mcl8 <- read.table("out.seq.mci.I20", fill = TRUE, sep = "\t", stringsAsFactors = FALSE, 
		comment.char="", colClasses = "character", nrows = 8 - 7, skip = 7)
	mcl9 <- read.table("out.seq.mci.I20", fill = TRUE, sep = "\t", stringsAsFactors = FALSE, 
		comment.char="", colClasses = "character", nrows = 9 - 8, skip = 8)
	mcl10 <- read.table("out.seq.mci.I20", fill = TRUE, sep = "\t", stringsAsFactors = FALSE, 
		comment.char="", colClasses = "character", nrows = 10 - 9, skip = 9)
	mcl <- read.table("out.seq.mci.I20", fill = TRUE, sep = "\t", stringsAsFactors = FALSE, 
		comment.char="", colClasses = "character", nrows = n - 10, skip = 10)
	unassigned <- read.table("out.seq.mci.I20", fill = TRUE, sep = "\t", stringsAsFactors = FALSE, 
		comment.char="", colClasses = "character", skip = n)
}

ids_1 <- ids[, 1]
ids_2 <- ids[, 2]

## match mcl1
inds <- match(mcl1[1, ], ids_1)
mcl1[1, ] <- ids_2[inds] 

## match mcl2
inds <- match(mcl2[1, ], ids_1)
mcl2[1, ] <- ids_2[inds]

## match mcl3
inds <- match(mcl3[1, ], ids_1)
mcl3[1, ] <- ids_2[inds]

## match mcl4
inds <- match(mcl4[1, ], ids_1)
mcl4[1, ] <- ids_2[inds]

## match mcl5
inds <- match(mcl5[1, ], ids_1)
mcl5[1, ] <- ids_2[inds]

## match mcl6
inds <- match(mcl6[1, ], ids_1)
mcl6[1, ] <- ids_2[inds]

## match mcl7
inds <- match(mcl7[1, ], ids_1)
mcl7[1, ] <- ids_2[inds]

## match mcl8
inds <- match(mcl8[1, ], ids_1)
mcl8[1, ] <- ids_2[inds]

## match mcl9
inds <- match(mcl9[1, ], ids_1)
mcl9[1, ] <- ids_2[inds]

## match mcl10
inds <- match(mcl10[1, ], ids_1)
mcl10[1, ] <- ids_2[inds]

if (I == 1.5) {
	write.table(mcl1, file = "mcl_families.processed15.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
	write.table(mcl2, file = "mcl_families.processed15.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
	write.table(mcl3, file = "mcl_families.processed15.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
	write.table(mcl4, file = "mcl_families.processed15.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
	write.table(mcl5, file = "mcl_families.processed15.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
	write.table(mcl6, file = "mcl_families.processed15.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
	write.table(mcl7, file = "mcl_families.processed15.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
	write.table(mcl8, file = "mcl_families.processed15.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
	write.table(mcl9, file = "mcl_families.processed15.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
	write.table(mcl10, file = "mcl_families.processed15.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
}

if (I == 2.0) {
	write.table(mcl1, file = "mcl_families.processed20.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
	write.table(mcl2, file = "mcl_families.processed20.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
	write.table(mcl3, file = "mcl_families.processed20.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
	write.table(mcl4, file = "mcl_families.processed20.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
	write.table(mcl5, file = "mcl_families.processed20.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
	write.table(mcl6, file = "mcl_families.processed20.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
	write.table(mcl7, file = "mcl_families.processed20.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
	write.table(mcl8, file = "mcl_families.processed20.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
	write.table(mcl9, file = "mcl_families.processed20.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
	write.table(mcl10, file = "mcl_families.processed20.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
}

rm(mcl1)
rm(mcl2)
rm(mcl3)
rm(mcl4)
rm(mcl5)
rm(mcl6)
rm(mcl7)
rm(mcl8)
rm(mcl9)
rm(mcl10)

if (I == 1.5) {
	nrow_i3 <- 27047 ## when 4 to 3 orthologous pairs (ind of last row with 4)
	nrow_i2 <- 36094 ## when 3 to 2
	
	for (i in 1:nrow_i3) { 
		print(i)
		## match mcl[i, ] with ids_1, write corresponding ids_2 to mcl[i, ] 
		if (i >= 1000) {
			mcl_i <- mcl[i, 1:720]
		} else {
			if (i < 10) {
				mcl_i <- mcl[i, ]
			}
			if (i >= 10 & i < 100) {
				mcl_i <- mcl[i, 1:10200]
			}
			if (i >= 100 & i < 1000) {
				mcl_i <- mcl[i, 1:3700]
			}
		}	
		mcl_i <- mcl_i[mcl_i != ""]
		inds <- match(mcl_i, ids_1)
		mcl[i, 1:length(mcl_i)] <- ids_2[inds]
	}
}
	
if (I == 2.0) { 
	nrow_i3 <- 32555 ## when 4 to 3
	nrow_i2 <- 44001 ## when 3 to 2

	for (i in 1:nrow_i3) { 
		print(i)
		## match mcl[i, ] with ids_1, write corresponding ids_2 to mcl[i, ] 
		if (i >= 1000) {
			mcl_i <- mcl[i, 1:670]
		} else {
			if (i < 10) {
				mcl_i <- mcl[i, ]
			}
			if (i >= 10 & i < 100) {
				mcl_i <- mcl[i, 1:9300]
			}
			if (i >= 100 & i < 1000) {
				mcl_i <- mcl[i, 1:3600]
			}
		}	
		mcl_i <- mcl_i[mcl_i != ""]
		inds <- match(mcl_i, ids_1)
		mcl[i, 1:length(mcl_i)] <- ids_2[inds]
	}
}
## match orthologous triplets (third column)
inds_row <- ((nrow_i3 + 1):nrow_i2)
inds <- match(mcl[inds_row, 3], ids_1)
mcl[inds_row, 3] <- ids_2[inds]

## match orthologous pairs and triplets (first two columns
inds_row <- (nrow_i3 + 1):nrow(mcl)
inds <- match(mcl[inds_row, 1], ids_1)
mcl[inds_row, 1] <- ids_2[inds]
inds <- match(mcl[inds_row, 2], ids_1)
mcl[inds_row, 2] <- ids_2[inds]

if (I == 1.5) {
	write.table(mcl, file = "mcl_families.processed15.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
}

if (I == 2.0) {
	write.table(mcl, file = "mcl_families.processed20.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
}

## match unassigned
inds <- match(unassigned[, 1], ids_1)
unassigned[, 1] <- ids_2[inds]

#mcl <- rbind(mcl, unassigned)

if (I == 1.5) {
	write.table(unassigned, file = "mcl_families.processed15.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
}
if (I == 2.0) {
	write.table(unassigned, file = "mcl_families.processed20.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
}