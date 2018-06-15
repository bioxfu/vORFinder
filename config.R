setwd("/srv/shiny-server/")
virus_id     <- paste0(getwd(),"/database/virus.id")
virus_db     <- paste0(getwd(),"/database/virus_db.fasta")
virus_msa_db <- paste0(getwd(),"/database/virus_msa.fasta")
orf_db       <- paste0(getwd(),"/database/orf_db.fasta")

mafft_path   <- "mafft"
blastdb_path <- paste0(getwd(),"/software/makeblastdb")
blastp_path  <- paste0(getwd(),"/software/blastp")
orf_path     <- paste0(getwd(),"/software/ORFfinder")

