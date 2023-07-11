library(here)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot(font_size = 12))
library(Biostrings)
library(ggseqlogo)
library(jsonlite)
library(stringr)

generate_summary_df <- function(JSON_df) {
    #extract columns and make tidy df
    summary_df <- data.frame(geneID = character(), 
                        protein_type = character(), 
                        other_p = numeric(), 
                        SP_p = numeric(),
                        site = character())
    for (i in seq_len(ncol(df))){
        geneID <- as.character(df[[1,i]])
        protein_type <- as.character(df[[2,i]])
        p_values <- df[[3,i]][[1]]
        other_p <- p_values[1]
        SP_p <- p_values[2]
        site <- as.character(df[[4,i]])
        summary_df <- rbind(summary_df, data.frame(geneID, protein_type, other_p, SP_p, site))
    }

    return(summary_df)
}

json <- jsonlite::fromJSON(paste(here::here(), 
                        "motif_conservation",
                        "data", 
                        "proteins",
                        "Asp_N.json",
                        sep= "/"))
df <- as_tibble(json$SEQUENCES)
df <- df[c(1,5,6,8),]
Asp_N_summary_df <- generate_summary_df(df)

json <- jsonlite::fromJSON(paste(here::here(), 
                        "motif_conservation",
                        "data", 
                        "proteins",
                        "Asp_F.json",
                        sep= "/"))
df <- as_tibble(json$SEQUENCES)
df <- df[c(1,5,6,8),]
Asp_F_summary_df <- generate_summary_df(df)


Asp_N_summary_df <- read.table(paste(here::here(), 
                        "motif_conservation",
                        "data", 
                        "proteins",
                        "results",
                        "prediction_results.txt",
                        sep= "/"), header = TRUE, sep = "\t")
