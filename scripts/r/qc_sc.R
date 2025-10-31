# Install packages not available on conda
library(devtools)
devtools::install_github("JLSteenwyk/ggpubfigs")

# Load libraries

library(dplyr)
library(ggplot2)
library(ggpubfigs)
library(glue)
library(tibble)
library(tidyr)

# Setup
sample_ids <- unlist(snakemake@params[["sample_ids"]])

# Initialize list to contain dataframes
on_aggr_bin_list <- vector(mode="list", length=length(sample_ids))
names(on_aggr_bin_list) <- sample_ids
off_aggr_bin_list <- vector(mode="list", length=length(sample_ids))
names(off_aggr_bin_list) <- sample_ids

on_aggr_cell_list <- vector(mode="list", length=length(sample_ids))
names(on_aggr_cell_list) <- sample_ids
off_aggr_cell_list <- vector(mode="list", length=length(sample_ids))
names(off_aggr_cell_list) <- sample_ids

# Fill lists
loop_list <- list(unlist(snakemake@input[["norm_depth_aggr_bin"]]), unlist(snakemake@input[["norm_depth_aggr_cell"]]))
idx_vect <- c(1:length(loop_list[[1]]))
on_ct <- 0
off_ct <- 0
for (i in idx_vect) {
    aggr_bin <- read.csv(loop_list[[1]][i],
                         sep="\t") %>%
      tibble::column_to_rownames("bin")

    if(grepl("on_target", loop_list[[1]][i])) {
        on_ct <- on_ct+1
        aggr_cell <- read.csv(loop_list[[2]][i],
                              sep="\t") %>%
            dplyr::mutate(Sample=sample_ids[on_ct])
        on_aggr_bin_list[[sample_ids[on_ct]]] <- aggr_bin
        on_aggr_cell_list[[sample_ids[on_ct]]] <- aggr_cell
    } else {
        off_ct <- off_ct+1
        aggr_cell <- read.csv(loop_list[[2]][i],
                              sep="\t") %>%
            dplyr::mutate(Sample=sample_ids[off_ct])
        off_aggr_bin_list[[sample_ids[off_ct]]] <- aggr_bin
        off_aggr_cell_list[[sample_ids[off_ct]]] <- aggr_cell
    }
}

# Create on and off target dataframes
## Create bin dataframes
aggr_bin <- list(
    "on_target" = Reduce(cbind, on_aggr_bin_list),
    "off_target" = Reduce(cbind, off_aggr_bin_list)
)
colnames(aggr_bin[["on_target"]]) <- sample_ids
colnames(aggr_bin[["off_target"]]) <- sample_ids

## Create chr dataframes and bring bin dataframes to long format
aggr_chr <- vector(mode="list", length=length(unlist(snakemake@params[["chr_regions"]])))
names(aggr_chr) <- unlist(snakemake@params[["chr_regions"]])
for(region in unlist(snakemake@params[["chr_regions"]])) {
    aggr_chr[[region]] <- aggr_bin[[region]] %>%
        tibble::rownames_to_column("Bin") %>%
        dplyr::mutate(Chr = factor(sub("^((chr[0-9XYM]+)).*", "\\1", Bin), levels = c(paste0("chr", 1:22), "chrX", "chrY", "chrM"))) %>%
        group_by(Chr) %>%
        summarize(across(where(is.numeric), sum), .groups = "drop") %>%
        tidyr::pivot_longer(
            cols = -Chr,
            names_to = "Sample",
            values_to = "Depth"
        )
    
    aggr_bin[[region]] <- aggr_bin[[region]] %>%
    tibble::rownames_to_column("Bin") %>%
    tidyr::pivot_longer(
        cols = -Bin,
        names_to = "Sample",
        values_to = "Depth"
    )
}

## Create cell dataframes
aggr_cell <- list(
    "on_target" = Reduce(rbind, on_aggr_cell_list),
    "off_target" = Reduce(rbind, off_aggr_cell_list)
)

# Plot QC metrics
plotdir <- snakemake@output[["plotdir"]]
if(!dir.exists(plotdir)) dir.create(plotdir, recursive=TRUE)

for(region in names(aggr_cell)) {
    dpcell <- ggplot(aggr_cell[[region]], aes(x = "", y = norm_depth, fill = Sample)) +
        geom_violin(color = "black") +
        labs(title = "Distribution of normalized depth per cell", y = "", x = "", fill = "") +
        ggpubfigs::theme_big_simple() +
        theme(axis.text.x = element_blank(),
                axis.title.x = element_blank())
    save_path <- file.path(plotdir, glue("{region}_cell_depth_distribution.pdf"))
    ggsave(save_path, dpcell, dpi=300, width=40, units="cm")
}

for(region in names(aggr_bin)) {
    dpb <- ggplot(aggr_bin[[region]], aes(x = 1:nrow(aggr_bin[[region]]), y = log1p(Depth))) +
        geom_line(linewidth = 0.7) +
        labs(x = "Bin", y = "Log1p-normalized depth") +
        ggpubfigs::theme_big_simple() +
        facet_wrap(~Sample, nrow = length(sample_ids)) +
        theme(axis.text.x = element_blank(),
                axis.title.x = element_blank())

    save_path <- file.path(plotdir, glue("{region}_sample_depth_per_bin.pdf"))
    ggsave(save_path, dpb, dpi=300, width=40, units="cm")

    dpb_1000 <- ggplot(aggr_bin[[region]][1:1000,], aes(x = 1:1000, y = Depth)) +
        geom_line(linewidth = 0.7) +
        labs(x = "Bin", y = "Normalized depth") +
        ggpubfigs::theme_big_simple() +
        facet_wrap(~Sample, nrow = length(sample_ids)) +
        theme(axis.text.x = element_blank(),
                axis.title.x = element_blank())

    save_path <- file.path(plotdir, glue("{region}_sample_depth_per_bin.zoomed_first_1000_bins.pdf"))
    ggsave(save_path, dpb_1000, dpi=300, width=40, units="cm")
}

for(region in names(aggr_chr)) {
    dpchr <- ggplot(aggr_chr[[region]], aes(x = Chr, y = log1p(Depth), fill = Sample)) +
        geom_bar(stat = "identity") +
        labs(x = "Chromosome", y = "Log1p-normalized depth", fill = "") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

    save_path <- file.path(plotdir, glue("{region}_sample_depth_per_chr.pdf"))
    ggsave(save_path, dpchr, dpi=300, width=40, units="cm")
}


